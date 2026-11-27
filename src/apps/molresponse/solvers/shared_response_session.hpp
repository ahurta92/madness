#pragma once
#include "shared_response_checkpoint.hpp"
#include "shared_es_checkpoint.hpp"
#include "response_metadata.hpp"
#include <madness/chem/IterationHistory.h>
#include <madness/chem/TimingRecords.h>
#include <madness/chem/RunFiles.h>
#include "../ResponseProtocol.hpp"

namespace molresponse_v3::response_shared {
// Root-only JSON work must finish collectively before any subsequent archive or solver call.
template<class Update>
void update_response_metadata(madness::World& world,const std::string& path,Update&& update) {
  cp::detail::on_root(world,[&] {
    auto metadata=ResponseMetadata::load_or_create(path);
    update(metadata);metadata.save();return std::string{};
  });
}
struct TaskSpec {
  std::string id;
  json settings;
  TaskDefinition definition=driven_task();
  std::string checkpoint_id;
  std::optional<std::string> start_from; // another driven node this one may start from
};
// A driven state's node name: its perturbation and frequency, joined by a
// hyphen, as vbc_node_name joins its parts. A perturbation description never
// holds a hyphen and the frequency part starts with 'f', so two different
// states never share a name.
inline std::string driven_node_name(const json& settings) {
  return settings.at("perturbation").get<std::string>()+"-"+
         ResponseMetadata::freq_key(settings.at("frequency").get<double>());
}
// One saved driven node as the choice of a starting state sees it.
struct StartCandidate {
  std::string name;
  double frequency=0.0;
  json ground;        // the ground state's identity in its settings
  json model;         // start_model of its settings, without the ground identity
  bool settled=true;  // completely saved, and written by no process group again in this run
};
// The node a state at `frequency` starts from: of the settled candidates with
// this ground identity and model, the one at the nearest frequency; of two
// equally near, the lower frequency. Nothing when none can serve.
inline std::optional<std::string> nearest_start(double frequency,const json& ground,const json& model,
                                                const std::vector<StartCandidate>& candidates) {
  const StartCandidate* best=nullptr;
  for(const auto& candidate:candidates) {
    if(!candidate.settled || candidate.ground!=ground || candidate.model!=model) continue;
    if(best==nullptr) {best=&candidate;continue;}
    const double distance=std::abs(candidate.frequency-frequency),best_distance=std::abs(best->frequency-frequency);
    const double tie=1e-12*std::max({1.0,std::abs(candidate.frequency),std::abs(best->frequency)});
    const bool nearer=distance<best_distance-tie;
    const bool equally_near=!nearer && distance<=best_distance+tie;
    if(nearer || (equally_near && (candidate.frequency<best->frequency ||
                                   (candidate.frequency==best->frequency && candidate.name<best->name))))
      best=&candidate;
  }
  if(best==nullptr) return std::nullopt;
  return best->name;
}
// The two-photon nodes of excited-state set `bundle` that no current root
// needs: driven nodes (name, settings) whose settings name the set as their
// origin and that are not in `needed`. An ordinary driven state names no
// origin, so it is never on the list.
inline std::vector<std::string> two_photon_removals(const std::string& bundle,const std::set<std::string>& needed,
                                                    const std::vector<std::pair<std::string,json>>& nodes) {
  std::vector<std::string> result;
  for(const auto& [name,settings]:nodes) {
    const auto origin=settings.find("derived_from");
    if(origin==settings.end() || !origin->is_object() || origin->value("bundle",std::string())!=bundle) continue;
    if(!needed.contains(name)) result.push_back(name);
  }
  return result;
}
// The ground state's identity in a driven state's settings.
inline json ground_identity(const json& settings) {
  return {{"ground_settings",settings.value("ground_settings",json())},{"ground_saved",settings.value("ground_saved",json())}};
}
// Every completely saved driven node under <root>/states/response, by name.
// Run on the root process only; unreadable nodes are skipped.
inline std::vector<std::pair<std::string,it::NodeDescription>> saved_driven_nodes(const std::filesystem::path& root) {
  std::vector<std::pair<std::string,it::NodeDescription>> result;
  const auto directory=root/"states"/"response";
  if(!std::filesystem::is_directory(directory)) return result;
  for(const auto& entry:std::filesystem::directory_iterator(directory)) {
    if(!completely_saved(entry.path())) continue;
    try {
      std::ifstream file(entry.path()/"description.json");
      auto d=it::read_description(json::parse(file));
      if(d.kind=="driven") result.emplace_back(entry.path().filename().string(),std::move(d));
    } catch(const std::exception&) {}
  }
  std::sort(result.begin(),result.end(),[](const auto& a,const auto& b){return a.first<b.first;});
  return result;
}
// A session of one process group (for_world) has a group number; the main
// session has none.
struct SessionScope {
  std::string metadata_shard;
  int group=-1;
};
// A group's timing record, which the main session merges and then removes.
inline std::filesystem::path group_timing_file(const std::filesystem::path& root,int group) {
  return root/"analysis"/"response"/("group-"+std::to_string(group)+".timings.json");
}
// One response request, many independently planned/checkpointed states.
class Session {
  madness::World& world_;
  std::filesystem::path root_;
  GroundSource ground_;
  it::Request request_;
  std::vector<it::Protocol> protocols_;
  std::map<std::string,TaskSpec> specifications_;
  SessionScope scope_;
  std::string metadata_path() const {
    return analysis_directory_+"/response_metadata"+
        (scope_.metadata_shard.empty()?"":".group"+scope_.metadata_shard)+".json";
  }
  std::string analysis_directory_;
  it::Protocol final_;
  it::TimingCollector<> timing_;
  std::optional<std::size_t> timing_part_; // this run's part of the launch's timing entry
  it::Verbosity verbosity_;
  json properties_=json::object();
  json run_=json::object();
public:
  std::map<std::string,std::unique_ptr<DrivenState>> states;
  std::map<std::string,std::unique_ptr<SecondOrderNode>> second_order;
  // Excited-state sets are nodes, solved on the main World only.
  std::map<std::string,std::unique_ptr<ExcitedSet>> excited;
  Session(madness::World& world,const std::filesystem::path& root,const GroundSource& ground,const it::Request& request,
          const std::vector<it::Protocol>& protocols,
          const std::vector<TaskSpec>& definitions,
          std::string analysis_directory, it::TimingCollector<> timing,it::Verbosity verbosity,
          SessionScope scope={})
      : world_(world),root_(std::filesystem::absolute(root)),ground_(ground),request_(request),protocols_(protocols),scope_(std::move(scope)),
        analysis_directory_(std::move(analysis_directory)),final_(protocols.back()),
        timing_(std::move(timing)),verbosity_(verbosity) {
    timing_.end(); // reference/input setup began before loading the ground state
    timing_.begin("checkpoint_planning");
    // Group timing files left by an earlier launch that failed are not this launch's.
    if(scope_.group<0) cp::detail::on_root(world_,[&] {
      const auto directory=root_/"analysis"/"response";
      if(std::filesystem::is_directory(directory))
        for(const auto& entry:std::filesystem::directory_iterator(directory)) {
          const auto name=entry.path().filename().string();
          if(name.starts_with("group-") && name.ends_with(".timings.json")) std::filesystem::remove(entry.path());
        }
      return std::string{};
    });
    for (const auto& task : definitions) {
      cp::records::require(!states.contains(task.id) && !excited.contains(task.id),"ambiguous response task identity");
      if(task.settings.contains("excited_model")) {
        cp::records::require(scope_.group<0,"excited-state sets are not solved in subworlds");
        excited.emplace(task.id,std::make_unique<ExcitedSet>(world,root_,ground,
            task.checkpoint_id.empty()?task.id:task.checkpoint_id,task.settings,request,protocols,task.definition,verbosity_));
        continue;
      }
      specifications_[task.id]=task;
      states.emplace(task.id,bind(task));
    }
    for (const auto& [id,state] : states) if (state->reusable()) record_summary(*state,0);
    for (const auto& [id,set] : excited) if (set->reusable()) record_summary(*set,0);
    timing_.end();
  }
  // A process group bound by for_world only solves or reads; the launch that
  // made the group already wrote each node's launch row.
  std::unique_ptr<DrivenState> bind(const TaskSpec& task,bool launch_header=true) const {
    return std::make_unique<DrivenState>(world_,root_,ground_,task.checkpoint_id.empty()?task.id:task.checkpoint_id,
        driven_node_name(task.settings),task.settings,request_,protocols_,task.definition,task.start_from,
        launch_header && scope_.group<0,verbosity_);
  }
  // Whether a saved driven state meets the current request at the final
  // level, so that its owner reuses it and no process group writes it again.
  [[nodiscard]] bool meets_request(const it::NodeDescription& d) const {
    if(d.protocol!=final_) return false;
    const auto definition=driven_task();
    const auto settings=it::resolve(request_,definition.model(final_),final_,true,definition.application);
    const auto context=json{{"settings",d.settings},{"inputs",json::array()}}.dump();
    const std::string task="nearby";
    const auto subjects=cp::records::subjects(definition.record);
    const it::ExpectedTask expected{task,context,{subjects.begin(),subjects.end()}};
    try {
      return it::saved_completion(settings,expected,it::stored_metrics(d,{task,d.protocol,0,context}),true).has_value();
    } catch(const std::exception&) {return false;}
  }
  // Whether the saved node `name` (description `d`, read now) meets the
  // current request, as the restart would judge it. `nodes` are all the saved
  // nodes, by folder name. A state or set of this request is judged by its
  // binding (meets_request). A driven node this request has not bound (a
  // two-photon state whose root is not expanded yet) is judged against its own
  // settings at the final level, with the current ground state. A VBC source
  // meets it when it was built at the final level from the current saves of
  // its two inputs, and both meet it. Nothing else does. Root process only.
  [[nodiscard]] bool node_meets_request(const std::string& name,const it::NodeDescription& d,
                                        const std::map<std::string,it::NodeDescription>& nodes,int depth=0) const {
    if(d.kind=="driven") {
      for(const auto& [id,state]:states) if(state->name()==name) return state->meets_request(d);
      const json ground{{"ground_settings",ground_.identity.at("settings")},
                        {"ground_saved",{{"protocol",ground_.identity.at("protocol")},{"saved_at",ground_.identity.at("saved_at")}}}};
      return d.settings.is_object() && ground_identity(d.settings)==ground && meets_request(d);
    }
    if(d.kind=="excited_set") {
      for(const auto& [id,set]:excited)
        if(std::filesystem::path(set->folder()).filename().string()==name) return set->meets_request(d);
      return false;
    }
    if(d.kind!="second_order" || depth>0 || d.settings!=vbc_source_settings(final_) || d.protocol!=final_ ||
       !d.lineage.is_array() || d.lineage.empty())
      return false;
    try {
      for(const auto& input:d.lineage) {
        const auto found=nodes.find(input.at("name").get<std::string>());
        if(found==nodes.end()) return false;
        const auto& saved=found->second;
        if(input.at("protocol")!=cp::records::protocol_record(saved.protocol) ||
           input.at("identity").at("saved_at")!=saved.saved_at ||
           !node_meets_request(found->first,saved,nodes,depth+1))
          return false;
      }
      return true;
    } catch(const std::exception&) {return false;}
  }
  // An ordinary driven state with no saved state of its own starts from the
  // saved state of the same perturbation, ground state and model at the
  // nearest frequency, chosen when it starts; that node is only read. The
  // main session solves one state at a time, so any completely saved node
  // can serve. A process group may start only from a state of its own or one
  // that already meets the request: no other group writes those again, so
  // none is replaced while it is read. Collective.
  DrivenState& start_nearby(const std::string& id) {
    auto& state=*states.at(id);
    const auto& task=specifications_.at(id);
    if(state.has_saved() || task.settings.contains("derived_from")) return state;
    const auto& settings=state.scientific_settings();
    const auto ground=ground_identity(settings);
    auto model=start_model(settings);model.erase("ground_settings");model.erase("ground_saved");
    std::set<std::string> own;
    for(const auto& [other,spec]:specifications_) own.insert(driven_node_name(spec.settings));
    const auto chosen=cp::detail::on_root(world_,[&] {
      std::vector<StartCandidate> candidates;
      for(const auto& [name,d]:saved_driven_nodes(root_)) {
        if(name==state.name()) continue;
        auto other=start_model(d.settings);other.erase("ground_settings");other.erase("ground_saved");
        candidates.push_back({name,d.settings.value("frequency",0.0),ground_identity(d.settings),std::move(other),
                              scope_.group<0 || own.contains(name) || meets_request(d)});
      }
      return nearest_start(settings.at("frequency").get<double>(),ground,model,candidates).value_or(std::string{});
    });
    if(chosen.empty()) return state;
    auto from=task;from.start_from=chosen;
    states[id]=bind(from,false);
    return *states.at(id);
  }
  // A state derived from an excited root keeps its root while the root's
  // energy, so its frequency and node name, moves. With no node of its own it
  // starts from the node of the same origin and perturbation at the nearest
  // frequency (of two equally near, the lower); that node is only read.
  DrivenState& bind_derived(const std::string& id,json settings) {
    const auto present=states.find(id);
    if(present!=states.end()) {
      cp::records::require(present->second->scientific_settings().at("frequency")==settings.at("frequency") &&
          present->second->scientific_settings().at("perturbation")==settings.at("perturbation"),
          "distinct requested frequencies share a response task identity");
      return *present->second;
    }
    const auto name=driven_node_name(settings);
    // The ground identity is not compared: a moved root follows a ground
    // state saved again.
    const json origin{{"derived_from",settings.at("derived_from")},{"perturbation",settings.at("perturbation")}};
    const auto nearby=cp::detail::on_root(world_,[&] {
      std::vector<StartCandidate> candidates;
      for(const auto& [other,d]:saved_driven_nodes(root_)) {
        if(other==name) continue;
        candidates.push_back({other,d.settings.value("frequency",0.0),json(),
            {{"derived_from",d.settings.value("derived_from",json())},{"perturbation",d.settings.value("perturbation",json())}}});
      }
      return nearest_start(settings.at("frequency").get<double>(),json(),origin,candidates).value_or(std::string{});
    });
    TaskSpec task{id,settings,derived_task(),{},nearby.empty()?std::nullopt:std::optional<std::string>(nearby)};
    auto binding=bind(task);
    specifications_[id]=task;
    if(binding->reusable()) record_summary(*binding,0);
    auto& result=*binding;states[id]=std::move(binding);return result;
  }
  // Once every two-photon state the current roots need has met the request,
  // the two-photon nodes of the same excited-state set that no current root
  // needs are removed. Ordinary driven states are never removed. Main session
  // only, after every process group has finished. Collective. Returns the
  // names removed.
  std::vector<std::string> prune_two_photon() {
    cp::records::require(scope_.group<0,"only the main session removes two-photon states");
    std::map<std::string,std::set<std::string>> needed;
    std::set<std::string> unfinished;
    for(const auto& [id,state]:states) {
      const auto& settings=state->scientific_settings();
      if(!settings.contains("derived_from")) continue;
      const auto bundle=settings.at("derived_from").at("bundle").get<std::string>();
      needed[bundle].insert(state->name());
      if(!state->complete() || !state->has_node()) unfinished.insert(bundle);
    }
    const auto removed=json::parse(cp::detail::on_root(world_,[&] {
      json result=json::array();
      std::vector<std::pair<std::string,json>> nodes;
      for(const auto& [name,d]:saved_driven_nodes(root_)) nodes.emplace_back(name,d.settings);
      std::vector<std::string> names_removed;
      for(const auto& [bundle,names]:needed) {
        if(unfinished.contains(bundle)) continue;
        for(const auto& name:two_photon_removals(bundle,names,nodes)) names_removed.push_back(name);
      }
      if(names_removed.empty()) return result.dump();
      // The metadata stops naming a state before its folder goes, so a run
      // stopped in between leaves only a folder, which the next run removes.
      auto meta=ResponseMetadata::load_or_create(metadata_path());
      for(const auto& name:names_removed) meta.drop_fd_node(driven_folder(name).generic_string());
      meta.save();
      for(const auto& name:names_removed) {
        std::filesystem::remove_all(root_/driven_folder(name));
        result.push_back(name);
      }
      return result.dump();
    }));
    return removed.get<std::vector<std::string>>();
  }
  // The lineage names response states, whose settings carry the ground identity.
  SecondOrderNode& bind_second_order(const std::string& id,const std::string& name,const json& settings,const json& lineage) {
    auto binding=std::make_unique<SecondOrderNode>(world_,root_,id,name,settings,lineage);
    auto& result=*binding;second_order[id]=std::move(binding);return result;
  }
  void record_vbc(const std::string& id,const SecondOrderNode& binding) const {
    cp::detail::on_root(world_,[&] {
      auto meta=ResponseMetadata::load_or_create(metadata_path());
      const auto key=protocol_key(final_.thresh,final_.k);
      meta.set_protocol(key,final_.thresh,final_.k,0);
      meta.set_vbc_state(id,key,{{"converged",true},{"diverged",false},
          {"node",binding.folder()}});
      meta.save();return std::string{};
    });
  }
  it::TimingCollector<>& timing() {return timing_;}
  void record_iteration(const it::IterationInfo& info,const it::IterationTiming& measured) {
    json children=json::array();
    for (const auto& [name,seconds]:measured.steps) children.push_back({{"name",name},{"seconds",seconds}});
    timing_.append({{"name","iteration"},{"iteration",info.iteration},
        {"seconds",measured.seconds},{"children",children}});
  }
  void publish_timings(const std::string& failure={}) {
    const it::Execution execution{static_cast<std::size_t>(world_.size()),
                                  static_cast<std::size_t>(madness::ThreadPool::size()+1)};
    json record=it::timing_entry(it::launch_started_at(),execution,timing_.record(),failure);
    if(scope_.group<0) {timing_part_=it::write_timing_entry(world_,root_,record,timing_part_);return;}
    record["second_order"]=json::array();
    for(const auto& [id,binding]:second_order) if(binding->reusable())
      record["second_order"].push_back({{"id",id},{"name",binding->name()},{"settings",binding->settings()},{"lineage",binding->lineage()}});
    // Each group writes only its own file.
    cp::detail::on_root(world_,[&] {
      const auto path=group_timing_file(root_,scope_.group);
      std::filesystem::create_directories(path.parent_path());
      cp::detail::replace_json(path,record);return std::string{};
    });
  }
  // A new binding owns the subgroup World. Only serializable specifications
  // and exact native source references cross the workflow boundary.
  std::shared_ptr<Session> for_world(madness::World& world,const std::set<std::string>& ids,int group) const {
    std::vector<TaskSpec> definitions;
    definitions.reserve(ids.size());
    for(const auto& id:ids) definitions.push_back(specifications_.at(id));
    it::TimingCollector<> timing;
    timing.begin("subworld_work",{{"clock_scope","subworld_driver_rank"},{"subworld_id",group},
        {"mpi_ranks",world.size()},{"threads_per_rank",madness::ThreadPool::size()+1}});
    timing.begin("response");timing.begin("binding_setup");
    return std::make_shared<Session>(world,root_,ground_,request_,protocols_,definitions,analysis_directory_,std::move(timing),verbosity_,
        SessionScope{std::to_string(group),group});
  }
  void collect_subworlds() {
    for(auto& [id,binding]:states) binding->refresh();
    const auto records=json::parse(cp::detail::on_root(world_,[&] {
      json result=json::array();const auto directory=root_/"analysis"/"response";
      if(std::filesystem::is_directory(directory)) for(const auto& entry:std::filesystem::directory_iterator(directory)) {
        const auto name=entry.path().filename().string();
        if(name.starts_with("group-") && name.ends_with(".timings.json")) {
          result.push_back({{"name",name},{"record",cp::detail::read_json(entry.path())}});
          std::filesystem::remove(entry.path()); // merged once, here
        }
      }
      return result.dump();
    }));
    for(const auto& entry:records) {
      const auto& record=entry.at("record");
      timing_.append(record.at("timing"));
      for(const auto& node:record.at("second_order")) {
        const auto id=node.at("id").get<std::string>();
        second_order[id]=std::make_unique<SecondOrderNode>(world_,root_,id,node.at("name"),node.at("settings"),node.at("lineage"));
      }
    }
  }
  void set_properties(const json& properties) {properties_.update(properties);}
  // The scheduler's end of this run, for the response entry of results.json.
  void set_run(const json& run) {run_=run;}
  // A task's node name: its folder name under states/response/. results.json
  // refers to every saved state by it.
  [[nodiscard]] std::string node_name(const std::string& id) const {
    if(const auto state=states.find(id);state!=states.end()) return state->second->name();
    if(const auto set=excited.find(id);set!=excited.end())
      return std::filesystem::path(set->second->folder()).filename().string();
    if(const auto node=second_order.find(id);node!=second_order.end()) return node->second->name();
    throw std::runtime_error("no response node for task "+id);
  }
  [[nodiscard]] bool knows(const std::string& id) const {
    return states.contains(id) || excited.contains(id) || second_order.contains(id);
  }
  // Add a task's node, or one subject of it (an excited root), to the states
  // a value was computed from, once, if this request binds the task. The
  // results file turns them into the value's convergence.
  void name_state(json& value,const std::string& id,const std::string& subject={}) const {
    if(!knows(id)) return;
    it::add_results_source(value,node_name(id),subject);
  }
  // The final-level target of each metric, per kind of state this request
  // solves. One kind, or kinds that agree, give one object; otherwise each
  // kind is named.
  [[nodiscard]] json targets() const {
    json kinds=json::object();
    if(!states.empty()) {
      const auto definition=driven_task();
      kinds["driven_states"]=it::results_targets(it::resolve(request_,definition.model(final_),final_,true,definition.application));
    }
    if(!excited.empty()) {
      const auto& definition=excited.begin()->second->definition();
      kinds["excited_states"]=it::results_targets(it::resolve(request_,definition.model(final_),final_,true,definition.application));
    }
    if(kinds.empty()) return json::object();
    const auto& first=kinds.begin().value();
    for(const auto& [kind,targets]:kinds.items()) if(targets!=first) return kinds;
    return first;
  }
  const json& properties() const {return properties_;}
  template<class Work> decltype(auto) with_failure_timing(Work&& work) {
    try {return work();}
    catch(const std::exception& error) {
      try {collect_subworlds();} catch(...) {}
      try {publish_timings(error.what());} catch(...) {}
      throw; // preserve the computation/publication failure even if timing I/O also fails
    }
  }
  DrivenState& at(const std::string& id) const {return *states.at(id);}
  ExcitedSet& excited_at(const std::string& id) const {return *excited.at(id);}
  // The run folder: states/, seeds/ and analysis/ are under it.
  const std::filesystem::path& root() const {return root_;}
  // The ground state every state of the session is solved on.
  [[nodiscard]] const GroundSource& ground() const {return ground_;}
  // Whether the named driven state or excited-state set has met the request.
  bool complete(const std::string& id) const {
    if(const auto set=excited.find(id);set!=excited.end()) return set->second->complete();
    const auto state=states.find(id);
    return state!=states.end() && state->second->complete();
  }
  bool complete() const {
    return std::all_of(states.begin(),states.end(),[](const auto& entry){return entry.second->complete();}) &&
           std::all_of(excited.begin(),excited.end(),[](const auto& entry){return entry.second->complete();});
  }
  // The metadata copy of an excited set, with its seed guard when the saved
  // roots carry one (`guard_seconds` when this launch measured it).
  void record_summary(const ExcitedSet& set,std::size_t evaluations,
                      std::optional<double> guard_seconds={}) const {
    const auto& d=set.description();
    const auto& physical=set.scientific_settings();
    cp::detail::on_root(world_,[&] {
      auto meta=ResponseMetadata::load_or_create(metadata_path());
      const auto key=protocol_key(d.protocol.thresh,d.protocol.k);
      meta.set_protocol(key,d.protocol.thresh,d.protocol.k,0);
      const auto& roots=d.details.at("roots");
      meta.set_es_bundle(key,{{"converged",set.complete()},{"diverged",false},{"stalled",false},
        {"type",physical.at("excited_model")},{"shell","closed_shell"},{"n_roots",roots.size()},
        {"roots",roots},{"iter",d.updates_at_protocol},{"node",set.folder()},
        {"shared_evaluations",evaluations},{"shared_reused",set.reusable()}});
      if(d.details.contains("seed_guard")) {
        auto guard=d.details.at("seed_guard");
        if(guard_seconds) guard["wall_s"]=*guard_seconds;
        meta.set_es_seed_guard(key,guard);
      }
      meta.save();return std::string{};
    });
  }
  void record_summary(const DrivenState& state,std::size_t evaluations) const {
    const auto& d=state.description();
    const auto frame=state.saved_metrics();
    const auto& physical=state.scientific_settings();
    cp::detail::on_root(world_,[&] {
      auto meta=ResponseMetadata::load_or_create(metadata_path());
      const auto key=protocol_key(d.protocol.thresh,d.protocol.k);
      meta.set_protocol(key,d.protocol.thresh,d.protocol.k,0);
      meta.set_fd_state(physical.at("perturbation"),key,ResponseMetadata::freq_key(physical.at("frequency")),
          {{"converged",state.complete()},{"frequency",physical.at("frequency")},{"iter",d.updates_at_protocol},
           {"bsh_residual",frame.states.front().metrics.at(it::Metric::bsh_residual).value.value_or(0.0)},
           {"diverged",false},{"stalled",false},{"node",state.folder()},
           {"shared_evaluations",evaluations},{"shared_reused",state.reusable()}});
      meta.save();return std::string{};
    });
  }
  void report() {
    cp::records::require(scope_.group<0,"only the main session writes the results file");
    // Each saved node's measured values, keyed by its folder name: the table
    // each property's convergence is taken from, with the number of updates
    // the node took at its saved level (one number for an excited set).
    json metrics=json::object();
    const auto add=[&](const std::string& name,const it::IterationMetrics& frame,std::size_t updates) {
      json subjects=json::object();
      for (const auto& subject:frame.states) subjects[subject.id]=measured_values(frame,subject.id);
      metrics[name]={{"protocol",cp::records::protocol_record(frame.info.protocol)},{"updates",updates},
                     {"subjects",subjects}};
    };
    for (const auto& [id,state] : states)
      if (state->has_node()) add(state->name(),state->saved_metrics(),state->description().updates_at_protocol);
    for (const auto& [id,set] : excited)
      if (set->has_saved()) add(node_name(id),set->saved_metrics(),set->description().updates_at_protocol);
    // A second-order source is built, not iterated: it has no measured values.
    for (const auto& [id,node] : second_order) if (node->reusable())
      metrics[node->name()]={{"protocol",cp::records::protocol_record(node->description().protocol)},
                             {"subjects",json::object()}};
    // Only the main session reports; it replaces only response's entry.
    it::write_results_entry(world_,root_,{{"calculation_id","response"},{"application","molresponse"},
       {"protocol",cp::records::protocol_record(final_)},
       {"targets",targets()},
       {"properties",properties_},
       {"run",run_with_unconverged()},
       {"state_metrics",metrics}});
  }
  // The scheduler's end of the run, plus `unconverged`: the folders of the
  // states that ended without meeting the request. A run that left any is not
  // `complete`; stalled and dropped work keep their own lists and stop reasons.
  [[nodiscard]] json run_with_unconverged() const {
    json run=run_;
    const auto stalled=run.value("stalled",json::array());
    const auto listed=[&](const std::string& id,const std::string& name) {
      return std::find(stalled.begin(),stalled.end(),json(id))!=stalled.end() ||
             std::find(stalled.begin(),stalled.end(),json(name))!=stalled.end();
    };
    json unconverged=json::array();
    for (const auto& [id,state] : states)
      if (!state->complete() && !listed(id,node_name(id))) unconverged.push_back(node_name(id));
    for (const auto& [id,set] : excited)
      if (!set->complete() && !listed(id,node_name(id))) unconverged.push_back(node_name(id));
    run["unconverged"]=unconverged;
    if (!unconverged.empty() && run.value("stop_reason",std::string())=="complete") run["stop_reason"]="unconverged";
    return run;
  }
};
} // namespace molresponse_v3::response_shared
