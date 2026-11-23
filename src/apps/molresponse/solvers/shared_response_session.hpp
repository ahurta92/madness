#pragma once
#include "shared_response_checkpoint.hpp"
#include "shared_es_checkpoint.hpp"
#include "response_metadata.hpp"
#include <madness/chem/IterationHistory.h>
#include <madness/chem/TimingRecords.h>
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
  std::optional<cp::GenerationRef> guess;
};
struct SessionScope {
  std::string writer_group;
  std::string metadata_shard;
  std::string parent_attempt;
  int subworld_id=-1;
};
// One response request, many independently planned/checkpointed states.
class Session {
  madness::World& world_;
  std::filesystem::path root_;
  cp::CheckpointStore store_;
  GroundSource ground_;
  it::Request request_;
  std::vector<it::Protocol> protocols_;
  std::map<std::string,std::string> replaced_tasks_;
  std::map<std::string,TaskSpec> specifications_;
  std::set<std::string> merged_groups_;
  SessionScope scope_;
  std::string metadata_path() const {
    return analysis_directory_+"/response_metadata"+
        (scope_.metadata_shard.empty()?"":".group"+scope_.metadata_shard)+".json";
  }
  std::string owner_, analysis_directory_;
  it::Protocol final_;
  it::TimingCollector<> timing_;
  it::Verbosity verbosity_;
  json properties_=json::object();
public:
  std::map<std::string,std::unique_ptr<Checkpoint>> states;
  std::map<std::string,std::unique_ptr<SecondOrderNode>> second_order;
  // Excited-state sets are nodes, solved on the main World only.
  std::map<std::string,std::unique_ptr<ExcitedSet>> excited;
  Session(madness::World& world,const std::filesystem::path& root,const GroundSource& ground,const it::Request& request,
          const std::vector<it::Protocol>& protocols,
          const std::vector<TaskSpec>& definitions,
          std::string analysis_directory, it::TimingCollector<> timing,it::Verbosity verbosity,
          SessionScope scope={})
      : world_(world),root_(std::filesystem::absolute(root)),store_(world,root_,{},scope.writer_group),ground_(ground),request_(request),protocols_(protocols),scope_(std::move(scope)),
        analysis_directory_(std::move(analysis_directory)),final_(protocols.back()),
        timing_(std::move(timing)),verbosity_(verbosity) {
    timing_.end(); // reference/input setup began before loading the ground state
    timing_.begin("checkpoint_planning");
    owner_=cp::detail::on_root(world_,[] {
      static std::size_t sequence=0;
      return "ra"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+"-"+std::to_string(++sequence);
    });
    if(!scope_.parent_attempt.empty()) owner_=scope_.parent_attempt+"-g"+std::to_string(scope_.subworld_id)+"-"+owner_;
    for (const auto& task : definitions) {
      cp::records::require(!states.contains(task.id) && !excited.contains(task.id),"ambiguous response task identity");
      if(task.settings.contains("excited_model")) {
        cp::records::require(scope_.parent_attempt.empty(),"excited-state sets are not solved in subworlds");
        excited.emplace(task.id,std::make_unique<ExcitedSet>(world,root_,ground,
            task.checkpoint_id.empty()?task.id:task.checkpoint_id,task.settings,request,protocols,owner_,task.definition,verbosity_));
        continue;
      }
      specifications_[task.id]=task;
      states.emplace(task.id,std::make_unique<Checkpoint>(world,root_,ground,
          task.checkpoint_id.empty()?task.id:task.checkpoint_id,task.settings,request,protocols,owner_,task.definition,task.guess,scope_.writer_group,verbosity_));
    }
    for (const auto& [id,state] : states) if (state->reusable()) record_summary(*state,0);
    for (const auto& [id,set] : excited) if (set->reusable()) record_summary(*set,0);
    timing_.end();
  }
  Checkpoint& bind_derived(const std::string& id,json settings) {
    const auto present=states.find(id);
    if(present!=states.end()) {
      cp::records::require(present->second->scientific_settings().at("frequency")==settings.at("frequency") &&
          present->second->scientific_settings().at("perturbation")==settings.at("perturbation"),
          "distinct requested frequencies share a response task identity");
      return *present->second;
    }
    std::optional<cp::GenerationRef> guess;
    // A replacement may have been published just before interruption. Prefer
    // that exact task; only its explicit lineage identifies an older task to retire.
    const auto current=store_.read_latest("response",id,true);
    if(current) {
      const auto& definition=current->record.at("task").at("definition");
      if(definition.contains("replaces")) {
        const auto previous=definition.at("replaces").at("task_id").get<std::string>();
        if(previous!=id) replaced_tasks_[id]=previous;
      }
    } else for(const auto& entry:store_.task_inventory("response")) {
      if(entry.value("retired",false) || entry.at("kind")!="iterative") continue;
      const auto old_id=entry.at("task_id").get<std::string>();
      const auto saved=store_.read_latest("response",old_id);
      if(!saved) continue;
      const auto& previous=saved->record.at("task").at("definition").at("settings");
      if(previous.value("derived_from",json())==settings.at("derived_from") &&
          previous.at("perturbation")==settings.at("perturbation")) {
        cp::records::require(!guess,"ambiguous saved derived response origin");
        guess=cp::GenerationRef{"response",old_id,saved->generation};
      }
    }
    auto binding=std::make_unique<Checkpoint>(world_,root_,ground_,id,settings,request_,protocols_,owner_,derived_task(),guess,scope_.writer_group,verbosity_);
    specifications_[id]={id,settings,derived_task(),{},guess};
    if(guess && guess->task_id!=id) replaced_tasks_[id]=guess->task_id;
    if(binding->reusable()) record_summary(*binding,0);
    auto& result=*binding;states[id]=std::move(binding);return result;
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
  void record_evaluation(const it::IterationMetrics& metrics,const it::Decision& decision,
                         const std::map<it::Metric,std::string>& metric_ids,
                         const std::vector<it::Diagnostic>& diagnostics={}) {
    cp::detail::on_root(world_,[&] {
      it::append_history(root_/"analysis"/"response"/(metrics.info.attempt_id+".jsonl"),
          it::history_row(metrics,decision,metric_ids,diagnostics));
      return std::string{};
    });
  }
  void publish_timings(const std::string& failure={}) {
    const it::Execution execution{static_cast<std::size_t>(world_.size()),
                                  static_cast<std::size_t>(madness::ThreadPool::size()+1)};
    json record=it::timing_entry(owner_,execution,timing_.record(),failure);
    if(scope_.parent_attempt.empty()) {store_.publish_timings(record);return;}
    record["parent_attempt"]=scope_.parent_attempt;
    record["second_order"]=json::array();
    for(const auto& [id,binding]:second_order) if(binding->reusable())
      record["second_order"].push_back({{"id",id},{"name",binding->name()},{"settings",binding->settings()},{"lineage",binding->lineage()}});
    cp::detail::on_root(world_,[&] {
      const auto directory=root_/"analysis"/"response";
      std::filesystem::create_directories(directory);
      const auto path=directory/(owner_+".timings.json");
      const auto temp=directory/("."+owner_+".timings.tmp");
      if(std::filesystem::is_regular_file(std::filesystem::symlink_status(temp))) std::filesystem::remove(temp);
      cp::detail::write_exclusive(temp,record.dump(2)+"\n");
      std::filesystem::rename(temp,path);return std::string{};
    });
  }
  // A new binding owns the subgroup World. Only serializable specifications
  // and exact native source references cross the workflow boundary.
  std::shared_ptr<Session> for_world(madness::World& world,const std::set<std::string>& ids,int group) const {
    std::vector<TaskSpec> definitions;
    for(const auto& id:ids) {
      auto task=specifications_.at(id);
      const auto& binding=at(id);
      if(binding.has_saved()) task.guess=cp::GenerationRef{"response",
          binding.manifest().record.at("task").at("task_id").get<std::string>(),binding.manifest().generation};
      definitions.push_back(std::move(task));
    }
    it::TimingCollector<> timing;
    timing.begin("subworld_work",{{"clock_scope","subworld_driver_rank"},{"subworld_id",group},
        {"mpi_ranks",world.size()},{"threads_per_rank",madness::ThreadPool::size()+1}});
    timing.begin("response");timing.begin("binding_setup");
    return std::make_shared<Session>(world,root_,ground_,request_,protocols_,definitions,analysis_directory_,std::move(timing),verbosity_,
        SessionScope{owner_,std::to_string(group),owner_,group});
  }
  void collect_subworlds() {
    for(auto& [id,binding]:states) binding->refresh();
    const auto records=json::parse(cp::detail::on_root(world_,[&] {
      json result=json::array();const auto directory=root_/"analysis"/"response";
      if(std::filesystem::is_directory(directory)) for(const auto& entry:std::filesystem::directory_iterator(directory)) {
        const auto name=entry.path().filename().string();
        if(name.starts_with(owner_+"-g") && name.ends_with(".timings.json") && !merged_groups_.contains(name))
          result.push_back({{"name",name},{"record",cp::detail::read_json(entry.path())}});
      }
      return result.dump();
    }));
    for(const auto& entry:records) {
      const auto& record=entry.at("record");
      cp::records::require(record.at("parent_attempt")==owner_,"subworld timing belongs to another invocation");
      timing_.append(record.at("timing"));
      for(const auto& node:record.at("second_order")) {
        const auto id=node.at("id").get<std::string>();
        second_order[id]=std::make_unique<SecondOrderNode>(world_,root_,id,node.at("name"),node.at("settings"),node.at("lineage"));
      }
      merged_groups_.insert(entry.at("name").get<std::string>());
    }
  }
  void set_properties(const json& properties) {properties_.update(properties);}
  const json& properties() const {return properties_;}
  template<class Work> decltype(auto) with_failure_timing(Work&& work) {
    try {return work();}
    catch(const std::exception& error) {
      try {collect_subworlds();} catch(...) {}
      try {publish_timings(error.what());} catch(...) {}
      throw; // preserve the computation/publication failure even if timing I/O also fails
    }
  }
  Checkpoint& at(const std::string& id) const {return *states.at(id);}
  ExcitedSet& excited_at(const std::string& id) const {return *excited.at(id);}
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
  void record_summary(const ExcitedSet& set,std::size_t evaluations) const {
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
      meta.save();return std::string{};
    });
  }
  void record_summary(const Checkpoint& state,std::size_t evaluations) const {
    const auto frame=cp::iteration_metrics_from_manifest(state.manifest());
    const auto& physical=state.scientific_settings();
    cp::detail::on_root(world_,[&] {
      auto meta=ResponseMetadata::load_or_create(metadata_path());
      const auto protocol=frame.info.protocol;
      const auto key=protocol_key(protocol.thresh,protocol.k);
      meta.set_protocol(key,protocol.thresh,protocol.k,0);
      meta.set_fd_state(physical.at("perturbation"),key,ResponseMetadata::freq_key(physical.at("frequency")),
          {{"converged",state.complete()},{"frequency",physical.at("frequency")},{"iter",frame.info.iteration},
           {"bsh_residual",frame.states.front().metrics.at(it::Metric::bsh_residual).value.value_or(0.0)},
           {"diverged",false},{"stalled",false},{"shared_generation",state.manifest().generation},
           {"shared_evaluations",evaluations},{"shared_reused",state.reusable()}});
      meta.save();return std::string{};
    });
  }
  void report() {
    json sources=json::object();
    std::map<std::string,json> groups;
    for (const auto& [id,state] : states) {
      if (!state->has_saved()) continue;
      const auto frame=cp::iteration_metrics_from_manifest(state->manifest());
      const auto protocol=frame.info.protocol;
      auto& group=groups[protocol_key(protocol.thresh,protocol.k)];
      group["protocol"]=cp::records::protocol_record(protocol);
      for (const auto& subject:frame.states) {
        const auto public_id=frame.states.size()==1 ? id : id+":"+subject.id;
        sources[public_id]=state->source_record(subject.id);
        group["states"][public_id]=state->measured_values(subject.id);
      }
    }
    // Excited-state sets are nodes: their measured values are reported from
    // the node description, and the source map lists no source for them.
    for (const auto& [id,set] : excited) {
      if (!set->has_saved()) continue;
      const auto frame=set->saved_metrics();
      auto& group=groups[protocol_key(frame.info.protocol.thresh,frame.info.protocol.k)];
      group["protocol"]=cp::records::protocol_record(frame.info.protocol);
      for (const auto& subject:frame.states) group["states"][id+":"+subject.id]=measured_values(frame,subject.id);
    }
    // Second-order intermediates are nodes, not store tasks: the source map
    // lists none (removed in part 4, with the source map).
    store_.publish_source_map({{"version",2},{"calculation_id","response"},{"attempt_id",owner_},
                               {"states",sources},{"tasks",json::array()}});
    auto result=json::parse(cp::detail::on_root(world_,[&]{return cp::detail::read_json(root_/"results.json").dump();}));
    const auto prior_attempt=result.at("attempt_id");
    json calculations=json::array(),metrics=json::array();
    for (auto entry:result.at("calculations")) if(entry.at("calculation_id")!="response") {
      if (!entry.contains("attempt_id")) entry["attempt_id"]=prior_attempt;
      calculations.push_back(std::move(entry));
    }
    for (const auto& [key,group]:groups) metrics.push_back(group);
    json second_order_nodes=json::array();
    for(const auto& [id,binding]:second_order) if(binding->reusable())
      second_order_nodes.push_back({{"id",id},{"node",binding->folder()}});
    json excited_nodes=json::array();
    for(const auto& [id,set]:excited) if(set->has_saved())
      excited_nodes.push_back({{"id",id},{"node",set->folder()}});
    calculations.push_back({{"calculation_id","response"},{"attempt_id",owner_},{"application","molresponse"},
       {"protocol",cp::records::protocol_record(final_)},
       {"properties",properties_},
       {"state_metrics",metrics},
       {"second_order",second_order_nodes},
       {"excited",excited_nodes}});
    result["attempt_id"]=owner_;result["calculations"]=std::move(calculations);
    const auto warning=store_.publish_results(result);
    if (!warning.empty() && world_.rank()==0) madness::print(warning);
    for(const auto& [current,previous]:replaced_tasks_)
      if(!states.contains(previous) && states.at(current)->has_saved() &&
          states.at(current)->manifest().record.at("task").at("task_id")==current)
        store_.retire_task("response",previous);
    const auto cleanup=store_.cleanup();
    if(!cleanup.empty() && world_.rank()==0) madness::print(cleanup);
  }
};
} // namespace molresponse_v3::response_shared
