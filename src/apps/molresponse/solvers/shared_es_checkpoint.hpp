#pragma once
#include "shared_response_checkpoint.hpp"
#include "es_solver.hpp"

namespace molresponse_v3::response_shared {
inline std::map<it::Metric,std::string> excited_metric_set() {
  return {{it::Metric::bsh_residual,"es-bsh-v1"},{it::Metric::density_change,"es-density-v1"},
          {it::Metric::energy_change,"es-energy-v1"},{it::Metric::kinetic_change,"es-kinetic-v1"}};
}
inline std::map<it::Metric,std::string> excited_metric_labels() {
  auto labels=metric_labels();
  labels[it::Metric::energy_change]="Excitation-energy change (Ha)";
  return labels;
}
inline TaskDefinition excited_task(std::size_t roots) {
  auto ids=excited_metric_set();
  json definitions_es=json::object();
  const auto fd=definitions();
  for(const auto& [metric,id]:ids) {
    auto definition=fd.at(metric_set().at(metric));
    if(metric==it::Metric::energy_change) {
      definition["meaning"]="absolute change of the prepared root excitation energy";
      definition["units"]="hartree";definition["fields"]="prepared root Ritz energy";
    }
    if(metric!=it::Metric::bsh_residual)
      definition["comparison"]="preceding prepared root with unambiguous positive-overlap match; density phase aligned";
    definitions_es[id]=std::move(definition);
  }
  json subjects=json::array();
  for(std::size_t root=0;root<roots;++root) subjects.push_back(make_root_id(static_cast<int>(root)));
  return {{{"subjects",subjects}},ids,definitions_es,
      [](it::Protocol protocol) {
        // Excited states also judge the change in excitation energy, and keep
        // the plain rule for every metric.
        auto metrics=model(protocol);
        for(auto& [metric,definition]:metrics) definition.rule=it::LimitRule::plain;
        metrics[it::Metric::energy_change]={it::Support::supported,true,5e-4,it::LimitRule::plain,{}};
        return metrics;
      },
      [](const json& previous,const json& current) {
        cp::records::require(previous.at("layout")==current.at("layout"),
            "unsupported excited checkpoint block layout");
      },"excited states",excited_metric_labels()};
}

// The set of excited-state roots of one model, solved together, is one node:
// <run>/states/response/excited-<model>/, with description.json, history.jsonl
// and archive/root_<slot> (one native archive per root). Every save replaces
// the whole set.
inline std::filesystem::path excited_folder(const std::string& model) {
  cp::records::require(!model.empty() && model.find_first_not_of(
      "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")==std::string::npos,
      "invalid excited-state model name: "+model);
  return std::filesystem::path("states")/"response"/("excited-"+model);
}

// The restart decision follows SCFCheckpoint: the stored metrics become
// measured values (stored_metrics) and saved_completion decides reuse. Saved
// roots with other settings (a ground state saved again among them) are only
// a starting guess. More requested roots than saved keep every saved root and
// resume from them; fewer are refused, as the store refused them.
class ExcitedSet {
  madness::World& world_;
  std::filesystem::path folder_;
  std::string task_,context_;
  json settings_,lineage_;
  TaskDefinition definition_;
  it::Request request_;
  std::vector<it::Protocol> protocols_;
  it::RestartPlan restart_plan_;
  it::ResolvedSettings settings;
  bool final_level_=true;
  std::optional<it::NodeDescription> saved_;
  bool reusable_=false,complete_=false;

  [[nodiscard]] it::IterationMetrics stored(const it::NodeDescription& d) const {
    return it::stored_metrics(d,{task_,d.protocol,0,context_});
  }
  [[nodiscard]] static std::vector<std::string> saved_subjects(const it::NodeDescription& d) {
    std::vector<std::string> ids;
    for(const auto& root:d.details.at("roots")) ids.push_back(root.at("root_id"));
    return ids;
  }
public:
  ExcitedSet(madness::World& world,const std::filesystem::path& root,const GroundSource& ground,
             std::string task,json scientific_settings,const it::Request& request,
             std::vector<it::Protocol> protocols,TaskDefinition definition,
             it::Verbosity verbosity=it::Verbosity::detail)
      : world_(world),task_(std::move(task)),settings_(std::move(scientific_settings)),
        definition_(std::move(definition)),request_(request),protocols_(std::move(protocols)),
        settings(it::resolve(request,definition_.model(protocols_.at(0)),protocols_.at(0),
                             protocols_.size()==1,definition_.application)) {
    // The number of roots is not a setting: a saved set with fewer roots is
    // still the same set.
    settings_.erase("n_roots");
    settings_["ground_settings"]=ground.identity.at("settings");
    settings_["ground_saved"]={{"protocol",ground.identity.at("protocol")},{"saved_at",ground.identity.at("saved_at")}};
    context_=json{{"settings",settings_},{"inputs",json::array()}}.dump();
    lineage_={{"ground_state",ground.identity}};
    folder_=std::filesystem::absolute(root)/excited_folder(settings_.at("excited_model"));
    madness::write_launch_header(world_,folder_);
    saved_=madness::read_node(world_,folder_); // collective; a half-present node is an error
    std::optional<it::SavedTask> saved;
    auto restart_expected=expected();
    if(saved_) {
      const auto& d=*saved_;
      cp::records::require(d.kind=="excited_set" && d.details.at("excited_node_version")==1,
          "unsupported excited-state node: "+folder_.string());
      cp::records::require(madness::scf_model_settings(d.settings.at("ground_settings"))==
                           madness::scf_model_settings(settings_.at("ground_settings")),
          "response ground-state model differs from its saved checkpoint");
      definition_.validate_settings(d.settings,settings_);
      const auto previous=saved_subjects(d);
      const auto wanted=restart_expected.states;
      const std::set<std::string> have(previous.begin(),previous.end()),want(wanted.begin(),wanted.end());
      cp::records::require(!previous.empty() && have.size()==previous.size() &&
          std::includes(want.begin(),want.end(),have.begin(),have.end()),
          "the saved excited-state set holds roots that are not requested; use a new output prefix");
      saved=it::SavedTask{};
      saved->metrics=stored(d);
      saved->protocol_complete=d.level_finished;
      saved->updates_at_protocol=d.updates_at_protocol;
      saved->conversion_supported=true;
      saved->bundle_expanded=previous.size()<wanted.size();
      if(d.settings!=settings_) saved->use=it::SavedEvidenceUse::guess_only;
      // Saved roots that started from imported vectors keep naming that source.
      else if(d.lineage.contains("guess")) lineage_["guess"]=d.lineage.at("guess");
    }
    for(const auto& protocol:protocols_) {
      const bool final_level=protocol==protocols_.back();
      const auto resolved=it::resolve(request_,definition_.model(protocol),protocol,final_level,definition_.application);
      if(final_level && world_.rank()==0 && it::prints(verbosity,it::needed::targets_below_floor))
        it::write_targets_below_floor(std::cout,"response",task_,resolved,definition_.labels);
    }
    restart_plan_=it::plan_restart(protocols_,request_,definition_.model(protocols_.back()),restart_expected,saved);
    activate(protocols_.at(restart_plan_.start_protocol));
  }
  [[nodiscard]] const it::RestartPlan& plan() const {return restart_plan_;}
  [[nodiscard]] const json& scientific_settings() const {return settings_;}
  void activate(it::Protocol protocol) {
    cp::records::require(std::find(protocols_.begin(),protocols_.end(),protocol)!=protocols_.end(),"unrequested response protocol");
    final_level_=protocol==protocols_.back();
    settings=it::resolve(request_,definition_.model(protocol),protocol,final_level_,definition_.application);
    reusable_=saved_ && final_level_ && saved_->settings==settings_ && saved_->protocol==protocol &&
        it::saved_completion(settings,expected(),stored(*saved_),true).has_value();
    complete_=reusable_;
  }
  [[nodiscard]] const it::ResolvedSettings& controls() const {return settings;}
  [[nodiscard]] bool final_level() const {return final_level_;}
  [[nodiscard]] it::ExpectedTask expected() const {
    const auto subjects=cp::records::subjects(definition_.record);
    return {task_,context_,{subjects.begin(),subjects.end()}};
  }
  [[nodiscard]] it::IterationInfo info() const {return {task_,settings.protocol,0,context_};}
  [[nodiscard]] bool has_saved() const {return saved_.has_value();}
  [[nodiscard]] bool reusable() const {return reusable_;}
  [[nodiscard]] bool complete() const {return complete_;}
  [[nodiscard]] const TaskDefinition& definition() const {return definition_;}
  // The node folder relative to the run folder, for metadata and results.
  [[nodiscard]] std::string folder() const {return excited_folder(settings_.at("excited_model")).generic_string();}
  [[nodiscard]] const it::NodeDescription& description() const {
    cp::records::require(saved_.has_value(),"no saved excited-state set");
    return *saved_;
  }
  // The saved measurements as this launch's measured values.
  [[nodiscard]] it::IterationMetrics saved_metrics() const {return stored(description());}
  [[nodiscard]] std::size_t previous_updates() const {
    return saved_ && saved_->protocol==settings.protocol ? saved_->updates_at_protocol : 0;
  }
  void note_guess(const json& source) {lineage_["guess"]=source;}
  // reader(archive folder, description), on every process.
  template<class Reader> auto load(const Reader& reader,bool require_complete=false) {
    cp::records::require(!require_complete || complete_,"property input has not met the current request");
    return reader(it::Node(folder_).archive(),description());
  }
  // One row per decision in the node's history, before any save.
  void record_evaluation(const it::IterationMetrics& metrics,const it::Decision& decision,
                         const std::vector<it::Diagnostic>& diagnostics={}) const {
    cp::detail::on_root(world_,[&] {
      it::append_history(it::Node(folder_).history(),
          it::history_row(metrics,decision,definition_.metric_ids,diagnostics));
      return std::string{};
    });
  }
  // writer(folder) saves root_<slot> for every slot collectively into an empty
  // staging folder; rank 0's one Node::save then replaces the whole set.
  template<class Writer>
  void publish(const it::IterationMetrics& metrics,const it::Decision& decision,
               const it::Progress& progress,const std::vector<int>& stable_index,const Writer& writer) {
    cp::records::require(metrics.info==it::IterationInfo{task_,settings.protocol,metrics.info.iteration,context_},
                         "response metrics belong to a different binding or protocol");
    complete_=decision.action==it::Action::complete;
    if(!decision.checkpoint) return;
    cp::records::require(stable_index.size()==metrics.states.size(),"excited bundle identity count differs");
    json roots=json::array();
    for(std::size_t slot=0;slot<stable_index.size();++slot) {
      const auto id=make_root_id(stable_index[slot]);
      const auto root=std::find_if(metrics.states.begin(),metrics.states.end(),[&](const auto& s){return s.id==id;});
      cp::records::require(root!=metrics.states.end(),"excited root has no metrics");
      const auto omega=root->metrics.at(it::Metric::energy_change).evidence.current;
      cp::records::require(omega && std::isfinite(*omega),"excited root has no finite excitation energy");
      roots.push_back({{"slot",slot},{"root_id",id},{"stable_index",stable_index[slot]},{"omega",*omega}});
    }
    const auto staging=folder_/"archive.staging";
    cp::detail::on_root(world_,[&] {
      std::filesystem::remove_all(staging);
      std::filesystem::create_directories(staging);
      return std::string{};
    });
    writer(staging);
    world_.gop.fence();
    cp::detail::on_root(world_,[&] {
      // The native writer count, from the header of the first shard.
      madness::archive::BinaryFstreamInputArchive header((staging/"root_0.00000").string());
      int writers=0;
      header & writers;
      it::NodeDescription d;
      d.kind="excited_set";
      d.protocol=metrics.info.protocol;
      d.updates_at_protocol=progress.updates_at_protocol;
      d.settings=settings_;
      d.metrics=it::node_metrics(metrics,definition_.metric_ids);
      d.lineage=lineage_;
      d.code_version=madness::info::git_commit();
      d.level_finished=decision.action==it::Action::finish_protocol || decision.action==it::Action::complete;
      d.details={{"excited_node_version",1},{"n_roots",roots.size()},{"roots",roots},{"writers",writers}};
      it::Node(folder_).save(d,[&](const std::filesystem::path& archive) {
        for(const auto& file:std::filesystem::directory_iterator(staging))
          std::filesystem::rename(file.path(),archive/file.path().filename());
      });
      std::filesystem::remove_all(staging);
      return std::string{};
    });
    saved_=madness::read_node(world_,folder_);
  }
};

// Use the existing per-root native serializer, one archive per root slot.
template<class State>
void save_excited_functions(madness::World& world,const State& state,const std::filesystem::path& directory) {
  for(std::size_t slot=0;slot<state.roots.size();++slot)
    state.roots[slot].save(world,madness::native_archive_name((directory/("root_"+std::to_string(slot))).string()),IoBackend::Native);
}

template<class Type,class Shell>
typename ESSolver<Type,Shell>::State load_excited_functions(
    madness::World& world,const std::filesystem::path& directory,
    const it::NodeDescription& description,it::Protocol target) {
  using Solver=ESSolver<Type,Shell>;
  const auto stored=description.protocol;
  const auto& roots=description.details.at("roots");
  cp::records::require(!roots.empty() && roots.size()==description.details.at("n_roots").get<std::size_t>() &&
      description.details.at("writers").get<int>()<=world.size(),"unsupported excited native payload");
  cp::records::require(it::Protocol{madness::FunctionDefaults<3>::get_thresh(),madness::FunctionDefaults<3>::get_k()}==target,
      "prepare excited representation before loading");
  typename Solver::State result;
  result.omega=madness::Tensor<double>(static_cast<long>(roots.size()));
  std::set<int> identities;
  for(std::size_t slot=0;slot<roots.size();++slot) {
    const auto& root=roots[slot];
    const int index=root.at("stable_index");
    cp::records::require(root.at("slot")==slot && index>=0 && root.at("root_id")==make_root_id(index) &&
        identities.insert(index).second,"unsupported excited native payload");
    const double omega=root.at("omega");
    cp::records::require(std::isfinite(omega),"excited checkpoint has no finite root energy");
    typename Solver::Storage functions;
    {
      struct Restore {
        it::Protocol target;
        ~Restore(){madness::FunctionDefaults<3>::set_k(target.k);madness::FunctionDefaults<3>::set_thresh(target.thresh);}
      } restore{target};
      madness::FunctionDefaults<3>::set_k(stored.k);madness::FunctionDefaults<3>::set_thresh(stored.thresh);
      functions=Solver::Storage::load(world,madness::native_archive_name((directory/("root_"+std::to_string(slot))).string()),IoBackend::Native);
    }
    for(auto* block:functions.blocks()) {
      for(auto& function:*block) {
        cp::records::require(function.k()==stored.k,"excited archive polynomial order differs");
        if(stored.k!=target.k) function=madness::project(function,target.k,target.thresh,false);
      }
      madness::set_thresh(world,*block,target.thresh);
      if(stored!=target) madness::truncate(world,*block,target.thresh);
    }
    result.roots.push_back(std::move(functions));result.stable_index.push_back(index);
    result.omega(static_cast<long>(slot))=omega;
  }
  world.gop.fence();
  return result;
}
} // namespace molresponse_v3::response_shared
