#pragma once
#include <madness/chem/SCFScientificSettings.h>

#include "response_state.hpp"
#include <chrono>
#include <madness/chem/RootCoordination.h>
#include <madness/chem/ConvergenceSettings.h>
#include <madness/chem/IterationOutput.h>
#include <madness/chem/Restart.h>
#include <madness/chem/LadderPlanner.h>
#include <madness/chem/SCFCheckpoint.h>
#include <madness/misc/info.h>

namespace molresponse_v3::response_shared {
namespace cp = madness::checkpoint;
namespace it = madness::iteration;
using json = nlohmann::json;
// The ground state as response sees it: where its orbitals are, the level they
// were saved at, and an identity (its settings, that level and its save time)
// that changes whenever the ground state is saved again. Response records the
// identity in its settings; it never names how the ground state is stored.
struct GroundSource {
  std::filesystem::path archive;
  it::Protocol protocol;
  json identity;
};
inline std::map<it::Metric, std::string> metric_set() {
  return {{it::Metric::bsh_residual, "fd-bsh-v1"},
          {it::Metric::density_change, "fd-density-v1"},
          {it::Metric::energy_change, "fd-energy-v1"},
          {it::Metric::kinetic_change, "fd-kinetic-v1"}};
}
inline std::map<it::Metric, std::string> metric_labels() {
  return {{it::Metric::bsh_residual, "BSH residual, maximum"},
          {it::Metric::density_change, "Induced-density change"},
          {it::Metric::energy_change, "Energy change (fixed frequency)"},
          {it::Metric::kinetic_change, "Response kinetic change (Ha * amplitude^2)"}};
}
inline json definitions() {
  json result;
  for (const auto &[metric, id] : metric_set()) {
    const bool residual = metric == it::Metric::bsh_residual;
    const bool density = metric == it::Metric::density_change;
    result[id] = {
        {"metric", it::metric_name(metric)},
        {"meaning",
         residual
             ? "raw response fixed-point difference before KAIN and restriction"
         : density ? "L2 change of kernel-induced response density"
         : metric == it::Metric::kinetic_change
             ? "absolute change in positive gradient kinetic sum"
             : "inapplicable: driven frequency is fixed input"},
        {"units", residual  ? "response_amplitude_L2_atomic_units"
                  : density ? "kernel_induced_density_L2_atomic_units"
                            : "hartree_times_response_amplitude_scale_squared"},
        {"fields", residual  ? "all stored spin and X/Y functions"
                   : density ? "kernel-induced density"
                             : "actual stored response blocks"},
        {"weights",
         density ? "kernel physical spin and X/Y factors" : "unit weights"},
        {"normalization",
         "existing response amplitude convention; no added rescaling"},
        {"reduction", residual  ? "maximum"
                      : density ? "field L2 difference"
                                : "absolute scalar difference"},
        {"comparison",
         residual ? "none" : "preceding comparable prepared iteration"},
        {"scale",
         residual
             ? "maximum current per-function L2 norm"
             : "maximum current/previous quantity magnitude or field norm"}};
  }
  return result;
}
// Driven response. Default targets are the branch's response defaults; the
// shared floor (5 * threshold) applies at intermediate levels. BSH and density
// use the scaled rule, T * (1 + scale), with the scale this solver measures.
// The frequency is fixed input, so an energy change is not applicable.
inline it::ModelMetrics model(it::Protocol = {1e-4, 6}) {
  return {{it::Metric::bsh_residual, {it::Support::supported, true, 5e-4, it::LimitRule::scaled, {}}},
          {it::Metric::density_change, {it::Support::supported, true, 5e-4, it::LimitRule::scaled, {}}},
          {it::Metric::energy_change, {it::Support::not_applicable, false, {}, it::LimitRule::plain, {}}},
          {it::Metric::kinetic_change, {it::Support::supported, false, {}, it::LimitRule::plain, {}}}};
}

// What a saved driven state must share with this one to be its starting
// state: every setting but the frequency and the static or full part of the
// layout. The ground state's identity, when the settings hold it, is part of it.
inline json start_model(json settings) {
  settings.erase("frequency");
  if(settings.contains("layout")) {
    const auto layout=settings.at("layout").get<std::string>();
    settings["layout"]=layout.substr(0,layout.rfind('_')+1); // "closed_shell_" or "open_shell_"
  }
  return settings;
}
// Applications supply their metric definitions and restart compatibility check.
// The binding owns records, never numerical state or comparison history.
struct TaskDefinition {
  json record;
  std::map<it::Metric,std::string> metric_ids;
  json metric_definitions;
  std::function<it::ModelMetrics(it::Protocol)> model;
  std::function<void(const json&,const json&)> validate_settings;
  std::string application="driven response"; // names the solver in input errors
  it::MetricLabels labels=metric_labels();
  // The check on another node's saved settings before this state starts from
  // it; validate_settings when unset.
  std::function<void(const json&,const json&)> validate_start;
};
inline TaskDefinition driven_task() {
  return {{{"subject","response"},{"fd_checkpoint_version",1}}, metric_set(), definitions(),
      [](it::Protocol protocol) {return model(protocol);},
      [](const json& previous,const json& current) {
        cp::records::require(previous.at("frequency")==current.at("frequency") &&
            previous.at("perturbation")==current.at("perturbation"),
            "response frequency conflicts with saved task identity; use a new output prefix");
        cp::records::require(previous.at("layout")==current.at("layout"),
            "unsupported driven checkpoint block layout");
      },
      "driven response",metric_labels(),
      [](const json& other,const json& current) {
        cp::records::require(start_model(other)==start_model(current),
            "a saved driven state of another perturbation or model cannot start this one");
      }};
}

inline TaskDefinition derived_task() {
  auto result=driven_task();
  result.validate_settings=[](const json& previous,const json& current) {
    cp::records::require(previous.at("perturbation")==current.at("perturbation") &&
        previous.at("layout")==current.at("layout") &&
        previous.at("derived_from")==current.at("derived_from"),
        "derived response origin or layout differs from its saved functions");
  };
  result.validate_start=result.validate_settings;
  return result;
}

// One subject's measured values and the reasons for the rest, for results.
inline json measured_values(const it::IterationMetrics& frame,std::string subject = {}) {
  if (subject.empty()) {
    cp::records::require(frame.states.size()==1,"bundle metrics need a subject");
    subject=frame.states.front().id;
  }
  const auto state=std::find_if(frame.states.begin(),frame.states.end(),[&](const auto& s){return s.id==subject;});
  cp::records::require(state!=frame.states.end(),"unknown response subject");
  json values=json::object(),availability=json::object();
  for (const auto& [metric,measurement] : state->metrics) {
    if (measurement.status==it::Availability::measured) values[it::metric_name(metric)]=*measurement.value;
    else availability[it::metric_name(metric)]={{"status",measurement.status==it::Availability::not_applicable ? "not_applicable" : "unavailable"},{"reason",measurement.reason}};
  }
  return {{"values",values},{"availability",availability}};
}
// A response node's folder, <run>/states/response/<name>/. The caller names it
// from the state's parts, joined by hyphens; a name holds only letters,
// digits, '.', '_' and '-'.
inline std::filesystem::path response_node_folder(const std::string& name,const std::string& what) {
  cp::records::require(!name.empty() && name.find_first_not_of(
      "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")==std::string::npos &&
      name.find_first_not_of('.')!=std::string::npos,"invalid "+what+" node name: "+name);
  return std::filesystem::path("states")/"response"/name;
}
inline std::filesystem::path driven_folder(const std::string& name) {return response_node_folder(name,"driven-state");}
// Another node's description, read without finishing or discarding anything
// in its folder: a node this launch does not own is never changed. Nothing
// when it has no readable description and archive. Collective.
inline std::optional<it::NodeDescription> read_other_node(madness::World& world,const std::filesystem::path& folder) {
  const auto text=cp::detail::on_root(world,[&] {
    try {
      if(!std::filesystem::is_directory(it::Node(folder).archive())) return std::string{};
      std::ifstream file(folder/"description.json");
      return it::description_record(it::read_description(json::parse(file))).dump();
    } catch(const std::exception&) {return std::string{};}
  });
  if(text.empty()) return std::nullopt;
  return it::read_description(json::parse(text));
}

// The functions of a saved driven state (archive/response in `directory`) at
// level `target`: read at the saved level, then projected and truncated to the
// target. FunctionDefaults must already be at the target. Collective.
template<class Storage>
Storage load_driven_functions(madness::World& world,const std::filesystem::path& directory,
                              const it::NodeDescription& d,it::Protocol target) {
  cp::records::require(d.details.at("writers").get<int>()<=world.size(),"unsupported native response payload");
  cp::records::require(it::Protocol{madness::FunctionDefaults<3>::get_thresh(),madness::FunctionDefaults<3>::get_k()}==
                       target,"prepare response representation before native loading");
  const auto stored=d.protocol;
  Storage functions;
  {
    struct Restore {
      it::Protocol target;
      ~Restore() {madness::FunctionDefaults<3>::set_k(target.k);madness::FunctionDefaults<3>::set_thresh(target.thresh);}
    } restore{target};
    madness::FunctionDefaults<3>::set_k(stored.k);
    madness::FunctionDefaults<3>::set_thresh(stored.thresh);
    functions=Storage::load(world,madness::native_archive_name((directory/"response").string()),IoBackend::Native);
  }
  for (auto* block : functions.blocks()) {
    for (auto& function : *block) {
      cp::records::require(function.k()==stored.k,"response payload disagrees with saved polynomial order");
      if (stored.k!=target.k) function=madness::project(function,target.k,target.thresh,false);
    }
    madness::set_thresh(world,*block,target.thresh);
    if (stored!=target) madness::truncate(world,*block,target.thresh);
  }
  world.gop.fence();
  return functions;
}

// A driven state, the response to one perturbation at one frequency, is one
// node: <run>/states/response/<perturbation>-f<frequency>/, with
// description.json, history.jsonl and archive/response (one native archive).
// The restart decision follows ExcitedSet: the stored metrics become measured
// values and saved_completion decides reuse. A saved state with other
// settings (a ground state saved again among them) is only a starting guess,
// and so is the saved state of another node this one starts from (a nearby
// frequency); that node is only read. Each node is written by the one process
// group that solves the state, with no lock.
class DrivenState {
  madness::World& world_;
  std::string name_;
  std::filesystem::path folder_,source_;
  std::string task_,context_;
  json settings_,lineage_;
  TaskDefinition definition_;
  it::Request request_;
  std::vector<it::Protocol> protocols_;
  it::RestartPlan restart_plan_;
  it::ResolvedSettings settings;
  bool final_level_=true;
  std::optional<it::NodeDescription> saved_; // this node's, or the node it starts from
  bool own_=false,reusable_=false,complete_=false;

  [[nodiscard]] it::IterationMetrics stored(const it::NodeDescription& d) const {
    return it::stored_metrics(d,{task_,d.protocol,0,context_});
  }
public:
  DrivenState(madness::World& world,const std::filesystem::path& root,const GroundSource& ground,
              std::string task,std::string name,json scientific_settings,const it::Request& request,
              std::vector<it::Protocol> protocols,TaskDefinition definition=driven_task(),
              const std::optional<std::string>& start_from={},bool launch_header=true,
              it::Verbosity verbosity=it::Verbosity::detail)
      : world_(world),name_(std::move(name)),task_(std::move(task)),
        settings_(std::move(scientific_settings)),definition_(std::move(definition)),request_(request),
        protocols_(std::move(protocols)),
        settings(it::resolve(request,definition_.model(protocols_.at(0)),protocols_.at(0),
                             protocols_.size()==1,definition_.application)) {
    // The ground state's identity in the settings turns saved response
    // states into guesses whenever the ground state is saved again.
    settings_["ground_settings"]=ground.identity.at("settings");
    settings_["ground_saved"]={{"protocol",ground.identity.at("protocol")},{"saved_at",ground.identity.at("saved_at")}};
    context_=json{{"settings",settings_},{"inputs",json::array()}}.dump();
    lineage_={{"ground_state",ground.identity}};
    folder_=std::filesystem::absolute(root)/driven_folder(name_);
    if(launch_header) madness::write_launch_header(world_,folder_);
    saved_=madness::read_node(world_,folder_); // collective; a half-present node is an error
    own_=saved_.has_value();source_=folder_;
    if(!saved_ && start_from && *start_from!=name_) {
      const auto other=std::filesystem::absolute(root)/driven_folder(*start_from);
      saved_=read_other_node(world_,other);
      if(saved_) {source_=other;lineage_["started_from"]=*start_from;}
    }
    std::optional<it::SavedTask> saved;
    if(saved_) {
      const auto& d=*saved_;
      cp::records::require(d.kind=="driven" && d.details.at("driven_node_version")==1,
          "unsupported driven-state node: "+source_.string());
      cp::records::require(madness::scf_model_settings(d.settings.at("ground_settings"))==
                           madness::scf_model_settings(settings_.at("ground_settings")),
          "response ground-state model differs from its saved checkpoint");
      if(own_) definition_.validate_settings(d.settings,settings_);
      else (definition_.validate_start ? definition_.validate_start : definition_.validate_settings)(d.settings,settings_);
      saved=it::SavedTask{};
      saved->metrics=stored(d);
      const auto subjects=cp::records::subjects(definition_.record);
      for(const auto& state:saved->metrics.states)
        cp::records::require(subjects.contains(state.id),"unsupported response task definition");
      saved->protocol_complete=d.level_finished;
      // Another node's update count is not this state's.
      saved->updates_at_protocol=own_ ? d.updates_at_protocol : 0;
      saved->conversion_supported=true;
      if(!own_ || d.settings!=settings_) saved->use=it::SavedEvidenceUse::guess_only;
      // A compatible restart keeps naming what the saved state started from.
      else for(const auto* key:{"guess","started_from"})
        if(d.lineage.contains(key)) lineage_[key]=d.lineage.at(key);
    }
    // Validate the request at every level before work; warn once about final
    // targets below the protocol floor.
    for(const auto& protocol:protocols_) {
      const bool final_level=protocol==protocols_.back();
      const auto resolved=it::resolve(request_,definition_.model(protocol),protocol,final_level,definition_.application);
      if(final_level && world_.rank()==0 && it::prints(verbosity,it::needed::targets_below_floor))
        it::write_targets_below_floor(std::cout,"response",task_,resolved,definition_.labels);
    }
    restart_plan_=it::plan_restart(protocols_,request_,definition_.model(protocols_.back()),expected(),saved);
    activate(protocols_.at(restart_plan_.start_protocol));
  }
  // Re-read this node after another process group saved it.
  void refresh() {
    auto latest=madness::read_node(world_,folder_);
    if(!latest) return;
    saved_=std::move(latest);own_=true;source_=folder_;
    if(std::find(protocols_.begin(),protocols_.end(),saved_->protocol)!=protocols_.end()) activate(saved_->protocol);
  }
  [[nodiscard]] const it::RestartPlan& plan() const {return restart_plan_;}
  [[nodiscard]] const json& scientific_settings() const {return settings_;}
  void activate(it::Protocol protocol) {
    cp::records::require(std::find(protocols_.begin(),protocols_.end(),protocol)!=protocols_.end(),"unrequested response protocol");
    final_level_=protocol==protocols_.back();
    settings=it::resolve(request_,definition_.model(protocol),protocol,final_level_,definition_.application);
    reusable_=own_ && final_level_ && saved_->settings==settings_ && saved_->protocol==protocol &&
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
  // A saved state to start from: this node's, or a nearby frequency's.
  [[nodiscard]] bool has_saved() const {return saved_.has_value();}
  // This node has its own saved state.
  [[nodiscard]] bool has_node() const {return own_;}
  // The node this state starts from, while it has no saved state of its own.
  [[nodiscard]] std::optional<std::string> started_from() const {
    if(own_ || !saved_) return std::nullopt;
    return source_.filename().string();
  }
  // The settings of the saved state it starts from: this node's or another's.
  [[nodiscard]] const json& source_settings() const {
    cp::records::require(saved_.has_value(),"no saved driven state");
    return saved_->settings;
  }
  [[nodiscard]] bool reusable() const {return reusable_;}
  [[nodiscard]] bool complete() const {return complete_;}
  [[nodiscard]] const TaskDefinition& definition() const {return definition_;}
  [[nodiscard]] const std::string& name() const {return name_;}
  // The node folder relative to the run folder, for metadata and results.
  [[nodiscard]] std::string folder() const {return driven_folder(name_).generic_string();}
  // This node's own saved description.
  [[nodiscard]] const it::NodeDescription& description() const {
    cp::records::require(own_,"no saved driven state");
    return *saved_;
  }
  // The saved measurements as this launch's measured values.
  [[nodiscard]] it::IterationMetrics saved_metrics() const {
    cp::records::require(saved_.has_value(),"no saved driven state");
    return stored(*saved_);
  }
  [[nodiscard]] json measured_values(std::string subject = {}) const {
    return response_shared::measured_values(stored(description()),std::move(subject));
  }
  [[nodiscard]] std::size_t previous_updates() const {
    return own_ && saved_->protocol==settings.protocol ? saved_->updates_at_protocol : 0;
  }
  void note_guess(const json& source) {lineage_["guess"]=source;}
  // reader(archive folder, description), on every process: this node's saved
  // state, or the one it starts from.
  template<class Reader> auto load(const Reader& reader,bool require_complete=false) {
    cp::records::require(!require_complete || complete_,"property input has not met the current request");
    cp::records::require(saved_.has_value(),"no saved driven state");
    return reader(it::Node(source_).archive(),*saved_);
  }
  // Driven storage's native conversion to the current protocol level.
  template<class Storage> Storage load(bool require_complete=false) {
    return load([&](const std::filesystem::path& directory,const it::NodeDescription& d) {
      return load_driven_functions<Storage>(world_,directory,d,settings.protocol);
    },require_complete);
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
  // writer(folder) saves the state's native archive "response" collectively
  // into an empty staging folder; rank 0's one Node::save then replaces it.
  template<class Writer>
  void publish(const it::IterationMetrics& metrics,const it::Decision& decision,
               const it::Progress& progress,const Writer& writer) {
    cp::records::require(metrics.info==it::IterationInfo{task_,settings.protocol,metrics.info.iteration,context_},
                         "response metrics belong to a different binding or protocol");
    complete_=decision.action==it::Action::complete;
    if(!decision.checkpoint) return;
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
      madness::archive::BinaryFstreamInputArchive header((staging/"response.00000").string());
      int writers=0;
      header & writers;
      it::NodeDescription d;
      d.kind="driven";
      d.protocol=metrics.info.protocol;
      d.updates_at_protocol=progress.updates_at_protocol;
      d.settings=settings_;
      d.metrics=it::node_metrics(metrics,definition_.metric_ids);
      d.lineage=lineage_;
      d.code_version=madness::info::git_commit();
      d.level_finished=decision.action==it::Action::finish_protocol || decision.action==it::Action::complete;
      d.details={{"driven_node_version",1},{"writers",writers}};
      it::Node(folder_).save(d,[&](const std::filesystem::path& archive) {
        for(const auto& file:std::filesystem::directory_iterator(staging))
          std::filesystem::rename(file.path(),archive/file.path().filename());
      });
      std::filesystem::remove_all(staging);
      return std::string{};
    });
    saved_=madness::read_node(world_,folder_);own_=true;source_=folder_;
  }
};
// A second-order intermediate (VBC sources, two-photon sources) is built once
// from converged response states and saved in its own node folder,
// <run>/states/response/<name>/, with no history: it is never iterated.
// The caller names the folder from the intermediate's parts, joined by
// hyphens; a name holds only letters, digits, '.', '_' and '-'.
inline std::filesystem::path second_order_folder(const std::string& name) {
  return response_node_folder(name,"second-order");
}
// One lineage entry: which input state (its node name), at which level, and
// its identity, the save time from the input's own node description.
inline json lineage_entry(const std::string& role,const DrivenState& input) {
  cp::records::require(input.complete(),"response input has not met the current final request");
  const auto& d=input.description();
  return {{"role",role},{"name",input.name()},{"protocol",cp::records::protocol_record(d.protocol)},
          {"identity",{{"saved_at",d.saved_at}}}};
}
// Reused when the saved settings and every lineage entry equal the current
// ones; otherwise built again and saved over the old node. Never a guess.
// Each node is written by the one process group that builds it.
class SecondOrderNode {
  madness::World& world_;
  std::string id_,name_;
  std::filesystem::path folder_;
  json settings_,lineage_;
  std::optional<it::NodeDescription> saved_;
  bool reusable_=false;
public:
  SecondOrderNode(madness::World& world,const std::filesystem::path& root,
      std::string id,std::string name,json settings,json lineage)
      : world_(world),id_(std::move(id)),name_(std::move(name)),
        folder_(std::filesystem::absolute(root)/second_order_folder(name_)),
        settings_(std::move(settings)),lineage_(std::move(lineage)) {
    cp::records::require(lineage_.is_array() && !lineage_.empty(),"second-order node needs its inputs");
    saved_=madness::read_node(world_,folder_); // collective; a half-present node is an error
    if (saved_) {
      cp::records::require(saved_->kind=="second_order","response node changed kind: "+folder_.string());
      reusable_=saved_->settings==settings_ && saved_->lineage==lineage_;
    }
  }
  [[nodiscard]] const std::string& id() const {return id_;}
  [[nodiscard]] const std::string& name() const {return name_;}
  [[nodiscard]] bool reusable() const {return reusable_;}
  [[nodiscard]] const json& settings() const {return settings_;}
  [[nodiscard]] const json& lineage() const {return lineage_;}
  // The node folder relative to the run folder, for metadata and results.
  [[nodiscard]] std::string folder() const {return second_order_folder(name_).generic_string();}
  [[nodiscard]] const it::NodeDescription& description() const {
    cp::records::require(saved_.has_value(),"second-order node has no saved state");
    return *saved_;
  }
  // reader(archive folder, description), on every process.
  template<class Reader> auto load(const Reader& reader) {
    cp::records::require(reusable_,"second-order node does not match its current inputs");
    return reader(it::Node(folder_).archive(),*saved_);
  }
  // writer(folder) writes the named archives collectively into an empty
  // staging folder; rank 0's one save then moves them into the node.
  template<class Writer> void save(const std::vector<std::string>& archives,const Writer& writer) {
    const auto staging=folder_/"archive.staging";
    cp::detail::on_root(world_,[&] {
      std::filesystem::remove_all(staging);
      std::filesystem::create_directories(staging);
      return std::string{};
    });
    writer(staging);
    world_.gop.fence();
    cp::detail::on_root(world_,[&] {
      it::NodeDescription d;
      d.kind="second_order";
      d.protocol=cp::records::read_protocol(settings_.at("protocol"));
      d.settings=settings_;
      d.metrics=json::array();
      d.lineage=lineage_;
      d.code_version=madness::info::git_commit();
      d.level_finished=true;
      d.details={{"second_order_node_version",1},{"archives",archives}};
      it::Node(folder_).save(d,[&](const std::filesystem::path& archive) {
        for (const auto& file:std::filesystem::directory_iterator(staging))
          std::filesystem::rename(file.path(),archive/file.path().filename());
      });
      std::filesystem::remove_all(staging);
      return std::string{};
    });
    saved_=madness::read_node(world_,folder_);
    reusable_=true;
  }
};
} // namespace molresponse_v3::response_shared
