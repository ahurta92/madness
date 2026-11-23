#pragma once
#include <madness/chem/SCFScientificSettings.h>

#include "response_state.hpp"
#include <chrono>
#include <madness/chem/CheckpointStore.h>
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

// Applications supply their metric definitions and restart compatibility check.
// The binding owns records, never numerical state or comparison history.
struct TaskDefinition {
  json record;
  std::map<it::Metric,std::string> metric_ids;
  json metric_definitions;
  std::function<it::ModelMetrics(it::Protocol)> model;
  std::function<void(const json&,const json&)> validate_settings;
  bool allow_subject_expansion=false;
  std::string application="driven response"; // names the solver in input errors
  it::MetricLabels labels=metric_labels();
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
  return result;
}

class Checkpoint {
  madness::World &world_;
  std::filesystem::path root_;
  cp::CheckpointStore store_;
  std::string task_, owner_, context_;
  json settings_, inputs_;
  TaskDefinition definition_;
  it::Request request_;
  std::vector<it::Protocol> protocols_;
  it::RestartPlan restart_plan_;
  it::ResolvedSettings settings;
  bool final_level_ = true;
  std::optional<cp::Manifest> selected_;
  bool reusable_ = false;
  bool complete_ = false;
  std::size_t publications_ = 0;
  std::optional<cp::GenerationRef> replaces_;
  std::filesystem::path ground_archive_;

public:
  Checkpoint(madness::World &world, const std::filesystem::path& root, const GroundSource &ground,
             std::string task, json scientific_settings,
             const it::Request &request, std::vector<it::Protocol> protocols,
             const std::string& attempt, TaskDefinition definition = driven_task(),
             std::optional<cp::GenerationRef> guess = {}, std::string writer_group = {},
             it::Verbosity verbosity = it::Verbosity::detail)
      : world_(world), root_(std::filesystem::absolute(root)),
        store_(world, root_,{},std::move(writer_group)), task_(std::move(task)),
        settings_(std::move(scientific_settings)),
        definition_(std::move(definition)), request_(request), protocols_(std::move(protocols)),
        ground_archive_(ground.archive),
        settings(it::resolve(request, definition_.model(protocols_.at(0)), protocols_.at(0),
                             protocols_.size() == 1, definition_.application)) {
    // The ground state is not a store input: its identity in the settings
    // makes a changed ground state turn saved response states into guesses.
    inputs_ = json::array();
    settings_["ground_settings"] = ground.identity.at("settings");
    settings_["ground_saved"] = {{"protocol", ground.identity.at("protocol")},
                                 {"saved_at", ground.identity.at("saved_at")}};
    context_ = json{{"settings", settings_}, {"inputs", inputs_}}.dump();
    owner_ = attempt;
    store_.initialize(owner_);
    selected_ = store_.read_latest("response", task_, true);
    if(!selected_ && guess) selected_=store_.read(*guess);
    std::optional<it::SavedTask> saved;
    if (selected_) {
      const auto &old = selected_->record;
      cp::records::require(old.at("metric_definitions") == definition_.metric_definitions,
          "unsupported response metric definitions");
      bool expanded=false;
      const auto& old_definition=old.at("task").at("definition");
      if(definition_.allow_subject_expansion) {
        const auto previous=cp::records::subjects(old_definition);
        const auto wanted=cp::records::subjects(definition_.record);
        expanded=previous.size()<wanted.size() &&
            std::includes(wanted.begin(),wanted.end(),previous.begin(),previous.end());
      }
      for (const auto& [key,value] : definition_.record.items()) {
        if(key=="subjects" && expanded) continue;
        cp::records::require(old_definition.at(key)==value,"unsupported response task definition");
      }
      const auto& previous_settings=old.at("task").at("definition").at("settings");
      cp::records::require(madness::scf_model_settings(previous_settings.at("ground_settings")) ==
                              madness::scf_model_settings(settings_.at("ground_settings")),
                          "response ground-state model differs from its saved checkpoint");
      definition_.validate_settings(previous_settings,settings_);
      const auto metrics = cp::iteration_metrics_from_manifest(*selected_);
      saved = it::SavedTask{};
      saved->metrics=metrics;
      saved->protocol_complete=old.at("task").at("iteration").at("protocol_complete");
      saved->updates_at_protocol=cp::records::count(old.at("task").at("iteration").at("updates_at_protocol"));
      saved->conversion_supported=true;
      saved->bundle_expanded=expanded;
      if (metrics.info.context_id!=context_ || metrics.info.task_id!=task_)
        saved->use=it::SavedEvidenceUse::guess_only;
      if(old_definition.contains("replaces")) {
        const auto& previous=old_definition.at("replaces");
        replaces_=cp::GenerationRef{"response",previous.at("task_id"),previous.at("generation")};
      }
      if(settings_.contains("derived_from") &&
          (previous_settings.at("frequency")!=settings_.at("frequency") || metrics.info.task_id!=task_))
        replaces_=cp::GenerationRef{"response",metrics.info.task_id,selected_->generation};
    }
    // Change provenance only after validating the saved model and native layout.
    // A refined same-model ground source keeps these functions as guesses; the
    // changed context prevents old measurements from qualifying the new solve.
    store_.register_task("response", "molresponse", task_, "iterative", true, inputs_);
    auto restart_expected=expected();
    // The application explicitly selected this old task as a guess. Validate
    // its original evidence identity; never relabel its saved measurements.
    if(saved && saved->use==it::SavedEvidenceUse::guess_only)
      restart_expected.task_id=saved->metrics.info.task_id;
    // Validate the request at every level before work; warn once about final
    // targets below the protocol floor.
    for(const auto& protocol:protocols_) {
      const bool final_level=protocol==protocols_.back();
      const auto resolved=it::resolve(request_,definition_.model(protocol),protocol,final_level,definition_.application);
      if(final_level && world_.rank()==0 && it::prints(verbosity,it::needed::targets_below_floor))
        it::write_targets_below_floor(std::cout,"response",task_,resolved,definition_.labels);
    }
    restart_plan_=it::plan_restart(protocols_,request_,definition_.model(protocols_.back()),restart_expected,saved);
    activate(protocols_.at(restart_plan_.start_protocol));
  }
  void refresh() {
    const auto latest=store_.read_latest("response",task_);
    if(!latest) return;
    selected_=latest;
    const auto protocol=cp::iteration_metrics_from_manifest(*latest).info.protocol;
    if(std::find(protocols_.begin(),protocols_.end(),protocol)!=protocols_.end()) activate(protocol);
  }
  const it::RestartPlan& plan() const {return restart_plan_;}
  const json& scientific_settings() const {return settings_;}
  it::Protocol final_protocol() const {return protocols_.back();}
  void activate(it::Protocol protocol) {
    cp::records::require(std::find(protocols_.begin(),protocols_.end(),protocol)!=protocols_.end(),"unrequested response protocol");
    final_level_=protocol==protocols_.back();
    settings=it::resolve(request_,definition_.model(protocol),protocol,final_level_,definition_.application);
    reusable_=false;
    if (selected_ && final_level_) {
      const auto metrics=cp::iteration_metrics_from_manifest(*selected_);
      reusable_=metrics.info.task_id==task_ && metrics.info.context_id==context_ && metrics.info.protocol==protocol &&
          it::saved_completion(settings,expected(),metrics,true).has_value();
    }
    complete_=reusable_;
  }
  const it::ResolvedSettings &controls() const { return settings; }
  bool final_level() const { return final_level_; }
  it::ExpectedTask expected() const {
    const auto subjects=cp::records::subjects(definition_.record);
    return {task_,context_,{subjects.begin(),subjects.end()}};
  }
  it::IterationInfo info() const {
    return {task_, owner_, settings.protocol, 0, context_};
  }
  bool has_saved() const { return selected_.has_value(); }
  bool reusable() const { return reusable_; }
  bool complete() const { return complete_; }
  const cp::Manifest &manifest() const {
    cp::records::require(selected_.has_value(), "no published response state");
    return *selected_;
  }
  std::size_t previous_updates() const {
    return selected_ && cp::iteration_metrics_from_manifest(*selected_).info.protocol==settings.protocol ? cp::records::count(selected_->record.at("task")
                                              .at("iteration")
                                              .at("updates_at_protocol"))
                     : 0;
  }
  template<class Reader> auto load(const Reader& reader,bool require_complete=false) {
    cp::records::require(!require_complete || complete_,
                         "property input has not met the current request");
    const auto source_task=manifest().record.at("task").at("task_id").get<std::string>();
    const auto verified=store_.read({"response",source_task,manifest().generation});
    return reader(store_.task_directory("response",source_task)/verified.generation,verified);
  }
  // Driven storage's native conversion stays in this application adapter.
  template<class Storage> Storage load(bool require_complete=false) {
    return load([&](const std::filesystem::path& directory,const cp::Manifest& verified) {
      cp::records::require(verified.archives.size()==1 &&
          verified.archives.front().basename=="response" &&
          verified.archives.front().role=="functions" &&
          verified.archives.front().writers<=world_.size(),"unsupported native response payload");
      cp::records::require(
          it::Protocol{madness::FunctionDefaults<3>::get_thresh(),
                       madness::FunctionDefaults<3>::get_k()} ==
              settings.protocol,
          "prepare response representation before native loading");
      const auto stored=cp::iteration_metrics_from_manifest(verified).info.protocol;
      const auto target=settings.protocol;
      Storage functions;
      {
        struct Restore {
          it::Protocol target;
          ~Restore() {madness::FunctionDefaults<3>::set_k(target.k);madness::FunctionDefaults<3>::set_thresh(target.thresh);}
        } restore{target};
        madness::FunctionDefaults<3>::set_k(stored.k);
        madness::FunctionDefaults<3>::set_thresh(stored.thresh);
        functions=Storage::load(world_,madness::native_archive_name((directory/"response").string()),IoBackend::Native);
      }
      for (auto* block : functions.blocks()) {
        for (auto& function : *block) {
          cp::records::require(function.k()==stored.k,"response payload disagrees with saved polynomial order");
          if (stored.k!=target.k) function=madness::project(function,target.k,target.thresh,false);
        }
        madness::set_thresh(world_,*block,target.thresh);
        if (stored!=target) madness::truncate(world_,*block,target.thresh);
      }
      world_.gop.fence();
      return functions;
    },require_complete);
  }

  template <class Writer>
  void publish(const it::IterationMetrics &metrics, const it::Decision &decision,
               const it::Progress &progress, const Writer& writer) {
    cp::records::require(metrics.info == info_with_iteration(metrics.info.iteration),
                         "response metrics belong to a different binding or protocol");
    complete_ = decision.action == it::Action::complete;
    if (!decision.checkpoint)
      return;
    const auto generation = owner_ + "-" + std::to_string(++publications_);
    const json stage{
        {"application", "molresponse"},
        {"settings", settings_},
        {"requested_properties", json::array()},
        {"inputs", inputs_},
        {"ground", {{"archive", ground_archive_.string()}}}, // for inspection tools only
        {"outcome", {{"complete", complete_}, {"reason", decision.reason}}}};
    auto definition=definition_.record;
    definition["settings"]=settings_;
    if(replaces_) definition["replaces"]={{"task_id",replaces_->task_id},{"generation",replaces_->generation}};
    auto saved = cp::iterative_manifest("response", generation, stage,
        definition, inputs_, definition_.metric_definitions, definition_.metric_ids,
        metrics, decision, progress);
    const auto warning = store_.publish(saved,writer);
    selected_ = store_.read({"response", task_, generation});
    if (!warning.empty() && world_.rank() == 0)
      madness::print(warning);
  }
  it::IterationInfo info_with_iteration(std::size_t iteration) const {
    auto current=info(); current.iteration=iteration; return current;
  }
  json qualified_input(const std::string& role) const {
    cp::records::require(complete_,"response input has not met the current final request");
    const auto& saved=manifest();
    return {{"role",role},{"calculation_id","response"},{"task_id",task_},{"generation",saved.generation}};
  }
  void note_guess(const json& source) {definition_.record["guess"]=source;}
  const TaskDefinition& definition() const {return definition_;}
  json source_record(std::string subject = {}) const {
    const auto& saved=manifest();const auto frame=cp::iteration_metrics_from_manifest(saved);
    if (subject.empty()) {
      cp::records::require(frame.states.size()==1,"bundle source needs a subject");
      subject=frame.states.front().id;
    }
    cp::records::require(cp::records::subjects(definition_.record).contains(subject),"unknown response subject");
    return {{"task_id",frame.info.task_id},{"subject",subject},{"attempt_id",frame.info.attempt_id},
            {"generation",saved.generation},{"protocol",cp::records::protocol_record(frame.info.protocol)},
            {"iteration",frame.info.iteration}};
  }
  json measured_values(std::string subject = {}) const {
    const auto frame=cp::iteration_metrics_from_manifest(manifest());
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
};
// A second-order intermediate (VBC sources, two-photon sources) is built once
// from converged response states and saved in its own node folder,
// <run>/states/response/<name>/, with no history: it is never iterated.
// The caller names the folder from the intermediate's parts, joined by
// hyphens; a name holds only letters, digits, '.', '_' and '-'.
inline std::filesystem::path second_order_folder(const std::string& name) {
  cp::records::require(!name.empty() && name.find_first_not_of(
      "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")==std::string::npos &&
      name.find_first_not_of('.')!=std::string::npos,"invalid second-order node name: "+name);
  return std::filesystem::path("states")/"response"/name;
}
// One lineage entry: which input state, at which level, and its identity.
// 4c part 1: the inputs are still store states, whose records hold no save
// time, so the generation stands in for it. Parts 2 and 3 replace the
// identity with saved_at from the input's own node description; only this
// function changes.
inline json lineage_entry(const std::string& role,const Checkpoint& input) {
  const auto ref=input.qualified_input(role); // the input has met the final request
  const auto protocol=cp::iteration_metrics_from_manifest(input.manifest()).info.protocol;
  return {{"role",role},{"name",ref.at("task_id")},{"protocol",cp::records::protocol_record(protocol)},
          {"identity",{{"generation",ref.at("generation")}}}};
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
  const std::string& id() const {return id_;}
  const std::string& name() const {return name_;}
  bool reusable() const {return reusable_;}
  const json& settings() const {return settings_;}
  const json& lineage() const {return lineage_;}
  // The node folder relative to the run folder, for metadata and results.
  std::string folder() const {return second_order_folder(name_).generic_string();}
  const it::NodeDescription& description() const {
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
