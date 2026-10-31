#pragma once
#include "ConvergenceReport.hpp"
#include "MoldftCheckpoint.hpp"
#include "NamedStages.hpp"
#include "ResponseCheckpoint.hpp"
#include "ResponseProperties.hpp"
#include <algorithm>
#include <apps/molresponse/madqc_adapter.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <madness/chem/Applications.hpp>
#include <optional>
#include <unistd.h>

// Execute the validated stages from NamedStages.hpp using the existing engines.
// Flow: read/validate -> create output directories -> run stages -> write
// reports. Application adapters publish immutable ground and response
// checkpoints.
namespace named_stages {
namespace fs = std::filesystem;
using json = nlohmann::json;

// All ranks call this helper. Rank 0 performs the operation, then broadcasts
// any standard exception so every rank takes the same error path.
template <class F>
void on_root(madness::World& world, F operation) {
  std::string error;
  if (world.rank() == 0) {
    try {
      operation();
    } catch (const std::exception& e) {
      error = e.what();
    }
  }
  world.gop.broadcast_serializable(error, 0);
  if (!error.empty())
    throw std::runtime_error(error);
}
// Write an ordinary result report. This is not atomic checkpoint publication.
inline void write_result(madness::World& world, const fs::path& file,
                         const json& value) {
  on_root(world, [&] {
    std::ofstream output(file);
    output.exceptions(std::ios::failbit | std::ios::badbit);
    output << value.dump(2) << '\n';
    output.close();
  });
}
// DESIGN A8: what a stage hands to the stages that consume it. Everything a
// consumer may use travels here, by value; the runner keeps one per stage
// label.
struct StageOutcome {
  std::string state_id; // the published generation (A1 identity)
  bool converged = false;
  json result; // this stage's results.json record
};

// DESIGN A1: one declared input of a stage, read only from its producer's
// published generation on disk; no engine object passes between stages.
struct StageInput {
  std::string label;                   // the producing stage
  madness::checkpoint::Manifest saved; // its published generation
  fs::path saved_prefix;               // the archive prefix a consumer reads
  bool accepted = false;               // the producer's completion verdict
};

// Read a declared input from its producer's published generation. A producer
// in this deck has just run (or been reused), and the generation must be the
// one it reported. A producer outside the deck is read from the run's output
// directory as it stands, and only a converged generation is accepted.
inline StageInput
read_input(madness::World& world, const fs::path& root,
           const std::string& label,
           const std::map<std::string, StageOutcome>& outcomes) {
  const fs::path store = root / label / "checkpoints";
  const auto produced = outcomes.find(label);
  on_root(world, [&] {
    if (!fs::exists(store / "latest"))
      throw std::runtime_error("reference " + label
                               + " has no published generation in "
                               + store.string());
  });
  StageInput input{
      label, madness::checkpoint::read_published(world, store), {}, false};
  input.saved_prefix = store / input.saved.generation / "state";
  if (produced != outcomes.end()) {
    if (produced->second.state_id != input.saved.generation)
      throw std::runtime_error("reference " + label + " published "
                               + input.saved.generation + ", not "
                               + produced->second.state_id);
    input.accepted = produced->second.converged;
  } else {
    input.accepted = input.saved.stop_reason == "converged";
    if (!input.accepted)
      throw std::runtime_error("reference " + label + " generation "
                               + input.saved.generation + " stopped as "
                               + input.saved.stop_reason + ", not converged");
  }
  return input;
}

// What a consumer records about its reference ground state, from the
// published generation alone.
[[nodiscard]] inline json reference_evidence(const StageInput& input) {
  const auto report = ground_convergence_report(
      input.saved.metrics, input.saved.generation, input.accepted);
  json evidence = {{"stage", input.label},
                   {"state_id", input.saved.generation},
                   {"engine_converged", input.accepted},
                   {"energy_change", report["measurements"]["energy_change"]}};
  if (report.contains("saved_state_validation"))
    evidence["saved_state_validation"] = report["saved_state_validation"];
  return evidence;
}

// DESIGN A8: everything a run function gets. The stage directory comes in
// here, never through an engine object (A7).
struct StageContext {
  madness::World& world;
  Stage& stage;
  fs::path stage_directory;   // <root>/<stage>
  fs::path attempt_directory; // <stage_directory>/attempt-<run_id>
  std::string run_id;
  bool continuing = false;
  std::map<std::string, StageInput> inputs; // by role, e.g. "reference"
};

using RunStage = StageOutcome (*)(const StageContext&);

// The record every stage starts its results.json from.
[[nodiscard]] inline json stage_record(const StageContext& ctx) {
  const auto& cp = ctx.stage.params.get<CalculationParameters>();
  return {{"stage", ctx.stage.name},
          {"application", ctx.stage.kind},
          {"state_id", ctx.run_id + ":" + ctx.stage.name},
          {"identity_kind", "fresh_run_state"},
          {"checkpoint_generation", nullptr},
          {"requested_properties", ctx.stage.properties},
          {"protocol", cp.protocol()},
          {"output_directory", ctx.attempt_directory.string()},
          {"properties", json::object()}};
}

// The application's part of each step of run_ladder_stage (S1, D-041):
// closures over its engine and checkpoint adapter. The working form of the
// application's AppSpec (S3, D-036).
struct LadderSteps {
  // The atom count the AppSpec's fields take (see AppSpec).
  int natom = 0;
  // 2. The stage's convergence settings (B9), cadence and retention (A3).
  std::function<madness::convergence::ConvergenceSettings()> settings;
  // 3. The application's set-up before any generation is read.
  std::function<void()> prepare;
  // 4. What the plan needs that the stage does not hold, for the generation
  // read; nothing when the application uses it without planning (response:
  // another reference, or the stage loads as published).
  std::function<std::optional<PlanContext>(
      const madness::checkpoint::Manifest&)>
      plan_context;
  // 4. The load of the generation read, with its state records and their
  // plans (the AppSpec's state_records and plan) when they were planned:
  // each plan as the engine's input.
  std::function<void(const madness::checkpoint::Manifest&, const json& records,
                     const std::optional<StagePlans>& plans)>
      load;
  // 5. Run the engine with the checkpoint observer.
  std::function<void()> run_engine;
  // 6. Validate, publish, retention; the final publish takes the skeleton's
  // verdict (D-077).
  std::function<void()> finish;
  // The published generation, the work done, the iterations of this run.
  std::function<const madness::checkpoint::Manifest&()> latest;
  std::function<std::string()> work;
  std::function<std::size_t()> iterations;
  // The application's fields of the generation record.
  std::function<void(json&)> generation_record;
  // A3: the deletions retention could not make in this stage's publishes
  // (extra files only), the same on every rank.
  std::function<std::vector<std::string>()> checkpoint_warnings;
  // D-094: the record of a stage loaded whole, with no per-state plans (the
  // response's met request: request equal, converged, properties present);
  // null otherwise. Unset when the application has no such load.
  std::function<json()> whole_load;
  // A3, D-090: the archive shards the application copied into the attempt
  // (a load's copies), by attempt-relative path, deleted at the end with the
  // latest generation's; unset or empty when it copies none beyond those.
  std::function<std::vector<std::string>()> attempt_copies;
  // 9. The convergence report on the skeleton's verdict (B8), with each
  // needed state's evaluation of the written targets (`states`, the shared
  // verdict's), the saved-state outcomes of the table, and the label the
  // table ends with.
  std::function<void(json&, bool, const json& states)> report;
  std::function<std::vector<madness::convergence::MetricOutcome>()> outcomes;
  std::function<std::string(bool)> outcome_label;
  // 8. The properties and the stage's own fields.
  std::function<void(json&)> properties;
};

// DESIGN C4, skeleton spec 3.2 (S3a, D-036): an application's row, the
// functions the skeleton needs from it. Each is a function of its arguments
// (no engine object, no files) and table-tested in test_named_stages
// through <app>_fields(); <app>_spec() adds the run binding.
// `natom` is the engine fact both applications need: the atom count of the
// stage's molecule (moldft) or of its reference (response). The state
// records and the plan are fields since S4b (D-048), called by the skeleton.
struct AppSpec {
  std::string name; // the deck's block kind
  // The deck check (B3, B7); empty until S6 (D-047): parse() checks every
  // kind today.
  std::function<void(const Stage&)> validate;
  // B9: the stage's convergence settings.
  std::function<madness::convergence::ConvergenceSettings(const Stage&,
                                                          int natom)>
      settings;
  // A2: the files a generation must hold, from its engine metadata.
  std::function<std::vector<std::string>(const json& engine_metadata)> members;
  // A2: the science checks on a complete generation; throws why not.
  std::function<void(const Stage&, int natom,
                     const madness::checkpoint::Manifest&)>
      validate_generation;
  // B6: the generation's per-state records (madqc.state.v1), for the plan;
  // for a generation written before them, the records it gives.
  std::function<json(const Stage&, int natom,
                     const madness::checkpoint::Manifest&)>
      state_records;
  // B6, D-064: one plan per state the stage needs, from the records and
  // today's request, and what --continue did and why.
  std::function<StagePlans(const Stage&, int natom, const PlanContext&,
                           const json& records)>
      plan;
  // The requested property records from the engine's output (D-024:
  // unavailable values are records, not errors).
  std::function<json(const Stage&, int natom, const json& engine_output)>
      properties;
  // D-045: one engine sample (as JSON) as its record: moldft's metrics
  // record, response's work item. `carried` is what the observer keeps
  // across the samples of one run. The run binding's observer calls it
  // before the publish decision.
  std::function<json(const Stage&, int natom, const json& sample,
                     json& carried)>
      observe;
  // A3 (S4c2, D-054): whether this save boundary publishes, from the same
  // sample and `carried` (moldft: the planner's resumed rung, resumed_rung,
  // seeded by the binding), with a second output, the sample's rung record
  // (PublishDecision::rung, the manifest's position.rung; moldft only). The
  // run binding's observer calls it after observe.
  std::function<PublishDecision(const Stage&, int natom, const json& sample,
                                const json& carried)>
      publish_now;
  // A3 (S4c3, D-054, D-067): the generation's manifest, every field but the
  // shared ones (set_shared_fields) and the archives the binding writes,
  // from the application's publish sample (JSON, D-045's pattern).
  // publish_generation calls it on every rank.
  std::function<madness::checkpoint::Manifest(const Stage&, int natom,
                                              const json& publish_sample)>
      manifest;
  // B8, D-046 (S4d1): the states the stage's request needs, each with the
  // rung it must reach (the request's final rung): [{key, rung: {thresh,
  // k}, measurements}], keyed as the state records are. `measurements`
  // (D-082) names the written targets the state's kind has rows for; absent,
  // every written target.
  std::function<json(const Stage&, int natom)> requested_states;
  // B1, B8 (S4d1): one state record's saved-state measurements, keyed as the
  // written targets are (convergence::requested_targets), for the skeleton's
  // one call of evaluate_explicit_criteria.
  std::function<json(const Stage&, int natom, const json& record)>
      saved_measurements;
  // A1: what a consuming stage reads about this one besides its archives,
  // from a description of the engine; empty for an application nothing
  // consumes (response).
  std::function<json(const json& engine)> consumer_context;
  // B5: may a dependent stage bind to this outcome.
  std::function<bool(const StageOutcome&)> release;
  // A8: the run binding, which runs the stage through run_ladder_stage.
  RunStage run = nullptr;
};

// B5, the shared release: only a converged outcome is a reference.
[[nodiscard]] inline bool released(const StageOutcome& outcome) {
  return outcome.converged;
}

// B6, B8, D-046 (S4d1; the stage's verdict since S4d2, D-077): the shared
// verdict (skeleton spec 3.3 step 7), from a generation's state records
// (AppSpec::state_records, which also rebuilds those of a generation written
// before S4a): converged when every state the request needs
// (AppSpec::requested_states) has a record at its rung (or finer, R3,
// D-085), not diverged, and was either loaded under a met request
// (decided_by planner, D-016) or stopped converged and, when the written
// targets decide (B8), its saved measurements meet them. Each state is judged
// on the written targets its kind has rows for (D-082; requested_states'
// `measurements`), by one call of evaluate_explicit_criteria per state (B1);
// its result is the state's requested_accuracy, the unmet targets B8 lists
// separately. A state none of whose kind's rows is a written target is judged
// on its stop alone (its requested_accuracy not_requested), as the planner
// judges it. A written target no needed state has rows for is listed once, in
// `unmeasured_targets`. Returns {converged, states: [...],
// unmeasured_targets: [...]}.
inline json
shared_verdict(const AppSpec& spec, const Stage& stage, int natom,
               const madness::convergence::ConvergenceSettings& settings,
               const madness::checkpoint::Manifest& published) {
  namespace convergence = madness::convergence;
  const auto records = spec.state_records(stage, natom, published);
  const auto targets = convergence::requested_targets(settings);
  const bool enforced
      = settings.decider == convergence::Decider::requested_criteria;
  bool converged = true;
  json states = json::array();
  std::set<std::string> measured;
  for (const auto& needed : spec.requested_states(stage, natom)) {
    const auto& key = needed.at("key");
    json judged = json::object();
    for (const auto& [name, target] : targets.items())
      if (!needed.contains("measurements")
          || std::find(needed.at("measurements").begin(),
                       needed.at("measurements").end(), name)
                 != needed.at("measurements").end()) {
        judged[name] = target;
        measured.insert(name);
      }
    const madness::scf_continuation::Rung rung{
        needed.at("rung").at("thresh").get<double>(),
        needed.at("rung").at("k").get<int>()};
    json state = {{"key", key}, {"rung", needed.at("rung")}};
    const auto record
        = std::find_if(records.begin(), records.end(),
                       [&](const json& r) { return r.at("key") == key; });
    bool met = false;
    if (record == records.end()) {
      state["found"] = false;
    } else {
      const madness::scf_continuation::Rung at{
          record->at("rung").at("thresh").get<double>(),
          record->at("rung").at("k").get<int>()};
      const auto accuracy = madness::evaluate_explicit_criteria(
          spec.saved_measurements(stage, natom, *record), judged);
      const bool loaded = record->at("decided_by") == "planner";
      const bool stopped = record->at("stop") == "converged";
      const bool targets_met = !enforced || accuracy.at("status") == "satisfied"
                               || accuracy.at("status") == "not_requested";
      state["found"] = true;
      state["at_rung"] = madness::scf_continuation::same_rung(at, rung);
      // D-085, R3: a finished state finer than the final rung counts, as the
      // planner loads it (LadderPlanner's at_or_beyond_final); the record
      // keeps its true rung.
      state["beyond_rung"]
          = !state.at("at_rung").get<bool>() && at.thresh < rung.thresh;
      state["stop"] = record->at("stop");
      state["diverged"] = record->at("diverged");
      state["decided_by"] = record->at("decided_by");
      state["requested_accuracy"] = accuracy;
      state["enforced"] = enforced;
      met = (state.at("at_rung").get<bool>()
             || state.at("beyond_rung").get<bool>())
            && !record->at("diverged").get<bool>()
            && (loaded || (stopped && targets_met));
    }
    state["converged"] = met;
    converged = converged && met;
    states.push_back(std::move(state));
  }
  json unmeasured = json::array();
  for (const auto& [name, target] : targets.items())
    if (!measured.contains(name))
      unmeasured.push_back(name);
  return {{"converged", converged},
          {"states", std::move(states)},
          {"unmeasured_targets", std::move(unmeasured)}};
}

// A1: the generations of stage `stage` that another stage's manifest under
// the same run root names as its reference, and (D-087) the other stages'
// generations whose manifest exists but cannot be read, each with the error:
// what such a manifest references is unknown. A generation directory with no
// manifest is not a published generation and references nothing. Rank 0
// reads the other stages' manifests; every rank gets the same answer.
// Collective.
struct References {
  std::set<std::string> generations;
  std::vector<std::string> unreadable;
};
inline References referenced_generations(madness::World& world,
                                         const fs::path& store,
                                         const std::string& stage) {
  namespace checkpoint = madness::checkpoint;
  const auto text = checkpoint::detail::on_root(world, [&] {
    json names = json::array();
    json unreadable = json::array();
    const auto stage_directory = store.parent_path();
    for (const auto& other_stage :
         fs::directory_iterator(stage_directory.parent_path())) {
      const auto other = other_stage.path() / "checkpoints";
      if (other_stage.path() == stage_directory || !fs::is_directory(other))
        continue;
      for (const auto& entry : fs::directory_iterator(other)) {
        if (!checkpoint::detail::generation_name(
                entry.path().filename().string()))
          continue;
        if (!fs::exists(fs::symlink_status(entry.path() / "manifest.json")))
          continue;
        try {
          const auto reference
              = checkpoint::read_manifest(entry.path()).reference;
          if (reference.is_object() && reference.value("stage", "") == stage)
            names.push_back(reference.value("generation", ""));
        } catch (const std::exception& e) {
          unreadable.push_back(entry.path().string() + " (" + e.what() + ")");
        }
      }
    }
    return json{{"names", names}, {"unreadable", unreadable}}.dump();
  });
  const auto read = json::parse(text);
  References references;
  for (const auto& name : read.at("names"))
    references.generations.insert(name.get<std::string>());
  for (const auto& path : read.at("unreadable"))
    references.unreadable.push_back(path.get<std::string>());
  return references;
}

// A3, D-070 (a): one copy of each saved function. After a publish of
// `latest` succeeds, every other generation of the stage is deleted except
// those another stage's manifest references (A1); checkpoint.keep_all is the
// request for more (the caller's). Other stages are read only when there is
// something to delete, which after a run's first publish is always: a handful
// of manifest reads, one generation per stage being left. A store written
// before D-070 loses every older generation at its first publish, old-form
// names included. Collective; returns the deletions it could not make (extra
// files only), the same on every rank, empty on success. D-087 (a): when
// another stage's manifest exists but cannot be read, what it references is
// unknown, so nothing is deleted this time and the unread manifests are
// returned instead.
inline std::string retain_generations(madness::World& world,
                                      const fs::path& store,
                                      const std::string& stage,
                                      const std::string& latest) {
  namespace checkpoint = madness::checkpoint;
  std::vector<std::string> generations;
  for (const auto& manifest :
       checkpoint::read_generation_manifests(world, store))
    if (manifest.stage == stage)
      generations.push_back(manifest.generation);
  if (checkpoint::generations_to_delete(generations, latest).empty())
    return {};
  const auto references = referenced_generations(world, store, stage);
  if (!references.unreadable.empty()) {
    std::string unread;
    for (const auto& path : references.unreadable)
      unread += (unread.empty() ? "" : "; ") + path;
    return "nothing deleted: another stage's manifest could not be read, so "
           "what it references is unknown: "
           + unread;
  }
  const auto failed = checkpoint::remove_generations(
      world, store,
      checkpoint::generations_to_delete(generations, latest,
                                        references.generations));
  return failed.empty() ? failed
                        : "could not delete every older generation (extra "
                          "files only): "
                              + failed;
}

// A3, D-072 (a): after the stage ends, the attempt directory's copies of the
// stage's latest generation's archive shards are deleted (the engine's
// working files and the copies a load makes), so the generation is the only
// copy at rest. The published inventory names the shards by the same paths
// the attempt uses, for every application; the attempt keeps its text output
// and JSON records. D-090 (a): `copied` adds the shards a load copied in that
// the latest generation does not list (a root or state dropped from the
// request). Collective (rank 0 deletes); returns the deletions it could not
// make, the same on every rank, empty on success.
inline std::string
delete_attempt_copies(madness::World& world, const fs::path& attempt,
                      const madness::checkpoint::Manifest& latest,
                      const std::vector<std::string>& copied = {}) {
  std::set<std::string> paths(copied.begin(), copied.end());
  for (const auto& archive : latest.archives)
    for (const auto& shard : archive.shards)
      paths.insert(shard.path);
  return madness::checkpoint::detail::on_root(world, [&] {
    std::string warning;
    for (const auto& path : paths) {
      const auto copy = attempt / path;
      std::error_code error;
      if (fs::is_regular_file(fs::symlink_status(copy, error))
          && !fs::remove(copy, error) && error)
        warning += (warning.empty() ? "" : "; ") + copy.string() + ": "
                   + error.message();
    }
    return warning;
  });
}

// A3, step 6 (skeleton spec 3.3; S3c2): publish one generation of `stage`
// under `store`, in one order for every application. The application
// decided to publish and gathered its publish sample (D-051, D-054); here:
//   0. every rank: for the stage's final generation, the verdict
//      (shared_verdict, D-077); the manifest, assembled by the AppSpec's
//      manifest field from the publish sample (S4c3); then (S4c1) the
//      shared fields
//      (set_shared_fields), the name
//      the next of `names` (D-058), before anything below reads the
//      manifest: rank 0's check before writing and validate_generation both
//      compare the property records' generation with the manifest's;
//   1. rank 0: its check before writing, then the generation directory;
//   2. every rank: its archives (write_archives);
//   3. rank 0 (D-052): the archives are the generation's members, the
//      manifest has its state records (S4a; the read side accepts a
//      generation written before them) and the finished manifest passes
//      validate_generation; rank 0's manifest is the one published, so the
//      checks read it there;
//   4. checkpoint::publish, whose collective_error carries any failure from
//      2-3 to every rank (the generation then stays unpublished);
//   5. read it back; then retention (retain_generations, D-070) unless the
//      deck writes checkpoint.keep_all, its warnings added to `warnings`.
inline madness::checkpoint::Manifest
publish_generation(madness::World& world, const fs::path& store,
                   const AppSpec& spec, const Stage& stage, int natom,
                   GenerationSequence& names,
                   std::vector<std::string>& warnings, GenerationWrite write) {
  namespace checkpoint = madness::checkpoint;
  if (write.final) {
    // D-077 (S4d2): the stage's verdict, from the state records of the
    // generation about to be published (they do not depend on the verdict);
    // the final publish takes it.
    write.publish_sample["converged"]
        = shared_verdict(spec, stage, natom, spec.settings(stage, natom),
                         spec.manifest(stage, natom, write.publish_sample))
              .at("converged");
  }
  auto m = spec.manifest(stage, natom, write.publish_sample);
  set_shared_fields(m, stage.name, names.next(), write.iterations_this_run);
  const auto directory = store / m.generation;
  checkpoint::detail::on_root(world, [&] {
    if (write.check_before_write)
      write.check_before_write(m);
    fs::create_directories(store);
    if (!fs::create_directory(directory))
      throw std::runtime_error(spec.name + " checkpoint: generation "
                               + m.generation + " already exists");
    return std::string{};
  });
  auto error = write.write_archives(directory, m);
  if (world.rank() == 0 && error.empty()) {
    try {
      std::vector<std::string> written;
      written.reserve(m.archives.size());
      for (const auto& archive : m.archives)
        written.push_back(archive.basename);
      // The engine metadata the generation records, which the read side's
      // check reads too (response); moldft records none.
      auto expected
          = spec.members(m.numerical.value("engine_metadata", json::object()));
      std::sort(written.begin(), written.end());
      std::sort(expected.begin(), expected.end());
      if (written != expected)
        throw std::runtime_error(spec.name + " checkpoint: generation "
                                 + m.generation
                                 + " does not hold its members' archives");
      if (!m.numerical.contains("states"))
        throw std::runtime_error(spec.name + " checkpoint: generation "
                                 + m.generation + " has no state records");
      spec.validate_generation(stage, natom, m);
    } catch (const std::exception& e) {
      error = e.what();
    }
  }
  checkpoint::publish(world, store, m, error);
  auto published = checkpoint::read_published(world, store);
  if (!stage.checkpoint_keep_all) {
    const auto warning
        = retain_generations(world, store, stage.name, published.generation);
    if (!warning.empty()) {
      warnings.push_back(warning);
      if (world.rank() == 0)
        std::cout << "stage " << stage.name
                  << ": retention: " << warning << '\n'
                  << std::flush;
    }
  }
  return published;
}

// D-049: --continue reads a generation only when the stage has one; a stage
// with none (added to the deck after the first run, or stopped before its
// first save) starts fresh, for every application. Collective.
inline bool stage_has_generation(madness::World& world, const fs::path& store) {
  bool exists = false;
  madness::checkpoint::detail::on_root(world, [&] {
    exists = fs::exists(store / "latest");
    return std::string{};
  });
  world.gop.broadcast(exists, 0);
  return exists;
}

// The two rows, defined below the run functions they name.
inline const AppSpec& moldft_spec();
inline const AppSpec& response_spec();

// The sequence every ladder stage runs (skeleton spec 2.4 and 3.3; moldft
// and response, D-043: one order). The shared part of each step is here; the
// application's part is a LadderSteps closure. finish() decides the verdict
// before the properties are written; the continuation prints before the
// outcome table.
// D-094: the plan lines printed before a --continue stage's work, one per
// planned state (moldft's one record has no `state`: the ground state), or
// the whole load's one line.
[[nodiscard]] inline std::vector<std::string>
continuation_plan_lines(const std::string& stage, const json& continuation) {
  const auto line = [&](const std::string& state, const json& plan) {
    if (!plan.contains("rung"))
      return "stage " + stage + " plan: " + state + " - "
             + plan.at("summary").get<std::string>() + ": "
             + plan.at("why").get<std::string>();
    const auto& rung = plan.at("rung");
    return "stage " + stage + " plan: " + state + " - "
           + plan.at("action").get<std::string>() + " at rung "
           + std::to_string(rung.at("index").get<int>()) + " (thresh "
           + madness::convergence::detail::scientific(
               rung.at("thresh").get<double>())
           + ", k " + std::to_string(rung.at("k").get<int>())
           + "): " + plan.at("why").get<std::string>();
  };
  std::vector<std::string> lines;
  if (continuation.contains("states"))
    for (const auto& state : continuation.at("states"))
      lines.push_back(line(state.at("state").get<std::string>(), state));
  else
    lines.push_back(
        line(continuation.contains("rung") ? "ground state" : "the stage",
             continuation));
  return lines;
}

// D-094: the one line after a --continue stage's work: how many states were
// loaded, resumed, seeded and started fresh (moldft's one state: its
// summary), and the iterations this run.
[[nodiscard]] inline std::string
continuation_after_line(const std::string& stage, const json& continuation,
                        std::size_t iterations) {
  std::string what;
  if (continuation.contains("states")) {
    std::map<std::string, int> counts;
    for (const auto& state : continuation.at("states"))
      ++counts[state.at("action").get<std::string>()];
    for (const char* action : {"loaded", "resumed", "seeded", "fresh"})
      if (counts[action] > 0)
        what += (what.empty() ? "" : ", ") + std::to_string(counts[action])
                + " " + action;
  } else {
    what = continuation.at("summary").get<std::string>();
  }
  return "stage " + stage + " continuation: " + what + "; "
         + std::to_string(iterations) + " iterations this run";
}

inline StageOutcome run_ladder_stage(const StageContext& ctx,
                                     const AppSpec& spec,
                                     const LadderSteps& steps) {
  auto& world = ctx.world;
  auto& stage = ctx.stage;
  // 1. The record starts.
  json result = stage_record(ctx);
  // 2. B9: the settings, built once and shown before the solve; the same
  // struct feeds the record and the stopping decision.
  const auto settings = steps.settings();
  if (world.rank() == 0)
    std::cout << '\n'
              << madness::convergence::render_settings_table(settings) << '\n'
              << std::flush;
  result["convergence_settings"] = madness::convergence::to_json(settings);
  // 3-4. Read the latest generation, validate it (A2), plan its states from
  // their records (B6), load what the plans use.
  steps.prepare();
  std::optional<StagePlans> plans;
  if (ctx.continuing
      && stage_has_generation(world, ctx.stage_directory / "checkpoints")) {
    const auto latest = madness::checkpoint::read_published(
        world, ctx.stage_directory / "checkpoints");
    spec.validate_generation(stage, steps.natom, latest);
    json records = json::array();
    if (const auto context = steps.plan_context(latest)) {
      records = spec.state_records(stage, steps.natom, latest);
      plans = spec.plan(stage, steps.natom, *context, records);
    }
    steps.load(latest, records, plans);
  }
  // D-094: what --continue will do, with the recorded reasons, before the
  // work; a stage loaded whole has one record instead of per-state plans.
  json continuation;
  if (plans && !plans->states.empty())
    continuation = plans->continuation;
  else if (steps.whole_load)
    continuation = steps.whole_load();
  if (world.rank() == 0 && !continuation.is_null())
    for (const auto& line : continuation_plan_lines(stage.name, continuation))
      std::cout << line << '\n' << std::flush;
  // 5-7. Run, publish, judge.
  steps.run_engine();
  steps.finish();
  // The published generation is the stage's identity (A1).
  const auto& saved = steps.latest();
  // 7. The verdict (D-046, D-057, D-077): the shared verdict of the stage's
  // latest generation. A generation published in this run was published
  // with it; one used as published (loaded) is judged here.
  const auto verdict
      = shared_verdict(spec, stage, steps.natom, settings, saved);
  const bool converged = verdict.at("converged").get<bool>();
  const std::string state_id = saved.generation;
  result["state_id"] = state_id;
  result["identity_kind"] = "checkpoint_generation";
  result["checkpoint_generation"] = saved.generation;
  result["work"] = steps.work();
  result["iterations_this_run"] = steps.iterations();
  steps.generation_record(result);
  result["acceleration_history"] = ctx.continuing ? "reset" : "fresh";
  if (!continuation.is_null()) {
    // B6, D-094: what --continue did and why; one line for the stage.
    result["continuation"] = continuation;
    if (world.rank() == 0)
      std::cout << continuation_after_line(stage.name, continuation,
                                           steps.iterations())
                << '\n'
                << std::flush;
  }
  // 9. The outcome table and the record.
  steps.report(result, converged, verdict.at("states"));
  // D-082: a written target no needed state has rows for, reported once.
  result["convergence_report"]["unmeasured_targets"]
      = verdict.at("unmeasured_targets");
  const auto outcomes = steps.outcomes();
  result["convergence_report"]["convergence_outcome"]
      = madness::convergence::outcome_record(settings, outcomes);
  if (world.rank() == 0)
    std::cout << madness::convergence::render_outcome_table(
        settings, outcomes, steps.outcome_label(converged))
              << '\n'
              << std::flush;
  // 8. The properties.
  steps.properties(result);
  // A3, D-072, D-083 (a): the stage has ended with its latest generation,
  // published in this run or used as published; the attempt's copies of its
  // shards go, and so do retention's leftovers' warnings into the record.
  auto warnings = steps.checkpoint_warnings();
  if (const auto warning = delete_attempt_copies(
          world, ctx.attempt_directory, saved,
          steps.attempt_copies ? steps.attempt_copies()
                               : std::vector<std::string>{});
      !warning.empty())
    warnings.push_back("attempt copies: " + warning);
  for (const auto& name : verdict.at("unmeasured_targets"))
    warnings.push_back("written target " + name.get<std::string>()
                       + ": no state this stage needs has a row for it, so "
                         "it judged nothing");
  result["checkpoint_warnings"] = warnings;
  // 10. What the stage hands on (A8).
  StageOutcome outcome;
  outcome.state_id = state_id;
  outcome.converged = converged;
  outcome.result = std::move(result);
  return outcome;
}

// A named moldft stage: the SCF engine with explicit paths, published through
// MoldftCheckpoint, run by run_ladder_stage.
inline StageOutcome run_moldft_stage(const StageContext& ctx) {
  auto& world = ctx.world;
  auto& stage = ctx.stage;
  const auto& directory = ctx.attempt_directory;
  auto& cp = stage.params.get<CalculationParameters>();
  // Run the normal SCF engine with explicit paths; never change cwd.
  moldft_lib library((directory / "resolved.in").string());
  cp.set_user_defined_value("save",
                            false); // immutable adapter owns persistence
  auto scf = library.calc(world, stage.params);
  std::optional<madness::convergence::ConvergenceSettings> settings;
  std::optional<MoldftCheckpoint> checkpoint;
  std::optional<decltype(library.run(world, stage.params, NextAction::Redo))>
      values;
  LadderSteps steps;
  const auto natom = static_cast<int>(scf->molecule.natom());
  steps.natom = natom;
  steps.settings = [&] {
    settings = moldft_spec().settings(stage, natom);
    return *settings;
  };
  GenerationSequence names{stage.name, ctx.run_id};
  std::vector<std::string> warnings;
  steps.checkpoint_warnings = [&] { return warnings; };
  steps.prepare = [&] {
    checkpoint.emplace(world, stage, ctx.stage_directory, ctx.run_id, *settings,
                       [&](GenerationWrite write) {
                         return publish_generation(
                             world, ctx.stage_directory / "checkpoints",
                             moldft_spec(), stage, natom, names, warnings,
                             std::move(write));
                       });
    checkpoint->prepare(*scf);
  };
  // moldft plans every generation it reads; its problem is the live SCF's.
  steps.plan_context = [&](const madness::checkpoint::Manifest&) {
    return std::optional<PlanContext>(PlanContext{checkpoint->identity(), {}});
  };
  steps.load = [&](const madness::checkpoint::Manifest& m, const json&,
                   const std::optional<StagePlans>& plans) {
    checkpoint->load(*scf, m, plans.value());
  };
  steps.run_engine = [&] {
    // D-045: what observe carries across the samples of this run; S4c2:
    // the rung the plan resumed, for publish_now.
    json carried = json::object();
    if (auto resumed = checkpoint->resumed_rung(); !resumed.is_null())
      carried["resumed_rung"] = std::move(resumed);
    scf->checkpoint_observer
        = [&](madness::World&, const madness::SCF& state,
              const madness::SCF::IterationCheckpoint& step) {
            const auto sample = checkpoint->sample(state, step);
            auto metrics = moldft_spec().observe(stage, natom, sample, carried);
            checkpoint->observe(
                state, step, metrics,
                moldft_spec().publish_now(stage, natom, sample, carried));
          };
    values.emplace(library.run(world, stage.params, NextAction::Redo));
    scf->checkpoint_observer = {};
  };
  steps.finish = [&] {
    const auto& evidence = std::get<2>(*values);
    // The engine's own verdict at the final rung (it records its claim only
    // when its check passed, Q6); finish() lets the selected group decide.
    const bool claimed = evidence.converged_for_thresh <= cp.protocol().back()
                         && evidence.converged_for_dconv
                                <= std::max(cp.dconv(), cp.protocol().back());
    checkpoint->finish(*scf, claimed);
  };
  steps.latest = [&]() -> const madness::checkpoint::Manifest& {
    return checkpoint->latest();
  };
  steps.work = [&] { return checkpoint->work; };
  steps.iterations = [&] { return checkpoint->iterations; };
  steps.generation_record = [](json&) {};
  const auto saved_measurements = [&] {
    return checkpoint->latest()
        .metrics.value("convergence_report", json::object())
        .value("saved_state_validation", json::object())
        .value("measurements", json::object());
  };
  steps.report = [&](json& result, bool converged, const json& states) {
    const auto& saved = checkpoint->latest();
    result["checkpoint_metrics"] = saved.metrics;
    result["convergence_report"]
        = ground_convergence_report(saved.metrics, saved.generation, converged);
    // B8: the current request judged on the saved state, also when the
    // generation was reused from an earlier request: the shared verdict's
    // evaluation of the ground state.
    result["convergence_report"]["requested_accuracy"]
        = requested_criteria_report(states.at(0).at("requested_accuracy"),
                                    *settings);
  };
  steps.outcomes = [&] { return saved_outcomes(saved_measurements()); };
  steps.outcome_label = [&](bool) { return checkpoint->latest().stop_reason; };
  steps.properties = [&](json& result) {
    const auto& saved = checkpoint->latest();
    const auto& evidence = std::get<2>(*values);
    result["method"] = "HF";
    result["convergence"] = evidence.to_json();
    result["targets"] = {{"dconv", cp.dconv()}};
    if (const auto econv = requested_econv(cp)) // D8: user-written only
      result["targets"]["econv"] = *econv;
    result["parameters"] = cp.to_json_if_precedence("defined");
    result["geometry"] = scf->molecule.to_json();
    result["representation"] = {{"L", scf->param.L()},
                                {"k", scf->amo.at(0).k()},
                                {"threshold", saved.numerical.at("threshold")}};
    // These properties were evaluated on the exact archived orbital state.
    auto selected = moldft_spec().properties(
        stage, natom, {{"properties", saved.properties}});
    for (auto& [name, record] : selected.items()) {
      record["source_state_id"] = saved.generation;
      result["properties"][name] = std::move(record);
    }
  };
  return run_ladder_stage(ctx, moldft_spec(), steps);
}

// A named molresponse stage: the response engine against its declared
// reference, published through ResponseCheckpoint, run by run_ladder_stage.
inline StageOutcome run_response_stage(const StageContext& ctx) {
  auto& world = ctx.world;
  auto& stage = ctx.stage;
  const auto& directory = ctx.attempt_directory;
  auto& cp = stage.params.get<CalculationParameters>();
  // DESIGN A1: the reference is its published generation: the orbital
  // archive, and the consumer context moldft wrote beside it (the molecule
  // as the SCF oriented it, the cell, the final thresh and k).
  const StageInput& reference = ctx.inputs.at("reference");
  if (!reference.saved.numerical.contains("consumer_context"))
    throw std::runtime_error(
        "reference " + reference.label + " generation "
        + reference.saved.generation
        + " has no consumer context (published before DESIGN A1); "
          "run that stage again");
  const auto& context = reference.saved.numerical.at("consumer_context");
  if (!context.at("spin_restricted").get<bool>()
      || context.at("multiplicity").get<int>() != 1)
    throw std::runtime_error(
        "named stages currently require a closed-shell reference");
  // Reject downstream tightening instead of silently refining the source.
  if (cp.protocol().back() < context.at("thresh").get<double>()
      || protocol_k(cp.protocol().back()) > context.at("k").get<int>())
    throw std::runtime_error("response protocol is tighter than reference "
                             + stage.reference);
  madness::Molecule molecule;
  molecule.from_json(context.at("molecule"));
  cp.set_user_defined_value("l", context.at("L").get<double>());
  std::optional<madness::convergence::ConvergenceSettings> settings;
  std::optional<ResponseCheckpoint> checkpoint;
  molresponse_v3_lib::Results response;
  bool loaded = false;
  // The selected property records, chosen by finish() on root and written
  // by properties() once the generation is known.
  json selected = json::object();
  LadderSteps steps;
  // B9: the solver stops on the same targets.
  const auto natom = static_cast<int>(molecule.natom());
  steps.natom = natom;
  steps.settings = [&] {
    settings = response_spec().settings(stage, natom);
    return *settings;
  };
  GenerationSequence names{stage.name, ctx.run_id};
  std::vector<std::string> warnings;
  steps.checkpoint_warnings = [&] { return warnings; };
  steps.attempt_copies = [&] { return checkpoint->attempt_copies(); };
  steps.whole_load = [&] { return checkpoint->whole_load(); };
  steps.prepare = [&] {
    checkpoint.emplace(
        world, stage, ctx.stage_directory, directory,
        json{{"stage", stage.reference},
             {"generation", reference.saved.generation}},
        static_cast<int>(molecule.natom()), [&](GenerationWrite write) {
          return publish_generation(world, ctx.stage_directory / "checkpoints",
                                    response_spec(), stage, natom, names,
                                    warnings, std::move(write));
        });
  };
  steps.plan_context = [&](const madness::checkpoint::Manifest& m) {
    if (!checkpoint->accept(m))
      return std::optional<PlanContext>();
    return std::optional<PlanContext>(
        PlanContext{{}, reference.saved.generation});
  };
  steps.load = [&](const madness::checkpoint::Manifest&, const json& records,
                   const std::optional<StagePlans>& plans) {
    loaded = checkpoint->load(records, plans);
  };
  steps.run_engine = [&] {
    if (loaded) {
      response.metadata = checkpoint->metadata();
      response.properties = checkpoint->properties();
      response.vibrational_analysis
          = checkpoint->analysis("vibrational_analysis");
      response.raman_spectra = checkpoint->analysis("raman_spectra");
    } else {
      response = molresponse_v3_lib::run_response_on_reference(
          world, stage.params, molecule, directory, reference.saved_prefix,
          [&](madness::World& w, const json& event) {
            json carried;
            auto work_item
                = response_spec().observe(stage, natom, event, carried);
            checkpoint->observe(
                w, event, work_item,
                response_spec().publish_now(stage, natom, event, carried));
          },
          checkpoint->state_plans());
    }
  };
  steps.finish = [&] {
    // The engine's property/evidence JSON exists on root only.
    on_root(world, [&] {
      // D-025: the bundle's verdict says who decided it, and its transition
      // analysis is the engine's file of this attempt, or the generation's
      // on a load.
      checkpoint->annotate_bundles(response.metadata);
      json transition = loaded ? checkpoint->transition() : json();
      const auto analysis
          = directory / ("es_analysis__" + final_protocol_key(stage) + ".json");
      if (!loaded && fs::exists(analysis))
        transition = madness::checkpoint::detail::read_json(analysis);
      selected = response_spec().properties(
          stage, natom,
          {{"properties", response.properties},
           {"metadata", response.metadata},
           {"vibrational_analysis", response.vibrational_analysis},
           {"raman_spectra", response.raman_spectra},
           {"transition", transition}});
    });
    checkpoint->finish(response.metadata);
  };
  steps.latest = [&]() -> const madness::checkpoint::Manifest& {
    return checkpoint->latest();
  };
  steps.work = [&] { return checkpoint->work; };
  steps.iterations = [&] { return checkpoint->iterations; };
  steps.generation_record = [&](json& result) {
    result["work_items"] = checkpoint->work_items;
    result["reference_changed"] = checkpoint->reference_changed;
  };
  steps.report = [&](json& result, bool, const json& states) {
    result["convergence_report"] = response_convergence_report(
        checkpoint->metadata(), checkpoint->latest().generation,
        reference_evidence(reference));
    result["convergence_report"]["decided_by"]
        = madness::convergence::to_string(settings->decider);
    // B8: each needed state's written targets judged on its saved state,
    // the shared verdict's evaluation, by state key.
    auto& accuracy = result["convergence_report"]["requested_accuracy"];
    accuracy = json::object();
    for (const auto& state : states)
      if (state.contains("requested_accuracy"))
        accuracy[state.at("key").get<std::string>()]
            = requested_criteria_report(state.at("requested_accuracy"),
                                        *settings);
  };
  steps.outcomes = [&] {
    const auto final_protocol = molresponse_v3::protocol_key(
        cp.protocol().back(), protocol_k(cp.protocol().back()));
    return response_saved_outcomes(checkpoint->metadata(), final_protocol);
  };
  steps.outcome_label = [](bool converged) -> std::string {
    return converged ? "converged" : "unconverged";
  };
  steps.properties = [&](json& result) {
    result["reference_representation"]
        = {{"L", context.at("L")}, {"k", context.at("k")}};
    result["method"] = "TDHF";
    result["reference"] = {{"stage", stage.reference},
                           {"state_id", reference.saved.generation}};
    result["targets"]
        = {{"dconv", stage.params.get<ResponseParameters>().dconv()}};
    result["parameters"]
        = stage.params.get<ResponseParameters>().to_json_if_precedence(
            "defined");
    on_root(world, [&] {
      for (auto& [name, record] : selected.items()) {
        record["source_state_id"] = checkpoint->latest().generation;
        record["reference_state_id"] = reference.saved.generation;
        result["properties"][name] = std::move(record);
      }
      result["engine_metadata"] = response.metadata;
    });
  };
  return run_ladder_stage(ctx, response_spec(), steps);
}

// moldft's row (D-040): the SCF stage. The fields alone, without the run
// binding, are what test_named_stages tables: it links no engine.
inline AppSpec moldft_fields() {
  return {.name = "moldft",
          .validate = {},
          .settings = moldft_stage_settings,
          .members = moldft_members,
          .validate_generation =
              [](const Stage& stage, int,
                 const madness::checkpoint::Manifest& m) {
                MoldftCheckpoint::validate_generation(stage.name, m);
              },
          .state_records = moldft_state_records,
          .plan = plan_moldft_state,
          .properties = moldft_properties,
          .observe = moldft_observe,
          .publish_now = moldft_publish_now,
          .manifest = moldft_manifest,
          .requested_states = moldft_requested_states,
          .saved_measurements = moldft_saved_measurements,
          .consumer_context = moldft_consumer_context,
          .release = released};
}
inline const AppSpec& moldft_spec() {
  static const AppSpec spec = [] {
    auto row = moldft_fields();
    row.run = run_moldft_stage;
    return row;
  }();
  return spec;
}

// molresponse's row (D-040): response properties against a reference.
inline AppSpec response_fields() {
  return {
      .name = "molresponse",
      .validate = {},
      .settings =
          [](const Stage& stage, int) {
            return response_convergence_settings(
                stage,
                stage.params.get<CalculationParameters>().protocol().back());
          },
      .members = ResponseCheckpoint::members,
      .validate_generation =
          [](const Stage& stage, int natom,
             const madness::checkpoint::Manifest& m) {
            ResponseCheckpoint::validate_generation(stage.name, natom, m);
          },
      .state_records = response_generation_records,
      .plan =
          [](const Stage& stage, int natom, const PlanContext& context,
             const json& records) {
            return plan_response_states(stage, response_request(stage), natom,
                                        context.reference_generation, records);
          },
      .properties =
          [](const Stage& stage, int natom, const json& output) {
            return select_response_properties(
                stage, output.at("properties"), output.at("metadata"),
                output.at("vibrational_analysis"), output.at("raman_spectra"),
                natom, output.at("transition"));
          },
      .observe = response_observe,
      // Every save event publishes; response's position has no rung.
      .publish_now = [](const Stage&, int, const json&,
                        const json&) { return PublishDecision{true, nullptr}; },
      .manifest = ResponseCheckpoint::manifest,
      .requested_states = ResponseCheckpoint::requested_states,
      .saved_measurements = ResponseCheckpoint::saved_measurements,
      .consumer_context = {},
      .release = released};
}
inline const AppSpec& response_spec() {
  static const AppSpec spec = [] {
    auto row = response_fields();
    row.run = run_response_stage;
    return row;
  }();
  return spec;
}

/// Returns false for the existing unnamed input surface. Missing default input
/// is left to that path, which also accepts geometry supplied on the command
/// line.
inline bool run_if_named(madness::World& world,
                         const madness::commandlineparser& parser) {
  // 1. Read once and give every rank the same input text.
  std::string text;
  on_root(world, [&] {
    std::ifstream input(parser.value("input"));
    if (!input)
      return;
    std::ostringstream contents;
    contents << input.rdbuf();
    if (input.bad())
      throw std::runtime_error("reading named input failed");
    text = contents.str();
  });
  world.gop.broadcast_serializable(text, 0);
  if (!named_input(text))
    return false;

  // 2. Settings belong to named blocks, so global scientific overrides are
  // rejected.
  for (const auto& option : parser.keyval) {
    if (option.first != "input" && option.first != "prefix"
        && option.first != "continue"
        && option.first != "user_defined_input_file"
        && option.first != "user_defined_prefix")
      throw std::runtime_error("named stages: put settings inside their stage; "
                               "unsupported option --"
                               + option.first);
  }
  auto stages = parse(
      world, text); // validate the complete deck before creating outputs

  const bool continuing = parser.key_exists("continue");

  // 3. Reopening an output root requires the explicit --continue option.
  // Each invocation gets separate working files; checkpoints remain immutable.
  fs::path root;
  std::string root_string, run_id;
  on_root(world, [&] {
    root = fs::absolute(parser.key_exists("user_defined_prefix")
                            ? fs::path(parser.value_raw("prefix"))
                            : fs::path(parser.value("input")).stem())
               .lexically_normal();
    // A1: a reference outside the deck is read from a published generation
    // under this output root; refuse one that has none before any stage runs
    // (and before a fresh root is created).
    std::string deck_stages;
    for (const auto& stage : stages)
      deck_stages += (deck_stages.empty() ? "" : ", ") + stage.name;
    for (const auto& stage : stages)
      if (!stage.reference.empty()
          && std::none_of(stages.begin(), stages.end(),
                          [&](const Stage& s) {
                            return s.name == stage.reference;
                          })
          && !fs::exists(root / stage.reference / "checkpoints" / "latest"))
        throw std::runtime_error(
            "stage " + stage.name + ": reference " + stage.reference
            + " names no stage of this deck (" + deck_stages
            + ") and no published generation under " + root.string()
            + "; check the name");
    fs::create_directories(root.parent_path());
    if (continuing ? !fs::is_directory(root) : !fs::create_directory(root))
      throw std::runtime_error(
          "fresh named-stage run requires a new output directory: "
          + root.string());
    root_string = root.string();
    run_id = std::to_string(
                 std::chrono::system_clock::now().time_since_epoch().count())
             + "-" + std::to_string(getpid());
    for (const auto& stage : stages)
      if (!fs::create_directory(root / stage.name) && !continuing)
        throw std::runtime_error("stage output directory collision: "
                                 + stage.name);
  });
  world.gop.broadcast_serializable(root_string, 0);
  world.gop.broadcast_serializable(run_id, 0);
  root = root_string;

  // A8, D-036: one AppSpec per application, looked up by its name; each
  // stage's outcome is kept by its label for the stages that consume it.
  const std::map<std::string, AppSpec> applications
      = {{"moldft", moldft_spec()}, {"molresponse", response_spec()}};
  std::map<std::string, StageOutcome> outcomes;
  json index = {{"schema", "madqc.named-stages.v1"},
                {"run_id", run_id},
                {"mode", continuing ? "continue" : "fresh"},
                {"mpi_ranks", world.size()},
                {"stages", json::array()}};

  // 4. Execute in input order; parse() has already checked the references.
  for (auto& stage : stages) {
    const fs::path stage_directory = root / stage.name;
    const fs::path directory = stage_directory / ("attempt-" + run_id);
    on_root(world, [&] {
      if (!fs::create_directory(directory))
        throw std::runtime_error("stage attempt directory collision");
    });
    stage.params.get<CalculationParameters>().set_user_defined_value(
        "prefix", (directory / "state").string());
    StageContext ctx{world,      stage, stage_directory, directory, run_id,
                     continuing, {}};
    // The inputs this stage's block declares (today: a response's reference),
    // each read from its producer's published generation (A1).
    if (!stage.reference.empty())
      ctx.inputs.emplace("reference",
                         read_input(world, root, stage.reference, outcomes));
    const auto& application = applications.at(stage.kind);
    auto outcome = application.run(ctx);

    // 5. Record the outcome before stopping an unconverged sequence. Later
    // stages must not run against a reference that failed its completion check.
    auto& result = outcome.result;
    result["status"] = outcome.converged ? "completed" : "unconverged";
    write_result(world, stage_directory / "results.json", result);
    index["stages"].push_back(result);
    write_result(world, root / "results.json", index);
    if (!application.release(outcome))
      throw std::runtime_error("stage " + stage.name
                               + " did not satisfy its requested convergence "
                                 "criteria; see results.json");
    outcomes.emplace(stage.name, std::move(outcome));
  }
  if (world.rank() == 0)
    madness::print("Named-stage results:", (root / "results.json").string());
  return true;
}
} // namespace named_stages
