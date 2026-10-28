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
      throw std::runtime_error(
          "reference " + label + " has no published generation in "
          + store.string()
          + (produced == outcomes.end()
                 ? "; a reference outside the deck is read from this run's "
                   "output directory, so continue that run (--continue)"
                 : ""));
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
  // 2. The stage's convergence settings (B9), cadence and retention (A3).
  std::function<madness::convergence::ConvergenceSettings()> settings;
  // 3. The application's set-up before any generation is read.
  std::function<void()> prepare;
  // 3. A2: the checks on the generation read (the AppSpec's
  // validate_generation); throws why not.
  std::function<void(const madness::checkpoint::Manifest&)> validate_generation;
  // 4. The state records and the plan (B6), from the generation read, and
  // the load of what the plan uses.
  std::function<void(const madness::checkpoint::Manifest&)> load;
  // 5. Run the engine with the checkpoint observer.
  std::function<void()> run_engine;
  // 6-7. Validate, publish, retention; returns the stage's verdict.
  std::function<bool()> finish;
  // The published generation, the work done, the iterations of this run.
  std::function<const madness::checkpoint::Manifest&()> latest;
  std::function<std::string()> work;
  std::function<std::size_t()> iterations;
  // The application's fields of the generation record.
  std::function<void(json&)> generation_record;
  // What --continue did and why (B6), when there is a record of it.
  std::function<std::optional<json>()> continuation;
  // 9. The convergence report on the verdict finish returned (B8), the
  // saved-state outcomes of the table, and the label the table ends with.
  std::function<void(json&, bool)> report;
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
// records and the plan are functions of their arguments since S3b, called
// by the run bindings; their fields come at S4 (D-048).
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

// A3, step 6 (skeleton spec 3.3; S3c2): publish one generation of `stage`
// under `store`, in one order for every application. The application
// assembled and named the manifest (D-051, D-054); here:
//   1. rank 0: its check before writing, then the generation directory;
//   2. every rank: its archives (write_archives);
//   3. rank 0 (D-052): the archives are the generation's members and the
//      finished manifest passes validate_generation; rank 0's manifest is
//      the one published, so the checks read it there;
//   4. checkpoint::publish, whose collective_error carries any failure from
//      2-3 to every rank (the generation then stays unpublished);
//   5. read it back; then `retain`, the application's retention (D-053).
inline madness::checkpoint::Manifest publish_generation(
    madness::World& world, const fs::path& store, const AppSpec& spec,
    const Stage& stage, int natom, GenerationWrite write,
    const std::function<void(const madness::checkpoint::Manifest&)>& retain) {
  namespace checkpoint = madness::checkpoint;
  auto& m = write.manifest;
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
      auto expected = spec.members(write.engine_metadata);
      std::sort(written.begin(), written.end());
      std::sort(expected.begin(), expected.end());
      if (written != expected)
        throw std::runtime_error(spec.name + " checkpoint: generation "
                                 + m.generation
                                 + " does not hold its members' archives");
      spec.validate_generation(stage, natom, m);
    } catch (const std::exception& e) {
      error = e.what();
    }
  }
  checkpoint::publish(world, store, m, error);
  auto published = checkpoint::read_published(world, store);
  if (retain)
    retain(published);
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
inline StageOutcome run_ladder_stage(const StageContext& ctx,
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
  // 3-4. Read the latest generation, validate it, plan (B6), load if met.
  steps.prepare();
  if (ctx.continuing
      && stage_has_generation(world, ctx.stage_directory / "checkpoints")) {
    const auto latest = madness::checkpoint::read_published(
        world, ctx.stage_directory / "checkpoints");
    steps.validate_generation(latest);
    steps.load(latest);
  }
  // 5-7. Run, publish, judge.
  steps.run_engine();
  const bool converged = steps.finish();
  // The published generation is the stage's identity (A1).
  const auto& saved = steps.latest();
  const std::string state_id = saved.generation;
  result["state_id"] = state_id;
  result["identity_kind"] = "checkpoint_generation";
  result["checkpoint_generation"] = saved.generation;
  result["work"] = steps.work();
  result["iterations_this_run"] = steps.iterations();
  steps.generation_record(result);
  result["acceleration_history"] = ctx.continuing ? "reset" : "fresh";
  if (const auto continuation = steps.continuation()) {
    // B6: what --continue did and why.
    result["continuation"] = *continuation;
    // One summary for the stage, or one per state.
    if (world.rank() == 0) {
      if (continuation->contains("states"))
        for (const auto& state : continuation->at("states"))
          madness::print("stage", stage.name,
                         "continuation:", state.at("state").get<std::string>(),
                         "-", state.at("summary").get<std::string>());
      else
        madness::print("stage", stage.name, "continuation:",
                       continuation->at("summary").get<std::string>());
    }
  }
  // 9. The outcome table and the record.
  steps.report(result, converged);
  const auto outcomes = steps.outcomes();
  result["convergence_report"]["convergence_outcome"]
      = madness::convergence::outcome_to_json(settings, outcomes);
  if (world.rank() == 0)
    std::cout << madness::convergence::render_outcome_table(
        settings, outcomes, steps.outcome_label(converged))
              << '\n'
              << std::flush;
  // 8. The properties.
  steps.properties(result);
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
  steps.settings = [&] {
    settings = moldft_spec().settings(stage, natom);
    return *settings;
  };
  steps.prepare = [&] {
    checkpoint.emplace(
        world, stage, ctx.stage_directory, ctx.run_id, *settings,
        [&](GenerationWrite write) {
          return publish_generation(
              world, ctx.stage_directory / "checkpoints", moldft_spec(), stage,
              natom, std::move(write),
              [&](const madness::checkpoint::Manifest& published) {
                // A3: retention unless the deck keeps every generation.
                if (!settings->checkpoints.keep_all)
                  checkpoint->apply_retention(published.generation);
              });
        });
    checkpoint->prepare(*scf);
  };
  steps.validate_generation = [&](const madness::checkpoint::Manifest& m) {
    moldft_spec().validate_generation(stage, natom, m);
  };
  bool read_generation = false;
  steps.load = [&](const madness::checkpoint::Manifest& m) {
    read_generation = true;
    checkpoint->load(*scf, m);
  };
  steps.run_engine = [&] {
    // D-045: what observe carries across the samples of this run.
    json carried = json::object();
    scf->checkpoint_observer
        = [&](madness::World&, const madness::SCF& state,
              const madness::SCF::IterationCheckpoint& step) {
            checkpoint->observe(
                state, step,
                moldft_spec().observe(
                    stage, natom, checkpoint->sample(state, step), carried));
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
    // A loaded state was judged against this request by prepare() and met
    // it (B6); a solved state is converged when its stop reason says so.
    return checkpoint->work == "loaded"
           || checkpoint->latest().stop_reason == "converged";
  };
  steps.latest = [&]() -> const madness::checkpoint::Manifest& {
    return checkpoint->latest();
  };
  steps.work = [&] { return checkpoint->work; };
  steps.iterations = [&] { return checkpoint->iterations; };
  steps.generation_record = [&](json& result) {
    result["checkpoint_warnings"] = checkpoint->retention_warnings();
  };
  steps.continuation = [&]() -> std::optional<json> {
    if (!read_generation)
      return std::nullopt;
    return checkpoint->continuation_record(cp);
  };
  const auto saved_measurements = [&] {
    return checkpoint->latest()
        .metrics.value("convergence_report", json::object())
        .value("saved_state_validation", json::object())
        .value("measurements", json::object());
  };
  steps.report = [&](json& result, bool converged) {
    const auto& saved = checkpoint->latest();
    result["checkpoint_metrics"] = saved.metrics;
    result["convergence_report"]
        = ground_convergence_report(saved.metrics, saved.generation, converged);
    // The current request, judged on the saved state (B6, B8), also when the
    // generation was reused from an earlier request.
    result["convergence_report"]["requested_accuracy"]
        = requested_criteria_report(saved_measurements(), *settings);
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
  return run_ladder_stage(ctx, steps);
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
  steps.settings = [&] {
    settings = response_spec().settings(stage, natom);
    return *settings;
  };
  steps.prepare = [&] {
    checkpoint.emplace(
        world, stage, ctx.stage_directory, directory, ctx.run_id,
        json{{"stage", stage.reference},
             {"generation", reference.saved.generation}},
        static_cast<int>(molecule.natom()), [&](GenerationWrite write) {
          // Response keeps every generation (D-053).
          return publish_generation(world, ctx.stage_directory / "checkpoints",
                                    response_spec(), stage, natom,
                                    std::move(write), {});
        });
  };
  steps.validate_generation = [&](const madness::checkpoint::Manifest& m) {
    response_spec().validate_generation(stage, natom, m);
  };
  steps.load = [&](const madness::checkpoint::Manifest& m) {
    loaded = checkpoint->load(m);
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
            checkpoint->observe(
                w, event,
                response_spec().observe(stage, natom, event, carried));
          },
          checkpoint->state_plans());
    }
  };
  steps.finish = [&] {
    bool converged = false;
    // The engine's property/evidence JSON exists on root only. The verdict
    // is whether the requested properties converged.
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
      converged = all_requested_converged(selected);
    });
    world.gop.broadcast(converged, 0);
    checkpoint->finish(response.metadata, converged);
    return converged;
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
  // B6: what --continue did per state and why.
  steps.continuation = [&]() -> std::optional<json> {
    if (checkpoint->state_plans().empty())
      return std::nullopt;
    return checkpoint->continuation();
  };
  steps.report = [&](json& result, bool) {
    result["convergence_report"] = response_convergence_report(
        checkpoint->metadata(), checkpoint->latest().generation,
        reference_evidence(reference));
    result["convergence_report"]["decided_by"]
        = madness::convergence::to_string(settings->decider);
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
  return run_ladder_stage(ctx, steps);
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
          .properties = moldft_properties,
          .observe = moldft_observe,
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
      .properties =
          [](const Stage& stage, int natom, const json& output) {
            return select_response_properties(
                stage, output.at("properties"), output.at("metadata"),
                output.at("vibrational_analysis"), output.at("raman_spectra"),
                natom, output.at("transition"));
          },
      .observe = response_observe,
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
