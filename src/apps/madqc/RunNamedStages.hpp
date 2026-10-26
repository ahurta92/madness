#pragma once
#include "MoldftCheckpoint.hpp"
#include "ConvergenceReport.hpp"
#include "ResponseCheckpoint.hpp"
#include "NamedStages.hpp"
#include <apps/molresponse/madqc_adapter.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <madness/chem/Applications.hpp>
#include <unistd.h>

// Execute the validated stages from NamedStages.hpp using the existing engines.
// Flow: read/validate -> create output directories -> run stages -> write
// reports. Application adapters publish immutable ground and response checkpoints.
namespace named_stages {
namespace fs = std::filesystem;
using json = nlohmann::json;

// All ranks call this helper. Rank 0 performs the operation, then broadcasts
// any standard exception so every rank takes the same error path.
template <class F> void on_root(madness::World &world, F operation) {
  std::string error;
  if (world.rank() == 0) {
    try {
      operation();
    } catch (const std::exception &e) {
      error = e.what();
    }
  }
  world.gop.broadcast_serializable(error, 0);
  if (!error.empty())
    throw std::runtime_error(error);
}
// Write an ordinary result report. This is not atomic checkpoint publication.
inline void write_result(madness::World &world, const fs::path &file,
                         const json &value) {
  on_root(world, [&] {
    std::ofstream output(file);
    output.exceptions(std::ios::failbit | std::ios::badbit);
    output << value.dump(2) << '\n';
    output.close();
  });
}
// Select requested frequencies/components from the engine's final-rung tensor.
// The response calculation has already run; this function only formats results.
inline json selected_alpha(const Stage &stage, const json &properties) {
  const auto &cp = stage.params.get<CalculationParameters>();
  const auto pkey = molresponse_v3::protocol_key(
      cp.protocol().back(), protocol_k(cp.protocol().back()));
  const auto &rows = properties.at("alpha").at(pkey);
  json selected = json::array();
  for (double omega : stage.frequencies) {
    const json *found = nullptr;
    for (const auto &row : rows) {
      if (std::abs(row.at("omega").get<double>() - omega) < 1e-12) {
        if (found)
          throw std::runtime_error(
              "duplicate alpha frequency in engine results");
        found = &row;
      }
    }
    if (!found)
      throw std::runtime_error("missing requested polarizability frequency");
    const std::string dirs = found->at("directions");
    json values = json::object();
    for (const auto &component : stage.components) {
      const auto i = dirs.find(component[0]), j = dirs.find(component[1]);
      if (i == std::string::npos || j == std::string::npos)
        throw std::runtime_error("missing requested alpha component");
      const double value = found->at("alpha").at(i).at(j).get<double>();
      if (!std::isfinite(value))
        throw std::runtime_error("nonfinite polarizability");
      values[component] = value;
    }
    selected.push_back({{"frequency", omega},
                        {"frequency_units", "au"},
                        {"units", "a0^3"},
                        {"components", values},
                        {"converged", found->value("converged", false)},
                        {"engine_evidence", *found}});
  }
  return selected;
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

// A named moldft stage: the SCF engine with explicit paths, published through
// MoldftCheckpoint.
inline StageOutcome run_moldft_stage(const StageContext& ctx) {
  // The names below are those of the former loop body in run_if_named, so
  // the code reads as it did there (A8: one run function per application).
  auto& world = ctx.world;
  auto& stage = ctx.stage;
  const auto& stage_directory = ctx.stage_directory;
  const auto& directory = ctx.attempt_directory;
  const auto& run_id = ctx.run_id;
  const bool continuing = ctx.continuing;
  auto& cp = stage.params.get<CalculationParameters>();
  json result = stage_record(ctx);
  std::string state_id = result.at("state_id").get<std::string>();
  bool converged = false;
  StageOutcome outcome;
  // Run the normal SCF engine with explicit paths; never change cwd.
  moldft_lib library((directory / "resolved.in").string());
  cp.set_user_defined_value("save",
                            false); // immutable adapter owns persistence
  auto scf = library.calc(world, stage.params);
  // B9: the stage's convergence settings, built once and shown before the
  // solve; the same struct feeds the record and the stopping decision.
  auto settings
      = moldft_convergence_settings(stage.name, cp, scf->molecule.natom());
  // A3: the observer's cadence and retention, shown with who set them.
  const auto set_by = [&](const char* key) {
    return stage.written.contains(key)
               ? madness::convergence::SetBy::user
               : madness::convergence::SetBy::engine_default;
  };
  settings.checkpoints = {
      stage.checkpoint_every_iterations, set_by("checkpoint.every_iterations"),
      stage.checkpoint_keep_all, set_by("checkpoint.keep_all"), ""};
  if (world.rank() == 0)
    std::cout << '\n'
              << madness::convergence::render_settings_table(settings) << '\n'
              << std::flush;
  result["convergence_settings"] = madness::convergence::to_json(settings);
  MoldftCheckpoint checkpoint(world, stage, stage_directory, run_id, settings);
  checkpoint.prepare(*scf, continuing);
  scf->checkpoint_observer
      = [&](madness::World&, const madness::SCF& state,
            const madness::SCF::IterationCheckpoint& step) {
          checkpoint.observe(state, step);
        };
  auto values = library.run(world, stage.params, NextAction::Redo);
  scf->checkpoint_observer = {};
  const auto& evidence = std::get<2>(values);

  // The engine's own verdict at the final rung (it records its claim only
  // when its check passed, Q6); finish() lets the selected group decide.
  converged = evidence.converged_for_thresh <= cp.protocol().back()
              && evidence.converged_for_dconv
                     <= std::max(cp.dconv(), cp.protocol().back());
  checkpoint.finish(*scf, converged);
  const auto& saved = checkpoint.latest();
  // A loaded state was judged against this request by prepare() and met it
  // (B6); a solved state is converged when its stop reason says so.
  converged = checkpoint.work == "loaded" || saved.stop_reason == "converged";
  state_id = saved.generation;
  result["state_id"] = state_id;
  result["identity_kind"] = "checkpoint_generation";
  result["checkpoint_generation"] = saved.generation;
  result["work"] = checkpoint.work;
  result["iterations_this_run"] = checkpoint.iterations;
  result["checkpoint_warnings"] = checkpoint.retention_warnings();
  result["acceleration_history"] = continuing ? "reset" : "fresh";
  if (continuing) {
    // B6: what --continue did and why.
    result["continuation"] = checkpoint.continuation_record(cp);
    if (world.rank() == 0)
      madness::print("stage", stage.name, "continuation:",
                     result["continuation"]["summary"].get<std::string>());
  }
  result["checkpoint_metrics"] = saved.metrics;
  result["convergence_report"]
      = ground_convergence_report(saved.metrics, state_id, converged);
  const auto saved_measurements
      = saved.metrics.value("convergence_report", json::object())
            .value("saved_state_validation", json::object())
            .value("measurements", json::object());
  // The current request, judged on the saved state (B6, B8), also when the
  // generation was reused from an earlier request.
  result["convergence_report"]["requested_accuracy"]
      = requested_criteria_report(saved_measurements, settings);
  result["convergence_report"]["convergence_outcome"]
      = madness::convergence::outcome_to_json(
          settings, saved_outcomes(saved_measurements));
  if (world.rank() == 0)
    std::cout << madness::convergence::render_outcome_table(
        settings, saved_outcomes(saved_measurements), saved.stop_reason)
              << '\n'
              << std::flush;
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
  for (const auto& property : saved.properties) {
    const auto name = property.at("name").get<std::string>();
    if (std::find(stage.properties.begin(), stage.properties.end(), name)
        != stage.properties.end()) {
      auto value = property;
      value.erase("name");
      value["source_state_id"] = state_id;
      result["properties"][name] = value;
    }
  }
  outcome.state_id = state_id;
  outcome.converged = converged;
  outcome.result = std::move(result);
  return outcome;
}

// A named molresponse stage: the response engine against its declared
// reference, published through ResponseCheckpoint.
inline StageOutcome run_response_stage(const StageContext& ctx) {
  // The names below are those of the former loop body in run_if_named, so
  // the code reads as it did there (A8: one run function per application).
  auto& world = ctx.world;
  auto& stage = ctx.stage;
  const auto& stage_directory = ctx.stage_directory;
  const auto& directory = ctx.attempt_directory;
  const auto& run_id = ctx.run_id;
  const bool continuing = ctx.continuing;
  auto& cp = stage.params.get<CalculationParameters>();
  json result = stage_record(ctx);
  std::string state_id = result.at("state_id").get<std::string>();
  bool converged = false;
  StageOutcome outcome;
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
  result["reference_representation"]
      = {{"L", context.at("L")}, {"k", context.at("k")}};
  // B9: the stage's convergence settings, shown before the solve and
  // recorded; the solver stops on the same targets.
  const auto settings
      = response_convergence_settings(stage, cp.protocol().back());
  if (world.rank() == 0)
    std::cout << '\n'
              << madness::convergence::render_settings_table(settings) << '\n'
              << std::flush;
  result["convergence_settings"] = madness::convergence::to_json(settings);
  ResponseCheckpoint checkpoint(
      world, stage, stage_directory, directory, run_id,
      {{"stage", stage.reference}, {"generation", reference.saved.generation}});
  molresponse_v3_lib::Results response;
  if (checkpoint.prepare(continuing)) {
    response.metadata = checkpoint.metadata();
    response.properties = checkpoint.properties();
  } else {
    response = molresponse_v3_lib::run_response_on_reference(
        world, stage.params, molecule, directory, reference.saved_prefix,
        [&](madness::World& w, const json& event) {
          checkpoint.observe(w, event);
        },
        checkpoint.state_plans());
  }
  result["method"] = "TDHF";
  result["reference"]
      = {{"stage", stage.reference}, {"state_id", reference.saved.generation}};
  result["targets"]
      = {{"dconv", stage.params.get<ResponseParameters>().dconv()}};
  result["parameters"]
      = stage.params.get<ResponseParameters>().to_json_if_precedence("defined");
  // The engine's property/evidence JSON exists on root only.
  on_root(world, [&] {
    auto alpha = selected_alpha(stage, response.properties);
    converged = std::all_of(alpha.begin(), alpha.end(), [](const json& row) {
      return row.at("converged").get<bool>();
    });
    result["properties"]["polarizability"]
        = {{"values", alpha},
           {"source_state_id", state_id},
           {"reference_state_id", reference.saved.generation}};
    result["engine_metadata"] = response.metadata;
  });
  world.gop.broadcast(converged, 0);
  checkpoint.finish(response.metadata, converged);
  state_id = checkpoint.latest().generation;
  result["state_id"] = state_id;
  result["identity_kind"] = "checkpoint_generation";
  result["checkpoint_generation"] = state_id;
  result["work"] = checkpoint.work;
  result["iterations_this_run"] = checkpoint.iterations;
  result["work_items"] = checkpoint.work_items;
  result["convergence_report"] = response_convergence_report(
      checkpoint.metadata(), state_id, reference_evidence(reference));
  result["convergence_report"]["decided_by"]
      = madness::convergence::to_string(settings.decider);
  {
    const auto final_protocol = molresponse_v3::protocol_key(
        cp.protocol().back(), protocol_k(cp.protocol().back()));
    const auto outcomes
        = response_saved_outcomes(checkpoint.metadata(), final_protocol);
    result["convergence_report"]["convergence_outcome"]
        = madness::convergence::outcome_to_json(settings, outcomes);
    if (world.rank() == 0)
      std::cout << madness::convergence::render_outcome_table(
          settings, outcomes, converged ? "converged" : "unconverged")
                << '\n'
                << std::flush;
  }
  result["reference_changed"] = checkpoint.reference_changed;
  result["acceleration_history"] = continuing ? "reset" : "fresh";
  if (!checkpoint.state_plans().empty()) {
    // B6: what --continue did per state and why.
    result["continuation"] = checkpoint.continuation();
    if (world.rank() == 0)
      for (const auto& state : result["continuation"]["states"])
        madness::print("stage", stage.name,
                       "continuation:", state.at("state").get<std::string>(),
                       "-", state.at("summary").get<std::string>());
  }
  result["properties"]["polarizability"]["source_state_id"] = state_id;
  outcome.state_id = state_id;
  outcome.converged = converged;
  outcome.result = std::move(result);
  return outcome;
}

/// Returns false for the existing unnamed input surface. Missing default input
/// is left to that path, which also accepts geometry supplied on the command
/// line.
inline bool run_if_named(madness::World &world,
                         const madness::commandlineparser &parser) {
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
  for (const auto &option : parser.keyval) {
    if (option.first != "input" && option.first != "prefix" &&
        option.first != "continue" &&
        option.first != "user_defined_input_file" &&
        option.first != "user_defined_prefix")
      throw std::runtime_error("named stages: put settings inside their stage; "
                               "unsupported option --" +
                               option.first);
  }
  auto stages =
      parse(world, text); // validate the complete deck before creating outputs

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
          "fresh named-stage run requires a new output directory: " +
          root.string());
    root_string = root.string();
    run_id = std::to_string(
                 std::chrono::system_clock::now().time_since_epoch().count()) +
             "-" + std::to_string(getpid());
    for (const auto &stage : stages)
      if (!fs::create_directory(root / stage.name) && !continuing)
        throw std::runtime_error("stage output directory collision: " +
                                 stage.name);
  });
  world.gop.broadcast_serializable(root_string, 0);
  world.gop.broadcast_serializable(run_id, 0);
  root = root_string;

  // A8: one run function per application, looked up by its name; each
  // stage's outcome is kept by its label for the stages that consume it.
  const std::map<std::string, RunStage> run_functions
      = {{"moldft", run_moldft_stage}, {"molresponse", run_response_stage}};
  std::map<std::string, StageOutcome> outcomes;
  json index = {{"schema", "madqc.named-stages.v1"},
                {"run_id", run_id},
                {"mode", continuing ? "continue" : "fresh"},
                {"mpi_ranks", world.size()},
                {"stages", json::array()}};

  // 4. Execute in input order; parse() has already checked the references.
  for (auto &stage : stages) {
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
    auto outcome = run_functions.at(stage.kind)(ctx);

    // 5. Record the outcome before stopping an unconverged sequence. Later
    // stages must not run against a reference that failed its completion check.
    auto& result = outcome.result;
    result["status"] = outcome.converged ? "completed" : "unconverged";
    write_result(world, stage_directory / "results.json", result);
    index["stages"].push_back(result);
    write_result(world, root / "results.json", index);
    if (!outcome.converged)
      throw std::runtime_error("stage " + stage.name +
                               " did not satisfy its requested convergence "
                               "criteria; see results.json");
    outcomes.emplace(stage.name, std::move(outcome));
  }
  if (world.rank() == 0)
    madness::print("Named-stage results:", (root / "results.json").string());
  return true;
}
} // namespace named_stages
