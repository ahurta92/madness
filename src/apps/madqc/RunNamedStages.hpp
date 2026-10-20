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

  // Keep completed SCF engines and their selected generation paths alive for
  // downstream response. Response states still have attempt-scoped identities.
  std::map<std::string, std::shared_ptr<madness::SCF>> references;
  std::map<std::string, std::string> state_ids;
  std::map<std::string, fs::path> reference_prefixes;
  std::map<std::string, json> reference_evidence;
  json index = {{"schema", "madqc.named-stages.v1"},
                {"run_id", run_id},
                {"mode", continuing ? "continue" : "fresh"},
                {"mpi_ranks", world.size()},
                {"stages", json::array()}};

  // 4. Execute in input order; parse() has already resolved earlier references.
  for (auto &stage : stages) {
    const fs::path stage_directory = root / stage.name;
    const fs::path directory = stage_directory / ("attempt-" + run_id);
    on_root(world, [&] {
      if (!fs::create_directory(directory))
        throw std::runtime_error("stage attempt directory collision");
    });
    auto &cp = stage.params.get<CalculationParameters>();
    cp.set_user_defined_value("prefix", (directory / "state").string());
    std::string state_id = run_id + ":" + stage.name;
    json result = {{"stage", stage.name},
                   {"application", stage.kind},
                   {"state_id", state_id},
                   {"identity_kind", "fresh_run_state"},
                   {"checkpoint_generation", nullptr},
                   {"requested_properties", stage.properties},
                   {"protocol", cp.protocol()},
                   {"output_directory", directory.string()},
                   {"properties", json::object()}};
    bool converged = false;
    if (stage.kind == "moldft") {
      // Run the normal SCF engine with explicit paths; never change cwd.
      moldft_lib library((directory / "resolved.in").string());
      cp.set_user_defined_value("save",
                                false); // immutable adapter owns persistence
      auto scf = library.calc(world, stage.params);
      MoldftCheckpoint checkpoint(world, stage, stage_directory, run_id);
      checkpoint.prepare(*scf, continuing);
      scf->checkpoint_observer =
          [&](madness::World &, const madness::SCF &state,
              const madness::SCF::IterationCheckpoint &step) {
            checkpoint.observe(state, step);
          };
      auto values = library.run(world, stage.params, NextAction::Redo);
      scf->checkpoint_observer = {};
      const auto &evidence = std::get<2>(values);
      scf->work_dir = directory;
      references.emplace(stage.name, scf);

      // Use the engine's existing convergence evidence. This adds no new
      // residual definitions or stopping rules.
      converged = evidence.converged_for_thresh <= cp.protocol().back() &&
                  evidence.converged_for_dconv <=
                      std::max(cp.dconv(), cp.protocol().back());
      checkpoint.finish(*scf, converged);
      const auto &saved = checkpoint.latest();
      state_id = saved.generation;
      state_ids.emplace(stage.name, state_id);
      reference_prefixes.emplace(stage.name, checkpoint.prefix());
      result["state_id"] = state_id;
      result["identity_kind"] = "checkpoint_generation";
      result["checkpoint_generation"] = saved.generation;
      result["work"] = checkpoint.work;
      result["iterations_this_run"] = checkpoint.iterations;
      result["acceleration_history"] = continuing ? "reset" : "fresh";
      result["checkpoint_metrics"] = saved.metrics;
      result["convergence_report"] = ground_convergence_report(saved.metrics, state_id, converged);
      reference_evidence[stage.name] = {{"stage", stage.name}, {"state_id", state_id},
          {"engine_converged", converged},
          {"energy_change", result["convergence_report"]["measurements"]["energy_change"]}};
      if (result["convergence_report"].contains("saved_state_validation"))
        reference_evidence[stage.name]["saved_state_validation"] =
            result["convergence_report"]["saved_state_validation"];
      result["method"] = "HF";
      result["convergence"] = evidence.to_json();
      result["targets"] = {{"dconv", cp.dconv()}, {"econv", cp.econv()}};
      result["parameters"] = cp.to_json_if_precedence("defined");
      result["geometry"] = scf->molecule.to_json();
      result["representation"] = {
          {"L", scf->param.L()},
          {"k", scf->amo.at(0).k()},
          {"threshold", saved.numerical.at("threshold")}};
      // These properties were evaluated on the exact archived orbital state.
      for (const auto &property : saved.properties) {
        const auto name = property.at("name").get<std::string>();
        if (std::find(stage.properties.begin(), stage.properties.end(), name) !=
            stage.properties.end()) {
          auto value = property;
          value.erase("name");
          value["source_state_id"] = state_id;
          result["properties"][name] = value;
        }
      }
    } else {
      // parse() guarantees that this branch is a molresponse stage.
      // The adapter reads the source SCF archive and runs the response engine.
      // SCF may orient/center the input geometry and rederive the cell.
      // Bind to the completed reference, not the pre-SCF input estimate.
      const auto &source = references.at(stage.reference);
      cp.set_user_defined_value("l", source->param.L());
      result["reference_representation"] = {{"L", source->param.L()},
                                            {"k", source->amo.at(0).k()}};
      ResponseCheckpoint checkpoint(world, stage, stage_directory, directory, run_id,
          {{"stage", stage.reference}, {"generation", state_ids.at(stage.reference)}});
      molresponse_v3_lib::Results response;
      if (checkpoint.prepare(continuing)) {
        response.metadata = checkpoint.metadata();
        response.properties = checkpoint.properties();
      } else {
        response = molresponse_v3_lib::run_response(
            world, stage.params, source, directory, reference_prefixes.at(stage.reference),
            [&](madness::World& w, const json& event) { checkpoint.observe(w, event); });
      }
      result["method"] = "TDHF";
      result["reference"] = {{"stage", stage.reference},
                             {"state_id", state_ids.at(stage.reference)}};
      result["targets"] = {
          {"dconv", stage.params.get<ResponseParameters>().dconv()}};
      result["parameters"] =
          stage.params.get<ResponseParameters>().to_json_if_precedence(
              "defined");
      // The engine's property/evidence JSON exists on root only.
      on_root(world, [&] {
        auto alpha = selected_alpha(stage, response.properties);
        converged =
            std::all_of(alpha.begin(), alpha.end(), [](const json &row) {
              return row.at("converged").get<bool>();
            });
        result["properties"]["polarizability"] = {
            {"values", alpha},
            {"source_state_id", state_id},
            {"reference_state_id", state_ids.at(stage.reference)}};
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
      result["convergence_report"] = response_convergence_report(checkpoint.metadata(), state_id,
          reference_evidence.at(stage.reference));
      result["reference_changed"] = checkpoint.reference_changed;
      result["acceleration_history"] = continuing ? "reset" : "fresh";
      result["properties"]["polarizability"]["source_state_id"] = state_id;
    }

    // 5. Record the outcome before stopping an unconverged sequence. Later
    // stages must not run against a reference that failed its completion check.
    result["status"] = converged ? "completed" : "unconverged";
    write_result(world, stage_directory / "results.json", result);
    index["stages"].push_back(result);
    write_result(world, root / "results.json", index);
    if (!converged)
      throw std::runtime_error("stage " + stage.name +
                               " did not satisfy the existing engine "
                               "convergence test; see results.json");
  }
  if (world.rank() == 0)
    madness::print("Named-stage results:", (root / "results.json").string());
  return true;
}
} // namespace named_stages
