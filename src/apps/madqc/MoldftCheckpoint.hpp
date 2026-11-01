#pragma once

#include <madness/chem/DensityChanges.h>
#include <madness/chem/FunctionMetricSummary.h>
#include <madness/chem/ExplicitCriteria.h>

#include "NamedStages.hpp"
#include <filesystem>
#include <optional>
#include <madness/chem/StageCheckpoint.h>
#include <madness/chem/FunctionResidualReport.h>

#include <madness/chem/CheckpointRetention.h>
#include <madness/chem/LadderPlanner.h>
#include <map>
#include <set>

namespace named_stages {
namespace fs = std::filesystem;
using json = nlohmann::json;

// D8: econv is the target of the energy-change metric, but only when the user
// wrote it. The engine's default econv is not a request (B3).
[[nodiscard]] inline std::optional<double>
requested_econv(const CalculationParameters& parameters) {
  if (!parameters.is_user_defined("econv"))
    return std::nullopt;
  return parameters.econv();
}

namespace convergence = madness::convergence;

// B9: one stage's convergence settings, built once from its parameters. The
// engine-check rows use the engine's own limit function, and the written
// targets come from SCF::requested_targets, the same source the engine
// stops on, so the table shows exactly what is compared.
[[nodiscard]] inline convergence::ConvergenceSettings
moldft_convergence_settings(const std::string& stage,
                            const CalculationParameters& cp,
                            std::size_t natom) {
  using convergence::Group;
  using convergence::Relation;
  using convergence::SetBy;
  const auto text
      = [](double value) { return convergence::detail::scientific(value); };
  convergence::ConvergenceSettings settings;
  settings.stage = stage;
  settings.application = "moldft";
  settings.decider = cp.get<bool>(convergence::use_default_criteria_key)
                         ? convergence::Decider::engine_check
                         : convergence::Decider::requested_criteria;
  settings.decider_set_by_user
      = cp.is_user_defined(convergence::use_default_criteria_key);
  const double thresh = cp.protocol().back();
  const auto limits = madness::SCF::engine_check_limits(
      cp.dconv(), thresh, natom, cp.get<bool>("conv_only_dens"));
  const auto dconv_by
      = cp.is_user_defined("dconv") ? SetBy::user : SetBy::engine_default;
  const std::string effective
      = "max(dconv " + text(cp.dconv()) + ", thresh " + text(thresh) + ")";
  settings.metrics.push_back(
      {Group::engine_check, "spin density change (max)", "dconv",
       "density_change", "", "a0^(-3/2)", cp.dconv(), dconv_by,
       Relation::less_than, limits.density,
       effective + " x max(5, natom " + std::to_string(natom) + ")", "", ""});
  if (limits.bsh_checked)
    settings.metrics.push_back({Group::engine_check, "BSH residual (max)",
                                "dconv", "bsh_residual", "", "dimensionless",
                                cp.dconv(), dconv_by, Relation::less_than,
                                limits.bsh, "5 x " + effective, "", ""});
  else
    settings.metrics.push_back({Group::engine_check, "BSH residual (max)",
                                "conv_only_dens", "bsh_residual", "",
                                "dimensionless", std::nullopt, SetBy::user,
                                Relation::less_than, std::nullopt, "",
                                "not checked: conv_only_dens", ""});

  // Written targets (B8): one row each, reported or enforced by the switch.
  const std::vector<convergence::RequestedMetric> requested = {
      {"charge density change", "dconv", "charge_density_change", "",
       "a0^(-3/2)",
       "compares spin-density change against " + text(limits.density), ""},
      {"BSH residual (max)", "accuracy.bsh_residual", "bsh_residual", "",
       "dimensionless",
       limits.bsh_checked ? "compares BSH residual against " + text(limits.bsh)
                          : std::string("does not check the BSH residual"),
       ""},
      {"energy change", "econv", "energy_change", "", "hartree",
       "has no energy term",
       "engine default econv " + text(cp.econv()) + " is not a target"},
      {"kinetic trace change", "accuracy.kinetic_trace_change",
       "absolute_kinetic_trace_change", "absolute_kinetic_trace_change",
       "hartree", "has no kinetic term", ""}};
  const auto targets = madness::SCF::requested_targets(cp);
  for (const auto& row : requested)
    settings.metrics.push_back(convergence::requested_metric(row, targets));

  const auto per_rung = cp.get<std::vector<int>>("maxiter_per_rung");
  settings.iteration_caps
      = {per_rung.empty() ? std::vector<int>{cp.maxiter()} : per_rung,
         cp.is_user_defined("maxiter") || !per_rung.empty()
             ? SetBy::user
             : SetBy::engine_default};
  settings.warnings = convergence::unenforced_target_warnings(settings);
  for (auto& warning : convergence::incomplete_request_warnings(
           settings, {"charge_density_change", "bsh_residual"}))
    settings.warnings.push_back(std::move(warning));
  return settings;
}

// B8: every written target judged on the saved state (`evaluation`, the
// evaluate_explicit_criteria result the shared verdict formed), and whether
// this mode enforces it. In engine mode an unmet target is reported, not
// enforced, and does not make --continue resume: resuming would stop on the
// engine check again at once (Adrian, 2026-10-21).
[[nodiscard]] inline json
requested_criteria_report(json evaluation,
                          const convergence::ConvergenceSettings& settings) {
  auto report = std::move(evaluation);
  const bool enforced
      = settings.decider == convergence::Decider::requested_criteria;
  report["enforced"] = enforced;
  for (auto& criterion : report["criteria"].items()) {
    const auto status = criterion.value().at("status").get<std::string>();
    if (!enforced && status != "satisfied")
      criterion.value()["outcome"]
          = status + ", not enforceable in engine mode";
  }
  return report;
}

// The saved-state measurements as end-of-stage outcomes (B9).
[[nodiscard]] inline std::vector<convergence::MetricOutcome>
saved_outcomes(const json& saved_measurements) {
  std::vector<convergence::MetricOutcome> outcomes;
  for (const auto& item : saved_measurements.items()) {
    const auto& record = item.value();
    const auto availability = record.value("availability", "unavailable");
    outcomes.push_back(
        {item.key(),
         availability == "measured" && record.at("value").is_number()
             ? std::optional<double>(record.at("value").get<double>())
             : std::nullopt,
         availability,
         std::nullopt,
         record.value("units", ""),
         record.value("phase", json()).is_string()
             ? record.at("phase").get<std::string>()
             : std::string{},
         {}});
  }
  return outcomes;
}

// B9, A3: the stage's convergence settings with the observer's cadence and
// retention, each shown with who set it. The moldft AppSpec's settings.
[[nodiscard]] inline convergence::ConvergenceSettings
moldft_stage_settings(const Stage& stage, int natom) {
  auto settings = moldft_convergence_settings(
      stage.name, stage.params.get<CalculationParameters>(), natom);
  const auto set_by = [&](const char* key) {
    return stage.written.contains(key) ? convergence::SetBy::user
                                       : convergence::SetBy::engine_default;
  };
  settings.checkpoints = {
      stage.checkpoint_every_iterations, set_by("checkpoint.every_iterations"),
      stage.checkpoint_keep_all, set_by("checkpoint.keep_all"), ""};
  const auto& cp = stage.params.get<CalculationParameters>();
  settings.scope = convergence::protocol_scope(
      cp.protocol(), cp.k() == -1 ? protocol_k(cp.protocol().back()) : cp.k());
  return settings;
}

// A2: the files a moldft generation holds: one orbital archive (with its
// shards). The moldft AppSpec's members.
[[nodiscard]] inline std::vector<std::string> moldft_members(const json&) {
  return {"state.restartdata"};
}

// The requested properties among those evaluated on the archived orbitals
// (`engine_output`: {"properties": the generation's property list}), keyed by
// name. The moldft AppSpec's properties.
[[nodiscard]] inline json moldft_properties(const Stage& stage, int,
                                            const json& engine_output) {
  json selected = json::object();
  for (const auto& property : engine_output.at("properties")) {
    const auto name = property.at("name").get<std::string>();
    if (std::find(stage.properties.begin(), stage.properties.end(), name)
        != stage.properties.end()) {
      auto value = property;
      value.erase("name");
      selected[name] = value;
    }
  }
  return selected;
}

namespace continuation = madness::scf_continuation;

// B6: the requested protocol ladder, each rung identified by value.
[[nodiscard]] inline std::vector<continuation::Rung>
requested_ladder(const CalculationParameters& parameters) {
  std::vector<continuation::Rung> ladder;
  for (const double thresh : parameters.protocol())
    ladder.push_back(
        {thresh, parameters.k() == -1 ? protocol_k(thresh) : parameters.k()});
  return ladder;
}

// What a saved generation says about its rung (B6). Generations written
// before the rung record existed are read from their representation and stop
// reason, with no iterations counted.
[[nodiscard]] inline continuation::SavedRung
saved_rung(const madness::checkpoint::Manifest& manifest) {
  // D-106: position.protocol_step, or a generation's before it,
  // position.rung with its v1 keys.
  const auto& position = manifest.position;
  if (position.contains("protocol_step") || position.contains("rung")) {
    const auto& rung = position.contains("protocol_step")
                           ? position.at("protocol_step")
                           : position.at("rung");
    return {{rung.at("thresh").get<double>(), rung.at("k").get<int>()},
            rung.contains("iterations_at_protocol_step")
                ? rung.at("iterations_at_protocol_step").get<int>()
                : rung.at("iterations_at_rung").get<int>(),
            continuation::rung_stop_from_string(
                rung.value("stop_reason", rung.value("stop", "")))};
  }
  const auto phase = manifest.position.value("phase", "");
  const bool finished = manifest.stop_reason == "converged"
                        || phase == "after_final_rotation"
                        || phase == "final_properties_on_saved_orbitals";
  const auto stop = manifest.stop_reason == "converged"
                        ? continuation::RungStop::converged
                    : finished ? continuation::RungStop::iteration_limit
                               : continuation::RungStop::in_progress;
  return {{manifest.numerical.at("threshold").get<double>(),
           manifest.numerical.at("k").get<int>()},
          0,
          stop};
}

namespace planning = madness::planning;

// B6: what the planner needs from a published moldft generation: its problem
// (the physical block), its rung and how it ended, and its raw measurements:
// the engine check's deciding sample and the values on the saved orbitals.
// Generations without a convergence report have no values to re-judge, so
// they resume (B6: never a saved verdict).
[[nodiscard]] inline planning::SavedState
saved_state(const madness::checkpoint::Manifest& manifest) {
  const auto rung = saved_rung(manifest);
  const auto report
      = manifest.metrics.value("convergence_report", json::object());
  planning::SavedState state;
  state.identity = manifest.physical.dump();
  state.rung = rung.rung;
  state.iterations_at_rung = rung.iterations_at_rung;
  state.stop = rung.stop;
  state.engine_measurements = report.value("measurements", json::object());
  state.saved_measurements
      = report.value("saved_state_validation", json::object())
            .value("measurements", json::object());
  // D-106: a generation written before the rename keeps the kinetic trace
  // change under its earlier name; it is read under the current one.
  for (auto* measurements :
       {&state.engine_measurements, &state.saved_measurements})
    if (measurements->contains("delta_kinetic_energy")) {
      (*measurements)["absolute_kinetic_trace_change"]
          = measurements->at("delta_kinetic_energy");
      measurements->erase("delta_kinetic_energy");
    }
  return state;
}

// S4a (D-056): the generation's per-state record, in the one form every
// ladder application writes (numerical.states): the ground state, read from
// the manifest's rung record and convergence report as saved_state reads
// them. `decided_by`: engine_check, or planner for a generation published
// again with the planner's verdict (D-016, as response's loaded states).
[[nodiscard]] inline json
moldft_state_record(const madness::checkpoint::Manifest& manifest,
                    const std::string& decided_by) {
  const auto state = saved_state(manifest);
  return {
      {"schema", state_schema},
      {"key", "ground_state"},
      {"kind", "ground_state"},
      {"identity", state.identity},
      {"protocol_step", {{"thresh", state.rung.thresh}, {"k", state.rung.k}}},
      {"iterations_at_protocol_step", state.iterations_at_rung},
      {"stop_reason", convergence::to_string(state.stop)},
      {"diverged", false}, // moldft has no divergence exit
      {"decided_by", decided_by},
      {"engine_measurements", state.engine_measurements},
      {"saved_measurements", state.saved_measurements}};
}

// B6, the moldft AppSpec's state_records (S4b): the generation's records,
// or, for one written before S4a, the record its rung record and report give
// (decided by the planner when it was published again with its verdict).
[[nodiscard]] inline json
moldft_state_records(const Stage&, int,
                     const madness::checkpoint::Manifest& manifest) {
  if (manifest.numerical.contains("states"))
    return upgrade_state_records(manifest.numerical.at("states"));
  return json::array({moldft_state_record(
      manifest,
      manifest.position.contains("reused_from") ? "planner" : "engine_check")});
}

// B6: today's request for the stage's one state: the problem (`identity`,
// the physical block of the live engine), the ladder, the settings (mode and
// written targets) and the engine check's limits at the final rung, from the
// engine's own limit function, all from the stage's parameters.
[[nodiscard]] inline planning::LadderRequest
moldft_ladder_request(const Stage& stage, int natom,
                      const PlanContext& context) {
  const auto& cp = stage.params.get<CalculationParameters>();
  planning::LadderRequest request;
  request.identity = context.identity;
  request.ladder = requested_ladder(cp);
  request.settings = moldft_stage_settings(stage, natom);
  const auto limits = madness::SCF::engine_check_limits(
      cp.dconv(), cp.protocol().back(), natom, cp.get<bool>("conv_only_dens"));
  request.engine_limits.push_back(
      {"density_change", limits.density, 0.0, false});
  if (limits.bsh_checked)
    request.engine_limits.push_back({"bsh_residual", limits.bsh, 0.0, false});
  request.forced_iterations = cp.get<int>("diagnostic.iterations_per_rung") > 0;
  return request;
}

// B6, the moldft AppSpec's plan (S4b): the ground state's plan from its
// record and today's request (D7: one planner), and what --continue did and
// why, the results.json record.
[[nodiscard]] inline StagePlans plan_moldft_state(const Stage& stage, int natom,
                                                  const PlanContext& context,
                                                  const json& records) {
  const auto request = moldft_ladder_request(stage, natom, context);
  const auto saved = saved_state_from_record(records.at(0));
  PlannedState state;
  state.identity = records.at(0).at("key").get<std::string>();
  state.plan = planning::plan_ladder(saved, std::nullopt, request);
  state.start = request.ladder.at(state.plan.start_rung);
  state.own_rung = saved.rung;
  StagePlans plans;
  plans.continuation = continuation_record(state, request.ladder);
  plans.continuation["planner"] = "plan_ladder";
  plans.states.push_back(std::move(state));
  return plans;
}

// Application policy lives here; StageCheckpoint only checks physical storage.
// DESIGN A1, the moldft AppSpec's consumer_context (S3c1): everything a
// consuming stage needs about this ground state, besides the orbital archive,
// so it never reads a live engine object; from the engine's description
// (MoldftCheckpoint::engine_description). The molecule is stored exactly as
// the engine used it (bohr, already oriented), so rebuilding it gives the
// same coordinates bit for bit.
[[nodiscard]] inline json moldft_consumer_context(const json& engine) {
  json molecule = {{"symbols", json::array()}, {"geometry", json::array()}};
  molecule["parameters"] = engine.at("molecule_parameters");
  molecule["parameters"]["units"] = "atomic";
  molecule["parameters"]["no_orient"] = true;
  // File provenance changes each attempt and is not molecular identity.
  molecule["parameters"].erase("source_name");
  molecule["parameters"].erase("source_type");
  for (const auto& atom : engine.at("atoms")) {
    molecule["symbols"].push_back(
        madness::get_atomic_data(atom.at("atomic_number").get<unsigned int>())
            .symbol);
    molecule["geometry"].push_back({atom.at("x"), atom.at("y"), atom.at("z")});
  }
  return {{"schema", "madqc.moldft.consumer_context.v1"},
          {"molecule", molecule},
          {"charge", engine.at("charge")},
          {"multiplicity",
           engine.at("nalpha").get<int>() - engine.at("nbeta").get<int>() + 1},
          {"L", engine.at("L")},
          {"functional", engine.at("functional")},
          {"spin_restricted", engine.at("spin_restricted")},
          {"thresh", engine.at("thresh")},
          {"k", engine.at("k")}};
}

// The unweighted kinetic trace of a function residual record's blocks.
[[nodiscard]] inline json moldft_kinetic_trace(const json& blocks) {
  auto result = madness::function_metric_trace(blocks);
  result["blocks"] = blocks;
  result["quantity"] = "unweighted_kinetic_trace";
  result["units"] = "hartree";
  result["collection_definition"]
      = "stored alpha and optional unrestricted beta orbitals; restricted "
        "alpha counted once";
  return result;
}

// D-045, the moldft AppSpec's observe (S3c1): one engine sample
// (MoldftCheckpoint::sample) as the generation's metrics record: the
// sample's values, each measurement against its limit, and the convergence
// report with the iteration history. `carried` is what the observer keeps
// across the samples of one run: the iteration history and the previous
// sample's kinetic blocks.
[[nodiscard]] inline json moldft_observe(const Stage& stage, int natom,
                                         const json& sample, json& carried) {
  const bool has_previous = sample.at("has_previous").get<bool>();
  const int iteration = sample.at("iteration").get<int>();
  const double threshold = sample.at("threshold").get<double>();
  const auto value
      = [&](const char* name) { return sample.at(name).get<double>(); };
  json metrics = {{"energy",
                   {{"value", value("input_energy")},
                    {"units", "hartree"},
                    {"phase", "iteration_input_before_rotation_and_update"}}}};
  if (has_previous) {
    metrics["density_change"]
        = {{"alpha", value("density_alpha")},
           {"beta", value("density_beta")},
           {"phase", "successive_iteration_input_densities"}};
    metrics["bsh_residual"] = {{"value", value("previous_bsh_residual")},
                               {"phase", "previous_iteration_update"}};
  }
  // Measurements describe iteration evaluations, not the archived output
  // after update/rotation. Keep that distinction even in a completed report.
  metrics["kinetic_energy"]
      = {{"value", value("input_kinetic")},
         {"units", "hartree"},
         {"phase", "iteration_kinetic_evaluation_before_update"}};
  metrics["occupation_weighted_physical_kinetic_energy"]
      = metrics.at("kinetic_energy");
  const auto current_blocks = sample.at("input_kinetic_diagonal").at("blocks");
  metrics["kinetic_trace"] = moldft_kinetic_trace(current_blocks);
  metrics["kinetic_trace"]["phase"]
      = "iteration_input_before_rotation_and_update";
  metrics["kinetic_trace"]["sample"] = {{"run_id", sample.at("run_id")},
                                        {"stage", stage.name},
                                        {"threshold", threshold},
                                        {"k", sample.at("k")},
                                        {"iteration", iteration}};
  const double effective = std::max(value("dconv"), threshold);
  auto measured = [&](double value, const std::string& units,
                      const std::string& phase, bool active, json target,
                      json limit, int sample_iteration) {
    const bool available = has_previous;
    const bool finite = std::isfinite(value);
    return json{{"value", available && finite ? json(value) : json(nullptr)},
                {"availability", !available ? "unavailable"
                                 : finite   ? "measured"
                                            : "invalid"},
                {"reason", available ? "" : "no_comparable_previous_sample"},
                {"units", units},
                {"phase", phase},
                {"active_gate", active},
                {"target", target},
                {"effective_limit", limit},
                {"satisfied", available && finite && !limit.is_null()
                                  ? json(value < limit.get<double>())
                                  : json(nullptr)},
                {"valid_for_saved_state", false},
                {"sample",
                 {{"run_id", sample.at("run_id")},
                  {"stage", stage.name},
                  {"threshold", threshold},
                  {"k", sample.at("k")},
                  {"iteration", sample_iteration}}}};
  };
  json evidence;
  evidence["density_change"]
      = measured(std::max(value("density_alpha"), value("density_beta")),
                 "a0^(-3/2)", "successive_iteration_input_densities", true,
                 value("dconv"), effective * std::max(5, natom), iteration);
  evidence["density_change"]["reduction"]
      = "maximum_over_spin_density_L2_norms";
  evidence["density_change"]["components"]
      = has_previous ? json{{"alpha", value("density_alpha")},
                            {"beta", value("density_beta")}}
                     : json(nullptr);
  // These candidate measures have no default or inferred threshold. Keep
  // the existing spin-resolved gate beside them for an explicit comparison.
  evidence["charge_density_change"]
      = measured(value("density_charge"), "a0^(-3/2)",
                 "successive_iteration_input_densities", false, nullptr,
                 nullptr, iteration);
  evidence["spin_density_change"]
      = measured(value("density_spin"), "a0^(-3/2)",
                 "successive_iteration_input_densities", false, nullptr,
                 nullptr, iteration);
  evidence["charge_density_change"]["density_definition"]
      = "rho_alpha + rho_beta";
  evidence["spin_density_change"]["density_definition"]
      = "rho_alpha - rho_beta";
  for (const auto* name : {"charge_density_change", "spin_density_change"})
    evidence[name]["reduction"] = "L2_norm_of_field_difference";
  evidence["bsh_residual"]
      = measured(value("previous_bsh_residual"), "dimensionless",
                 "raw_fixed_point_difference_before_previous_update",
                 !sample.at("conv_only_dens").get<bool>(), nullptr,
                 5.0 * effective, iteration - 1);
  evidence["bsh_residual"]["reduction"]
      = "maximum_over_orbital_and_spin_L2_norms";
  if (has_previous && !sample.at("previous_function_residuals").is_null())
    evidence["bsh_residual"]["per_function"]
        = sample.at("previous_function_residuals");
  // D8: reported against a user-written econv only, never enforced here; the
  // engine check alone still stops the run.
  const json econv_target = sample.at("requested_econv");
  const bool econv = econv_target.is_number();
  evidence["energy_change"]
      = measured(value("energy_change"), "hartree",
                 "successive_iteration_energy_evaluations", false, econv_target,
                 econv_target, iteration);
  auto& energy_change = evidence["energy_change"];
  energy_change["set_by"] = econv ? "user" : "engine_default_not_a_target";
  if (!econv)
    energy_change["engine_default"] = value("econv");
  if (econv) {
    // Requested criteria compare with <= (evaluate_explicit_criteria).
    energy_change["relation"] = "less_than_or_equal";
    energy_change["satisfied"]
        = energy_change.at("value").is_number()
              ? json(energy_change.at("value").get<double>()
                     <= econv_target.get<double>())
              : json(nullptr);
  }
  evidence["absolute_kinetic_trace_change"]
      = measured(value("kinetic_change"), "hartree",
                 "successive_iteration_kinetic_evaluations", false, nullptr,
                 nullptr, iteration);
  if (!has_previous)
    carried["previous_kinetic_blocks"] = json::object();
  evidence["absolute_kinetic_trace_change"].update(
      madness::function_metric_trace_change(
          current_blocks,
          carried.value("previous_kinetic_blocks", json::object())));
  evidence["absolute_kinetic_trace_change"]["quantity"]
      = "absolute_kinetic_trace_change";
  carried["previous_kinetic_blocks"] = current_blocks;
  // B2: each row keeps the raw total energy of its sample, so every change
  // in the history can be re-derived from the record alone.
  if (!carried.contains("iteration_history"))
    carried["iteration_history"] = json::array();
  carried["iteration_history"].push_back(
      {{"measurements", evidence},
       {"energy", metrics.at("energy")},
       {"kinetic_trace", metrics.at("kinetic_trace")},
       {"occupation_weighted_physical_kinetic_energy",
        metrics.at("occupation_weighted_physical_kinetic_energy")},
       {"explicit_accuracy", sample.at("explicit_accuracy")},
       {"configured_progress", sample.at("configured_progress")},
       {"engine_check_passed", sample.at("engine_check_passed")}});
  metrics["convergence_report"]
      = {{"schema", "madqc.convergence.v1"},
         {"decided_by", sample.at("decided_by")},
         {"engine_check_passed", sample.at("engine_check_passed")},
         {"explicit_accuracy", sample.at("explicit_accuracy")},
         {"configured_progress", sample.at("configured_progress")},
         {"diagnostic_iterations_per_rung",
          sample.at("diagnostic_iterations_per_rung")},
         {"measurements", evidence},
         {"iteration_history", carried.at("iteration_history")}};
  const auto& earlier = sample.at("earlier_runs_iteration_history");
  if (!earlier.empty())
    metrics["convergence_report"]["earlier_runs_iteration_history"] = earlier;
  return metrics;
}

// A3, the moldft AppSpec's publish_now (S4c2, D-054). Second output, the
// sample's rung record (B5/B6, the manifest's position.rung): the rung
// identified by value (the sample's thresh and k) on the stage's requested
// ladder, with the iterations spent there in this run (what maxiter caps)
// and in total, earlier --continue runs included (the planner's
// iterations_done, from carried's resumed_rung, when the plan resumed this
// rung; recorded only), and how it ended. The decision
// (checkpoint::publish_now): a rung stop, this run's iterations reaching the
// rung's cap, or every K-th iteration at the rung
// (checkpoint.every_iterations).
[[nodiscard]] inline PublishDecision moldft_publish_now(const Stage& stage, int,
                                                        const json& sample,
                                                        const json& carried) {
  const auto& cp = stage.params.get<CalculationParameters>();
  const continuation::Rung rung{sample.at("threshold").get<double>(),
                                sample.at("k").get<int>()};
  const auto ladder = requested_ladder(cp);
  std::size_t index = ladder.size();
  for (std::size_t i = 0; i < ladder.size() && index == ladder.size(); ++i)
    if (continuation::same_rung(ladder[i], rung))
      index = i;
  if (index == ladder.size())
    throw std::runtime_error(
        "moldft checkpoint: sample rung is not in the protocol");
  const auto resumed = carried.value("resumed_rung", json());
  const int earlier
      = !resumed.is_null()
                && resumed.at("start_rung").get<std::size_t>() == index
            ? resumed.at("iterations_done").get<int>()
            : 0;
  const int this_run = sample.at("iteration").get<int>() + 1;
  const auto stop = sample.at("rung_stop").get<std::string>();
  json record = {{"index", index + 1},
                 {"thresh", rung.thresh},
                 {"k", rung.k},
                 {"iterations_at_protocol_step", earlier + this_run},
                 {"iterations_this_run", this_run},
                 {"stop_reason", stop}};
  const bool publish = madness::checkpoint::publish_now(
      stop != continuation::to_string(continuation::RungStop::in_progress),
      earlier + this_run, this_run,
      planning::maxiter_at_rung(
          cp.get<std::vector<int>>("maxiter_per_rung"), cp.maxiter(),
          cp.get<int>("diagnostic.iterations_per_rung"), index),
      stage.checkpoint_every_iterations);
  return {publish, std::move(record)};
}

// B5: which exit ended the stage, from the verdict, whether the saved state
// was evaluated, and the final rung's record. "unmet_requested_accuracy" is
// the one extra case of criteria mode: the final rung stopped on the written
// targets, but the saved orbitals no longer meet them.
[[nodiscard]] inline std::string
moldft_stop_reason(bool converged, bool evaluated, const json& rung) {
  if (converged)
    return "converged";
  if (!evaluated)
    return "in_progress";
  const auto stop
      = rung.is_null()
            ? std::string("iteration_limit")
            : rung.value("stop_reason", rung.value("stop", "iteration_limit"));
  if (stop == "stalled")
    return "stalled";
  if (stop == "converged")
    return "unmet_requested_accuracy";
  return "iteration_limit";
}

// A3, the moldft AppSpec's manifest (S4c3, D-054): every field of a moldft
// generation but the shared ones (set_shared_fields) and the archives, from
// one of two kinds of publish sample (D-067):
//  - "engine", a generation of the live engine at a save boundary or with
//    the final properties: the engine's values (FunctionDefaults' threshold
//    and k, its parameters, which hold the values the SCF derives, and its
//    description for the consumer context), the problem (`physical`), the
//    metrics, properties and position, the verdict, whether the saved state
//    was evaluated, and the previous generation (D-066: moldft omits it
//    after a fresh plan);
//  - "reuse", a loaded generation that met today's request, published again
//    with that verdict (A1, B6; generations are immutable, A3): the loaded
//    manifest's record (read back by checkpoint::from_json_record), the
//    plan's why and unmet_not_enforceable. The orbital shards are the
//    source's, byte for byte, so the inventory is the source's.
// The written targets and the settings are the stage's (B9; D8: econv only
// when written).
[[nodiscard]] inline madness::checkpoint::Manifest
moldft_manifest(const Stage& stage, int natom, const json& sample) {
  const auto& cp = stage.params.get<CalculationParameters>();
  const auto settings = moldft_stage_settings(stage, natom);
  const auto targets = [&](double threshold) {
    json written = {{"dconv", cp.dconv()},
                    {"effective_dconv", std::max(cp.dconv(), threshold)}};
    if (const auto econv = requested_econv(cp))
      written["econv"] = *econv;
    written["explicit_accuracy"] = convergence::requested_targets(settings);
    written["decided_by"] = convergence::to_string(settings.decider);
    return written;
  };
  const auto kind = sample.at("kind").get<std::string>();
  if (kind == "reuse") {
    auto m = madness::checkpoint::from_json_record(sample.at("loaded"));
    const auto source = m.generation;
    const auto stop_before = m.stop_reason;
    // The skeleton's verdict, as the other final publishes take it (D-077):
    // a loaded state's planner record converges; otherwise the loaded
    // generation's stop reason stands.
    if (sample.at("converged").get<bool>())
      m.stop_reason = "converged";
    m.position["previous_generation"] = source;
    m.position["reused_from"] = source;
    // D-106: the new generation writes the current keys; a source written
    // before them carries its rung record as position.rung.
    if (m.position.contains("rung")) {
      auto step = m.position.at("rung");
      for (const auto& [from, to] :
           {std::pair{"iterations_at_rung", "iterations_at_protocol_step"},
            std::pair{"stop", "stop_reason"}})
        if (step.contains(from)) {
          step[to] = step.at(from);
          step.erase(from);
        }
      m.position["protocol_step"] = std::move(step);
      m.position.erase("rung");
    }
    m.targets = targets(m.numerical.at("threshold").get<double>());
    m.metrics["convergence_report"]["reuse_verdict"]
        = {{"work", "loaded"},
           {"from_generation", source},
           {"stop_reason_before", stop_before},
           {"decided_by", convergence::to_string(settings.decider)},
           {"why", sample.at("why")},
           {"unmet_not_enforceable", sample.at("unmet_not_enforceable")},
           {"convergence_settings", convergence::to_json(settings)}};
    m.numerical["states"] = json::array({moldft_state_record(m, "planner")});
    return m;
  }
  if (kind != "engine")
    throw std::runtime_error("moldft checkpoint: unknown publish sample kind "
                             + kind);
  const auto& engine = sample.at("engine");
  const auto threshold = engine.at("threshold").get<double>();
  madness::checkpoint::Manifest m;
  m.method = "HF";
  m.physical = sample.at("physical");
  m.numerical = {
      {"threshold", threshold},
      {"k", engine.at("k")},
      {"representation", "mo"},
      {"protocol", cp.protocol()},
      {"parameters", engine.at("parameters")},
      {"consumer_context", moldft_consumer_context(engine.at("description"))}};
  m.working_orbitals = {{"ordering", "archive_order"}}; // state_id: S4c1
  m.position = sample.at("position");
  if (!sample.at("previous_generation").is_null())
    m.position["previous_generation"] = sample.at("previous_generation");
  m.metrics = sample.at("metrics");
  m.targets = targets(threshold);
  const auto& position = m.position;
  m.stop_reason = moldft_stop_reason(
      sample.at("converged").get<bool>(), sample.at("evaluated").get<bool>(),
      position.value("protocol_step", position.value("rung", json())));
  m.properties = sample.at("properties");
  m.numerical["states"] = json::array({moldft_state_record(m, "engine_check")});
  return m;
}

// B8 (S4d1), the moldft AppSpec's requested_states: the ground state, at the
// requested ladder's final rung.
[[nodiscard]] inline json moldft_requested_states(const Stage& stage, int) {
  const auto final_rung
      = requested_ladder(stage.params.get<CalculationParameters>()).back();
  return json::array(
      {{{"key", "ground_state"},
        {"protocol_step",
         {{"thresh", final_rung.thresh}, {"k", final_rung.k}}}}});
}

// B1, B8 (S4d1), the moldft AppSpec's saved_measurements: the record's
// saved-state measurements are keyed as moldft's written targets already
// (requested_criteria_report reads the same values).
[[nodiscard]] inline json moldft_saved_measurements(const Stage&, int,
                                                    const json& record) {
  return record.at("saved_measurements");
}

// B6, D-091: whether a loaded generation is published again with the
// planner's verdict: only when the load met today's request (the plan's
// `unmet` empty) and the generation did not already stop converged. An R2
// load (a stalled final rung, the request unmet) publishes nothing: the
// stalled generation stays latest and the shared verdict judges its record.
[[nodiscard]] inline bool
republish_on_load(const std::string& stop_reason,
                  const std::vector<std::string>& unmet) {
  return stop_reason != "converged" && unmet.empty();
}

// D-106: the kinetic trace change of a measurement record, under its name
// or, in a generation written before D-106, the earlier one.
[[nodiscard]] inline const json&
kinetic_trace_change(const json& measurements) {
  return measurements.contains("absolute_kinetic_trace_change")
             ? measurements.at("absolute_kinetic_trace_change")
             : measurements.at("delta_kinetic_energy");
}

class MoldftCheckpoint {
  madness::World &world_;
  const Stage &stage_;
  fs::path store_;
  std::string run_id_;
  json physical_;
  madness::checkpoint::Manifest latest_;
  bool has_latest_ = false;
  std::optional<madness::SCF::IterationCheckpoint> comparison_;
  // A3: the per-iteration history of earlier --continue runs, carried into
  // every manifest of this run so it survives retention. Record only: it is
  // not re-validated against this run's settings.
  json earlier_history_ = json::array();
  // The metrics of the last observed sample, published or not (A3: within a
  // rung only some samples are published); finish() starts from them.
  json last_metrics_ = nullptr;
  // The stage's convergence settings (B9), built once by the runner and kept
  // as an immutable copy: observe() is an engine callback and needs them at
  // every iteration (A7: owned, explicit, never changed).
  const convergence::ConvergenceSettings settings_;
  // The reuse decision on --continue (B6, B8).

  // Where this run started (B6) and the rung of the latest sample (B5).
  // The planner's decision on --continue (B6, D7).
  planning::LadderPlan plan_;
  bool resumed_by_plan_ = false;
  json current_rung_ = nullptr;
  // The skeleton's publish (publish_generation), from the run binding.
  GenerationPublisher publish_;
  static json kinetic_blocks(const madness::FunctionResidualBlocks& values) {
    return madness::function_residual_record(values).at("blocks");
  }
  static json trace_record(const madness::FunctionResidualBlocks& values) {
    return moldft_kinetic_trace(kinetic_blocks(values));
  }

  static void require(bool condition, const std::string &message) {
    if (!condition)
      throw std::runtime_error("moldft checkpoint: " + message);
  }

  // A2: the science checks on a generation's measurements.
  static void validate_measurements(const json& metrics, const std::string& stage, const json& numerical) {
    const auto configured = numerical.value("parameters", json::object()).value("parameters", json::object());
    const bool has_progress_controls = configured.value("progress.window", 0) > 0;
    bool has_accuracy_controls = false;
    for (const char* name :
         {"accuracy.bsh_residual", "accuracy.kinetic_trace_change"})
      has_accuracy_controls = has_accuracy_controls || configured.value(name, -1.0) >= 0.0;
    if (!metrics.contains("convergence_report")) {
      require(!has_progress_controls && !has_accuracy_controls,
              "invalid ground-state measurement: missing configured accuracy/progress report");
      return;
    }
    const auto& report = metrics.at("convergence_report");
    require(!has_progress_controls || report.contains("configured_progress"),
            "invalid ground-state measurement: missing configured progress");
    auto validate = [&](const json& measurements) {
      if (measurements.contains("bsh_residual") && measurements.at("bsh_residual").contains("per_function"))
        require(madness::valid_function_residual_record(measurements.at("bsh_residual").at("per_function")),
                "invalid ground-state measurement: per-function BSH");
      for (const auto* name : {"charge_density_change", "spin_density_change"}) {
        if (!measurements.contains(name)) continue; // older reports lack these diagnostics
        const auto& m = measurements.at(name);
        const bool measured = m.value("availability", "") == "measured";
        const bool valid_value = measured
            ? m.at("value").is_number() && std::isfinite(m.at("value").get<double>()) &&
              m.at("value").get<double>() >= 0.0
            : m.value("availability", "") == "unavailable" && m.at("value").is_null();
        require(valid_value && !m.at("active_gate").get<bool>() &&
                m.at("target").is_null() && m.at("effective_limit").is_null() &&
                m.at("satisfied").is_null(),
                std::string("invalid ground-state measurement: ") + name);
      }
      for (const auto& metric : measurements.items())
        require(metric.value().value("availability", "") != "invalid",
                "invalid ground-state measurement: " + metric.key());
    };
    const auto check_accuracy = [&](const json& entry) {
      if (!entry.contains("explicit_accuracy")) return;
      const auto& decision = entry.at("explicit_accuracy");
      json targets = json::object();
      for (const auto& criterion : decision.at("criteria").items())
        targets[criterion.key()] = criterion.value().at("target");
      const auto expected = madness::evaluate_explicit_criteria(entry.at("measurements"), targets);
      require(expected.at("status") != "invalid" && expected == decision,
              "invalid ground-state measurement: explicit accuracy decision");
    };
    check_accuracy(report);
    const auto check_trace = [&](const json& trace, const std::string& expected_phase) {
      require(trace.at("phase") == expected_phase,
              "invalid ground-state measurement: kinetic trace phase");
      require(trace.at("quantity") == "unweighted_kinetic_trace" && trace.at("units") == "hartree",
              "invalid ground-state measurement: kinetic trace quantity or units");
      const auto expected = madness::function_metric_trace(trace.at("blocks"));
      require(expected.at("availability") == "measured", "invalid ground-state measurement: kinetic trace");
      for (const auto& field : expected.items())
        require(trace.at(field.key()) == field.value(), "invalid ground-state measurement: inconsistent kinetic trace");
    };
    const auto check_change = [&](const json& change, const json& current, const json& previous) {
      if (change.value("quantity", "") != "absolute_kinetic_trace_change") return;
      const auto expected = madness::function_metric_trace_change(current, previous);
      for (const auto& field : expected.items())
        require(change.at(field.key()) == field.value(), "invalid ground-state measurement: inconsistent kinetic trace change");
    };
    if (metrics.contains("kinetic_trace")) check_trace(metrics.at("kinetic_trace"), "iteration_input_before_rotation_and_update");
    if (report.contains("iteration_history")) {
      json previous = json::object(), previous_sample = nullptr, progress_history = json::array();
      for (const auto& iteration : report.at("iteration_history")) {
        validate(iteration.at("measurements"));
        check_accuracy(iteration);
        const auto& trace = iteration.at("kinetic_trace");
        check_trace(trace, "iteration_input_before_rotation_and_update");
        const auto& sample = trace.at("sample");
        require(sample.at("stage") == stage, "invalid ground-state measurement: history stage");
        const bool same_segment = !previous_sample.is_null() &&
            sample.at("threshold") == previous_sample.at("threshold") &&
            sample.at("k") == previous_sample.at("k") && sample.at("run_id") == previous_sample.at("run_id");
        require(sample.at("iteration").get<int>() == (same_segment ? previous_sample.at("iteration").get<int>() + 1 : 0),
                "invalid ground-state measurement: nonconsecutive history");
        if (!same_segment) {
          previous = json::object();
          progress_history = json::array();
        }
        require(!has_progress_controls || iteration.contains("configured_progress"),
                "invalid ground-state measurement: missing configured progress history");
        if (iteration.contains("configured_progress")) {
          progress_history.push_back(iteration.at("explicit_accuracy").at("criteria"));
          const auto& progress = iteration.at("configured_progress");
          const auto& controls = numerical.at("parameters").at("parameters");
          const int window = controls.at("progress.window").get<int>();
          const double improvement = controls.at("progress.minimum_relative_improvement").get<double>();
          madness::validate_progress_configuration(window, improvement, improvement != -1.0,
              !iteration.at("explicit_accuracy").at("criteria").empty());
          const int diagnostic = controls.at("diagnostic.iterations_per_rung").get<int>();
          require(report.at("diagnostic_iterations_per_rung") == diagnostic,
                  "invalid ground-state measurement: diagnostic configuration");
          bool has_next = false;
          const auto& protocol = numerical.at("protocol");
          for (std::size_t rung = 0; rung + 1 < protocol.size(); ++rung)
            has_next = has_next || protocol.at(rung) == sample.at("threshold");
          auto expected = window == 0
              ? json{{"status", "disabled"}, {"action", "continue"}}
              : madness::explicit_progress(progress_history, window, improvement, has_next, false);
          expected["has_next_supplied_rung"] = has_next;
          expected["action_applied"] = expected.at("status") == "stalled" &&
              diagnostic == 0;
          require(expected.at("action") != "stop_invalid" && expected == progress,
                  "invalid ground-state measurement: configured progress");
        }
        check_change(kinetic_trace_change(iteration.at("measurements")),
                     trace.at("blocks"), previous);
        previous = trace.at("blocks");
        previous_sample = sample;
      }
    }
    if (metrics.contains("kinetic_trace")) {
      const auto& trace = metrics.at("kinetic_trace");
      if (report.contains("configured_progress"))
        require(report.at("configured_progress") == report.at("iteration_history").back().at("configured_progress"),
                "invalid ground-state measurement: latest configured progress");
      require(!report.at("iteration_history").empty() && trace == report.at("iteration_history").back().at("kinetic_trace"),
              "invalid ground-state measurement: latest trace history");
      require(trace.at("sample").at("threshold") == numerical.at("threshold") &&
              trace.at("sample").at("k") == numerical.at("k"),
              "invalid ground-state measurement: latest trace representation");
    }
    if (report.contains("saved_state_validation") && report.at("saved_state_validation").contains("kinetic_trace")) {
      const auto& saved = report.at("saved_state_validation");
      check_accuracy(saved);
      check_trace(saved.at("kinetic_trace"), "saved_orbitals_before_archive");
      json previous = json::object();
      const auto& last_sample = metrics.at("kinetic_trace").at("sample");
      require(saved.contains("comparison_boundary") &&
              (saved.at("comparison_boundary") == "after_final_rotation" || saved.at("comparison_boundary") == "after_orbital_update"),
              "invalid ground-state measurement: comparison boundary");
      const int expected_iteration = last_sample.at("iteration").get<int>() -
          (saved.at("comparison_boundary") == "after_final_rotation" ? 1 : 0);
      require(saved.at("comparison").is_null() == (expected_iteration < 0),
              "invalid ground-state measurement: missing comparison");
      if (!saved.at("comparison").is_null()) {
        const auto& comparison = saved.at("comparison");
        auto expected_sample = last_sample;
        expected_sample["iteration"] = expected_iteration;
        for (const auto& field : expected_sample.items())
          require(comparison.at(field.key()) == field.value(), "invalid ground-state measurement: comparison identity");
        check_trace(comparison.at("kinetic_trace"), "iteration_input_before_rotation_and_update");
        const auto& history = report.at("iteration_history");
        const auto offset = saved.at("comparison_boundary") == "after_final_rotation" ? 2 : 1;
        require(history.size() >= std::size_t(offset), "invalid ground-state measurement: missing comparison history");
        const auto& source = history.at(history.size() - offset).at("kinetic_trace");
        require(source.at("sample") == expected_sample && source.at("blocks") == comparison.at("kinetic_trace").at("blocks"),
                "invalid ground-state measurement: comparison history mismatch");
        previous = comparison.at("kinetic_trace").at("blocks");
      }
      check_change(kinetic_trace_change(saved.at("measurements")),
                   saved.at("kinetic_trace").at("blocks"), previous);
    }
    validate(report.at("measurements"));
    if (report.contains("saved_state_validation"))
      validate(report.at("saved_state_validation").at("measurements"));
  }

public:
  // A2: the checks on a generation --continue reads, before the plan: its
  // stage and method, its orbital inventory and representation, and its
  // measurements. The moldft AppSpec's validate_generation (S3b).
  static void validate_generation(const std::string& stage,
                                  const madness::checkpoint::Manifest& m) {
    require(m.stage == stage && m.method == "HF", "stage/method mismatch");
    require(m.archives.size() == 1 && m.archives[0].role == "orbitals"
                && m.archives[0].basename == "state.restartdata",
            "required orbital inventory is missing");
    require(m.numerical.at("representation") == "mo",
            "unsupported representation");
    validate_measurements(m.metrics, m.stage, m.numerical);
    // S4a: the state record, when the generation has one (D-052: the same
    // check on both sides; publish_generation requires it): one current
    // record, the one this generation's own rung record and report give.
    if (m.numerical.contains("states")) {
      validate_state_records(m.numerical.at("states"));
      // D-106: a record of either schema, read in the current one.
      const auto states = upgrade_state_records(m.numerical.at("states"));
      require(states.size() == 1 && states[0].at("schema") == state_schema,
              "the ground state needs one " + std::string(state_schema)
                  + " record");
      const auto decided_by = states[0].at("decided_by").get<std::string>();
      require((decided_by == "engine_check" || decided_by == "planner")
                  && states[0] == moldft_state_record(m, decided_by),
              "state record disagrees with its generation");
    }
  }

private:
  // DESIGN A1: what moldft_consumer_context needs from the engine: the
  // molecule as the engine used it (bohr, already oriented), the electrons,
  // the cell, the functional, and the final thresh and k.
  static json engine_description(const madness::SCF& scf) {
    json atoms = json::array();
    for (std::size_t i = 0; i < scf.molecule.natom(); ++i) {
      const auto& atom = scf.molecule.get_atom(i);
      atoms.push_back({{"atomic_number", atom.atomic_number},
                       {"x", atom.x},
                       {"y", atom.y},
                       {"z", atom.z}});
    }
    return {{"molecule_parameters",
             scf.molecule.parameters.to_json()["parameters"]},
            {"atoms", atoms},
            {"charge", scf.param.charge()},
            {"nalpha", scf.param.nalpha()},
            {"nbeta", scf.param.nbeta()},
            {"L", scf.param.L()},
            {"functional", scf.param.xc()},
            {"spin_restricted", scf.param.spin_restricted()},
            {"thresh", madness::FunctionDefaults<3>::get_thresh()},
            {"k", madness::FunctionDefaults<3>::get_k()}};
  }

  static json physical(const madness::SCF &scf) {
    auto molecule = scf.molecule.to_json();
    // File provenance changes each attempt and is not molecular identity.
    molecule["parameters"].erase("source_name");
    molecule["parameters"].erase("source_type");
    return {{"schema", "madqc.moldft.v1"},
            {"molecule", molecule},
            {"xc", scf.param.xc()},
            {"charge", scf.param.charge()},
            {"spin_restricted", scf.param.spin_restricted()},
            {"nalpha", scf.param.nalpha()},
            {"nbeta", scf.param.nbeta()},
            {"localize", scf.param.localize_method()},
            {"L", scf.param.L()},
            {"representation", "mo"}};
  }

public:
  std::size_t iterations = 0;
  std::string work = "solved";

  // `publish`: the skeleton's publish (publish_generation, S3c2).
  MoldftCheckpoint(madness::World& world, const Stage& stage,
                   const fs::path& directory, std::string run_id,
                   convergence::ConvergenceSettings settings,
                   GenerationPublisher publish)
      : world_(world), stage_(stage), store_(directory / "checkpoints"),
        run_id_(std::move(run_id)), settings_(std::move(settings)),
        publish_(std::move(publish)) {}

  // The engine's progress controls and the problem, before any generation is
  // read.
  void prepare(madness::SCF& scf) {
    scf.validate_progress_controls();
    physical_ = physical(scf);
  }

  // The problem, as the plan's context reads it (PlanContext::identity).
  [[nodiscard]] std::string identity() const { return physical_.dump(); }

  // The generation --continue read, validated (validate_generation), and its
  // plan (the AppSpec's plan), applied before the engine opens an archive. A
  // private copy lets the existing restart planner use its ordinary prefix
  // without ever changing an immutable generation or using an AO sidecar.
  void load(madness::SCF& scf, const madness::checkpoint::Manifest& latest,
            const StagePlans& plans) {
    latest_ = latest;
    const auto &cp = scf.param;
    const auto report
        = latest_.metrics.value("convergence_report", json::object());
    plan_ = plans.states.at(0).plan;
    if (plan_.action == planning::Action::fresh) {
      // R1: a different problem starts fresh; the old generations stay.
      work = "solved";
      has_latest_ = false;
      return;
    }
    const bool load = plan_.action == planning::Action::load;
    work = load ? "loaded" : "resumed";
    madness::checkpoint::detail::on_root(world_, [&] {
      const auto directory = store_ / latest_.generation;
      for (const auto &shard : latest_.archives[0].shards) {
        const auto suffix = shard.path.substr(std::string("state").size());
        fs::copy_file(directory / shard.path, cp.prefix() + suffix,
                      fs::copy_options::none);
      }
      return std::string{};
    });
    // Explicit iteration cannot silently fall back to an atomic guess; a
    // loaded state is read without iterating, whatever its convergence claim.
    scf.param.set_user_defined_value("restart", load ? std::string("read_only")
                                                     : std::string("iterate"));
    if (!load) {
      scf.param.set_user_defined_value("restart.start_rung",
                                       static_cast<int>(plan_.start_rung));
      resumed_by_plan_ = true;
      // A3: the resumed-from generation may be deleted by retention, so its
      // history (and what it carried from earlier runs) moves forward.
      earlier_history_
          = report.value("earlier_runs_iteration_history", json::array());
      for (const auto& entry : report.value("iteration_history", json::array()))
        earlier_history_.push_back(entry);
    }
    has_latest_ = true;
  }

  const madness::checkpoint::Manifest &latest() const { return latest_; }

  // S4c2: the rung the plan resumed, as the AppSpec's publish_now reads it
  // from `carried` (resumed_rung): its index in the requested ladder and the
  // iterations spent there in earlier runs; null unless the plan resumed.
  [[nodiscard]] json resumed_rung() const {
    if (!resumed_by_plan_)
      return nullptr;
    return {{"start_rung", plan_.start_rung},
            {"iterations_done", plan_.iterations_done}};
  }

  // D-045: the engine's sample at a save boundary as JSON, for the AppSpec's
  // observe (moldft_observe): the values measured, with the engine settings
  // and the records of earlier runs it reads.
  [[nodiscard]] json
  sample(const madness::SCF& scf,
         const madness::SCF::IterationCheckpoint& step) const {
    const auto econv = requested_econv(scf.param);
    return {{"run_id", run_id_},
            {"iteration", step.iteration},
            {"threshold", step.threshold},
            {"k", scf.amo.at(0).k()},
            {"rung_stop", continuation::to_string(step.rung_stop)},
            {"has_previous", step.has_previous},
            {"input_energy", step.input_energy},
            {"density_alpha", step.density_alpha},
            {"density_beta", step.density_beta},
            {"density_charge", step.density_charge},
            {"density_spin", step.density_spin},
            {"previous_bsh_residual", step.previous_bsh_residual},
            {"previous_function_residuals",
             step.previous_function_residuals.empty()
                 ? json(nullptr)
                 : madness::function_residual_record(
                       step.previous_function_residuals)},
            {"input_kinetic", step.input_kinetic},
            {"input_kinetic_diagonal",
             madness::function_residual_record(step.input_kinetic_diagonal)},
            {"energy_change", step.energy_change},
            {"kinetic_change", step.kinetic_change},
            {"explicit_accuracy", step.explicit_accuracy},
            {"configured_progress", step.configured_progress},
            {"engine_check_passed", step.engine_check_passed},
            {"dconv", scf.param.dconv()},
            {"conv_only_dens", scf.param.get<bool>("conv_only_dens")},
            {"econv", scf.param.econv()},
            {"requested_econv", econv ? json(*econv) : json(nullptr)},
            {"diagnostic_iterations_per_rung",
             scf.param.get<int>("diagnostic.iterations_per_rung")},
            {"decided_by", convergence::to_string(settings_.decider)},
            {"earlier_runs_iteration_history", earlier_history_}};
  }

  // A3: one save boundary. `metrics` is the AppSpec's observe of this
  // sample; `decision`, its publish_now (S4c2): whether it publishes, and
  // the sample's rung record.
  void observe(const madness::SCF& scf,
               const madness::SCF::IterationCheckpoint& step,
               const json& metrics, const PublishDecision& decision) {
    ++iterations;
    comparison_ = step;
    current_rung_ = decision.rung;
    last_metrics_ = metrics;
    if (decision.publish)
      publish(scf, false, json::array(), metrics,
              {{"iteration", step.iteration},
               {"phase", step.final_rotation ? "after_final_rotation"
                                             : "after_orbital_update"},
               {"protocol_step", current_rung_}},
              false);
  }

  // Only these final properties have actually been evaluated on the saved
  // state.
  // `engine_converged`: the engine's own verdict at the final rung, recorded
  // in the report; the stage's verdict is the skeleton's (D-077).
  void finish(madness::SCF& scf, bool engine_converged) {
    if (work == "loaded" && iterations == 0) {
      // A loaded state (B6): when it met today's request and its generation
      // stopped as anything but converged, record the new verdict on disk,
      // so a consumer reading the generation agrees with this run (A1).
      if (republish_on_load(latest_.stop_reason, plan_.unmet))
        publish_reuse_verdict();
      return;
    }
    const auto evaluation = scf.evaluate_restricted_hf(world_);
    const double energy = evaluation.energy;
    require(std::isfinite(energy), "nonfinite final energy");
    const auto& dipole = evaluation.dipole;
    std::vector<double> dip(3);
    for (int i = 0; i < 3; ++i) {
      dip[i] = dipole(i);
      require(std::isfinite(dip[i]), "nonfinite final dipole");
    }
    scf.current_energy = energy;
    json properties = {
        {{"name", "energy"}, {"value", energy}, {"units", "hartree"}},
        {{"name", "dipole"}, {"value", dip}, {"units", "e*a0"}}};
    // The last observed sample, whether or not the cadence published it (A3).
    auto metrics = json::object();
    if (!last_metrics_.is_null())
      metrics = last_metrics_;
    else if (has_latest_)
      metrics = latest_.metrics;
    // History remains intact. These additional checks evaluate the saved
    // orbitals and use a same-rung iteration input as the comparison sample.
    json evidence = metrics.at("convergence_report").at("measurements");
    const bool comparable = comparison_ && comparison_->comparison_iteration >= 0;
    double density = 0.0;
    if (comparable)
      density = (evaluation.alpha_density - comparison_->comparison_alpha).norm2();
    // This adapter is restricted HF: the evaluated beta density equals alpha.
    // The comparison functions still identify the genuine preceding input.
    const auto density_changes = comparable
        ? madness::density_change_norms(world_, evaluation.alpha_density, evaluation.alpha_density,
              comparison_->comparison_alpha, comparison_->comparison_beta)
        : madness::DensityChanges{};
    // Without a same-rung comparison sample a change has no value at all (B2),
    // so it is absent rather than zero.
    const auto change = [&](double value) {
      return comparable ? std::optional<double>(value) : std::nullopt;
    };
    const std::map<std::string, std::optional<double>> values = {
        {"bsh_residual", evaluation.bsh_residual},
        {"density_change", change(density)},
        {"charge_density_change", change(density_changes.charge)},
        {"spin_density_change", change(density_changes.spin)},
        {"energy_change",
         comparable ? std::optional<double>(
                          std::abs(energy - comparison_->comparison_energy))
                    : std::nullopt},
        {"absolute_kinetic_trace_change",
         comparable ? std::optional<double>(std::abs(
                          evaluation.kinetic - comparison_->comparison_kinetic))
                    : std::nullopt}};
    for (auto& item : evidence.items()) {
      auto& m = item.value();
      const auto measured_value = values.at(item.key());
      const bool available = measured_value.has_value();
      const double value
          = measured_value.value_or(0.0); // read only when available
      const bool finite = std::isfinite(value);
      m["value"] = available && finite ? json(value) : json(nullptr);
      m["availability"] = !available ? "unavailable" : finite ? "measured" : "invalid";
      m["reason"] = available ? "" : "no_comparable_previous_sample";
      m["valid_for_saved_state"] = available && finite;
      m["phase"] = "saved_orbitals_before_archive";
      m["active_gate"] = false; // diagnostic; the engine verdict is unchanged
      const bool inclusive
          = m.value("relation", json()) == "less_than_or_equal";
      const auto& limit = m.at("effective_limit");
      m["satisfied"] = available && finite && !limit.is_null()
                           ? json(inclusive ? value <= limit.get<double>()
                                            : value < limit.get<double>())
                           : json(nullptr);
      m.erase("sample");
      m.erase("components");
    }
    auto saved_trace = trace_record(evaluation.kinetic_diagonal);
    saved_trace["phase"] = "saved_orbitals_before_archive";
    const auto comparison_blocks = comparable ? kinetic_blocks(comparison_->comparison_kinetic_diagonal) : json::object();
    evidence["absolute_kinetic_trace_change"].update(
        madness::function_metric_trace_change(saved_trace.at("blocks"),
                                              comparison_blocks));
    evidence["absolute_kinetic_trace_change"]["valid_for_saved_state"]
        = evidence["absolute_kinetic_trace_change"].at("availability")
          == "measured";
    // Replace historical components with evidence from the actual saved orbitals.
    evidence["bsh_residual"]["per_function"] = madness::function_residual_record(evaluation.function_residuals);
    auto comparison_trace = comparable ? trace_record(comparison_->comparison_kinetic_diagonal) : json(nullptr);
    if (comparable) comparison_trace["phase"] = "iteration_input_before_rotation_and_update";
    metrics["convergence_report"]["saved_state_validation"] = {
        {"comparison_boundary", comparison_->final_rotation ? "after_final_rotation" : "after_orbital_update"},
        {"affects_stopping", false}, {"measurements", evidence},
        {"energy", energy}, {"kinetic_energy", evaluation.kinetic},
        {"kinetic_trace", saved_trace},
        {"comparison", comparable ? json{{"run_id", run_id_}, {"stage", stage_.name},
            {"iteration", comparison_->comparison_iteration}, {"threshold", comparison_->threshold},
            {"k", scf.amo.at(0).k()}, {"energy", comparison_->comparison_energy},
            {"kinetic_energy", comparison_->comparison_kinetic},
            {"kinetic_trace", comparison_trace}} : json(nullptr)}};
    const auto requested_accuracy = madness::evaluate_explicit_criteria(
        evidence, convergence::requested_targets(settings_));
    auto& report = metrics["convergence_report"];
    report["saved_state_validation"]["explicit_accuracy"] = requested_accuracy;
    report["engine_check_passed"] = engine_converged;
    report["decided_by"] = convergence::to_string(settings_.decider);
    report["convergence_settings"] = convergence::to_json(settings_);
    report["convergence_outcome"]
        = convergence::outcome_record(settings_, saved_outcomes(evidence));
    json position = {{"phase", "final_properties_on_saved_orbitals"}};
    if (!current_rung_.is_null())
      position["protocol_step"] = current_rung_;
    // The skeleton forms the verdict from this generation's state record
    // and the publish takes it (D-077).
    publish(scf, false, properties, metrics, position, true, true);
  }

private:
  void publish(const madness::SCF& scf, bool converged, json properties,
               const json& metrics, const json& position, bool evaluated,
               bool final = false) {
    namespace checkpoint = madness::checkpoint;
    GenerationWrite write;
    // The engine kind of publish sample (moldft_manifest).
    write.publish_sample
        = {{"kind", "engine"},
           {"engine",
            {{"threshold", madness::FunctionDefaults<3>::get_thresh()},
             {"k", madness::FunctionDefaults<3>::get_k()},
             {"parameters", scf.param.to_json()},
             {"description", engine_description(scf)}}},
           {"physical", physical_},
           {"metrics", metrics},
           {"properties", std::move(properties)},
           {"position", position},
           {"previous_generation",
            has_latest_ ? json(latest_.generation) : json(nullptr)},
           {"converged", converged},
           {"evaluated", evaluated}};
    write.iterations_this_run = iterations;
    write.final = final;
    write.check_before_write = [](const checkpoint::Manifest& m) {
      validate_measurements(m.metrics, m.stage, m.numerical);
    };
    write.write_archives
        = [&](const fs::path& directory, checkpoint::Manifest& m) {
            auto meta = scf.orbital_metadata();
            if (!evaluated)
              meta.current_energy
                  = 1.e10; // input energy is not an output-state property
            // The verdict the generation was published with (D-077).
            if (m.stop_reason != "converged") {
              meta.converged_for_thresh = 1.e10;
              meta.converged_for_dconv = 1.e10;
            }
            std::string error = scf.write_mos(
                world_, (directory / "state.restartdata").string(), meta);
            // Inventory is built only after every writer has closed. The header
            // records the actual writer count, including the parallel archive's
            // clamp.
            checkpoint::detail::collective_error(world_, error);
            world_.gop.fence();
            if (world_.rank() == 0) {
              try {
                const int writers = checkpoint::detail::header(
                    directory / "state.restartdata.00000", true);
                m.archives.push_back(checkpoint::inspect_archive(
                    directory, "orbitals", "state.restartdata", {}, writers));
              } catch (const std::exception& e) {
                error = e.what();
              }
            }
            return error;
          };
    // Reading back also gives every rank the exact, validated published
    // inventory.
    latest_ = publish_(std::move(write));
    has_latest_ = true;
  }

  // A1, B6: a loaded generation re-judged as meeting today's request is
  // published again as a new generation with that verdict: the same orbital
  // shards, byte for byte, the same measurements, stop reason converged, and
  // a record of the request it met and the generation it came from.
  // Generations are immutable (A3), so the verdict cannot be written into the
  // old one.
  void publish_reuse_verdict() {
    namespace checkpoint = madness::checkpoint;
    const auto source = latest_.generation;
    GenerationWrite write;
    // The reuse kind of publish sample (moldft_manifest, D-067).
    write.publish_sample
        = {{"kind", "reuse"},
           {"converged", false}, // the skeleton's verdict replaces it (D-077)
           {"loaded", checkpoint::to_json_record(latest_)},
           {"why", plan_.why},
           {"unmet_not_enforceable", plan_.unmet_not_enforceable}};
    write.final = true;
    write.write_archives
        = [&](const fs::path& directory, const checkpoint::Manifest&) {
            return checkpoint::detail::on_root(world_, [&] {
              const auto& orbitals = latest_.archives.at(0);
              for (const auto& shard : orbitals.shards)
                fs::copy_file(store_ / source / shard.path,
                              directory / shard.path, fs::copy_options::none);
              const auto copied = checkpoint::inspect_archive(
                  directory, orbitals.role, orbitals.basename, orbitals.indices,
                  orbitals.writers);
              for (std::size_t i = 0; i < orbitals.shards.size(); ++i)
                require(copied.shards.at(i).checksum
                            == orbitals.shards.at(i).checksum,
                        "copied orbital shard differs from its source");
              return std::string{};
            });
          };
    latest_ = publish_(std::move(write));
  }
};
} // namespace named_stages
