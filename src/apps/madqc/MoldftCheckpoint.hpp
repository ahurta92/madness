#pragma once

#include <madness/chem/DensityChanges.h>
#include <madness/chem/FunctionMetricSummary.h>
#include <madness/chem/ExplicitCriteria.h>

#include "NamedStages.hpp"
#include <filesystem>
#include <optional>
#include <madness/chem/StageCheckpoint.h>
#include <madness/chem/FunctionResidualReport.h>

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

// The econv request in the form evaluate_explicit_criteria() takes; empty when
// the deck has no econv.
[[nodiscard]] inline json
energy_change_target(const CalculationParameters& parameters) {
  json targets = json::object();
  if (const auto econv = requested_econv(parameters))
    targets["energy_change"] = {{"value", *econv}, {"units", "hartree"}};
  return targets;
}

// D8 with B8 engine mode: a written econv is judged on the energy change
// measured on the saved orbitals (never on the econv a generation recorded),
// and reported. The engine check decides, so an unmet econv is not enforceable
// and does not make --continue resume: resuming would stop on the engine check
// again at once (Adrian, 2026-10-21).
[[nodiscard]] inline json
energy_change_request_report(const json& saved_measurements,
                             const CalculationParameters& parameters) {
  auto report = madness::evaluate_explicit_criteria(
      saved_measurements, energy_change_target(parameters));
  report["enforced"] = false;
  const auto status = report.at("status").get<std::string>();
  if (status != "not_requested" && status != "satisfied")
    report["outcome"] = status + ", not enforceable in engine mode";
  return report;
}

// Application policy lives here; StageCheckpoint only checks physical storage.
class MoldftCheckpoint {
  madness::World &world_;
  const Stage &stage_;
  fs::path store_;
  std::string run_id_;
  std::size_t sequence_ = 0;
  json physical_;
  madness::checkpoint::Manifest latest_;
  bool has_latest_ = false;
  std::optional<madness::SCF::IterationCheckpoint> comparison_;
  json previous_kinetic_blocks_ = json::object();
  json iteration_history_ = json::array();
  static json kinetic_blocks(const madness::FunctionResidualBlocks& values) {
    return madness::function_residual_record(values).at("blocks");
  }
  static json trace_record(const madness::FunctionResidualBlocks& values) {
    auto blocks = kinetic_blocks(values);
    auto result = madness::function_metric_trace(blocks);
    result["blocks"] = blocks;
    result["quantity"] = "unweighted_kinetic_trace";
    result["units"] = "hartree";
    result["collection_definition"] = "stored alpha and optional unrestricted beta orbitals; restricted alpha counted once";
    return result;
  }

  static void require(bool condition, const std::string &message) {
    if (!condition)
      throw std::runtime_error("moldft checkpoint: " + message);
  }

  static void validate_measurements(const json& metrics, const std::string& stage, const json& numerical) {
    const auto configured = numerical.value("parameters", json::object()).value("parameters", json::object());
    const bool has_progress_controls = configured.value("progress.window", 0) > 0;
    bool has_accuracy_controls = false;
    for (const char* name : {"accuracy.density_change", "accuracy.bsh_residual", "accuracy.energy_change",
                             "accuracy.kinetic_trace_change"})
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
        check_change(iteration.at("measurements").at("delta_kinetic_energy"), trace.at("blocks"), previous);
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
      check_change(saved.at("measurements").at("delta_kinetic_energy"), saved.at("kinetic_trace").at("blocks"), previous);
    }
    validate(report.at("measurements"));
    if (report.contains("saved_state_validation"))
      validate(report.at("saved_state_validation").at("measurements"));
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

  MoldftCheckpoint(madness::World &world, const Stage &stage,
                   const fs::path &directory, std::string run_id)
      : world_(world), stage_(stage), store_(directory / "checkpoints"),
        run_id_(std::move(run_id)) {}

  // Validate the complete selected generation before the engine opens an
  // archive. A private copy lets the existing restart planner use its ordinary
  // prefix without ever changing an immutable generation or using an AO
  // sidecar.
  void prepare(madness::SCF &scf, bool continuing) {
    scf.validate_progress_controls();
    physical_ = physical(scf);
    if (!continuing)
      return;
    latest_ = madness::checkpoint::read_published(world_, store_);
    require(latest_.stage == stage_.name && latest_.method == "HF",
            "stage/method mismatch");
    require(latest_.physical == physical_,
            "geometry, model or cell changed; use a new run");
    require(latest_.archives.size() == 1 &&
                latest_.archives[0].role == "orbitals" &&
                latest_.archives[0].basename == "state.restartdata",
            "required orbital inventory is missing");
    require(latest_.numerical.at("representation") == "mo",
            "unsupported representation");
    const auto &cp = scf.param;
    const auto accuracy_targets = madness::requested_accuracy_targets(cp, "a0^(-3/2)", "dimensionless", "hartree");
    const auto saved_measurements = latest_.metrics.value("convergence_report", json::object())
        .value("saved_state_validation", json::object()).value("measurements", json::object());
    const auto requested_accuracy = madness::evaluate_explicit_criteria(saved_measurements, accuracy_targets);
    const bool satisfied
        = (accuracy_targets.empty()
           || requested_accuracy.at("status") == "satisfied")
          && cp.get<int>("diagnostic.iterations_per_rung") == 0
          && latest_.stop_reason == "converged"
          && latest_.numerical.at("threshold").get<double>()
                 <= cp.protocol().back()
          && latest_.numerical.at("k").get<int>()
                 == (cp.k() == -1 ? protocol_k(cp.protocol().back()) : cp.k())
          && latest_.targets.at("effective_dconv").get<double>()
                 <= std::max(cp.dconv(), cp.protocol().back());
    work = satisfied ? "loaded" : "resumed";
    madness::checkpoint::detail::on_root(world_, [&] {
      validate_measurements(latest_.metrics, latest_.stage, latest_.numerical);
      const auto directory = store_ / latest_.generation;
      for (const auto &shard : latest_.archives[0].shards) {
        const auto suffix = shard.path.substr(std::string("state").size());
        fs::copy_file(directory / shard.path, cp.prefix() + suffix,
                      fs::copy_options::none);
      }
      return std::string{};
    });
    // Explicit iteration cannot silently fall back to an atomic guess.
    scf.param.set_user_defined_value(
        "restart", satisfied ? std::string("auto") : std::string("iterate"));
    has_latest_ = true;
  }

  const madness::checkpoint::Manifest &latest() const { return latest_; }
  fs::path prefix() const { return store_ / latest_.generation / "state"; }

  void observe(const madness::SCF &scf,
               const madness::SCF::IterationCheckpoint &step) {
    ++iterations;
    comparison_ = step;
    json metrics = {
        {"energy",
         {{"value", step.input_energy},
          {"units", "hartree"},
          {"phase", "iteration_input_before_rotation_and_update"}}}};
    if (step.has_previous) {
      metrics["density_change"] = {
          {"alpha", step.density_alpha},
          {"beta", step.density_beta},
          {"phase", "successive_iteration_input_densities"}};
      metrics["bsh_residual"] = {{"value", step.previous_bsh_residual},
                                 {"phase", "previous_iteration_update"}};
    }
    // Measurements describe iteration evaluations, not the archived output
    // after update/rotation. Keep that distinction even in a completed report.
    metrics["kinetic_energy"] = {{"value", step.input_kinetic}, {"units", "hartree"},
        {"phase", "iteration_kinetic_evaluation_before_update"}};
    metrics["occupation_weighted_physical_kinetic_energy"] = metrics.at("kinetic_energy");
    metrics["kinetic_trace"] = trace_record(step.input_kinetic_diagonal);
    metrics["kinetic_trace"]["phase"] = "iteration_input_before_rotation_and_update";
    metrics["kinetic_trace"]["sample"] = {{"run_id", run_id_}, {"stage", stage_.name},
        {"threshold", step.threshold}, {"k", scf.amo.at(0).k()}, {"iteration", step.iteration}};
    const double effective = std::max(scf.param.dconv(), step.threshold);
    auto measured = [&](double value, const std::string& units,
                        const std::string& phase, bool active, json target,
                        json limit, int sample_iteration) {
      const bool available = step.has_previous;
      const bool finite = std::isfinite(value);
      return json{{"value", available && finite ? json(value) : json(nullptr)},
                  {"availability", !available ? "unavailable" : finite ? "measured" : "invalid"},
                  {"reason", available ? "" : "no_comparable_previous_sample"},
                  {"units", units}, {"phase", phase}, {"active_gate", active},
                  {"target", target}, {"effective_limit", limit},
                  {"satisfied", available && finite && !limit.is_null()
                      ? json(value < limit.get<double>()) : json(nullptr)},
                  {"valid_for_saved_state", false},
                  {"sample", {{"run_id", run_id_}, {"stage", stage_.name},
                              {"threshold", step.threshold}, {"k", scf.amo.at(0).k()},
                              {"iteration", sample_iteration}}}};
    };
    json evidence;
    evidence["density_change"] = measured(std::max(step.density_alpha, step.density_beta),
        "a0^(-3/2)", "successive_iteration_input_densities", true, scf.param.dconv(),
        effective * std::max(std::size_t(5), scf.molecule.natom()), step.iteration);
    evidence["density_change"]["reduction"] = "maximum_over_spin_density_L2_norms";
    evidence["density_change"]["components"] = step.has_previous
        ? json{{"alpha", step.density_alpha}, {"beta", step.density_beta}} : json(nullptr);
    // These candidate measures have no default or inferred threshold. Keep
    // the existing spin-resolved gate beside them for an explicit comparison.
    evidence["charge_density_change"] = measured(step.density_charge, "a0^(-3/2)",
        "successive_iteration_input_densities", false, nullptr, nullptr, step.iteration);
    evidence["spin_density_change"] = measured(step.density_spin, "a0^(-3/2)",
        "successive_iteration_input_densities", false, nullptr, nullptr, step.iteration);
    evidence["charge_density_change"]["density_definition"] = "rho_alpha + rho_beta";
    evidence["spin_density_change"]["density_definition"] = "rho_alpha - rho_beta";
    for (const auto* name : {"charge_density_change", "spin_density_change"})
      evidence[name]["reduction"] = "L2_norm_of_field_difference";
    evidence["bsh_residual"] = measured(step.previous_bsh_residual, "dimensionless",
        "raw_fixed_point_difference_before_previous_update", !scf.param.get<bool>("conv_only_dens"),
        nullptr, 5.0 * effective, step.iteration - 1);
    evidence["bsh_residual"]["reduction"] = "maximum_over_orbital_and_spin_L2_norms";
    if (step.has_previous && !step.previous_function_residuals.empty())
      evidence["bsh_residual"]["per_function"] = madness::function_residual_record(step.previous_function_residuals);
    // D8: reported against a user-written econv only, never enforced here; the
    // engine check alone still stops the run.
    const auto econv = requested_econv(scf.param);
    const json econv_target = econv ? json(*econv) : json(nullptr);
    evidence["energy_change"]
        = measured(step.energy_change, "hartree",
                   "successive_iteration_energy_evaluations", false,
                   econv_target, econv_target, step.iteration);
    auto& energy_change = evidence["energy_change"];
    energy_change["set_by"] = econv ? "user" : "engine_default_not_a_target";
    if (!econv)
      energy_change["engine_default"] = scf.param.econv();
    if (econv) {
      // Requested criteria compare with <= (evaluate_explicit_criteria).
      energy_change["relation"] = "less_than_or_equal";
      energy_change["satisfied"]
          = energy_change.at("value").is_number()
                ? json(energy_change.at("value").get<double>() <= *econv)
                : json(nullptr);
    }
    evidence["delta_kinetic_energy"] = measured(step.kinetic_change, "hartree",
        "successive_iteration_kinetic_evaluations", false, nullptr, nullptr, step.iteration);
    const auto current_blocks = kinetic_blocks(step.input_kinetic_diagonal);
    if (!step.has_previous) previous_kinetic_blocks_ = json::object();
    evidence["delta_kinetic_energy"].update(madness::function_metric_trace_change(current_blocks, previous_kinetic_blocks_));
    evidence["delta_kinetic_energy"]["quantity"] = "absolute_kinetic_trace_change";
    previous_kinetic_blocks_ = current_blocks;
    iteration_history_.push_back({{"measurements", evidence}, {"kinetic_trace", metrics.at("kinetic_trace")},
        {"occupation_weighted_physical_kinetic_energy", metrics.at("occupation_weighted_physical_kinetic_energy")},
        {"explicit_accuracy", step.explicit_accuracy}, {"configured_progress", step.configured_progress}});
    metrics["convergence_report"] = {{"schema", "madqc.convergence.v1"},
        {"gate_policy", step.explicit_accuracy.value("status", "not_requested") == "not_requested"
            ? "existing_engine" : "existing_engine_and_explicit_accuracy"},
        {"explicit_accuracy", step.explicit_accuracy}, {"configured_progress", step.configured_progress}, {"diagnostic_iterations_per_rung", scf.param.get<int>("diagnostic.iterations_per_rung")}, {"measurements", evidence}, {"iteration_history", iteration_history_}};
    publish(scf, false, json::array(), metrics,
            {{"iteration", step.iteration},
             {"iterations_this_run", iterations},
             {"phase", step.final_rotation ? "after_final_rotation"
                                           : "after_orbital_update"}},
            false);
  }

  // Only these final properties have actually been evaluated on the saved
  // state.
  void finish(madness::SCF &scf, bool converged) {
    if (work == "loaded" && iterations == 0)
      return;
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
    auto metrics = has_latest_ ? latest_.metrics : json::object();
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
        {"delta_kinetic_energy",
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
    evidence["delta_kinetic_energy"].update(madness::function_metric_trace_change(
        saved_trace.at("blocks"), comparison_blocks));
    evidence["delta_kinetic_energy"]["valid_for_saved_state"] =
        evidence["delta_kinetic_energy"].at("availability") == "measured";
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
    const auto accuracy_targets = madness::requested_accuracy_targets(scf.param, "a0^(-3/2)", "dimensionless", "hartree");
    const auto requested_accuracy = madness::evaluate_explicit_criteria(evidence, accuracy_targets);
    metrics["convergence_report"]["saved_state_validation"]["explicit_accuracy"] = requested_accuracy;
    if (!accuracy_targets.empty() && requested_accuracy.at("status") != "satisfied") converged = false;
    publish(scf, converged, properties, metrics,
            {{"phase", "final_properties_on_saved_orbitals"},
             {"iterations_this_run", iterations}},
            true);
  }

private:
  void publish(const madness::SCF &scf, bool converged, json properties,
               const json &metrics, const json &position, bool evaluated) {
    namespace checkpoint = madness::checkpoint;
    checkpoint::Manifest m;
    m.stage = stage_.name;
    m.method = "HF";
    m.generation = "g-" + run_id_ + "-" + std::to_string(++sequence_);
    m.physical = physical_;
    m.numerical = {{"threshold", madness::FunctionDefaults<3>::get_thresh()},
                   {"k", madness::FunctionDefaults<3>::get_k()},
                   {"representation", "mo"},
                   {"protocol", scf.param.protocol()},
                   {"parameters", scf.param.to_json()}};
    m.working_orbitals = {{"state_id", m.generation},
                          {"ordering", "archive_order"}};
    m.position = position;
    if (has_latest_)
      m.position["previous_generation"] = latest_.generation;
    m.metrics = metrics;
    m.targets = {{"dconv", scf.param.dconv()},
                 {"effective_dconv",
                  std::max(scf.param.dconv(),
                           madness::FunctionDefaults<3>::get_thresh())}};
    if (const auto econv = requested_econv(scf.param)) // D8: user-written only
      m.targets["econv"] = *econv;
    m.targets["explicit_accuracy"] = madness::requested_accuracy_targets(scf.param, "a0^(-3/2)", "dimensionless", "hartree");
    m.stop_reason = converged   ? "converged"
                    : evaluated ? "iteration_limit"
                                : "in_progress";
    if (evaluated && !converged && !m.targets.at("explicit_accuracy").empty() &&
        m.metrics.at("convergence_report").at("saved_state_validation")
            .at("explicit_accuracy").at("status") != "satisfied")
      m.stop_reason = "unmet_requested_accuracy";
    for (auto &p : properties)
      p["generation"] = m.generation;
    m.properties = properties;
    const auto directory = store_ / m.generation;
    checkpoint::detail::on_root(world_, [&] {
      validate_measurements(m.metrics, m.stage, m.numerical);
      fs::create_directories(store_);
      require(fs::create_directory(directory), "generation already exists");
      return std::string{};
    });
    auto meta = scf.orbital_metadata();
    if (!evaluated)
      meta.current_energy =
          1.e10; // input energy is not an output-state property
    if (!converged) {
      meta.converged_for_thresh = 1.e10;
      meta.converged_for_dconv = 1.e10;
    }
    std::string error =
        scf.write_mos(world_, (directory / "state.restartdata").string(), meta);
    // Inventory is built only after every writer has closed. The header records
    // the actual writer count, including the parallel archive's clamp.
    checkpoint::detail::collective_error(world_, error);
    world_.gop.fence();
    if (world_.rank() == 0) {
      try {
        const int writers = checkpoint::detail::header(
            directory / "state.restartdata.00000", true);
        m.archives.push_back(checkpoint::inspect_archive(
            directory, "orbitals", "state.restartdata", {}, writers));
      } catch (const std::exception &e) {
        error = e.what();
      }
    }
    checkpoint::publish(world_, store_, m, error);
    // Reading back also gives every rank the exact, validated published
    // inventory.
    latest_ = checkpoint::read_published(world_, store_);
    has_latest_ = true;
  }
};
} // namespace named_stages
