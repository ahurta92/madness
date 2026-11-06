#pragma once

#include <madness/chem/ExplicitCriteria.h>
#include <madness/chem/FunctionMetricSummary.h>
#include <vector>

namespace molresponse_v3 {
// One channel/frequency only. These are opt-in criteria in addition to the
// existing engine gate; the inherited quadrature gate is not relabeled a max.
template<class Parameters>
nlohmann::json fd_requested_accuracy_targets(const Parameters& parameters) {
  // Written targets (DESIGN B8): dconv targets the change of the
  // kernel-induced density, compared absolutely; a driven response has no
  // local energy, so there is no energy target.
  auto targets = madness::requested_accuracy_targets(
      parameters, "kernel_induced_density_L2_atomic_units",
      "response_amplitude_L2_atomic_units",
      "hartree_times_response_amplitude_scale_squared", "dconv", "");
  if (targets.contains("bsh_residual"))
    targets["bsh_residual"]["quantity"] = "maximum_raw_fixed_point_function_L2";
  if (targets.contains("charge_density_change"))
    targets["charge_density_change"]["quantity"] = "kernel_induced_charge_density_change_L2";
  return targets;
}

// An absent policy adds no stall threshold. Validate before cached reuse too.
template<class Parameters>
nlohmann::json fd_requested_progress(const Parameters& parameters) {
  const int window = parameters.template get<int>("progress.window");
  const double improvement = parameters.template get<double>("progress.minimum_relative_improvement");
  madness::validate_progress_configuration(window, improvement,
      parameters.is_user_defined("progress.minimum_relative_improvement"),
      !fd_requested_accuracy_targets(parameters).empty());
  if (window == 0) return nlohmann::json::object();
  return {{"window", window}, {"minimum_relative_improvement", improvement}};
}

inline nlohmann::json fd_progress_decision(const nlohmann::json& history,
    const nlohmann::json& configuration, double threshold, bool invalid, int diagnostic,
    bool fixed_intermediate_work = false) {
  using json = nlohmann::json;
  if (!configuration.at("window").is_number_integer())
    throw std::invalid_argument("configured FD progress window must be an integer");
  const auto protocol = configuration.at("protocol").get<std::vector<double>>();
  auto current = protocol.end();
  for (auto rung = protocol.begin(); rung != protocol.end(); ++rung) {
    if (!std::isfinite(*rung) || *rung <= 0.0 || (rung != protocol.begin() && *rung >= *(rung - 1)))
      throw std::invalid_argument("invalid configured FD progress protocol");
    if (*rung == threshold) current = rung;
  }
  if (current == protocol.end()) throw std::invalid_argument("FD progress rung absent from supplied protocol");
  const bool has_next = current + 1 != protocol.end();
  auto decision = madness::explicit_progress(history, configuration.at("window"),
      configuration.at("minimum_relative_improvement"), has_next, invalid);
  decision["has_next_supplied_rung"] = has_next;
  decision["action_applied"] = decision.at("status") == "stalled" && diagnostic == 0 && !invalid;
  if (fixed_intermediate_work && has_next) {
    decision["fixed_intermediate_work"] = true;
    decision["action_applied"] = false;
  }
  return decision;
}

// Input records retain their sample/phase identities. Missing components are
// unavailable, never replaced with the aggregate quadrature or a zero value.
inline nlohmann::json fd_accuracy_measurements(const nlohmann::json& measurements) {
  using json = nlohmann::json;
  json result = json::object();
  if (measurements.contains("density_change")) {
    result["charge_density_change"] = measurements.at("density_change");
    auto& density = result["charge_density_change"];
    const bool convention = density.value("reduction", "") == "kernel_induced_density_L2_norm" &&
        density.value("units", "") == "kernel_induced_density_L2_atomic_units" &&
        (!density.contains("quantity") || density.at("quantity") == "kernel_induced_charge_density_change_L2" ||
         density.at("quantity") == "kernel_induced_density_L2_change");
    density["quantity"] = "kernel_induced_charge_density_change_L2";
    if (!convention) {
      density["availability"] = "unavailable";
      density["value"] = nullptr;
      density["reason"] = "unestablished_density_definition";
    }
    if (measurements.at("density_change").value("availability", "") == "invalid") {
      density["availability"] = "invalid";
      density["value"] = nullptr;
      density["reason"] = "invalid_density_measurement";
    }
  }
  // D-106: the engine records the kinetic trace change as
  // delta_kinetic_energy; the target names it by its quantity.
  if (measurements.contains("delta_kinetic_energy"))
    result["absolute_kinetic_trace_change"] = measurements.at("delta_kinetic_energy");
  if (measurements.contains("bsh_residual")) {
    const auto& aggregate = measurements.at("bsh_residual");
    json metric = aggregate;
    metric["quantity"] = "maximum_raw_fixed_point_function_L2";
    metric["value"] = nullptr;
    metric["availability"] = "unavailable";
    metric["reason"] = "raw_per_function_residuals_missing";
    metric["reduction"] = "maximum_over_stored_response_functions";
    // Legacy gate fields do not define this explicitly requested criterion.
    for (const auto* field : {"target", "effective_limit", "absolute_limit", "relative_factor",
                             "normalization_scale", "satisfied", "maximum_comparison"})
      metric.erase(field);
    metric["active_gate"] = false;
    if (aggregate.contains("per_function")) {
      const auto& raw = aggregate.at("per_function");
      const bool convention = raw.value("quantity", "") == "raw_fixed_point_function_L2" &&
          raw.value("normalization", "") == "none" && raw.value("reduction", "") == "none";
      const auto maximum = madness::function_metric_maximum(raw.value("blocks", json()));
      metric["availability"] = convention ? maximum.at("availability") : json("invalid");
      metric["reason"] = convention ? maximum.at("reason") : json("incompatible_raw_residual_definition");
      metric["value"] = convention ? maximum.at("value") : json(nullptr);
      metric["limiting_function"] = maximum.at("limiting_function");
      if (raw.value("availability", "") != "measured") {
        metric["availability"] = "invalid";
        metric["reason"] = "invalid_raw_residuals";
        metric["value"] = nullptr;
      }
    }
    if (aggregate.value("availability", "") == "invalid") {
      metric["availability"] = "invalid";
      metric["value"] = nullptr;
      metric["reason"] = "invalid_aggregate_measurement";
    }
    result["bsh_residual"] = std::move(metric);
  }
  return result;
}
// A named checkpoint supplies the reference and channel identity. Within it,
// every history and saved comparison must belong to this exact solve/rung.
inline bool valid_fd_explicit_accuracy(const nlohmann::json& entry, double threshold, int k) {
  using json = nlohmann::json;
  try {
    if (!entry.contains("checkpoint_accuracy_targets"))
      return !entry.contains("explicit_accuracy") && !entry.contains("explicit_accuracy_history") &&
          !entry.contains("checkpoint_progress") && !entry.contains("configured_progress") &&
          !entry.contains("configured_progress_history") &&
          !entry.value("saved_state_validation", json::object()).contains("explicit_accuracy");
    const auto& targets = entry.at("checkpoint_accuracy_targets");
    if (!targets.is_object() || targets.empty()) return false;
    const json definitions = {
        {"bsh_residual", {{"quantity", "maximum_raw_fixed_point_function_L2"}, {"units", "response_amplitude_L2_atomic_units"}}},
        {"charge_density_change", {{"quantity", "kernel_induced_charge_density_change_L2"}, {"units", "kernel_induced_density_L2_atomic_units"}}},
        {"absolute_kinetic_trace_change", {{"quantity", "absolute_kinetic_trace_change"}, {"units", "hartree_times_response_amplitude_scale_squared"}}},
        // D-106: a request saved before the rename.
        {"delta_kinetic_energy", {{"quantity", "absolute_kinetic_trace_change"}, {"units", "hartree_times_response_amplitude_scale_squared"}}}};
    for (const auto& target : targets.items()) {
      if (!definitions.contains(target.key()) || target.value().at("quantity") != definitions.at(target.key()).at("quantity") ||
          target.value().at("units") != definitions.at(target.key()).at("units")) return false;
    }
    if (!entry.at("checkpoint_iteration_count").is_number_integer() ||
        !entry.at("checkpoint_diagnostic_iterations").is_number_integer() ||
        !entry.at("iter").is_number_integer()) return false;
    const int iterations = entry.at("checkpoint_iteration_count").get<int>();
    const bool pending = entry.value("continuation_pending", false);
    if (iterations < 1 || entry.at("iter") != (pending ? 0 : iterations)) return false;
    if (pending && (entry.at("iteration_at_checkpoint") != iterations ||
        entry.value("converged", false) || entry.value("stalled", false))) return false;
    const int diagnostic = entry.at("checkpoint_diagnostic_iterations").get<int>();
    if (diagnostic < 0 || (diagnostic > 0 && iterations != diagnostic) ||
        entry.at("rung_exhausted") != (!pending && diagnostic > 0)) return false;
    const auto& bsh = entry.at("function_bsh_history");
    const auto& density = entry.at("density_history");
    const auto& kinetic = entry.at("function_kinetic_history");
    const auto& criteria = entry.at("explicit_accuracy_history");
    if (bsh.size() != static_cast<std::size_t>(iterations) || density.size() != bsh.size() ||
        criteria.size() != bsh.size() || kinetic.size() != bsh.size() + 1) return false;
    const auto attempt = kinetic.at(0).at("sample").at("attempt_directory").get<std::string>();
    if (attempt.empty()) return false;
    const auto sample_valid = [&](const json& sample, int iteration) {
      return sample.at("iteration") == iteration && sample.at("threshold") == threshold &&
          sample.at("k") == k && sample.at("attempt_directory") == attempt;
    };
    for (int i = 0; i <= iterations; ++i) {
      if (!sample_valid(kinetic.at(i).at("sample"), i) ||
          kinetic.at(i).at("phase") != (i == 0 ? "initial_solver_iterate" : "solver_iterate_after_acceleration")) return false;
      if (i == iterations) break;
      if (!sample_valid(bsh.at(i).at("sample"), i) || !sample_valid(density.at(i).at("sample"), i) ||
          bsh.at(i).at("phase") != "solver_step_input_before_acceleration" ||
          density.at(i).at("phase") != "solver_step_input_before_acceleration") return false;
      const json measurements = {{"bsh_residual", bsh.at(i)}, {"density_change", density.at(i)},
          {"delta_kinetic_energy", kinetic.at(i).at("change_from_previous")}};
      const auto expected = madness::evaluate_explicit_criteria(fd_accuracy_measurements(measurements), targets);
      if (expected.at("status") == "invalid" || criteria.at(i) !=
          json{{"sample", bsh.at(i).at("sample")}, {"criteria", expected}}) return false;
    }
    const auto& measurements = entry.at("convergence_measurements");
    for (const char* name : {"bsh_residual", "density_change", "kinetic_energy"}) {
      const auto& measured = measurements.at(name);
      if (!sample_valid(measured.at("sample"), iterations - 1) ||
          measured.at("phase") != "solver_step_input_before_acceleration" ||
          measured.at("valid_for_saved_state") != false) return false;
    }
    const auto matches_last = [](const json& measured, const json& history) {
      for (const char* field : {"sample", "phase", "availability", "value", "units", "reduction",
                               "normalization_scale", "effective_limit", "absolute_limit",
                               "relative_factor", "normalization_phase"})
        if (measured.at(field) != history.at(field)) return false;
      return true;
    };
    if (!matches_last(measurements.at("bsh_residual"), bsh.back()) ||
        !matches_last(measurements.at("density_change"), density.back()) ||
        measurements.at("bsh_residual").at("per_function") != bsh.back().at("per_function")) return false;
    const auto input = madness::evaluate_explicit_criteria(fd_accuracy_measurements(measurements), targets);
    const auto& saved = entry.at("saved_state_validation");
    const auto output = madness::evaluate_explicit_criteria(fd_accuracy_measurements(saved.at("measurements")), targets);
    if (input.at("status") == "invalid" || output.at("status") == "invalid"
        || entry.at("explicit_accuracy") != input
        || criteria.back().at("criteria") != input
        || saved.at("explicit_accuracy") != output
        || saved.at("affects_stopping")
               != (entry.value("checkpoint_decided_by", "") != "engine_check")
        || saved.at("comparison") != kinetic.at(iterations - 1).at("sample"))
      return false;
    for (const char* name : {"bsh_residual", "density_change", "kinetic_energy"}) {
      const auto& measured = saved.at("measurements").at(name);
      if (!sample_valid(measured.at("sample"), iterations) ||
          measured.at("phase") != "saved_response_before_archive" || measured.at("valid_for_saved_state") != true) return false;
    }
    if (entry.contains("checkpoint_progress")) {
      const auto& configuration = entry.at("checkpoint_progress");
      const int window = configuration.at("window").get<int>();
      const double improvement = configuration.at("minimum_relative_improvement").get<double>();
      madness::validate_progress_configuration(window, improvement, true, !targets.empty());
      if (window == 0 || configuration.size() != 3) return false;
      const auto& decisions = entry.at("configured_progress_history");
      if (!decisions.is_array() || decisions.size() != criteria.size()) return false;
      json prefix = json::array();
      for (std::size_t i = 0; i < criteria.size(); ++i) {
        prefix.push_back(criteria.at(i).at("criteria").at("criteria"));
        const auto expected = fd_progress_decision(
            prefix, configuration, threshold, false, diagnostic,
            decisions.at(i).value("fixed_intermediate_work", false));
        if (expected.at("action") == "stop_invalid" || decisions.at(i) != expected) return false;
        if (expected.at("action_applied") == true && i + 1 != criteria.size()) return false;
      }
      if (entry.at("configured_progress") != decisions.back() ||
          entry.at("stalled") != (!pending && decisions.back().at("action_applied").get<bool>())) return false;
    } else if (entry.contains("configured_progress") || entry.contains("configured_progress_history")) {
      return false;
    }
    if (entry.value("accepted", false)) return false;
    // When the written targets decided (or, before the switch, both checks
    // did), a converged entry must meet them on input and saved state. When
    // the engine check decided, they are only reported (B8).
    if (entry.value("checkpoint_decided_by", "") != "engine_check"
        && entry.value("converged", false)
        && (input.at("status") != "satisfied"
            || output.at("status") != "satisfied"))
      return false;
    return true;
  } catch (const std::exception&) {
    return false;
  }
}
} // namespace molresponse_v3
