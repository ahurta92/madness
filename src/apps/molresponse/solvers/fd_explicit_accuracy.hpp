#pragma once

#include <madness/chem/ExplicitCriteria.h>
#include <madness/chem/FunctionMetricSummary.h>

namespace molresponse_v3 {
// One channel/frequency only. These are opt-in criteria in addition to the
// existing engine gate; the inherited quadrature gate is not relabeled a max.
template<class Parameters>
nlohmann::json fd_requested_accuracy_targets(const Parameters& parameters) {
  if (parameters.is_user_defined("accuracy.energy_change"))
    throw std::invalid_argument("response accuracy.energy_change is not applicable: no local driven-response energy");
  auto targets = madness::requested_accuracy_targets(parameters,
      "kernel_induced_density_L2_atomic_units", "response_amplitude_L2_atomic_units",
      "hartree_times_response_amplitude_scale_squared");
  if (targets.contains("bsh_residual"))
    targets["bsh_residual"]["quantity"] = "maximum_raw_fixed_point_function_L2";
  if (targets.contains("charge_density_change"))
    targets["charge_density_change"]["quantity"] = "kernel_induced_charge_density_change_L2";
  return targets;
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
  if (measurements.contains("delta_kinetic_energy"))
    result["delta_kinetic_energy"] = measurements.at("delta_kinetic_energy");
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
          !entry.value("saved_state_validation", json::object()).contains("explicit_accuracy");
    const auto& targets = entry.at("checkpoint_accuracy_targets");
    if (!targets.is_object() || targets.empty()) return false;
    const json definitions = {
        {"bsh_residual", {{"quantity", "maximum_raw_fixed_point_function_L2"}, {"units", "response_amplitude_L2_atomic_units"}}},
        {"charge_density_change", {{"quantity", "kernel_induced_charge_density_change_L2"}, {"units", "kernel_induced_density_L2_atomic_units"}}},
        {"delta_kinetic_energy", {{"quantity", "absolute_kinetic_trace_change"}, {"units", "hartree_times_response_amplitude_scale_squared"}}}};
    for (const auto& target : targets.items()) {
      if (!definitions.contains(target.key()) || target.value().at("quantity") != definitions.at(target.key()).at("quantity") ||
          target.value().at("units") != definitions.at(target.key()).at("units")) return false;
    }
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
    if (input.at("status") == "invalid" || output.at("status") == "invalid" ||
        entry.at("explicit_accuracy") != input || criteria.back().at("criteria") != input ||
        saved.at("explicit_accuracy") != output || saved.at("affects_stopping") != true ||
        saved.at("comparison") != kinetic.at(iterations - 1).at("sample")) return false;
    for (const char* name : {"bsh_residual", "density_change", "kinetic_energy"}) {
      const auto& measured = saved.at("measurements").at(name);
      if (!sample_valid(measured.at("sample"), iterations) ||
          measured.at("phase") != "saved_response_before_archive" || measured.at("valid_for_saved_state") != true) return false;
    }
    if (entry.value("accepted", false)) return false;
    if (entry.value("converged", false) && (input.at("status") != "satisfied" || output.at("status") != "satisfied")) return false;
    return true;
  } catch (const std::exception&) {
    return false;
  }
}
} // namespace molresponse_v3
