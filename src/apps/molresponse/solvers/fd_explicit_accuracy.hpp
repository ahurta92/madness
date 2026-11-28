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
} // namespace molresponse_v3
