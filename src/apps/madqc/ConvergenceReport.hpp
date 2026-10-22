#pragma once

#include <nlohmann/json.hpp>

namespace named_stages {
// Reporting only. Preserve the gate policy recorded by the application.
// Older generations may lack samples; a successful engine verdict must never
// manufacture missing measurement evidence.
inline nlohmann::json unavailable_measurement(const char* reason) {
  return {{"value", nullptr}, {"availability", "unavailable"},
          {"reason", reason}, {"satisfied", nullptr}};
}

inline nlohmann::json ground_convergence_report(const nlohmann::json& metrics,
                                                const std::string& state_id,
                                                bool converged) {
  auto report = metrics.value("convergence_report", nlohmann::json::object());
  report["schema"] = "madqc.convergence.v1";
  // Generations from before the switch were decided by the engine check.
  if (!report.contains("decided_by"))
    report["decided_by"] = "engine_check";
  report["state_id"] = state_id;
  report["outcome"] = converged ? "converged" : "unfinished";
  if (report.contains("saved_state_validation"))
    report["saved_state_validation"]["state_id"] = state_id;
  for (const char* name : {"density_change", "bsh_residual", "energy_change", "delta_kinetic_energy"})
    if (!report["measurements"].contains(name))
      report["measurements"][name] = unavailable_measurement("checkpoint_predates_measurement_reporting");
  return report;
}

inline nlohmann::json response_convergence_report(const nlohmann::json& metadata,
                                                  const std::string& state_id,
                                                  const nlohmann::json& reference) {
  using json = nlohmann::json;
  // The runner overwrites decided_by from the stage's settings; generations
  // from before the switch were decided by the engine check.
  json report = {{"schema", "madqc.convergence.v1"},
                 {"decided_by", "engine_check"},
                 {"state_id", state_id},
                 {"reference", reference},
                 {"channels", json::array()}};
  for (const auto& channel : metadata.at("fd_states").items()) {
    for (const auto& rung : channel.value().items()) {
      for (const auto& point : rung.value().items()) {
        const auto& entry = point.value();
        auto measurements = entry.value("convergence_measurements", json::object());
        for (const char* name : {"density_change", "bsh_residual", "energy_change", "delta_kinetic_energy"})
          if (!measurements.contains(name))
            measurements[name] = unavailable_measurement("checkpoint_predates_measurement_reporting");
        // A driven response has no local energy-change criterion, regardless
        // of whether its generation predates the reporting schema.
        measurements["energy_change"] = {{"value", nullptr}, {"availability", "not_applicable"},
            {"active_gate", false}, {"reason", "ground_reference_evidence_required"}};
        const bool converged = entry.value("converged", false) && !entry.value("accepted", false);
        const char* outcome = entry.value("diverged", false) ? "invalid"
            : converged ? "converged" : entry.value("stalled", false) ? "stalled" : "unfinished";
        report["channels"].push_back({{"channel", channel.key()}, {"frequency", entry.at("freq")},
            {"protocol", metadata.at("protocols").at(rung.key())}, {"archive", entry.at("archive")},
            {"state_id", state_id}, {"outcome", outcome}, {"iterations_at_save", entry.at("iter")},
            {"convention", entry.value("measurement_convention", json::object())},
            {"measurements", measurements}});
        if (entry.contains("density_history"))
          report["channels"].back()["density_history"] = entry.at("density_history");
        if (entry.contains("function_bsh_history"))
          report["channels"].back()["function_bsh_history"] = entry.at("function_bsh_history");
        if (entry.contains("function_kinetic_history"))
          report["channels"].back()["function_kinetic_history"] = entry.at("function_kinetic_history");
        if (entry.contains("saved_state_validation")) {
          auto validation = entry.at("saved_state_validation");
          validation["state_id"] = state_id;
          report["channels"].back()["saved_state_validation"] = validation;
        }
      }
    }
  }
  return report;
}
} // namespace named_stages
