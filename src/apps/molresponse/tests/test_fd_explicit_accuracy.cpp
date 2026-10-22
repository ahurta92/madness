#include "../solvers/fd_explicit_accuracy.hpp"
#include <iostream>

struct Controls {
  nlohmann::json values;
  [[nodiscard]] bool is_user_defined(const std::string& key) const { return values.contains(key); }
  template<class T> [[nodiscard]] T get(const std::string& key) const { return values.at(key).get<T>(); }
};

int main() {
  using nlohmann::json;
  const auto require = [](bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
  };
  try {
    json raw = {{"quantity", "raw_fixed_point_function_L2"}, {"normalization", "none"},
        {"reduction", "none"}, {"availability", "measured"},
        {"blocks", {{"x_alpha", json::array({
            {{"function_index", 0}, {"value", 3.0}, {"availability", "measured"}},
            {{"function_index", 1}, {"value", 4.0}, {"availability", "measured"}}})}}}};
    json measurements = {
        {"bsh_residual", {{"value", 5.0}, {"availability", "measured"},
            {"units", "response_amplitude_L2_atomic_units"}, {"per_function", raw},
            {"phase", "solver_step_input_before_acceleration"}, {"sample", {{"iteration", 2}}},
            {"effective_limit", 6.0}, {"satisfied", true}}},
        {"density_change", {{"value", nullptr}, {"availability", "unavailable"},
            {"units", "kernel_induced_density_L2_atomic_units"}, {"reduction", "kernel_induced_density_L2_norm"}}}};
    auto mapped = molresponse_v3::fd_accuracy_measurements(measurements);
    require(mapped.at("bsh_residual").at("value") == 4.0, "must use raw maximum, not quadrature");
    require(mapped.at("bsh_residual").at("sample") == measurements.at("bsh_residual").at("sample"), "sample changed");
    require(!mapped.at("bsh_residual").contains("effective_limit"), "legacy target inherited");
    require(mapped.at("charge_density_change").at("value").is_null(), "missing density fabricated");
    require(!mapped.contains("delta_kinetic_energy"), "missing trace fabricated");
    json target = {{"bsh_residual", {{"value", 4.5}, {"units", "response_amplitude_L2_atomic_units"},
        {"quantity", "maximum_raw_fixed_point_function_L2"}}}};
    require(madness::evaluate_explicit_criteria(mapped, target).at("status") == "satisfied", "maximum criterion");
    Controls controls{{{"accuracy.bsh_residual", 4.5}, {"dconv", 0.1}}};
    const auto configured = molresponse_v3::fd_requested_accuracy_targets(controls);
    require(configured.at("bsh_residual") == target.at("bsh_residual"), "target mapping");
    require(configured.at("charge_density_change").at("quantity") == "kernel_induced_charge_density_change_L2", "density target mapping");
    // A driven response has no energy target: econv is never read (B8).
    controls.values["econv"] = 1.0;
    require(!molresponse_v3::fd_requested_accuracy_targets(controls).contains(
                "energy_change"),
            "driven response energy target created");
    measurements["bsh_residual"]["availability"] = "invalid";
    require(madness::evaluate_explicit_criteria(molresponse_v3::fd_accuracy_measurements(measurements), target).at("status") == "invalid", "parent invalidity hidden");
    measurements["bsh_residual"]["availability"] = "measured";
    auto wrong_density = measurements;
    wrong_density["density_change"]["reduction"] = "legacy_max_spin";
    require(molresponse_v3::fd_accuracy_measurements(wrong_density).at("charge_density_change").at("availability") == "unavailable", "unknown density relabeled");
    measurements["bsh_residual"].erase("per_function");
    require(madness::evaluate_explicit_criteria(molresponse_v3::fd_accuracy_measurements(measurements), target).at("status") == "unavailable", "quadrature fallback");
    measurements["bsh_residual"]["per_function"] = raw;
    measurements["bsh_residual"]["per_function"]["blocks"] = json::object();
    require(madness::evaluate_explicit_criteria(molresponse_v3::fd_accuracy_measurements(measurements), target).at("status") == "unavailable", "empty maximum fabricated");
    auto kinetic_target = json{{"delta_kinetic_energy", {{"value", 1.0}, {"units", "hartree_times_response_amplitude_scale_squared"},
        {"quantity", "absolute_kinetic_trace_change"}}}};
    measurements["delta_kinetic_energy"] = {{"value", 0.0}, {"availability", "measured"},
        {"units", "hartree_times_response_amplitude_scale_squared"}, {"quantity", "max_diagonal_change"}};
    require(madness::evaluate_explicit_criteria(molresponse_v3::fd_accuracy_measurements(measurements), kinetic_target).at("status") == "unavailable", "legacy kinetic convention accepted");
    raw["blocks"]["x_alpha"][1]["value"] = -1.0;
    measurements["bsh_residual"]["per_function"] = raw;
    require(madness::evaluate_explicit_criteria(molresponse_v3::fd_accuracy_measurements(measurements), target).at("status") == "invalid", "invalid component hidden");
    Controls progress_controls{{{"accuracy.bsh_residual", 0.0}, {"progress.window", 1},
                                {"progress.minimum_relative_improvement", 0.99999999}}};
    auto configuration = molresponse_v3::fd_requested_progress(progress_controls);
    require(configuration.at("minimum_relative_improvement") == 0.99999999, "progress precision lost");
    configuration["protocol"] = {1e-4, 1e-6};
    auto progress_target = target;
    progress_target["bsh_residual"]["value"] = 0.0;
    json progress_history = json::array();
    for (double value : {5.0, 4.0}) {
      json metric = {{"bsh_residual", {{"quantity", "maximum_raw_fixed_point_function_L2"},
          {"units", "response_amplitude_L2_atomic_units"}, {"availability", "measured"}, {"value", value}}}};
      progress_history.push_back(madness::evaluate_explicit_criteria(metric, progress_target).at("criteria"));
    }
    const auto coarse = molresponse_v3::fd_progress_decision(progress_history, configuration, 1e-4, false, 0);
    require(coarse.at("action") == "advance_supplied_rung" && coarse.at("action_applied") == true, "coarse stall action");
    const auto fine = molresponse_v3::fd_progress_decision(progress_history, configuration, 1e-6, false, 0);
    require(fine.at("action") == "retain_unmet", "final stall invented a rung");
    require(molresponse_v3::fd_progress_decision(progress_history, configuration, 1e-4, false, 5).at("action_applied") == false,
            "forced count acted on stall");
    require(molresponse_v3::fd_progress_decision(progress_history, configuration, 1e-4, true, 0).at("action") == "stop_invalid",
            "invalid state promoted");
    bool missing_rung = false;
    try { molresponse_v3::fd_progress_decision(progress_history, configuration, 1e-8, false, 0); }
    catch (const std::invalid_argument&) { missing_rung = true; }
    require(missing_rung, "unlisted rung accepted");
    std::cout << "FD explicit accuracy scalar tests passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
