#pragma once

#include <madness/chem/FunctionMetricSummary.h>
#include <algorithm>
#include <limits>

namespace molresponse_v3 {
inline nlohmann::json density_history_record(double value, bool available,
    double absolute, double relative, double norm, const nlohmann::json& sample) {
    using json = nlohmann::json;
    const double limit = absolute + relative * norm;
    const bool valid = std::isfinite(value) && value >= 0.0 &&
        std::isfinite(absolute) && absolute >= 0.0 && std::isfinite(relative) && relative >= 0.0 &&
        std::isfinite(norm) && norm >= 0.0 && std::isfinite(limit);
    return {{"quantity", "kernel_induced_density_L2_change"},
        {"value", available && valid ? json(value) : json(nullptr)},
        {"availability", !valid ? "invalid" : available ? "measured" : "unavailable"},
        {"reason", !valid ? "invalid_density_evidence" : available ? "" : "no_comparable_previous_sample"},
        {"units", "kernel_induced_density_L2_atomic_units"},
        {"phase", "solver_step_input_before_acceleration"},
        {"normalization_phase", "solver_step_input_before_acceleration"},
        {"reduction", "kernel_induced_density_L2_norm"},
        {"absolute_limit", absolute}, {"relative_factor", relative},
        {"normalization_scale", norm}, {"effective_limit", valid ? json(limit) : json(nullptr)},
        {"satisfied", available && valid ? json(value <= limit) : json(nullptr)},
        {"active_gate", true}, {"valid_for_saved_state", false}, {"sample", sample}};
}

// The containing response entry supplies channel/equation/reference identity.
// The caller must keep comparisons inside that entry's unchanged solve segment;
// preparation clears histories on restart, reference change and each new rung.
inline nlohmann::json kinetic_change_record(const nlohmann::json& current,
                                            const nlohmann::json& previous = nullptr) {
    using json = nlohmann::json;
    json result = {{"quantity", "absolute_kinetic_trace_change"},
        {"units", current.at("units")}, {"reduction", "absolute_difference_of_traces"},
        {"weighting", "none"}, {"active_gate", false}, {"target", nullptr},
        {"satisfied", nullptr}, {"availability", "unavailable"},
        {"reason", "no_comparable_previous_sample"}, {"value", nullptr},
        {"blocks", json::object()}, {"limiting_function", nullptr},
        {"phase", current.at("phase")}, {"sample", current.at("sample")},
        {"comparison", previous.is_null() ? json(nullptr) : previous.at("sample")},
        {"valid_for_saved_state", false}};
    if (previous.is_null()) return result;
    const auto& now = current.at("sample");
    const auto& before = previous.at("sample");
    bool comparable = current.at("quantity") == previous.at("quantity") &&
                      current.at("units") == previous.at("units") &&
                      now.at("iteration").get<int>() == before.at("iteration").get<int>() + 1;
    for (const char* key : {"attempt_directory", "threshold", "k"})
        comparable = comparable && now.at(key) == before.at(key);
    // Prior samples are genuine iterates, never a second measurement of output.
    comparable = comparable && (previous.at("phase") == "initial_solver_iterate" ||
                                 previous.at("phase") == "solver_iterate_after_acceleration");
    comparable = comparable && (current.at("phase") == "solver_iterate_after_acceleration" ||
                                  current.at("phase") == "saved_response_before_archive");
    if (!comparable) return result;
    result.update(madness::function_metric_trace_change(current.at("blocks"), previous.at("blocks")));
    result["valid_for_saved_state"] = current.at("valid_for_saved_state").get<bool>() &&
                                      result.at("availability") == "measured";
    return result;
}

// Counterfactual for this one BSH sample: change only the residual reduction.
// Keep the existing output-state norm, absolute/relative limits and <= relation.
// A true result is NOT a complete convergence verdict or permission to stop.
inline nlohmann::json bsh_maximum_comparison(const nlohmann::json& metric) {
    using json = nlohmann::json;
    auto result = madness::function_metric_maximum(metric.at("per_function").at("blocks"));
    result.update({{"quantity", "raw_fixed_point_function_L2"},
        {"policy", "maximum_with_existing_channel_limit_v1"}, {"active_gate", false},
        {"units", metric.at("units")}, {"phase", metric.at("phase")},
        {"sample", metric.at("sample")}, {"effective_limit", nullptr},
        {"satisfied", nullptr}, {"comparison_relation", "less_than_or_equal"},
        {"valid_for_saved_state", metric.at("valid_for_saved_state")}});
    for (const char* key : {"absolute_limit", "relative_factor", "normalization_scale", "effective_limit"})
        if (!metric.contains(key) || metric.at(key).is_null()) return result;
    try {
        const double absolute = metric.at("absolute_limit").get<double>();
        const double relative = metric.at("relative_factor").get<double>();
        const double norm = metric.at("normalization_scale").get<double>();
        const double limit = metric.at("effective_limit").get<double>();
        const double reconstructed = absolute + relative * norm;
        // Readers may use separate multiply/add where the writer used FMA.
        // Allow only machine-rounding differences; keep the recorded gate limit.
        const double rounding = 4.0 * std::numeric_limits<double>::epsilon() *
            std::max({limit, reconstructed, std::numeric_limits<double>::min()});
        if (!std::isfinite(reconstructed) || !std::isfinite(absolute) || !std::isfinite(relative) || !std::isfinite(norm) ||
            !std::isfinite(limit) || absolute < 0.0 || relative < 0.0 || norm < 0.0 ||
            std::abs(limit - reconstructed) > rounding)
            throw std::runtime_error("invalid effective limit");
        const auto& phase = metric.at("phase");
        const auto& norm_phase = metric.at("normalization_phase");
        if (!((phase == "solver_step_input_before_acceleration" &&
               norm_phase == "solver_step_output_after_acceleration") ||
              (phase == "saved_response_before_archive" && norm_phase == phase)))
            throw std::runtime_error("invalid normalization phase");
        result["effective_limit"] = limit;
        result["normalization_phase"] = norm_phase;
        if (result.at("availability") == "measured")
            result["satisfied"] = result.at("value").get<double>() <= limit;
    } catch (const std::exception&) {
        result["availability"] = "invalid";
        result["reason"] = "invalid_effective_limit";
        result["effective_limit"] = nullptr;
    }
    return result;
}
// Optional additions on old checkpoints stay optional. When present, derived
// evidence must agree with its raw samples, not merely say "measured".
inline bool valid_metric_diagnostics(const nlohmann::json& entry) {
    using json = nlohmann::json;
    try {
        const auto check_bsh = [](const json& metric) {
            if (!metric.contains("maximum_comparison")) return true;
            const auto expected = bsh_maximum_comparison(metric);
            return expected.at("availability") != "invalid" &&
                   metric.at("maximum_comparison") == expected;
        };
        if (entry.contains("function_bsh_history"))
            for (const auto& sample : entry.at("function_bsh_history"))
                if (!check_bsh(sample)) return false;
        if (entry.contains("density_history")) {
            if (!entry.at("density_history").is_array()) return false;
            for (const auto& sample : entry.at("density_history")) {
                const bool available = sample.at("availability") == "measured";
                if (!available && sample.at("availability") != "unavailable") return false;
                const auto expected = density_history_record(available ? sample.at("value").get<double>() : 0.,
                    available, sample.at("absolute_limit").get<double>(), sample.at("relative_factor").get<double>(),
                    sample.at("normalization_scale").get<double>(), sample.at("sample"));
                if (expected.at("availability") == "invalid" || sample != expected) return false;
            }
        }
        const auto history = entry.value("function_kinetic_history", json::array());
        if (!history.is_array()) return false;
        json previous = nullptr;
        for (const auto& sample : history) {
            if (sample.contains("change_from_previous")) {
                const auto expected = kinetic_change_record(sample, previous);
                if (expected.at("availability") == "invalid" ||
                    sample.at("change_from_previous") != expected) return false;
            }
            previous = sample;
        }
        const auto check_measurements = [&](const json& measurements, bool saved) {
            if (measurements.contains("bsh_residual") &&
                !check_bsh(measurements.at("bsh_residual"))) return false;
            if (!measurements.contains("delta_kinetic_energy")) return true;
            const auto& delta = measurements.at("delta_kinetic_energy");
            if (delta.value("quantity", "") != "absolute_kinetic_trace_change")
                return true; // older reporting convention
            if (history.size() < 2) return false;
            const auto& input = history.at(history.size() - 2);
            if (!saved) {
                const auto& measured = measurements.at("kinetic_energy");
                for (const char* key : {"sample", "blocks", "quantity", "units", "availability"})
                    if (measured.at(key) != input.at(key)) return false;
                if (measured.at("phase") != "solver_step_input_before_acceleration") return false;
            }
            const auto expected = saved ? kinetic_change_record(measurements.at("kinetic_energy"), input) :
                kinetic_change_record(input, history.size() >= 3 ? history.at(history.size() - 3) : json(nullptr));
            return expected.at("availability") != "invalid" && delta == expected;
        };
        if (entry.contains("convergence_measurements") &&
            !check_measurements(entry.at("convergence_measurements"), false)) return false;
        if (entry.contains("saved_state_validation") &&
            !check_measurements(entry.at("saved_state_validation").at("measurements"), true)) return false;
        return true;
    } catch (const std::exception&) {
        return false;
    }
}
} // namespace molresponse_v3
