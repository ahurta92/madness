#pragma once

#include <nlohmann/json.hpp>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace madness {
// Evaluate only explicitly supplied targets, with no universal threshold or
// fallback to another quantity. Solver progress/validity are separate inputs
// to a policy; satisfying these scalar criteria alone is not solver success.
inline nlohmann::json evaluate_explicit_criteria(const nlohmann::json& measurements,
                                                const nlohmann::json& targets) {
    using json = nlohmann::json;
    json result = {{"status", "not_requested"}, {"criteria", json::object()},
                   {"comparison_relation", "less_than_or_equal"}};
    bool invalid = false, unavailable = false, unmet = false;
    if (!targets.is_object() || !measurements.is_object()) {
        result["status"] = "invalid";
        result["reason"] = "expected_metric_and_target_objects";
        return result;
    }
    for (const auto& target : targets.items()) {
        json criterion = {{"target", target.value()}, {"value", nullptr}, {"satisfied", nullptr}};
        try {
            const double limit = target.value().at("value").get<double>();
            const auto units = target.value().at("units").get<std::string>();
            if (!std::isfinite(limit) || limit < 0.0 || units.empty())
                throw std::runtime_error("invalid_target");
            if (!measurements.contains(target.key()) ||
                measurements.at(target.key()).value("availability", "") == "unavailable") {
                criterion["status"] = "unavailable";
                unavailable = true;
            } else {
                const auto& metric = measurements.at(target.key());
                if (metric.at("availability") != "measured" || metric.at("units") != units)
                    throw std::runtime_error("incompatible_or_invalid_measurement");
                const double value = metric.at("value").get<double>();
                if (!std::isfinite(value) || value < 0.0)
                    throw std::runtime_error("invalid_measurement_value");
                criterion["value"] = value;
                criterion["satisfied"] = value <= limit;
                criterion["status"] = value <= limit ? "satisfied" : "unmet";
                unmet = unmet || value > limit;
            }
        } catch (const std::exception& error) {
            criterion["status"] = "invalid";
            criterion["reason"] = error.what();
            invalid = true;
        }
        result["criteria"][target.key()] = std::move(criterion);
    }
    if (!targets.empty())
        result["status"] = invalid ? "invalid" : unavailable ? "unavailable" : unmet ? "unmet" : "satisfied";
    return result;
}
// History contains evaluate_explicit_criteria(...)["criteria"] snapshots
// from one comparable solve segment. The caller resets it at every change of
// representation/reference and supplies whether an immediate next rung exists.
// This helper never declares engine convergence or creates a protocol rung.
inline nlohmann::json explicit_progress(const nlohmann::json& history,
                                       int window, double minimum_relative_improvement,
                                       bool has_next_supplied_rung, bool invalid_state) {
    using json = nlohmann::json;
    json result = {{"status", "progressing"}, {"action", "continue"},
                   {"accuracy_satisfied", false}, {"window", window},
                   {"minimum_relative_improvement", minimum_relative_improvement}};
    const auto finish = [&](const char* status, const char* action) {
        auto output = result;
        output["status"] = status;
        output["action"] = action;
        return output;
    };
    if (invalid_state) return finish("invalid", "stop_invalid");
    if (window < 0 || !std::isfinite(minimum_relative_improvement) ||
        minimum_relative_improvement <= 0.0 || minimum_relative_improvement >= 1.0)
        return finish("invalid_configuration", "stop_invalid");
    if (window == 0) return finish("disabled", "continue");
    if (!history.is_array()) return finish("invalid_history", "stop_invalid");
    if (history.empty()) return finish("insufficient_history", "continue");
    try {
        const auto& current = history.back();
        if (!current.is_object() || current.empty()) return finish("insufficient_history", "continue");
        bool all_satisfied = true;
        for (const auto& metric : current.items()) {
            const auto status = metric.value().at("status").get<std::string>();
            if (status == "invalid") return finish("invalid", "stop_invalid");
            if (status != "satisfied" && status != "unmet")
                return finish("insufficient_history", "continue");
            const double value = metric.value().at("value").get<double>();
            const double target = metric.value().at("target").at("value").get<double>();
            if (!std::isfinite(value) || value < 0.0 || !std::isfinite(target) || target < 0.0 ||
                ((value <= target) != (status == "satisfied")))
                return finish("invalid", "stop_invalid");
            all_satisfied = all_satisfied && status == "satisfied";
        }
        if (all_satisfied) {
            result["accuracy_satisfied"] = true;
            result["status"] = "criteria_satisfied";
            return result; // engine eligibility remains a separate decision
        }
        if (history.size() <= static_cast<std::size_t>(window))
            return finish("insufficient_history", "continue");
        const auto& previous = history.at(history.size() - 1 - window);
        if (!previous.is_object() || previous.size() != current.size())
            return finish("insufficient_history", "continue");
        for (const auto& metric : current.items()) {
            const auto& now = metric.value();
            if (!previous.contains(metric.key())) return finish("insufficient_history", "continue");
            const auto& before = previous.at(metric.key());
            if (before.at("target") != now.at("target"))
                return finish("insufficient_history", "continue");
            const auto prior_status = before.at("status").get<std::string>();
            if (prior_status == "invalid") return finish("invalid", "stop_invalid");
            if (prior_status != "satisfied" && prior_status != "unmet")
                return finish("insufficient_history", "continue");
            const double value = now.at("value").get<double>();
            const double prior = before.at("value").get<double>();
            if (!std::isfinite(value) || !std::isfinite(prior) || value < 0.0 || prior < 0.0)
                return finish("invalid", "stop_invalid");
            if (now.at("status") == "unmet" && value <= (1.0 - minimum_relative_improvement) * prior)
                return result;
        }
        return finish("stalled", has_next_supplied_rung ? "advance_supplied_rung" : "retain_unmet");
    } catch (const std::exception&) {
        return finish("invalid_history", "stop_invalid");
    }
}
} // namespace madness
