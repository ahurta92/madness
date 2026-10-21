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
} // namespace madness
