#pragma once

#include <madness/chem/IterationMetrics.h>

#include <algorithm>
#include <cmath>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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
            if (measurements.contains(target.key()) &&
                measurements.at(target.key()).value("availability", "") == "invalid")
                throw std::runtime_error("invalid_measurement");
            if (!measurements.contains(target.key()) ||
                measurements.at(target.key()).value("availability", "") == "unavailable" ||
                (target.value().contains("quantity") &&
                 measurements.at(target.key()).value("quantity", "") != target.value().at("quantity").get<std::string>())) {
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
// Numeric controls are opt-in. Sentinel defaults are never accuracy targets.
// Units belong to the application's actual density/amplitude definitions.
// Only a value the user wrote is a target (B3, D8). The density and energy
// targets come from the keywords the application names: SCF uses its own
// dconv and econv (DESIGN B8, "written values are targets"); an empty
// energy keyword means the application has no energy-change target.
template <typename Parameters>
nlohmann::json requested_accuracy_targets(const Parameters& parameters,
                                          const std::string& density_units,
                                          const std::string& bsh_units,
                                          const std::string& kinetic_units,
                                          const std::string& density_keyword,
                                          const std::string& energy_keyword) {
    using json = nlohmann::json;
    struct Definition {
        std::string keyword, metric, units;
    };
    std::vector<Definition> definitions
        = {{density_keyword, "charge_density_change", density_units},
           {"accuracy.bsh_residual", "bsh_residual", bsh_units},
           {"accuracy.kinetic_trace_change", "absolute_kinetic_trace_change",
            kinetic_units}};
    if (!energy_keyword.empty())
        definitions.push_back({energy_keyword, "energy_change", "hartree"});
    json targets = json::object();
    for (const auto& definition : definitions) {
        if (!parameters.is_user_defined(definition.keyword))
            continue;
        const auto value = parameters.template get<double>(definition.keyword);
        if (!std::isfinite(value) || value < 0.0)
            throw std::invalid_argument(definition.keyword
                                        + " must be finite and nonnegative");
        targets[definition.metric]
            = {{"value", value}, {"units", definition.units}};
    }
    if (targets.contains("absolute_kinetic_trace_change"))
        targets["absolute_kinetic_trace_change"]["quantity"] = "absolute_kinetic_trace_change";
    return targets;
}

// Validate even when restart can reuse a solved state without entering a loop.
inline void validate_progress_configuration(int window, double improvement,
                                            bool improvement_supplied, bool has_targets) {
    if (window < 0) throw std::invalid_argument("progress.window must be nonnegative");
    if (improvement_supplied && window == 0)
        throw std::invalid_argument("progress improvement requires a positive window");
    if (window > 0 && (!has_targets || !improvement_supplied || !std::isfinite(improvement) ||
                      improvement <= 0.0 || improvement >= 1.0))
        throw std::invalid_argument("progress requires explicit accuracy targets and a supplied improvement fraction strictly between zero and one");
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

// Shared measurement-based assessment, used by live and restart controllers.
// Legacy JSON entry points above remain until their callers migrate.
namespace madness::iteration {

inline Scale measurement_scale(Metric metric, const MetricValue& measurement) {
    if (measurement.scale.status != Availability::unavailable)
        return measurement.scale;
    const auto& evidence = measurement.evidence;
    if (metric == Metric::energy_change || metric == Metric::kinetic_change) {
        if (evidence.current && evidence.previous && std::isfinite(*evidence.current)
            && std::isfinite(*evidence.previous))
            return Scale::measured(std::max(std::abs(*evidence.current), std::abs(*evidence.previous)));
    } else if (metric == Metric::density_change) {
        if (evidence.current_norm && evidence.previous_norm
            && finite_nonnegative(*evidence.current_norm) && finite_nonnegative(*evidence.previous_norm))
            return Scale::measured(std::max(*evidence.current_norm, *evidence.previous_norm));
    } else if (!evidence.components.empty()) {
        double maximum = 0.0;
        for (const auto& component : evidence.components) {
            if (!component.function_norm || !finite_nonnegative(*component.function_norm))
                return measurement.scale;
            maximum = std::max(maximum, *component.function_norm);
        }
        return Scale::measured(maximum);
    }
    return measurement.scale;
}

inline MetricCheck evaluate_metric(Metric metric, const MetricValue& measurement,
                                   const std::optional<Comparison>& predecessor,
                                   const IterationInfo& info, double absolute, double relative) {
    MetricCheck check;
    check.metric = metric; check.absolute = absolute; check.relative = relative;
    if (!finite_nonnegative(absolute) || !finite_nonnegative(relative))
        throw std::invalid_argument("invalid resolved convergence tolerance");
    if (measurement.status != Availability::measured) {
        check.result = measurement.status == Availability::invalid ? CheckResult::invalid : CheckResult::unavailable;
        check.reason = measurement.reason.empty() ? "measurement_unavailable" : measurement.reason;
        return check;
    }
    if (!measurement.value || !finite_nonnegative(*measurement.value)) {
        check.result = CheckResult::invalid; check.reason = "invalid_measurement"; return check;
    }
    if (metric != Metric::bsh_residual) {
        const auto& comparison = measurement.comparison ? measurement.comparison : predecessor;
        if (!comparison || comparison->attempt_id.empty() || comparison->subject.empty()
            || comparison->task_id != info.task_id
            || !std::isfinite(comparison->protocol.thresh) || comparison->protocol.thresh <= 0.0
            || comparison->protocol.k <= 0
            || (comparison->attempt_id == info.attempt_id && comparison->protocol == info.protocol
                && comparison->iteration >= info.iteration)) {
            check.reason = "comparison_unavailable"; return check;
        }
    }
    const auto scale = measurement_scale(metric, measurement);
    if (scale.status == Availability::measured && scale.value && finite_nonnegative(*scale.value))
        check.scale = scale.value;
    if (relative > 0.0 && !check.scale) {
        check.result = scale.status == Availability::invalid || scale.status == Availability::measured
                           ? CheckResult::invalid : CheckResult::unavailable;
        check.reason = "scale_unavailable"; return check;
    }
    const double limit = absolute + (relative > 0.0 ? relative * *check.scale : 0.0);
    if (!finite_nonnegative(limit)) {
        check.result = CheckResult::invalid; check.reason = "invalid_effective_limit"; return check;
    }
    check.limit = limit;
    check.result = *measurement.value <= limit ? CheckResult::met : CheckResult::unmet;
    return check;
}

} // namespace madness::iteration
