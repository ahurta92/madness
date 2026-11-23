#pragma once

// Convergence decision module: resolves targets, applies the protocol floor
// and decides the action after each iteration. Pure functions of plain values;
// no MADNESS runtime, MPI or files. Applications measure; this module judges.
#include <madness/chem/IterationMetrics.h>

#include <algorithm>
#include <set>

namespace madness::iteration {

// What the user wrote. A metric appears here only if a target was written.
struct Request {
    std::map<Metric, double> targets;
    std::optional<double> floor_factor;  // absent: model value, else shared default
    std::size_t max_updates = 0;
    std::size_t checkpoint_interval = 0; // zero: endpoints only
};

enum class Support : std::uint8_t { supported, not_applicable, unsupported };

// plain:  value <= T
// scaled: value <= T * (1 + scale), with the scale measured by the application
enum class LimitRule : std::uint8_t { plain, scaled };

// What one application says about one metric, at one protocol level.
struct ModelMetric {
    Support support = Support::unsupported;
    bool default_criterion = false;       // always a criterion for this solver
    std::optional<double> default_target; // used when the user wrote no target
    LimitRule rule = LimitRule::plain;
    std::optional<double> floor_factor;   // overrides the shared default
};
using ModelMetrics = std::map<Metric, ModelMetric>;

inline constexpr double shared_floor_factor = 5.0;

// The floor bounds only metrics that track the representation's resolution.
inline bool has_floor(Metric metric) {
    return metric == Metric::bsh_residual || metric == Metric::density_change;
}

struct ResolvedCriterion {
    Metric metric = Metric::bsh_residual;
    double requested = 0.0; // the user's target, or the default
    double floor = 0.0;     // factor times the level's threshold; zero without a floor
    double effective = 0.0; // intermediate: max(requested, floor); final: requested
    LimitRule rule = LimitRule::plain;
};

struct ResolvedSettings {
    Protocol protocol;
    std::size_t max_updates = 0;
    std::size_t checkpoint_interval = 0;
    std::vector<ResolvedCriterion> criteria;
    std::vector<Metric> not_applicable; // written targets this application cannot judge
};

// The floor raises targets at intermediate levels only; the final level judges
// the request as written. `application` only names the solver in error messages.
inline ResolvedSettings resolve(const Request& request, const ModelMetrics& model, Protocol protocol,
                                bool final_protocol, const std::string& application = "this application") {
    const auto fail = [&](const std::string& what) {
        throw std::invalid_argument(what + " (" + application + ")");
    };
    if (!std::isfinite(protocol.thresh) || protocol.thresh <= 0.0 || protocol.k <= 0)
        fail("invalid protocol representation");
    if (request.floor_factor && !finite_nonnegative(*request.floor_factor))
        fail("convergence.floor_factor must be finite and not negative");
    for (const auto& [metric, target] : request.targets) {
        const std::string key = std::string("convergence.") + metric_name(metric);
        if (!finite_nonnegative(target)) fail(key + " must be finite and not negative");
        const auto found = model.find(metric);
        const auto support = found == model.end() ? Support::unsupported : found->second.support;
        if (support == Support::unsupported) fail(key + " is not supported by " + application);
    }
    ResolvedSettings settings{protocol, request.max_updates, request.checkpoint_interval, {}, {}};
    for (const auto& [metric, target] : request.targets)
        if (model.at(metric).support == Support::not_applicable) settings.not_applicable.push_back(metric);
    for (const auto& [metric, definition] : model) {
        const auto written = request.targets.find(metric);
        if (definition.support == Support::not_applicable) continue;
        if (!definition.default_criterion && written == request.targets.end()) continue;
        const std::string name = metric_name(metric);
        if (definition.support != Support::supported)
            fail("default criterion " + name + " is not supported by its own model");
        ResolvedCriterion criterion;
        criterion.metric = metric;
        criterion.rule = definition.rule;
        if (written != request.targets.end()) criterion.requested = written->second;
        else if (definition.default_target && finite_nonnegative(*definition.default_target))
            criterion.requested = *definition.default_target;
        else fail("no target for convergence criterion " + name);
        if (definition.floor_factor && !finite_nonnegative(*definition.floor_factor))
            fail("invalid model floor factor for " + name);
        const double factor = request.floor_factor ? *request.floor_factor
                            : definition.floor_factor ? *definition.floor_factor
                                                      : shared_floor_factor;
        criterion.floor = has_floor(metric) ? factor * protocol.thresh : 0.0;
        criterion.effective = final_protocol ? criterion.requested
                                             : std::max(criterion.requested, criterion.floor);
        settings.criteria.push_back(criterion);
    }
    if (settings.criteria.empty()) fail("no convergence criteria");
    return settings;
}

// The scale an application measured for a metric, taken from the measurement
// or derived from its evidence. Unavailable when neither supplies one.
inline Scale measurement_scale(Metric metric, const MetricValue& measurement) {
    if (measurement.scale.status != Availability::unavailable) return measurement.scale;
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

// Unavailable or invalid evidence is never met; scaled never falls back to plain.
inline MetricCheck check_metric(const ResolvedCriterion& criterion, const MetricValue& measurement,
                                const std::optional<Comparison>& predecessor, const IterationInfo& info) {
    MetricCheck check;
    const bool scaled = criterion.rule == LimitRule::scaled;
    check.metric = criterion.metric;
    check.absolute = criterion.effective;
    check.relative = scaled ? criterion.effective : 0.0;
    if (!finite_nonnegative(criterion.effective)) throw std::invalid_argument("invalid resolved target");
    if (measurement.status != Availability::measured) {
        check.result = measurement.status == Availability::invalid ? CheckResult::invalid : CheckResult::unavailable;
        check.reason = measurement.reason.empty() ? "measurement_unavailable" : measurement.reason;
        return check;
    }
    if (!measurement.value || !finite_nonnegative(*measurement.value)) {
        check.result = CheckResult::invalid; check.reason = "invalid_measurement"; return check;
    }
    if (criterion.metric != Metric::bsh_residual) {
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
    const auto scale = measurement_scale(criterion.metric, measurement);
    if (scale.status == Availability::measured && scale.value && finite_nonnegative(*scale.value))
        check.scale = scale.value;
    if (scaled && !check.scale) {
        const bool bad = scale.status == Availability::invalid || scale.status == Availability::measured;
        check.result = bad ? CheckResult::invalid : CheckResult::unmet;
        check.reason = "scale_unavailable"; return check;
    }
    const double limit = scaled ? criterion.effective * (1.0 + *check.scale) : criterion.effective;
    if (!finite_nonnegative(limit)) {
        check.result = CheckResult::invalid; check.reason = "invalid_effective_limit"; return check;
    }
    check.limit = limit;
    check.result = *measurement.value <= limit ? CheckResult::met : CheckResult::unmet;
    return check;
}

inline Assessment assess(const ResolvedSettings& settings, const ExpectedTask& expected,
                         const IterationMetrics& metrics) {
    Assessment assessment;
    const auto& info = metrics.info;
    if (expected.task_id.empty() || expected.context_id.empty() || info.attempt_id.empty()
        || info.task_id != expected.task_id || info.context_id != expected.context_id
        || info.protocol != settings.protocol) {
        assessment.reason = "iteration_context_mismatch"; return assessment;
    }
    const std::set<std::string> wanted(expected.states.begin(), expected.states.end());
    std::set<std::string> present;
    for (const auto& subject : metrics.states)
        if (!present.insert(subject.id).second) { assessment.reason = "duplicate_subject"; return assessment; }
    if (wanted.empty() || wanted.contains("") || wanted.size() != expected.states.size() || wanted != present) {
        assessment.reason = "subject_membership_mismatch"; return assessment;
    }
    std::set<Metric> selected;
    for (const auto& criterion : settings.criteria)
        if (!selected.insert(criterion.metric).second)
            throw std::invalid_argument("duplicate resolved convergence criterion");
    if (selected.empty()) throw std::invalid_argument("no resolved convergence criteria");
    assessment.valid = true; assessment.all_met = true;
    const auto missing = MetricValue::unavailable("not_measured");
    for (const auto& subject : metrics.states)
        for (const auto& criterion : settings.criteria) {
            const auto found = subject.metrics.find(criterion.metric);
            auto check = check_metric(criterion, found == subject.metrics.end() ? missing : found->second,
                                      subject.comparison, info);
            assessment.all_met = assessment.all_met && check.result == CheckResult::met;
            check.subject = subject.id;
            assessment.checks.push_back(std::move(check));
        }
    for (const auto& subject : metrics.states)
        for (const auto metric : settings.not_applicable) {
            MetricCheck check;
            check.subject = subject.id; check.metric = metric;
            check.result = CheckResult::not_applicable; check.reason = "model_not_applicable";
            assessment.checks.push_back(std::move(check));
        }
    return assessment;
}

// Actions: stop on invalid state or context; complete (final) or finish the
// level (intermediate) when all are met, or at the cap (final: stop); else
// update. There is no stall exit; limits do not depend on the level's role.
inline Decision decide(const ResolvedSettings& settings, const ExpectedTask& expected,
                       const Progress& progress, const IterationMetrics& metrics, bool final_protocol) {
    auto assessment = assess(settings, expected, metrics);
    Decision decision;
    decision.info = metrics.info;
    decision.checks = std::move(assessment.checks);
    if (!progress.state_valid || !assessment.valid) {
        decision.reason = progress.state_valid ? assessment.reason : "invalid_state";
        return decision;
    }
    if (metrics.info.iteration != progress.main_updates || progress.main_updates > progress.updates_at_protocol) {
        decision.reason = "iteration_context_mismatch";
        return decision;
    }
    if (assessment.all_met) {
        decision.action = final_protocol ? Action::complete : Action::finish_protocol;
        decision.reason = final_protocol ? "final_criteria_met" : "intermediate_criteria_met";
    } else if (progress.main_updates >= settings.max_updates) {
        decision.action = final_protocol ? Action::stop : Action::finish_protocol;
        decision.reason = "iteration_limit";
    } else {
        decision.action = Action::update; decision.reason = "criteria_unmet";
    }
    decision.checkpoint = decision.action != Action::update
        || (settings.checkpoint_interval > 0 && progress.updates_at_protocol > 0
            && progress.updates_at_protocol % settings.checkpoint_interval == 0);
    return decision;
}

// A completion only at the final level, from valid saved measurements that
// meet the current request. Otherwise nothing, and the caller resumes.
inline std::optional<Decision> saved_completion(const ResolvedSettings& settings, const ExpectedTask& expected,
                                               const IterationMetrics& saved, bool final_protocol) {
    if (!final_protocol) return std::nullopt;
    auto assessment = assess(settings, expected, saved);
    if (!assessment.valid || !assessment.all_met) return std::nullopt;
    return Decision{saved.info, Action::complete, "saved_metrics_satisfy_request", false,
                    std::move(assessment.checks)};
}

} // namespace madness::iteration
