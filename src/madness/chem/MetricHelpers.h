#pragma once

#include <madness/chem/IterationMetrics.h>
#include <algorithm>
#include <set>

namespace madness::iteration {

/// Applications supply norms in their declared working representation. The
/// reduction never weights components or normalizes by their number.
inline MetricValue maximum_residual(std::vector<ComponentEvidence> components) {
    if (components.empty()) return MetricValue::invalid("empty_residual");
    double maximum = -1.0, scale = 0.0;
    bool have_scale = true;
    std::string limiting;
    std::set<std::string> ids;
    for (const auto& component : components) {
        if (component.id.empty() || !ids.insert(component.id).second ||
            !finite_nonnegative(component.residual_norm) ||
            (component.function_norm && !finite_nonnegative(*component.function_norm)))
            return MetricValue::invalid("invalid_residual_component");
        if (component.residual_norm > maximum) {
            maximum = component.residual_norm;
            limiting = component.id;
        }
        if (component.function_norm) scale = std::max(scale, *component.function_norm);
        else have_scale = false;
    }
    auto metric = MetricValue::measured(maximum);
    if (have_scale) metric.scale = Scale::measured(scale);
    metric.evidence.components = std::move(components);
    metric.evidence.limiting_id = std::move(limiting);
    return metric;
}

/// Physical scalar changes (energy or kinetic energy). A previous value is
/// usable only when the application also establishes a comparable predecessor.
inline MetricValue scalar_change(double current, std::optional<double> previous,
                                 const std::optional<Comparison>& comparison) {
    if (!std::isfinite(current) || (previous && !std::isfinite(*previous)))
        return MetricValue::invalid("invalid_scalar");
    auto metric = MetricValue::unavailable("comparison_unavailable");
    metric.evidence.current = current;
    if (previous && comparison) {
        const double change = std::abs(current - *previous);
        if (!std::isfinite(change)) return MetricValue::invalid("invalid_scalar_change");
        metric = MetricValue::measured(change, std::max(std::abs(current), std::abs(*previous)));
        metric.evidence.current = current;
        metric.evidence.previous = previous;
        metric.comparison = comparison;
    }
    return metric;
}

struct FieldNorms {
    double current = 0.0;
    std::optional<double> previous, change;
};

/// The application owns the physical fields and calculates their L2 norms and
/// differences. This helper takes the maximum across those fields, not their sum.
inline MetricValue maximum_field_change(const std::vector<FieldNorms>& fields,
                                        const std::optional<Comparison>& comparison) {
    if (fields.empty()) return MetricValue::invalid("empty_fields");
    double current = 0.0, previous = 0.0, change = 0.0;
    bool comparable = comparison.has_value();
    for (const auto& field : fields) {
        if (!finite_nonnegative(field.current) ||
            (field.previous && !finite_nonnegative(*field.previous)) ||
            (field.change && !finite_nonnegative(*field.change)))
            return MetricValue::invalid("invalid_field_norm");
        current = std::max(current, field.current);
        if (field.previous && field.change) {
            previous = std::max(previous, *field.previous);
            change = std::max(change, *field.change);
        } else comparable = false;
    }
    auto metric = MetricValue::unavailable("comparison_unavailable");
    metric.evidence.current_norm = current;
    if (comparable) {
        metric = MetricValue::measured(change, std::max(current, previous));
        metric.evidence.current_norm = current;
        metric.evidence.previous_norm = previous;
        metric.comparison = comparison;
    }
    return metric;
}

} // namespace madness::iteration
