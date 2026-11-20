#pragma once

// Value records for the shared iteration contract. Applications retain all
// numerical functions, comparison history, proposals and accelerator workspace.
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace madness::iteration {

enum class Metric : std::uint8_t { bsh_residual, density_change, energy_change, kinetic_change };

inline const char* metric_name(Metric metric) {
    switch (metric) {
    case Metric::bsh_residual: return "bsh_residual";
    case Metric::density_change: return "density_change";
    case Metric::energy_change: return "energy_change";
    case Metric::kinetic_change: return "kinetic_change";
    }
    throw std::invalid_argument("unknown convergence metric");
}
inline Metric metric_from_name(const std::string& name) {
    for (auto metric : {Metric::bsh_residual, Metric::density_change, Metric::energy_change, Metric::kinetic_change})
        if (name == metric_name(metric)) return metric;
    throw std::invalid_argument("unknown convergence metric: " + name);
}
inline bool finite_nonnegative(double value) { return std::isfinite(value) && value >= 0.0; }

struct Protocol {
    double thresh = 0.0;
    int k = 0;
    bool operator==(const Protocol&) const = default;
};

// Within the enclosing calculation, (task, attempt, protocol, iteration)
// identifies the prepared state. Zero precedes the first main update; N follows
// N completed main updates at this protocol in this solve session.
struct IterationInfo {
    std::string task_id;
    std::string attempt_id;
    Protocol protocol;
    std::size_t iteration = 0;
    // Application-validated definitions/reference/operator context. Tolerances
    // do not enter this context, and it is not an independently assigned label.
    std::string context_id;
    bool operator==(const IterationInfo&) const = default;
};

struct Comparison {
    std::string task_id;
    std::string attempt_id;
    Protocol protocol;
    std::size_t iteration = 0;
    std::string subject;
};

enum class Availability : std::uint8_t { measured, unavailable, not_applicable, invalid };
struct Scale {
    Availability status = Availability::unavailable;
    std::optional<double> value;
    std::string reason = "scale_unavailable";
    static Scale measured(double value) { return {Availability::measured, value, {}}; }
};
struct ComponentEvidence {
    std::string id;
    double residual_norm = 0.0;
    std::optional<double> function_norm;
};
struct ScalarEvidence {
    std::optional<double> current, previous;
    std::optional<double> current_norm, previous_norm;
    std::vector<ComponentEvidence> components;
    std::optional<std::string> limiting_id;
};

struct MetricValue {
    Availability status = Availability::unavailable;
    std::optional<double> value;
    std::string reason = "not_measured";
    Scale scale;
    std::optional<Comparison> comparison;
    ScalarEvidence evidence;

    static MetricValue measured(double value, std::optional<double> scale = std::nullopt) {
        MetricValue m;
        m.status = Availability::measured; m.value = value; m.reason.clear();
        if (scale) m.scale = Scale::measured(*scale);
        return m;
    }
    static MetricValue unavailable(std::string reason) {
        MetricValue m; m.reason = std::move(reason); return m;
    }
    static MetricValue not_applicable(std::string reason) {
        auto m = unavailable(std::move(reason)); m.status = Availability::not_applicable; return m;
    }
    static MetricValue invalid(std::string reason) {
        auto m = unavailable(std::move(reason)); m.status = Availability::invalid; return m;
    }
};
struct StateMetrics {
    std::string id;
    // The application supplies this only after physical/root correspondence is
    // established. Per-metric references override it when predecessors differ.
    std::optional<Comparison> comparison;
    std::map<Metric, MetricValue> metrics;
};
struct IterationMetrics {
    IterationInfo info;
    std::vector<StateMetrics> states;
};
struct ExpectedTask {
    std::string task_id;
    std::string context_id;
    std::vector<std::string> states;
};

enum class CheckResult : std::uint8_t { met, unmet, unavailable, not_applicable, invalid };
struct MetricCheck {
    std::string subject;
    Metric metric = Metric::bsh_residual;
    CheckResult result = CheckResult::unavailable;
    std::string reason;
    double absolute = 0.0;
    double relative = 0.0;
    std::optional<double> scale;
    std::optional<double> limit;
};
struct Assessment {
    bool valid = false; // state membership and scientific/iteration context
    bool all_met = false;
    std::string reason;
    std::vector<MetricCheck> checks;
};
enum class Action : std::uint8_t { update, finish_protocol, complete, stop };
struct Decision {
    IterationInfo info;
    Action action = Action::stop;
    std::string reason;
    bool checkpoint = false;
    std::vector<MetricCheck> checks;
};
struct Progress {
    std::size_t main_updates = 0;        // current invocation, current protocol
    std::size_t updates_at_protocol = 0; // includes prior attempts, for cadence
    bool state_valid = true;            // application's numerical validity check
};

} // namespace madness::iteration
