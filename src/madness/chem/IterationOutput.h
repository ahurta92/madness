#pragma once

// Shared presentation and driver-rank elapsed intervals. This adds no fences
// or numerical work; broad steps include their existing completion boundaries.
#include <madness/chem/IterationMetrics.h>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <ostream>
#include <sstream>

namespace madness::iteration {
template<class Clock = std::chrono::steady_clock>
class ElapsedTimer {
    typename Clock::time_point start_ = Clock::now();
    typename Clock::time_point previous_ = start_;
public:
    double seconds() const { return std::chrono::duration<double>(Clock::now() - start_).count(); }
    double lap() {
        const auto now = Clock::now();
        const double value = std::chrono::duration<double>(now - previous_).count();
        previous_ = now;
        return value;
    }
};

struct IterationTiming {
    std::vector<std::pair<std::string, double>> steps;
    double seconds = 0.0;
};

// Labels belong to the application: similarly named metrics can measure
// different model quantities. These strings affect presentation only.
using MetricLabels = std::map<Metric, std::string>;

namespace output_detail {
inline std::string number(double value, int precision, bool scientific) {
    std::ostringstream text;
    text << (scientific ? std::scientific : std::fixed) << std::setprecision(precision) << value;
    return text.str();
}
inline std::string words(std::string text) {
    std::replace(text.begin(), text.end(), '_', ' ');
    return text;
}
inline std::string label(Metric metric, const MetricLabels& labels) {
    if (const auto found = labels.find(metric); found != labels.end()) return found->second;
    switch (metric) {
    case Metric::bsh_residual: return "BSH residual, maximum";
    case Metric::density_change: return "Density change";
    case Metric::energy_change: return "Energy change";
    case Metric::kinetic_change: return "Kinetic change";
    }
    return metric_name(metric);
}
inline void table(std::ostream& os, const std::vector<std::vector<std::string>>& rows,
                  std::vector<std::size_t> widths, const std::vector<bool>& right) {
    for (const auto& row : rows)
        for (std::size_t i = 0; i < row.size(); ++i) widths.at(i) = std::max(widths.at(i), row[i].size());
    std::ostringstream text;
    for (const auto& row : rows) {
        text << "  ";
        for (std::size_t i = 0; i < row.size(); ++i) {
            if (i) text << "  ";
            text << (right.at(i) ? std::right : std::left)
                 << std::setw(static_cast<int>(widths.at(i))) << row[i];
        }
        text << '\n';
    }
    os << text.str();
}
inline void protocol(std::ostream& os, const Protocol& protocol) {
    os << "Threshold: " << number(protocol.thresh, 1, true)
       << "    Polynomial order: " << std::to_string(protocol.k);
}
inline std::string step_label(const std::string& name) {
    if (name == "prepare/evaluate/decide") return "Preparation, evaluation and checks";
    if (name == "evaluate") return "Evaluation";
    if (name == "decide") return "Convergence checks";
    if (name == "checkpoint/report") return "Checkpoint and reporting";
    if (name == "advance") return "Solver update";
    return name;
}
} // namespace output_detail

inline void write_value_table(std::ostream& os, const std::string& title,
                              const std::vector<std::pair<std::string, double>>& values,
                              int precision = 3) {
    os << title << '\n';
    std::vector<std::vector<std::string>> rows;
    rows.reserve(values.size());
    for (const auto& [name, value] : values)
        rows.push_back({name, output_detail::number(value, precision, false)});
    output_detail::table(os, rows, {38, 14}, {false, true});
}

inline void write_iteration_timing(std::ostream& os, const std::string& application,
                                   const IterationInfo& info, const IterationTiming& timing) {
    os << "Timing " << application << ": " << info.task_id
       << " (iteration " << std::to_string(info.iteration) << ")\n";
    std::vector<std::pair<std::string, double>> rows;
    rows.reserve(timing.steps.size() + 1);
    for (const auto& [name, seconds] : timing.steps) rows.emplace_back(output_detail::step_label(name), seconds);
    rows.emplace_back("Total", timing.seconds);
    write_value_table(os, "Iteration wall time (seconds, driver rank)", rows);
}

inline void write_checkpoint_timing(std::ostream& os, const std::string& application,
                                    const std::string& task, double seconds) {
    os << "Checkpoint " << application << ": " << task << '\n';
    write_value_table(os, "Wall time (seconds, driver rank)", {{"Load and checks", seconds}});
}

inline void write_metric_table(std::ostream& os, const IterationMetrics& metrics,
                               const std::vector<MetricCheck>& checks, const MetricLabels& labels = {}) {
    // One set of widths for the whole coupled bundle; state names are headings.
    std::vector<std::vector<std::vector<std::string>>> tables;
    std::vector<std::vector<std::string>> notes;
    std::size_t label_width = 30;
    for (const auto& state : metrics.states) {
        std::vector<std::vector<std::string>> rows{{"Metric", "Measured", "Limit", "Result"}};
        std::vector<std::string> state_notes;
        for (const auto metric : {Metric::bsh_residual, Metric::density_change,
                                  Metric::energy_change, Metric::kinetic_change}) {
            const auto check = std::find_if(checks.begin(), checks.end(), [&](const auto& value) {
                return value.subject == state.id && value.metric == metric;
            });
            const auto found = state.metrics.find(metric);
            const bool selected = check != checks.end();
            if (!selected && found == state.metrics.end()) continue;
            const MetricValue* value = found == state.metrics.end() ? nullptr : &found->second;
            const auto name = output_detail::label(metric, labels);
            label_width = std::max(label_width, name.size());
            std::string measured = "-", limit = "-", result, reason;
            if (value && value->status == Availability::measured && value->value &&
                finite_nonnegative(*value->value))
                measured = output_detail::number(*value->value, 4, true);
            if (selected) {
                if (check->limit) limit = output_detail::number(*check->limit, 4, true);
                else if (check->relative == 0.0 && check->result != CheckResult::not_applicable)
                    limit = output_detail::number(check->absolute, 4, true);
                reason = check->reason;
                switch (check->result) {
                case CheckResult::met: result = "met"; break;
                case CheckResult::unmet: result = "unmet"; break;
                case CheckResult::unavailable: result = "unavailable (required)"; break;
                case CheckResult::not_applicable: result = "not applicable"; measured = "-"; break;
                case CheckResult::invalid: result = "invalid"; measured = "-"; break;
                }
            } else if (measured != "-") result = "reported only";
            else {
                result = value && value->status == Availability::not_applicable ? "not applicable"
                       : value && value->status == Availability::invalid ? "invalid" : "unavailable";
                reason = value ? value->reason : "not_measured";
            }
            rows.push_back({name, measured, limit, result});
            if (!reason.empty() && result != "met" && result != "unmet")
                state_notes.push_back(name + ": " + output_detail::words(reason));
            if (selected && check->relative > 0.0)
                state_notes.push_back(name + ": limit = absolute + relative * scale (relative " +
                    output_detail::number(check->relative, 4, true) + ")");
        }
        tables.push_back(std::move(rows));
        notes.push_back(std::move(state_notes));
    }
    for (std::size_t i = 0; i < metrics.states.size(); ++i) {
        os << "State: " << metrics.states[i].id << '\n';
        output_detail::table(os, tables[i], {label_width, 13, 13, 14}, {false, true, true, false});
        for (const auto& note : notes[i]) os << "  " << note << '\n';
    }
}

inline void write_iteration_metrics(std::ostream& os, const std::string& application,
                                    const IterationMetrics& metrics, const Decision& decision,
                                    const MetricLabels& labels = {}, bool final_protocol = true) {
    os << application << ": " << metrics.info.task_id
       << (final_protocol ? " (final protocol)\n" : " (intermediate protocol)\n");
    output_detail::protocol(os, metrics.info.protocol);
    os << "    Iteration: " << std::to_string(metrics.info.iteration) << "\nConvergence\n";
    write_metric_table(os, metrics, decision.checks, labels);
    if (decision.reason == "final_criteria_met") os << "All applicable selected criteria met.\n";
    else if (decision.reason == "intermediate_criteria_met") os << "Applicable intermediate criteria met; continue to the next protocol.\n";
    else if (decision.reason == "iteration_limit")
        os << (decision.action == Action::finish_protocol
            ? "Iteration limit reached; continue to the next protocol.\n"
            : "Iteration limit reached; final criteria not met.\n");
    else if (decision.action == Action::update) os << "Selected criteria not yet met; continue iteration.\n";
    else os << "Stopped: " << output_detail::words(decision.reason) << ".\n";
}

// Call only after successful publication, never merely from a checkpoint request.
inline void write_checkpoint_saved(std::ostream& os, const Decision& decision) {
    os << "Checkpoint saved.";
    if (decision.action != Action::update) os << " No further update at this protocol.";
    os << '\n';
}

// The caller has already loaded/validated the payload. Assessment is scalar-only
// and may be invalid after conversion; those saved values cannot qualify reuse.
inline void write_restart_summary(std::ostream& os, const std::string& application,
                                   const IterationMetrics& saved, const Protocol& target,
                                   const Assessment& assessment, bool reused,
                                   const MetricLabels& labels = {},
                                   const std::string& disposition = {}) {
    os << application << " restart: " << saved.info.task_id << "\nLoaded checkpoint: ";
    output_detail::protocol(os, saved.info.protocol);
    os << '\n';
    if (!(target == saved.info.protocol)) {
        os << "Current protocol: "; output_detail::protocol(os, target); os << '\n';
    }
    if (assessment.valid) {
        os << "Saved metrics against current tolerances\n";
        write_metric_table(os, saved, assessment.checks, labels);
    }
    if (!disposition.empty()) os << disposition << '\n';
    else if (reused) os << "Saved metrics satisfy all applicable selected tolerances. Reusing saved functions; no solver iterations required.\n";
    else os << "Continuing from saved functions with fresh convergence checks.\n";
}
} // namespace madness::iteration
