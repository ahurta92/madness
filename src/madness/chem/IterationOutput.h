#pragma once

// Shared presentation and driver-rank elapsed intervals. This adds no fences
// or numerical work; broad steps include their existing completion boundaries.
#include <madness/chem/ConvergenceDecision.h>
#include <madness/chem/IterationHistory.h>
#include <madness/chem/IterationMetrics.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <ostream>
#include <sstream>

namespace madness::iteration {

/// How much the shared print functions show. The levels are defined by what
/// the ground state and response print today, at print_level 0, 1, 2 and 3+.
enum class Verbosity : std::uint8_t {
    silent = 0,     // warnings only
    progress = 1,   // one table per update, its timing, checkpoint and restart notices
    breakdown = 2,  // energy breakdown per update
    detail = 3      // timings of the final orbital preparation and spectrum
};

inline bool prints(Verbosity current, Verbosity needed) { return current >= needed; }

/// The one mapping from the input key print_level, for every application.
inline Verbosity verbosity(int print_level) {
    if (print_level <= 0) return Verbosity::silent;
    if (print_level == 1) return Verbosity::progress;
    if (print_level == 2) return Verbosity::breakdown;
    return Verbosity::detail;
}

/// The level each shared print call needs: the one place that says what each
/// level shows. write_value_table is listed by what it is used to print.
namespace needed {
// silent: printed at every level
inline constexpr Verbosity targets_below_floor = Verbosity::silent;     // write_targets_below_floor
inline constexpr Verbosity unconverged_diagnostics = Verbosity::silent; // write_unconverged_diagnostics
// progress
inline constexpr Verbosity iteration_metrics = Verbosity::progress;     // write_iteration_metrics
inline constexpr Verbosity iteration_timing = Verbosity::progress;      // write_iteration_timing
inline constexpr Verbosity checkpoint_saved = Verbosity::progress;      // write_checkpoint_saved
inline constexpr Verbosity restart_summary = Verbosity::progress;       // write_restart_summary
inline constexpr Verbosity checkpoint_timing = Verbosity::progress;     // write_checkpoint_timing
inline constexpr Verbosity start_source = Verbosity::progress;          // write_start_source, write_seed_start
inline constexpr Verbosity seed_refused = Verbosity::progress;          // write_seed_refused
inline constexpr Verbosity run_summary = Verbosity::progress;           // madqc's summary.txt at the end of the log
// breakdown
inline constexpr Verbosity energy_breakdown = Verbosity::breakdown;     // write_value_table, energy terms
// detail
inline constexpr Verbosity final_orbital_timing = Verbosity::detail;    // write_value_table, final orbital steps
} // namespace needed

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
                state_notes.push_back(name + ": limit = target * (1 + scale) (target " +
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

// Once per task, before its first iteration: final-level targets below the
// protocol floor (factor * threshold). They are used as written; silent if none.
inline void write_targets_below_floor(std::ostream& os, const std::string& application,
                                      const std::string& task_id, const ResolvedSettings& final_settings,
                                      const MetricLabels& labels = {}) {
    bool any = false;
    for (const auto& criterion : final_settings.criteria) {
        if (!(criterion.requested < criterion.floor)) continue;
        if (!any) os << "Warning: " << application << " " << task_id
                     << ": requested targets are below the final protocol floor\n";
        any = true;
        os << "  " << output_detail::label(criterion.metric, labels) << ": requested "
           << output_detail::number(criterion.requested, 4, true) << ", floor "
           << output_detail::number(criterion.floor, 4, true) << '\n';
    }
    if (any) os << "  The requested targets are used; this protocol may not reach them. A finer final level would.\n";
}

// A run that ends at its iteration limit with diagnostic values: each value is
// marked NOT CONVERGED, followed by every check of the state it came from.
inline void write_unconverged_diagnostics(std::ostream& os, const IterationMetrics& metrics,
                                          const Decision& decision, const std::vector<Diagnostic>& diagnostics,
                                          const MetricLabels& labels = {}) {
    for (const auto& diagnostic : diagnostics) {
        std::ostringstream value;
        value << std::setprecision(10) << diagnostic.value;
        os << "Diagnostic value, NOT CONVERGED: " << output_detail::words(diagnostic.name)
           << (diagnostic.component.empty() ? "" : " " + diagnostic.component) << " = " << value.str() << '\n'
           << "  task " << metrics.info.task_id << ", state " << diagnostic.subject << ", protocol "
           << output_detail::number(metrics.info.protocol.thresh, 0, true) << ", after "
           << metrics.info.iteration << " updates (" << output_detail::words(decision.reason) << ")\n";
        std::vector<std::vector<std::string>> rows;
        for (const auto& check : decision.checks) {
            if (check.subject != diagnostic.subject || check.result == CheckResult::not_applicable) continue;
            const auto state = std::find_if(metrics.states.begin(), metrics.states.end(),
                                            [&](const auto& s) { return s.id == check.subject; });
            std::string measured = "-";
            if (state != metrics.states.end())
                if (const auto found = state->metrics.find(check.metric);
                    found != state->metrics.end() && found->second.value)
                    measured = output_detail::number(*found->second.value, 1, true);
            // A check without a usable measurement names its reason, as the metric table does.
            std::string result = check.result == CheckResult::met ? "met"
                               : check.result == CheckResult::unmet ? "NOT MET"
                               : output_detail::words(check_result_name(check.result));
            if (check.result != CheckResult::met && check.result != CheckResult::unmet && !check.reason.empty())
                result += " (" + output_detail::words(check.reason) + ")";
            rows.push_back({"  " + output_detail::label(check.metric, labels), measured,
                            "limit " + (check.limit ? output_detail::number(*check.limit, 1, true) : std::string("-")),
                            result});
        }
        output_detail::table(os, rows, {16, 7, 13, 7}, {false, true, false, false});
    }
}

// Call only after successful publication, never merely from a checkpoint request.
inline void write_checkpoint_saved(std::ostream& os, const Decision& decision) {
    os << "Checkpoint saved.";
    if (decision.action != Action::update) os << " No further update at this protocol.";
    os << '\n';
}

// A state that starts from another saved state of the same subject at
// another frequency (a nearby frequency): one line naming both frequencies.
inline void write_start_source(std::ostream& os, const std::string& subject, double frequency,
                               double source_frequency) {
    std::ostringstream line;
    line << std::fixed << std::setprecision(5) << "Start: " << subject << " at frequency " << frequency
         << ", from the saved state at frequency " << source_frequency << '\n';
    os << line.str();
}

// A state that starts from a seed: a starting guess made elsewhere, kept in
// its own folder.
inline void write_seed_start(std::ostream& os, const std::string& subject, const std::string& seed_folder) {
    os << "Start: " << subject << ", from the seed in " << seed_folder << '\n';
}

// A seed whose settings differ from this calculation's is not used.
inline void write_seed_refused(std::ostream& os, const std::string& seed_folder) {
    os << "Seed ignored: " << seed_folder << " was made for another geometry or model\n";
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
