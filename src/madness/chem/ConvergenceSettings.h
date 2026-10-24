#pragma once

// Who decides convergence (DESIGN B8, decision D9) and what every convergence
// metric of a stage is set to (B9).
//
// A stage has two groups of checks:
//   - the engine check: the engine's built-in stopping test (for SCF, density
//     and BSH limits derived from dconv), and
//   - the requested criteria: targets the user wrote in the stage (econv,
//     dconv, accuracy.*).
// One switch, convergence.use_default_criteria, selects which group stops the
// run. Both groups are always measured and recorded; only the selected group
// is enforced. Each application builds a ConvergenceSettings once per stage;
// the printed table, the record and the stopping decision all read that one
// struct, so they cannot disagree (B9 "one source").
//
// Pure: no World, MPI, archives or solver state (B1).

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace madness::convergence {

/// The single spelling of the switch (Adrian, 2026-10-21). true (default): the
/// engine check decides; false: the requested criteria decide.
inline constexpr const char* use_default_criteria_key
    = "convergence.use_default_criteria";

enum class Decider : std::uint8_t { engine_check, requested_criteria };
enum class Group : std::uint8_t { engine_check, requested_criteria };
enum class SetBy : std::uint8_t { engine_default, user, not_set };
enum class Relation : std::uint8_t { less_than, less_than_or_equal };

/// One row of the settings table.
struct MetricSetting {
  Group group = Group::requested_criteria;
  std::string name;    ///< printed name, e.g. "BSH residual (max)"
  std::string keyword; ///< deck keyword that sets the target, e.g. "dconv"
  std::string
      measurement;      ///< key of the measurement record it is compared with
  std::string quantity; ///< required measurement quantity; empty: any
  std::string units;
  std::optional<double> target; ///< the number the user (or engine) set
  SetBy set_by = SetBy::not_set;
  Relation relation = Relation::less_than_or_equal;
  /// What is actually compared at the stage's final protocol rung. Absent when
  /// the limit is not one number (e.g. a relative response gate) or when
  /// there is no target.
  std::optional<double> effective_limit;
  /// How the effective limit follows from the target; empty when it is the
  /// target itself. Shows the hidden factors (B9).
  std::string limit_formula;
  std::string note; ///< e.g. "engine default econv is not a target (D8)"
  /// For a requested row: what the engine check does instead, completing
  /// "the engine check ..." in the engine-mode warning, e.g. "compares
  /// spin-density change against 5.0e-04" or "has no energy term".
  std::string engine_check_instead;
};

/// The iteration cap per protocol rung (B5): maxiter, one number or one per
/// rung, counted across --continue runs.
struct IterationCaps {
  std::vector<int> per_rung; ///< one entry: the same cap at every rung
  SetBy set_by = SetBy::engine_default;
};

struct ConvergenceSettings {
  std::string stage;
  std::string application;
  Decider decider = Decider::engine_check;
  bool decider_set_by_user = false;
  std::vector<MetricSetting> metrics;
  IterationCaps iteration_caps;
  std::vector<std::string> warnings;
};

/// A final measurement of one metric, for the end-of-stage table (B9).
struct MetricOutcome {
  std::string measurement; ///< matches MetricSetting::measurement
  std::optional<double> value;
  std::string availability
      = "unavailable"; ///< measured | unavailable | invalid
};

[[nodiscard]] inline bool enforced(const ConvergenceSettings& settings,
                                   const MetricSetting& metric) {
  if (!metric.target && !metric.effective_limit && metric.limit_formula.empty())
    return false; // nothing to compare against
  return settings.decider == Decider::engine_check
             ? metric.group == Group::engine_check
             : metric.group == Group::requested_criteria;
}

/// The single number a row is compared against: its effective limit, else
/// its target when no formula derives the limit from it. None when the limit
/// is only a formula (e.g. a relative, per-leg response gate).
[[nodiscard]] inline std::optional<double>
comparison_limit(const MetricSetting& metric) {
  if (metric.effective_limit)
    return metric.effective_limit;
  return metric.limit_formula.empty() ? metric.target : std::nullopt;
}

/// B8: the selected group alone decides; there is no AND mode.
[[nodiscard]] inline bool converged(Decider decider, bool engine_check_passed,
                                    bool requested_criteria_met) {
  return decider == Decider::engine_check ? engine_check_passed
                                          : requested_criteria_met;
}

/// How a protocol rung ended (B5): one of the three exits, or not yet
/// (`in_progress`, DESIGN's "none"). The strings are those saved in manifests.
enum class StopReason : std::uint8_t {
  in_progress,
  converged,
  stalled,
  iteration_limit
};

[[nodiscard]] inline const char* to_string(StopReason stop) {
  switch (stop) {
  case StopReason::converged:
    return "converged";
  case StopReason::stalled:
    return "stalled";
  case StopReason::iteration_limit:
    return "iteration_limit";
  case StopReason::in_progress:
    break;
  }
  return "in_progress";
}

[[nodiscard]] inline StopReason
stop_reason_from_string(const std::string& text) {
  if (text == "converged")
    return StopReason::converged;
  if (text == "stalled")
    return StopReason::stalled;
  if (text == "iteration_limit")
    return StopReason::iteration_limit;
  return StopReason::in_progress;
}

/// One sample at an engine's save boundary, as verdicts and counts (B5:
/// engines compare live state; everything outside takes verdicts).
struct RungSample {
  Decider decider = Decider::engine_check;
  bool engine_check_passed = false;
  bool criteria_met = false;
  bool stall_detected = false;
  int iterations_at_rung = 0; ///< in this run, this sample included
  int cap = 0;                ///< maxiter at this rung, per run
  int forced_iterations = 0;  ///< diagnostic.iterations_per_rung; 0: off
};

/// B5: the one place the exits of a rung are ordered. Converged (by the B8
/// decider), else the iteration cap, else a stall, else the rung goes on.
/// Forced iterations (the trajectory-study control) run exactly that many:
/// no exit before, and at the last one converged or iteration limit; the
/// stall is ignored. Divergence is the engine's validity check, made before.
[[nodiscard]] inline StopReason rung_exit(const RungSample& sample) {
  const bool decided = converged(sample.decider, sample.engine_check_passed,
                                 sample.criteria_met);
  if (sample.forced_iterations > 0) {
    if (sample.iterations_at_rung < sample.forced_iterations)
      return StopReason::in_progress;
    return decided ? StopReason::converged : StopReason::iteration_limit;
  }
  if (decided)
    return StopReason::converged;
  if (sample.iterations_at_rung >= sample.cap)
    return StopReason::iteration_limit;
  if (sample.stall_detected)
    return StopReason::stalled;
  return StopReason::in_progress;
}

[[nodiscard]] inline bool
has_requested_criteria(const ConvergenceSettings& settings) {
  return std::any_of(settings.metrics.begin(), settings.metrics.end(),
                     [](const MetricSetting& metric) {
                       return metric.group == Group::requested_criteria
                              && metric.target.has_value();
                     });
}

/// The requested targets in the form evaluate_explicit_criteria() takes.
[[nodiscard]] inline nlohmann::json
requested_targets(const ConvergenceSettings& settings) {
  nlohmann::json targets = nlohmann::json::object();
  for (const auto& metric : settings.metrics) {
    if (metric.group != Group::requested_criteria || !metric.target)
      continue;
    auto& target = targets[metric.measurement];
    target = {{"value", *metric.target}, {"units", metric.units}};
    if (!metric.quantity.empty())
      target["quantity"] = metric.quantity;
  }
  return targets;
}

/// B8: in criteria mode, a request that leaves out a quantity the engine check
/// covers is accepted with a warning. `engine_quantities` names the requested
/// measurements that stand for the engine's checks (e.g. density, BSH).
[[nodiscard]] inline std::vector<std::string>
incomplete_request_warnings(const ConvergenceSettings& settings,
                            const std::vector<std::string>& engine_quantities) {
  std::vector<std::string> warnings;
  if (settings.decider != Decider::requested_criteria)
    return warnings;
  for (const auto& quantity : engine_quantities) {
    bool requested = false;
    for (const auto& metric : settings.metrics)
      requested = requested
                  || (metric.group == Group::requested_criteria
                      && metric.measurement == quantity && metric.target);
    if (!requested)
      warnings.push_back("stage " + settings.stage
                         + ": your criteria decide convergence but set no "
                         + quantity
                         + " target, which the engine check would have tested");
  }
  return warnings;
}

/// B8 engine mode (Adrian, 2026-10-21): every target the user wrote that the
/// engine check doesn't enforce gets a warning at stage start, in the output
/// and the record.
[[nodiscard]] inline std::vector<std::string>
unenforced_target_warnings(const ConvergenceSettings& settings) {
  std::vector<std::string> warnings;
  if (settings.decider != Decider::engine_check)
    return warnings;
  for (const auto& metric : settings.metrics) {
    if (metric.group != Group::requested_criteria
        || metric.set_by != SetBy::user || !metric.target)
      continue;
    std::ostringstream target;
    target << std::setprecision(1) << std::scientific << *metric.target;
    warnings.push_back(metric.keyword + " " + target.str()
                       + " is reported, not enforced: the engine check "
                       + metric.engine_check_instead + ". Set "
                       + use_default_criteria_key
                       + " false to enforce your targets.");
  }
  return warnings;
}

enum class ReuseAction : std::uint8_t { load, resume };

/// What --continue does with a saved state (B6, B8; Adrian, 2026-10-21).
struct ReuseDecision {
  ReuseAction action = ReuseAction::resume;
  std::string why;
  /// Requested targets the saved state doesn't meet but that the deciding
  /// check doesn't enforce; reported, they don't force a resume.
  std::vector<std::string> unmet_not_enforceable;
};

/// Decide reuse from saved evidence, for the current request:
///   - criteria decide: load only if every requested target is met by the
///     saved measurements (`requested_status`, from
///     evaluate_explicit_criteria);
///   - engine decides: load if the saved state passes the engine check for the
///     current request; requested targets it doesn't meet are reported as
///     "unmet, not enforceable in engine mode", because resuming would stop on
///     the engine check again at once.
/// Where a resumed run starts (which rung) is the D7 planner's job.
[[nodiscard]] inline ReuseDecision
decide_reuse(const ConvergenceSettings& settings,
             bool saved_engine_check_passed,
             const nlohmann::json& requested_status) {
  ReuseDecision decision;
  std::vector<std::string> unmet;
  // Held in a local: iterating items() of the temporary value() returns
  // would dangle.
  const auto criteria
      = requested_status.value("criteria", nlohmann::json::object());
  for (const auto& criterion : criteria.items())
    if (criterion.value().value("status", "") != "satisfied")
      unmet.push_back(criterion.key());
  if (settings.decider == Decider::requested_criteria) {
    decision.action = unmet.empty() ? ReuseAction::load : ReuseAction::resume;
    decision.why = unmet.empty()
                       ? "saved state meets every requested target"
                       : "saved state does not meet every requested target";
    return decision;
  }
  decision.action
      = saved_engine_check_passed ? ReuseAction::load : ReuseAction::resume;
  decision.why = saved_engine_check_passed
                     ? "saved state passes the engine check for this request"
                     : "saved state does not pass the engine check for this "
                       "request";
  if (saved_engine_check_passed)
    decision.unmet_not_enforceable = unmet;
  return decision;
}

[[nodiscard]] inline const char* to_string(Decider decider) {
  return decider == Decider::engine_check ? "engine_check"
                                          : "requested_criteria";
}
[[nodiscard]] inline const char* to_string(Group group) {
  return group == Group::engine_check ? "engine_check" : "requested_criteria";
}
[[nodiscard]] inline const char* to_string(SetBy set_by) {
  switch (set_by) {
  case SetBy::engine_default:
    return "engine_default";
  case SetBy::user:
    return "user";
  case SetBy::not_set:
    break;
  }
  return "not_set";
}
[[nodiscard]] inline const char* to_string(Relation relation) {
  return relation == Relation::less_than ? "less_than" : "less_than_or_equal";
}

namespace detail {
inline std::string scientific(std::optional<double> value) {
  if (!value)
    return "-";
  std::ostringstream text;
  text << std::setprecision(1) << std::scientific << *value;
  return text.str();
}
/// A record spelling for print: "engine_default" -> "engine default".
inline std::string words(std::string text) {
  std::replace(text.begin(), text.end(), '_', ' ');
  return text;
}
inline std::string caps_text(const IterationCaps& caps) {
  std::string text;
  if (caps.per_rung.size() == 1) {
    text = std::to_string(caps.per_rung.front()) + " at every rung";
  } else {
    text = "[";
    for (std::size_t i = 0; i < caps.per_rung.size(); ++i)
      text += (i > 0 ? ", " : "") + std::to_string(caps.per_rung[i]);
    text += "] per rung";
  }
  return text;
}
/// Column widths of the printed tables.
inline constexpr std::size_t name_width = 32;
inline constexpr std::size_t header_name_width = name_width + 2; // row indent
inline constexpr std::size_t number_width = 8;
inline constexpr std::size_t set_by_width = 14;

inline std::string padded(std::string text, std::size_t width) {
  if (text.size() < width)
    text.append(width - text.size(), ' ');
  return text + "  ";
}
} // namespace detail

[[nodiscard]] inline nlohmann::json
to_json(const ConvergenceSettings& settings) {
  nlohmann::json metrics = nlohmann::json::array();
  for (const auto& metric : settings.metrics) {
    metrics.push_back(
        {{"group", to_string(metric.group)},
         {"name", metric.name},
         {"measurement", metric.measurement},
         {"quantity", metric.quantity.empty()
                          ? nlohmann::json(nullptr)
                          : nlohmann::json(metric.quantity)},
         {"units", metric.units},
         {"target", metric.target ? nlohmann::json(*metric.target)
                                  : nlohmann::json(nullptr)},
         {"set_by", to_string(metric.set_by)},
         {"enforced", enforced(settings, metric)},
         {"relation", to_string(metric.relation)},
         {"effective_limit", metric.effective_limit
                                 ? nlohmann::json(*metric.effective_limit)
                                 : nlohmann::json(nullptr)},
         {"limit_formula", metric.limit_formula},
         {"keyword", metric.keyword},
         {"note", metric.note}});
  }
  return {{"schema", "madqc.convergence_settings.v1"},
          {"stage", settings.stage},
          {"application", settings.application},
          {"decided_by", to_string(settings.decider)},
          {"decided_by_set_by",
           settings.decider_set_by_user ? "user" : "engine_default"},
          {"metrics", metrics},
          {"maxiter",
           {{"per_rung", settings.iteration_caps.per_rung},
            {"set_by", to_string(settings.iteration_caps.set_by)},
            {"counts", "every iteration at a rung, across --continue runs"}}},
          {"warnings", settings.warnings}};
}

/// The "Effective limit" cell of a settings row: the value compared, how it
/// follows from the target, and whether it is only reported.
[[nodiscard]] inline std::string limit_text(const ConvergenceSettings& settings,
                                            const MetricSetting& metric) {
  std::string limit;
  if (!comparison_limit(metric) && metric.limit_formula.empty()) {
    limit = "reported only";
  } else {
    if (comparison_limit(metric))
      limit = detail::scientific(comparison_limit(metric));
    if (!metric.limit_formula.empty())
      limit += (limit.empty() ? "" : " = ") + metric.limit_formula;
    if (!enforced(settings, metric))
      limit += ", reported only";
  }
  if (!metric.note.empty())
    limit += " (" + metric.note + ")";
  return limit;
}

/// The table printed at the start of a stage (B9). Built from the same struct
/// as to_json().
[[nodiscard]] inline std::string
render_settings_table(const ConvergenceSettings& settings) {
  using detail::padded;
  using detail::scientific;
  std::string table = "Convergence settings: stage " + settings.stage + " ("
                      + settings.application + ")\n";
  table += settings.decider == Decider::engine_check
               ? "  Decided by: engine check"
               : "  Decided by: your requested criteria";
  table += settings.decider_set_by_user ? " (set in the stage)\n"
                                        : " (default)\n";
  if (settings.decider == Decider::engine_check)
    table += std::string("  Set ") + use_default_criteria_key
             + " false to let your own targets decide.\n";
  table += "\n  " + padded("Metric", detail::header_name_width)
           + padded("Target", detail::number_width)
           + padded("Set by", detail::set_by_width)
           + padded("Enforced", detail::number_width)
           + padded("Relation", detail::number_width) + "Effective limit\n";
  for (const auto group : {Group::engine_check, Group::requested_criteria}) {
    table += group == Group::engine_check ? "  engine check\n"
                                          : "  requested criteria\n";
    for (const auto& metric : settings.metrics) {
      if (metric.group != group)
        continue;
      table += "    " + padded(metric.name, detail::name_width)
               + padded(scientific(metric.target), detail::number_width)
               + padded(detail::words(to_string(metric.set_by)),
                        detail::set_by_width)
               + padded(enforced(settings, metric) ? "yes" : "no",
                        detail::number_width)
               + padded(metric.relation == Relation::less_than ? "<" : "<=",
                        detail::number_width)
               + limit_text(settings, metric) + "\n";
    }
  }
  if (!settings.iteration_caps.per_rung.empty())
    table += "  Iteration cap (maxiter): "
             + detail::caps_text(settings.iteration_caps) + " ("
             + detail::words(to_string(settings.iteration_caps.set_by))
             + "), counted across --continue runs\n";
  for (const auto& warning : settings.warnings)
    table += "  WARNING: " + warning + "\n";
  return table;
}

/// One row's verdict at the end of a stage: met, unmet, unavailable, or
/// reported (no single number to compare with).
[[nodiscard]] inline std::string verdict(const MetricSetting& metric,
                                         const MetricOutcome& outcome) {
  const auto limit = comparison_limit(metric);
  if (!limit)
    return "reported";
  if (outcome.availability != "measured" || !outcome.value)
    return outcome.availability == "invalid" ? "invalid" : "unavailable";
  const bool met = metric.relation == Relation::less_than
                       ? *outcome.value < *limit
                       : *outcome.value <= *limit;
  return met ? "met" : "unmet";
}

namespace detail {
inline const MetricOutcome*
find_outcome(const std::vector<MetricOutcome>& outcomes,
             const std::string& measurement) {
  for (const auto& outcome : outcomes)
    if (outcome.measurement == measurement)
      return &outcome;
  return nullptr;
}
/// The outcome for a measurement, or an "unavailable" one when none was given.
inline MetricOutcome outcome_for(const std::vector<MetricOutcome>& outcomes,
                                 const std::string& measurement) {
  const auto* found = find_outcome(outcomes, measurement);
  return found != nullptr
             ? *found
             : MetricOutcome{measurement, std::nullopt, "unavailable"};
}
} // namespace detail

/// The enforced row furthest from its limit (largest value / limit); none when
/// nothing enforced was measured.
[[nodiscard]] inline std::optional<std::size_t>
limiting_metric(const ConvergenceSettings& settings,
                const std::vector<MetricOutcome>& outcomes) {
  std::optional<std::size_t> limiting;
  double worst = -1.0;
  for (std::size_t i = 0; i < settings.metrics.size(); ++i) {
    const auto& metric = settings.metrics[i];
    const auto limit = comparison_limit(metric);
    const auto outcome = detail::outcome_for(outcomes, metric.measurement);
    if (!enforced(settings, metric) || !limit || *limit <= 0.0
        || !outcome.value)
      continue;
    const double ratio = *outcome.value / *limit;
    if (ratio > worst) {
      worst = ratio;
      limiting = i;
    }
  }
  return limiting;
}

/// The end-of-stage record (B9): every row with its final value and verdict.
[[nodiscard]] inline nlohmann::json
outcome_to_json(const ConvergenceSettings& settings,
                const std::vector<MetricOutcome>& outcomes) {
  nlohmann::json rows = nlohmann::json::array();
  for (const auto& metric : settings.metrics) {
    const auto outcome = detail::outcome_for(outcomes, metric.measurement);
    rows.push_back({{"group", to_string(metric.group)},
                    {"name", metric.name},
                    {"measurement", metric.measurement},
                    {"value", outcome.value ? nlohmann::json(*outcome.value)
                                            : nlohmann::json(nullptr)},
                    {"availability", outcome.availability},
                    {"enforced", enforced(settings, metric)},
                    {"verdict", verdict(metric, outcome)}});
  }
  const auto limiting = limiting_metric(settings, outcomes);
  return {{"decided_by", to_string(settings.decider)},
          {"metrics", rows},
          {"limiting_metric",
           limiting ? nlohmann::json(settings.metrics[*limiting].name)
                    : nlohmann::json(nullptr)}};
}

/// The table printed at the end of a stage (B9), from the same struct.
[[nodiscard]] inline std::string
render_outcome_table(const ConvergenceSettings& settings,
                     const std::vector<MetricOutcome>& outcomes,
                     const std::string& stop_reason) {
  using detail::padded;
  using detail::scientific;
  std::string table = "Convergence outcome: stage " + settings.stage + " ("
                      + settings.application + "), stop reason " + stop_reason
                      + "\n";
  table += "  " + padded("Metric", detail::header_name_width)
           + padded("Value", detail::number_width)
           + padded("Limit", detail::number_width)
           + padded("Enforced", detail::number_width) + "Verdict\n";
  for (const auto& metric : settings.metrics) {
    const auto outcome = detail::outcome_for(outcomes, metric.measurement);
    table
        += "    " + padded(metric.name, detail::name_width)
           + padded(scientific(outcome.value), detail::number_width)
           + padded(scientific(comparison_limit(metric)), detail::number_width)
           + padded(enforced(settings, metric) ? "yes" : "no",
                    detail::number_width)
           + verdict(metric, outcome)
           + (enforced(settings, metric) || !comparison_limit(metric)
                  ? ""
                  : " (reported only)")
           + "\n";
  }
  const auto limiting = limiting_metric(settings, outcomes);
  table += "  Limiting metric: "
           + (limiting ? settings.metrics[*limiting].name : std::string("-"))
           + "\n";
  return table;
}

} // namespace madness::convergence
