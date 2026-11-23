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

#include <madness/chem/ConvergenceDecision.h>
#include <madness/chem/ExplicitCriteria.h>

#include <nlohmann/json.hpp>
#include <set>

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

/// B9 (D-107): an application's definition of one written-target row; the
/// rule that turns it into a settings row is requested_metric's.
struct RequestedMetric {
  const char *name, *keyword, *measurement, *quantity, *units;
  std::string engine_instead;
  std::string unset_note; ///< shown when no target was written; may be empty
};

/// B9 (D-107): a written-target row: the target when the deck wrote one (set
/// by the user), else none and the row's unset note; compared with <=.
[[nodiscard]] inline MetricSetting
requested_metric(const RequestedMetric& row, const nlohmann::json& targets) {
  const bool written = targets.contains(row.measurement);
  return {Group::requested_criteria,
          row.name,
          row.keyword,
          row.measurement,
          row.quantity,
          row.units,
          written ? std::optional<double>(
                        targets.at(row.measurement).at("value").get<double>())
                  : std::nullopt,
          written ? SetBy::user : SetBy::not_set,
          Relation::less_than_or_equal,
          std::nullopt,
          "",
          written ? std::string() : row.unset_note,
          row.engine_instead};
}

/// The iteration cap per protocol rung (B5): maxiter, one number or one per
/// rung, per run: a restart starts the count again (D-007).
struct IterationCaps {
  std::vector<int> per_rung; ///< one entry: the same cap at every rung
  SetBy set_by = SetBy::engine_default;
};

/// When a stage publishes checkpoints and what it keeps (A3). By default only
/// at decision points (every rung end and every stop), with every older
/// generation deleted after each publish except those another stage
/// references (D-070).
struct CheckpointCadence {
  int every_iterations
      = 0; ///< K within a rung (opt-in); 0: decision points only
  SetBy every_iterations_set_by = SetBy::engine_default;
  bool keep_all = false; ///< true: retention off (debugging)
  SetBy keep_all_set_by = SetBy::engine_default;
  /// Non-empty when the application cannot change when it publishes yet:
  /// what it does instead (e.g. molresponse publishes at solve ends only).
  /// What it keeps is keep_all's either way.
  std::string fixed;
};

struct ConvergenceSettings {
  std::string stage;
  std::string application;
  Decider decider = Decider::engine_check;
  bool decider_set_by_user = false;
  std::vector<MetricSetting> metrics;
  IterationCaps iteration_caps;
  CheckpointCadence checkpoints;
  /// A stage's starting-guess source, as written (e.g. a Dalton seed, S5);
  /// empty: none. Shown and recorded only when present.
  std::string seed;
  std::vector<std::string> warnings;
  /// D-095: which protocol step the limits are, e.g. "limits of the final
  /// protocol step (thresh 1e-06, k 8); the run starts at step 1 (thresh
  /// 1e-04)". Shown under both tables and in the outcome record; not in
  /// to_json, so the settings record and the manifests are unchanged.
  std::string scope;
};

/// A final measurement of one metric, for the end-of-stage table (B9).
struct MetricOutcome {
  std::string measurement; ///< matches MetricSetting::measurement
  std::optional<double> value;
  std::string availability
      = "unavailable"; ///< measured | unavailable | invalid
  /// D-095: the limit this value was judged against when it is its own (a
  /// per-leg effective limit); absent: the row's.
  std::optional<double> limit;
  std::string units; ///< the measurement's units; empty: the row's
  std::string phase; ///< when in the solve the value was taken
  /// What the value is about when a row has one per leg or root (e.g.
  /// "dipole_z f0.02000"); empty for a stage-wide value.
  std::string subject;
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
  bool final_protocol = true; ///< false: consume the full intermediate allowance
};

/// Intermediate protocols consume their full cap before reporting the
/// selected convergence verdict. Final protocols retain the B8 decider,
/// followed by the cap and optional stall exit.
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
  if (!sample.final_protocol) {
    if (sample.iterations_at_rung < sample.cap)
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
    // D-097: one warning per target; the instruction to let the targets
    // decide is the table's one line above them.
    warnings.push_back(metric.keyword + " " + target.str()
                       + " is reported, not enforced: the engine check "
                       + metric.engine_check_instead);
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
checkpoints_to_json(const CheckpointCadence& checkpoints) {
  if (!checkpoints.fixed.empty())
    return {{"fixed", checkpoints.fixed},
            {"keep_all", checkpoints.keep_all},
            {"keep_all_set_by", to_string(checkpoints.keep_all_set_by)}};
  return {
      {"every_iterations", checkpoints.every_iterations > 0
                               ? nlohmann::json(checkpoints.every_iterations)
                               : nlohmann::json(nullptr)},
      {"every_iterations_set_by",
       to_string(checkpoints.every_iterations_set_by)},
      {"at_decision_points", true},
      {"keep_all", checkpoints.keep_all},
      {"keep_all_set_by", to_string(checkpoints.keep_all_set_by)}};
}

/// The "Checkpoints" line of the settings table (A3, B9).
[[nodiscard]] inline std::string
checkpoints_text(const CheckpointCadence& checkpoints) {
  const std::string when
      = !checkpoints.fixed.empty() ? checkpoints.fixed
        : checkpoints.every_iterations > 0
            ? "every " + std::to_string(checkpoints.every_iterations)
                  + " iterations within a rung ("
                  + detail::words(
                      to_string(checkpoints.every_iterations_set_by))
                  + "), plus every rung end and stop"
            : "at every rung end and stop ("
                  + detail::words(
                      to_string(checkpoints.every_iterations_set_by))
                  + "; checkpoint.every_iterations K adds one every K "
                    "iterations)";
  const std::string kept
      = checkpoints.keep_all
            ? "every generation kept (checkpoint.keep_all, "
                  + detail::words(to_string(checkpoints.keep_all_set_by)) + ")"
            : "older generations deleted after each publish, except any "
              "another stage references ("
                  + detail::words(to_string(checkpoints.keep_all_set_by)) + ")";
  return when + "; " + kept;
}

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
  nlohmann::json json
      = {{"schema", "madqc.convergence_settings.v1"},
         {"stage", settings.stage},
         {"application", settings.application},
         {"decided_by", to_string(settings.decider)},
         {"decided_by_set_by",
          settings.decider_set_by_user ? "user" : "engine_default"},
         {"metrics", metrics},
         {"maxiter",
          {{"per_rung", settings.iteration_caps.per_rung},
           {"set_by", to_string(settings.iteration_caps.set_by)},
           {"counts", "iterations at a rung in this run"}}},
         {"checkpoints", checkpoints_to_json(settings.checkpoints)},
         {"warnings", settings.warnings}};
  if (!settings.seed.empty())
    json["seed"] = settings.seed;
  return json;
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
  if (!settings.scope.empty())
    table += "  Scope: " + settings.scope + "\n";
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
             + "), per run; a restart starts the count again\n";
  table += "  Checkpoints: " + checkpoints_text(settings.checkpoints) + "\n";
  if (!settings.seed.empty())
    table += "  Seed: " + settings.seed + "\n";
  for (const auto& warning : settings.warnings)
    table += "  WARNING: " + warning + "\n";
  return table;
}

/// D-095: the scope line of a stage's tables: which protocol step the
/// limits belong to, and where the run starts.
[[nodiscard]] inline std::string
protocol_scope(const std::vector<double>& protocol, int final_k) {
  if (protocol.empty())
    return {};
  const auto step = "(thresh " + detail::scientific(protocol.back()) + ", k "
                    + std::to_string(final_k) + ")";
  if (protocol.size() == 1)
    return "limits of the one protocol step " + step;
  return "limits of the final protocol step " + step
         + "; the run starts at step 1 (thresh "
         + detail::scientific(protocol.front()) + ")";
}

/// One row's verdict at the end of a stage: met, unmet, unavailable, not
/// applicable, or reported (no single number to compare with). An outcome
/// that carries its own limit (D-095, a per-leg limit) is judged against it.
[[nodiscard]] inline std::string verdict(const MetricSetting& metric,
                                         const MetricOutcome& outcome) {
  const auto limit = outcome.limit ? outcome.limit : comparison_limit(metric);
  if (!limit)
    return metric.note.starts_with("not applicable") ? "not applicable"
                                                       : "reported";
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
/// D-095: a row's per-leg (or per-root) outcomes, in the order given; empty
/// when the row has one stage-wide outcome or none.
inline std::vector<MetricOutcome>
leg_outcomes(const std::vector<MetricOutcome>& outcomes,
             const std::string& measurement) {
  std::vector<MetricOutcome> legs;
  for (const auto& outcome : outcomes)
    if (outcome.measurement == measurement && !outcome.subject.empty())
      legs.push_back(outcome);
  return legs;
}
} // namespace detail

/// D-095: a row's verdict over its legs when it has them (unmet when any
/// leg is, else the first leg that is not met, else met), else verdict().
[[nodiscard]] inline std::string
row_verdict(const MetricSetting& metric,
            const std::vector<MetricOutcome>& outcomes) {
  const auto legs = detail::leg_outcomes(outcomes, metric.measurement);
  if (legs.empty())
    return verdict(metric, detail::outcome_for(outcomes, metric.measurement));
  std::string result = "met";
  for (const auto& leg : legs) {
    const auto v = verdict(metric, leg);
    if (v == "unmet")
      return v;
    if (v != "met" && result == "met")
      result = v;
  }
  return result;
}

/// The enforced row furthest from its limit (largest value / limit); none when
/// nothing enforced was measured.
[[nodiscard]] inline std::optional<std::size_t>
limiting_metric(const ConvergenceSettings& settings,
                const std::vector<MetricOutcome>& outcomes) {
  std::optional<std::size_t> limiting;
  double worst = -1.0;
  for (std::size_t i = 0; i < settings.metrics.size(); ++i) {
    const auto& metric = settings.metrics[i];
    if (!enforced(settings, metric))
      continue;
    auto candidates = detail::leg_outcomes(outcomes, metric.measurement);
    if (candidates.empty())
      candidates.push_back(detail::outcome_for(outcomes, metric.measurement));
    for (const auto& outcome : candidates) {
      const auto limit
          = outcome.limit ? outcome.limit : comparison_limit(metric);
      if (!limit || *limit <= 0.0 || !outcome.value)
        continue;
      const double ratio = *outcome.value / *limit;
      if (ratio > worst) {
        worst = ratio;
        limiting = i;
      }
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

/// D-095: the end-of-stage record in results.json: outcome_to_json's rows,
/// each with its limit, units and phase, and the legs (or roots) of a row
/// judged per leg; the scope of the limits. The manifests keep
/// outcome_to_json.
[[nodiscard]] inline nlohmann::json
outcome_record(const ConvergenceSettings& settings,
               const std::vector<MetricOutcome>& outcomes) {
  auto record = outcome_to_json(settings, outcomes);
  auto& rows = record["metrics"];
  for (std::size_t i = 0; i < settings.metrics.size(); ++i) {
    const auto& metric = settings.metrics[i];
    auto& row = rows[i];
    const auto outcome = detail::outcome_for(outcomes, metric.measurement);
    const auto limit = comparison_limit(metric);
    row["limit"] = limit ? nlohmann::json(*limit) : nlohmann::json(nullptr);
    row["units"] = outcome.units.empty() ? metric.units : outcome.units;
    row["phase"] = outcome.phase.empty() ? nlohmann::json(nullptr)
                                         : nlohmann::json(outcome.phase);
    row["verdict"] = row_verdict(metric, outcomes);
    const auto legs = detail::leg_outcomes(outcomes, metric.measurement);
    if (legs.empty())
      continue;
    nlohmann::json each = nlohmann::json::array();
    std::optional<double> largest;
    bool measured = true;
    for (const auto& leg : legs) {
      each.push_back({{"subject", leg.subject},
                      {"value", leg.value ? nlohmann::json(*leg.value)
                                          : nlohmann::json(nullptr)},
                      {"limit", leg.limit ? nlohmann::json(*leg.limit)
                                          : nlohmann::json(nullptr)},
                      {"availability", leg.availability},
                      {"verdict", verdict(metric, leg)},
                      {"phase", leg.phase}});
      measured = measured && leg.availability == "measured";
      if (leg.value && (!largest || *leg.value > *largest))
        largest = leg.value;
    }
    row["legs"] = std::move(each);
    row["value"] = largest ? nlohmann::json(*largest) : nlohmann::json(nullptr);
    row["availability"] = measured ? "measured" : "unavailable";
    row["phase"] = legs.front().phase;
  }
  record["scope"] = settings.scope;
  return record;
}

/// The table printed at the end of a stage (B9), from the same struct: the
/// engine-check rows, then the requested ones, each group with the phase of
/// its values; a row judged per leg is followed by one line per leg (D-095).
[[nodiscard]] inline std::string
render_outcome_table(const ConvergenceSettings& settings,
                     const std::vector<MetricOutcome>& outcomes,
                     const std::string& stop_reason) {
  using detail::padded;
  using detail::scientific;
  std::string table = "Convergence outcome: stage " + settings.stage + " ("
                      + settings.application + "), stop reason " + stop_reason
                      + "\n";
  if (!settings.scope.empty())
    table += "  Scope: " + settings.scope + "\n";
  table += "  " + padded("Metric", detail::header_name_width)
           + padded("Value", detail::number_width)
           + padded("Limit", detail::number_width)
           + padded("Enforced", detail::number_width) + "Verdict\n";
  for (const auto group : {Group::engine_check, Group::requested_criteria}) {
    std::string rows;
    std::string phase;
    for (const auto& metric : settings.metrics) {
      if (metric.group != group)
        continue;
      const auto outcome = detail::outcome_for(outcomes, metric.measurement);
      const auto legs = detail::leg_outcomes(outcomes, metric.measurement);
      if (phase.empty())
        phase = outcome.phase;
      const auto units = outcome.units.empty() ? metric.units : outcome.units;
      const auto row = row_verdict(metric, outcomes);
      rows += "    " + padded(metric.name, detail::name_width)
              + padded(legs.empty() ? scientific(outcome.value) : "per leg",
                       detail::number_width)
              + padded(legs.empty() ? scientific(comparison_limit(metric))
                                    : "per leg",
                       detail::number_width)
              + padded(enforced(settings, metric) ? "yes" : "no",
                       detail::number_width)
              + row
              + (enforced(settings, metric) || row == "reported"
                         || row == "not applicable"
                     ? ""
                     : " (reported only)")
              + (units.empty() || row == "not applicable" ? ""
                                                          : "  [" + units + "]")
              + "\n";
      for (const auto& leg : legs)
        rows += "      " + padded(leg.subject, detail::name_width)
                + padded(scientific(leg.value), detail::number_width)
                + padded(scientific(leg.limit), detail::number_width)
                + padded("", detail::number_width) + verdict(metric, leg)
                + "\n";
    }
    if (rows.empty())
      continue;
    table += std::string(group == Group::engine_check ? "  Engine check"
                                                      : "  Requested targets")
             + (phase.empty() ? "" : " (values: " + detail::words(phase) + ")")
             + "\n" + rows;
  }
  const auto limiting = limiting_metric(settings, outcomes);
  table += "  Limiting metric: "
           + (limiting ? settings.metrics[*limiting].name : std::string("-"))
           + "\n";
  return table;
}

} // namespace madness::convergence

// The shared convergence decision (resolve, limits, decide) lives in
// ConvergenceDecision.h. Input keys it retired fail with a migration message.
namespace madness::iteration {
inline std::optional<std::string> retired_convergence_input(const std::string& key) {
    if (key == "convergence.criteria")
        return "convergence.criteria has been retired: the criteria are now the application's defaults "
               "plus every metric with a written target (convergence.<metric>)";
    if (key.starts_with("convergence.") && key.ends_with(".relative"))
        return key + " has been retired: the limit rule belongs to the application; "
               "write only the target, convergence.<metric>";
    return std::nullopt;
}
} // namespace madness::iteration
