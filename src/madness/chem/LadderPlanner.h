#pragma once

// The planner (DESIGN B6, decision D7): what --continue does with one saved
// state of a stage that climbs a supplied protocol ladder (B5).
//
// One shared function, plan_ladder, decides fresh, seed, load or resume for
// one state, from that state's saved record and the current request. moldft
// calls it once for its orbitals; molresponse once per (perturbation,
// frequency) state. Saved values are re-judged against today's limits, never
// trusted as a saved verdict (B6); the comparisons happen here, on saved
// records, in pure functions (B5: engines compare live state, everything else
// takes verdicts or judges saved records).
//
// Rules (B6, with the simplified-policy final-resume amendment): a finished
// rung is never redone; rungs are matched
// by value; R1 a changed reference or problem starts fresh (or seeded when the
// deck says seed); R2 an unmet final rung resumes with a fresh allowance;
// R3 a saved rung tighter than every requested rung is judged
// at the final one; R4 a new response state is seeded from the nearest saved
// frequency, R7 the lower on a tie; R5 removed states are not planned; R6 a
// diverged state is never resumed (seeded from its converged coarser rung,
// else fresh). maxiter is per run (B5, D-007): iterations already spent are
// recorded, never compared.
//
// Pure: no World, MPI, archives or solver state (B1).

#include <madness/chem/ConvergenceSettings.h>
#include <madness/chem/ExplicitCriteria.h>
#include <madness/chem/ScfContinuation.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace madness::planning {

using Rung = scf_continuation::Rung;
using StopReason = convergence::StopReason;

/// One row of the engine check, judged today on a saved value (B6). The limit
/// is absolute + relative * norm, where the norm is the scale the saved
/// measurement records (`normalization_scale`); relative = 0 for SCF.
struct EngineLimit {
  std::string measurement; ///< key in SavedState::engine_measurements
  double absolute = 0.0;
  double relative = 0.0;
  bool inclusive = false; ///< <= (response) instead of < (SCF)
};

/// What one saved state's record says.
struct SavedState {
  std::string identity;             ///< moldft: the physical block; response:
                                    ///< perturbation|frequency
  std::string reference_generation; ///< empty for moldft
  std::string perturbation;         ///< response only (R4)
  double frequency = 0.0;           ///< response only (R4, R7)
  Rung rung;                        ///< the rung it was produced at
  int iterations_at_rung = 0; ///< all runs; recorded, never compared (D-007)
  StopReason stop = StopReason::in_progress;
  bool diverged = false; ///< an invalid state (B5); R6
  /// The measurements of the sample the engine check decided on, raw.
  nlohmann::json engine_measurements = nlohmann::json::object();
  /// The measurements on the saved state itself, raw (the written targets
  /// are judged on these).
  nlohmann::json saved_measurements = nlohmann::json::object();
};

/// The current request for one state.
struct LadderRequest {
  std::string identity;
  std::string reference_generation;
  std::vector<Rung> ladder;
  convergence::ConvergenceSettings
      settings; ///< mode and written targets (B8, B9)
  std::vector<EngineLimit>
      engine_limits;              ///< today's engine check at the final rung
  bool forced_iterations = false; ///< diagnostic.iterations_per_rung > 0
};

enum class Action : std::uint8_t { fresh, seed, load, resume };

struct LadderPlan {
  Action action = Action::fresh;
  std::size_t start_rung = 0; ///< index into the requested ladder
  /// Spent at start_rung in earlier runs: recorded only (D-007).
  int iterations_done = 0;
  std::string seed_from;          ///< the seed, for Action::seed
  std::string why;                ///< one line, recorded (B6)
  std::vector<std::string> unmet; ///< e.g. a stalled final rung (R2)
  std::vector<std::string> unmet_not_enforceable; ///< engine mode (B8)
};

[[nodiscard]] inline const char* to_string(Action action) {
  switch (action) {
  case Action::seed:
    return "seed";
  case Action::load:
    return "load";
  case Action::resume:
    return "resume";
  case Action::fresh:
    break;
  }
  return "fresh";
}

/// B5: a rung's maxiter in this run, in one place: the forced count when
/// diagnostic.iterations_per_rung is set, else one entry of a per-rung list
/// (one entry: every rung), else maxiter.
[[nodiscard]] inline int maxiter_at_rung(const std::vector<int>& per_rung,
                                         int maxiter, int forced_iterations,
                                         std::size_t index) {
  if (forced_iterations > 0)
    return forced_iterations;
  if (per_rung.empty())
    return maxiter;
  return per_rung.size() == 1 ? per_rung.front() : per_rung.at(index);
}

/// The engine check, judged today on the saved deciding sample (B6). None
/// when it cannot be judged: a value or a needed norm is missing (D-005.4).
[[nodiscard]] inline std::optional<bool>
engine_check_today(const SavedState& saved,
                   const std::vector<EngineLimit>& limits) {
  if (limits.empty())
    return std::nullopt;
  for (const auto& limit : limits) {
    const auto& record = saved.engine_measurements.value(
        limit.measurement, nlohmann::json::object());
    const auto value = record.value("value", nlohmann::json());
    if (!value.is_number())
      return std::nullopt;
    double bound = limit.absolute;
    if (limit.relative != 0.0) {
      const auto norm = record.value("normalization_scale", nlohmann::json());
      if (!norm.is_number())
        return std::nullopt;
      bound += limit.relative * norm.get<double>();
    }
    const double saved_value = value.get<double>();
    if (!std::isfinite(saved_value)
        || !(limit.inclusive ? saved_value <= bound : saved_value < bound))
      return false;
  }
  return true;
}

namespace detail {
inline std::string rung_text(std::size_t index, const Rung& rung) {
  return scf_continuation::detail::describe(index, rung);
}
// "thresh 1e-04, k 6": a rung by value, without a position.
inline std::string rung_values(const Rung& rung) {
  std::ostringstream text;
  text << "thresh " << std::setprecision(0) << std::scientific << rung.thresh
       << ", k " << rung.k;
  return text.str();
}
// The first requested rung at least as tight as `rung` (a finished rung: the
// first strictly tighter); the last rung when there is none.
inline std::size_t start_on_ladder(const std::vector<Rung>& ladder,
                                   const Rung& rung, bool finished) {
  for (std::size_t i = 0; i < ladder.size(); ++i) {
    const bool usable = finished
                            ? scf_continuation::tighter(ladder[i], rung)
                            : (scf_continuation::same_rung(ladder[i], rung)
                               || scf_continuation::tighter(ladder[i], rung));
    if (usable)
      return i;
  }
  return ladder.size() - 1;
}
} // namespace detail

/// Plan one state. `own`: its saved record, if any. `seed`: a saved state
/// allowed as a starting guess (the deck's seed for R1, the nearest frequency
/// for R4/R7, the converged coarser rung for R6), if any.
[[nodiscard]] inline LadderPlan
plan_ladder(const std::optional<SavedState>& own,
            const std::optional<SavedState>& seed,
            const LadderRequest& request) {
  LadderPlan plan;
  if (request.ladder.empty()) {
    plan.why = "empty protocol";
    return plan;
  }
  const auto seeded = [&](std::size_t start, std::string why) {
    plan.start_rung = start;
    if (seed) {
      plan.action = Action::seed;
      plan.seed_from
          = seed->identity + " at " + detail::rung_values(seed->rung);
      plan.why = std::move(why) + "; seeded from " + plan.seed_from;
    } else {
      plan.action = Action::fresh;
      plan.why = std::move(why) + "; no seed, starting fresh";
    }
    return plan;
  };

  // R1: no record of this state, or a different problem or reference.
  if (!own)
    return seeded(0, "no saved state");
  if (own->identity != request.identity)
    return seeded(0, "the problem changed, so the state starts fresh");
  if (own->reference_generation != request.reference_generation)
    return seeded(0, "the reference changed, so the state starts fresh");

  const std::size_t last = request.ladder.size() - 1;
  const auto& final_rung = request.ladder[last];

  // R6: a diverged state is never resumed; its rung restarts from the seed.
  if (own->diverged)
    return seeded(detail::start_on_ladder(request.ladder, own->rung, false),
                  "the saved state diverged at "
                      + detail::rung_values(own->rung)
                      + "; that rung starts again from its seed");

  // B6, B8: the deciding check, judged today on the saved values.
  const auto engine = engine_check_today(*own, request.engine_limits);
  const auto requested = evaluate_explicit_criteria(
      own->saved_measurements,
      convergence::requested_targets(request.settings));
  const auto reuse = convergence::decide_reuse(
      request.settings, engine.value_or(false), requested);
  const bool engine_unjudgeable
      = request.settings.decider == convergence::Decider::engine_check
        && !engine;
  const bool meets = reuse.action == convergence::ReuseAction::load
                     && !request.forced_iterations;

  const bool finished = own->stop != StopReason::in_progress;
  const bool same_final = scf_continuation::same_rung(own->rung, final_rung);
  const bool at_or_beyond_final
      = same_final || own->rung.thresh < final_rung.thresh;

  if (finished && at_or_beyond_final) {
    plan.start_rung = last;
    if (meets) {
      plan.action = Action::load;
      plan.why = "saved state finished " + detail::rung_text(last, final_rung)
                 + " and meets the request (" + reuse.why + ")";
      plan.unmet_not_enforceable = reuse.unmet_not_enforceable;
      return plan;
    }
    // An earlier stall is evidence about that attempt, not a prohibition
    // on a fresh allowance in this run.
    plan.action = Action::resume;
    plan.iterations_done = same_final ? own->iterations_at_rung : 0;
    plan.why = engine_unjudgeable
                   ? "the saved values cannot be re-judged (a value or a "
                     "norm is missing); resuming the final rung"
                   : "saved state does not meet the request at the final rung";
    return plan;
  }

  // Otherwise continue on the ladder, matched by value; a finished rung is
  // never redone.
  plan.action = Action::resume;
  plan.start_rung
      = detail::start_on_ladder(request.ladder, own->rung, finished);
  plan.iterations_done
      = scf_continuation::same_rung(request.ladder[plan.start_rung], own->rung)
            ? own->iterations_at_rung
            : 0;
  plan.why = finished ? std::string("saved state finished its rung (")
                            + convergence::to_string(own->stop)
                            + "); a finished rung is never redone"
                      : std::string("saved state stopped inside its rung");
  return plan;
}

/// R4, R7: among `saved`, the state of the same perturbation and reference
/// nearest in frequency; on a tie the lower frequency. Diverged states are
/// never seeds.
[[nodiscard]] inline std::optional<SavedState>
nearest_seed(const std::vector<SavedState>& saved,
             const std::string& perturbation,
             const std::string& reference_generation, double frequency) {
  std::optional<SavedState> best;
  for (const auto& state : saved) {
    if (state.perturbation != perturbation
        || state.reference_generation != reference_generation || state.diverged)
      continue;
    const double distance = std::abs(state.frequency - frequency);
    if (!best || distance < std::abs(best->frequency - frequency)
        || (distance == std::abs(best->frequency - frequency)
            && state.frequency < best->frequency))
      best = state;
  }
  return best;
}

/// R6: the finest converged rung of the same state coarser than `rung`.
[[nodiscard]] inline std::optional<SavedState>
coarser_seed(const std::vector<SavedState>& same_state, const Rung& rung) {
  std::optional<SavedState> best;
  for (const auto& state : same_state)
    if (state.stop == StopReason::converged && !state.diverged
        && scf_continuation::tighter(rung, state.rung)
        && (!best || state.rung.thresh < best->rung.thresh))
      best = state;
  return best;
}

/// One-line record of what happened, e.g. "resumed at rung 2 (thresh 1e-06,
/// k 8), iteration 3 of that rung".
[[nodiscard]] inline std::string describe(const LadderPlan& plan,
                                          const std::vector<Rung>& ladder) {
  if (ladder.empty())
    return plan.why;
  const auto rung = detail::rung_text(plan.start_rung, ladder[plan.start_rung]);
  switch (plan.action) {
  case Action::load:
    return "loaded at " + rung
           + (plan.unmet.empty() ? "" : "; unmet: " + plan.unmet.front());
  case Action::resume:
    return "resumed at " + rung + ", iteration "
           + std::to_string(plan.iterations_done + 1) + " of that rung";
  case Action::seed:
    return "seeded from " + plan.seed_from + ", starting at " + rung;
  case Action::fresh:
    break;
  }
  return "started fresh at " + rung;
}

} // namespace madness::planning

namespace madness::iteration {

enum class SavedEvidenceUse : std::uint8_t { compatible, guess_only, invalid };
struct SavedTask {
    IterationMetrics metrics;
    bool protocol_complete = false;
    std::size_t updates_at_protocol = 0;
    // The application classifies archive integrity, numerical validity and
    // scientific compatibility before calling. guess_only permits loading
    // functions but discards metric/comparison/completed-stage claims.
    SavedEvidenceUse use = SavedEvidenceUse::compatible;
    bool conversion_supported = false;
    bool bundle_expanded = false;
};
enum class RestartAction : std::uint8_t { fresh, guess, resume, reuse };
struct RestartPlan {
    RestartAction action = RestartAction::fresh;
    std::size_t start_protocol = 0;
    std::size_t main_updates = 0; // deliberately never initialized from history
    std::size_t max_updates = 0;
    std::size_t updates_at_protocol = 0;
    bool reset_comparisons = true; // scalar evidence does not restore history
    std::optional<Decision> completion;
    std::string reason;
};

// One task or coupled bundle per call. Seeds for other tasks never confer
// progress. Restart always reconstructs raw equation work before advancement;
// neither this planner nor saved scalar measurements restores KAIN.
inline RestartPlan plan_restart(const std::vector<Protocol>& protocols,
                                const Request& request, const ModelMetrics& model,
                                const ExpectedTask& expected,
                                const std::optional<SavedTask>& saved) {
    if (protocols.empty()) throw std::invalid_argument("empty protocol sequence");
    for (auto position = protocols.begin(); position != protocols.end(); ++position) {
        const auto& protocol = *position;
        if (!std::isfinite(protocol.thresh) || protocol.thresh <= 0.0 || protocol.k <= 0)
            throw std::invalid_argument("invalid protocol representation");
        if (std::find(protocols.begin(), position, protocol) != position)
            throw std::invalid_argument("duplicate protocol representation");
    }
    const auto final_settings = resolve(request, model, protocols.back(), true);
    RestartPlan plan;
    plan.max_updates = request.max_updates;
    plan.reason = "new_task";
    if (!saved) return plan;
    if (saved->use == SavedEvidenceUse::invalid)
        throw std::invalid_argument("invalid selected checkpoint");

    const auto& metrics = saved->metrics;
    ExpectedTask stored_expected = expected;
    if (saved->use == SavedEvidenceUse::guess_only)
        stored_expected.context_id = metrics.info.context_id;
    if (saved->bundle_expanded) {
        const std::set<std::string> wanted(expected.states.begin(), expected.states.end());
        if (wanted.size() != expected.states.size() || wanted.contains("")
            || metrics.states.empty() || metrics.states.size() >= wanted.size())
            throw std::invalid_argument("invalid expanded bundle membership");
        stored_expected.states.clear();
        for (const auto& subject : metrics.states) {
            if (!wanted.contains(subject.id)) throw std::invalid_argument("saved subject absent from expanded bundle");
            stored_expected.states.push_back(subject.id);
        }
    }
    // Validate the saved frame independently of whether its metrics pass today's
    // limits. Final settings avoid applying a model floor for the wrong rung.
    auto saved_settings = final_settings;
    saved_settings.protocol = metrics.info.protocol;
    if (!std::isfinite(saved_settings.protocol.thresh) || saved_settings.protocol.thresh <= 0.0
        || saved_settings.protocol.k <= 0
        || !assess(saved_settings, stored_expected, metrics).valid)
        throw std::invalid_argument("invalid or mismatched saved metrics");

    const auto position = std::find(protocols.begin(), protocols.end(), metrics.info.protocol);
    if (position == protocols.end()) {
        if (!saved->conversion_supported)
            throw std::invalid_argument("checkpoint representation conversion is unsupported");
        plan.action = RestartAction::guess;
        plan.reason = "saved_protocol_not_requested";
        return plan;
    }
    plan.start_protocol = static_cast<std::size_t>(position - protocols.begin());
    plan.updates_at_protocol = saved->updates_at_protocol;
    plan.action = RestartAction::resume;
    plan.reason = "fresh_evaluation_required";
    if (saved->use == SavedEvidenceUse::guess_only || saved->bundle_expanded)
        return plan;
    if (plan.start_protocol + 1 == protocols.size()) {
        plan.completion = saved_completion(final_settings, expected, metrics);
        if (plan.completion) {
            plan.action = RestartAction::reuse;
            plan.reason = "saved_metrics_satisfy_request";
        }
    } else if (saved->protocol_complete) {
        ++plan.start_protocol;
        plan.updates_at_protocol = 0;
        plan.reason = "completed_intermediate_protocol";
    }
    return plan;
}

} // namespace madness::iteration
