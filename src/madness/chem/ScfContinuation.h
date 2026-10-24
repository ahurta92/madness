#pragma once

// Where a named moldft stage continues from on --continue (DESIGN B5, B6).
//
// This is the first piece of the D7 planner for SCF: plain data in, plain data
// out, no World, MPI, archives or solver state. The adapter reads the saved
// generation, calls plan_scf_continuation(), records the result, and passes
// the starting rung to the engine explicitly; the engine's own restart planner
// no longer chooses the rung on this path (Q6).
//
// Rules (Adrian, 2026-10-21):
//   - Rungs are matched by value (thresh, k), never by position, so a changed
//     protocol list still works.
//   - A finished rung is never redone: it has already given all the accuracy
//     it can. Continue at the first requested rung tighter than it.
//   - An unfinished rung continues at the first requested rung at least as
//     tight as it (normally the same rung).
//   - maxiter is per run (B5, D-007): a resumed rung gets maxiter iterations
//     again. The iterations already spent at the rung are recorded, not
//     compared against the cap.

#include <madness/chem/ConvergenceSettings.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace madness::scf_continuation {

/// One rung of a protocol ladder, identified by value.
struct Rung {
  double thresh = 0.0;
  int k = 0;
};

/// Thresholds come from decks and JSON; equal to this relative tolerance
/// means the same rung.
inline constexpr double same_thresh_tolerance = 1e-9;

[[nodiscard]] inline bool same_rung(const Rung& left, const Rung& right) {
  return left.k == right.k
         && std::abs(left.thresh - right.thresh)
                <= same_thresh_tolerance * std::abs(right.thresh);
}

/// `candidate` is tighter than `reference`: a smaller threshold (k follows
/// the threshold).
[[nodiscard]] inline bool tighter(const Rung& candidate,
                                  const Rung& reference) {
  return candidate.thresh < reference.thresh
         && !same_rung(candidate, reference);
}

/// How a rung ended (B5); defined beside rung_exit, which decides it.
using RungStop = convergence::StopReason;
using convergence::to_string;

[[nodiscard]] inline RungStop rung_stop_from_string(const std::string& text) {
  return convergence::stop_reason_from_string(text);
}

/// What the saved generation says about the rung it was produced at.
struct SavedRung {
  Rung rung;
  int iterations_at_rung = 0; ///< all runs, not only the last one
  RungStop stop = RungStop::in_progress;
};

enum class Action : std::uint8_t { load, resume };

struct Continuation {
  Action action = Action::resume;
  std::size_t start_rung = 0; ///< index into the requested ladder
  /// Already spent at start_rung in earlier runs: recorded, so the total at
  /// the rung stays known, never counted against maxiter (B5, D-007).
  int iterations_done = 0;
  std::string why;
};

namespace detail {
inline std::string describe(std::size_t index, const Rung& rung) {
  std::ostringstream text;
  text << "rung " << index + 1 << " (thresh " << std::setprecision(0)
       << std::scientific << rung.thresh << ", k " << rung.k << ")";
  return text.str();
}
} // namespace detail

/// Plan the continuation of a saved SCF state under the current request.
///
/// `saved_meets_request` is the reuse verdict of the deciding check from saved
/// evidence (ConvergenceSettings.h, decide_reuse): true when the saved state
/// may be loaded as far as the targets are concerned. This function adds the
/// rung logic.
[[nodiscard]] inline Continuation
plan_scf_continuation(const SavedRung& saved, const std::vector<Rung>& request,
                      bool saved_meets_request) {
  Continuation plan;
  if (request.empty()) {
    plan.why = "empty protocol";
    return plan;
  }
  const std::size_t last = request.size() - 1;
  const bool finished = saved.stop != RungStop::in_progress;
  const bool at_or_beyond_final = !tighter(request[last], saved.rung);

  // The saved state is at (or tighter than) the final requested rung and that
  // rung is finished: load it if it meets the request, otherwise resume there
  // with a fresh maxiter.
  if (finished && at_or_beyond_final) {
    const bool same_final = same_rung(saved.rung, request[last]);
    plan.start_rung = last;
    if (saved_meets_request) {
      plan.action = Action::load;
      plan.why = "saved state finished " + detail::describe(last, request[last])
                 + " and meets the request";
      return plan;
    }
    plan.action = Action::resume;
    plan.iterations_done = same_final ? saved.iterations_at_rung : 0;
    plan.why = "saved state does not meet the request at the final rung";
    return plan;
  }

  // Otherwise continue on the ladder, matched by value.
  std::size_t start = last;
  for (std::size_t i = 0; i < request.size(); ++i) {
    const auto& rung = request[i];
    const bool usable
        = finished ? tighter(rung, saved.rung)
                   : (same_rung(rung, saved.rung) || tighter(rung, saved.rung));
    if (usable) {
      start = i;
      break;
    }
  }
  plan.action = Action::resume;
  plan.start_rung = start;
  plan.iterations_done
      = same_rung(request[start], saved.rung) ? saved.iterations_at_rung : 0;
  plan.why = finished ? std::string("saved state finished its rung (")
                            + to_string(saved.stop)
                            + "); a finished rung is never redone"
                      : std::string("saved state stopped inside its rung");
  return plan;
}

/// One-line record of what happened, e.g. "resumed at rung 2 (thresh 1e-06,
/// k 8), iteration 3 of that rung".
[[nodiscard]] inline std::string describe(const Continuation& plan,
                                          const std::vector<Rung>& request) {
  if (request.empty())
    return plan.why;
  const auto rung = detail::describe(plan.start_rung, request[plan.start_rung]);
  if (plan.action == Action::load)
    return "loaded at " + rung;
  return "resumed at " + rung + ", iteration "
         + std::to_string(plan.iterations_done + 1) + " of that rung";
}

} // namespace madness::scf_continuation
