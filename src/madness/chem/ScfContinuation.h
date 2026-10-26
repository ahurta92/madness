#pragma once

// Protocol rungs identified by value, and what a saved generation says about
// its rung (DESIGN B5, B6). The planner that uses them is LadderPlanner.h.
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

namespace detail {
inline std::string describe(std::size_t index, const Rung& rung) {
  std::ostringstream text;
  text << "rung " << index + 1 << " (thresh " << std::setprecision(0)
       << std::scientific << rung.thresh << ", k " << rung.k << ")";
  return text.str();
}
} // namespace detail

} // namespace madness::scf_continuation
