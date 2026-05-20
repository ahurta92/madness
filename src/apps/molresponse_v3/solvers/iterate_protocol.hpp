#ifndef MOLRESPONSE_V3_SOLVERS_ITERATE_PROTOCOL_HPP
#define MOLRESPONSE_V3_SOLVERS_ITERATE_PROTOCOL_HPP

// =========================================================================
// solvers::iterate_protocol — wraps solvers::iterate with a protocol
// ladder. At each protocol step:
//
//   1. prepare(thresh, solver, state) — caller-supplied callback that
//      brings ALL protocol-sensitive data to the new (k, thresh):
//         * FunctionDefaults<3> via set_response_protocol(...)
//         * ground-state operators via gs.prepare(...)
//         * solver target via solver.set_target(...)
//         * state response-function basis via madness::project(...)
//      The callback is a free function so this header has zero
//      coupling to v3's ResponseProtocol / GroundState / kernels.
//   2. solver.refresh_convergence_targets() — re-derives BSH /
//      density-residual targets from the new thresh.
//   3. solvers::iterate(solver, state, IteratePolicy{max_iters}) —
//      drive to convergence at this protocol.
//
// The intermediate state propagates between protocol steps as the
// initial guess for the next one (post-projection inside prepare).
// ESSolver and the upcoming FDSolver share this driver — the FD/ES
// axis is entirely inside the solver's step() body.
//
// See docs/12_solver_architecture_sketch.md and
// solvers/convergence_policy.hpp.
// =========================================================================

#include "iterate.hpp"

#include <utility>
#include <vector>

namespace molresponse_v3::solvers {

struct IterateProtocolPolicy {
  int max_iters_per_step = 25;
};

/// SFINAE helper: call solver.sort_state_by_omega(state) if the solver
/// has one, otherwise no-op. ESSolver provides this for ascending-omega
/// canonicalization. FDSolver won't (slots are perturbation channels —
/// already canonical).
template <typename Solver, typename State>
auto maybe_canonicalize(Solver &solver, State &state, int)
    -> decltype(solver.sort_state_by_omega(state), void()) {
  solver.sort_state_by_omega(state);
}
template <typename Solver, typename State>
void maybe_canonicalize(Solver &, State &, long) { /* no-op */ }

template <typename Solver, typename PrepareFn>
auto iterate_protocol(Solver &solver,
                      typename Solver::State state,
                      const std::vector<double> &thresholds,
                      const PrepareFn &prepare,
                      const IterateProtocolPolicy &policy = {})
    -> typename Solver::State {
  for (double thresh : thresholds) {
    prepare(thresh, solver, state);
    solver.refresh_convergence_targets();
    state = solvers::iterate(
        solver, std::move(state),
        IteratePolicy{policy.max_iters_per_step});
    if (state.diverged) break;  // bail rather than tighten on garbage
  }
  // Canonicalize output (ES: sort by ascending omega; FD: no-op).
  // Caller's PASS/FAIL comparison can then compare by index against a
  // canonical reference set.
  maybe_canonicalize(solver, state, 0);
  return state;
}

} // namespace molresponse_v3::solvers

#endif // MOLRESPONSE_V3_SOLVERS_ITERATE_PROTOCOL_HPP
