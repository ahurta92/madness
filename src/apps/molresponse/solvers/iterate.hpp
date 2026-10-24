#ifndef MOLRESPONSE_V3_SOLVERS_ITERATE_HPP
#define MOLRESPONSE_V3_SOLVERS_ITERATE_HPP

// =========================================================================
// solvers::iterate — single driver template for every solver type.
//
// Any concrete solver (ClosedShellTdaSolver, ClosedShellFullSolver,
// OpenShellTdaSolver, ClosedShellFrequencyDependentSolver, ...) is
// driven by this one function. The step order is auditable here in
// one place; no solver class can re-order its own pipeline silently.
//
// The Solver concept (duck-typed, no inheritance required):
//
//   typename Solver::State
//
//   State Solver::step(State in);
//   bool  Solver::converged(const State&) const;
//
// Solvers that need optional acceleration (KAIN) or step restriction
// fold those into ::step. The driver itself stays minimal — one loop,
// one convergence check, no policy.
//
// See solvers/convergence_policy.hpp.
// =========================================================================

#include <stdexcept>

namespace molresponse_v3::solvers {

struct IteratePolicy {
  int max_iters = 25;
  // Explicit diagnostic only: zero preserves the ordinary convergence loop.
  int forced_iterations = 0;
};

template <typename State>
auto iteration_diverged(const State& state, int) -> decltype(bool(state.diverged)) { return state.diverged; }
template <typename State>
bool iteration_diverged(const State&, long) { return false; }
template <typename State>
auto iteration_stalled(const State& state, int) -> decltype(bool(state.stalled)) { return state.stalled; }
template <typename State>
bool iteration_stalled(const State&, long) { return false; }

template <typename Solver>
auto iterate(Solver &solver, typename Solver::State state,
             const IteratePolicy &policy) -> typename Solver::State {
  if (policy.forced_iterations < 0) throw std::invalid_argument("forced_iterations must be nonnegative");
  bool converged = false;
  const int count = policy.forced_iterations ? policy.forced_iterations : policy.max_iters;
  for (int k = 0; k < count; ++k) {
    state = solver.step(state);
    // DESIGN B5: a solver with stop_reason() (the FD solver) leaves the order
    // of the exits to the shared rung_exit; divergence, its validity check,
    // comes first. The others keep the loop below until their named ports.
    if constexpr (requires { solver.stop_reason(state, 0, 0, 0); }) {
      if (iteration_diverged(state, 0))
        break;
      state.stop = solver.stop_reason(state, k + 1, policy.max_iters,
                                      policy.forced_iterations);
      if (state.stop == decltype(state.stop)::in_progress)
        continue;
      converged = state.stop == decltype(state.stop)::converged;
      break;
    } else if (policy.forced_iterations) {
      if (iteration_diverged(state, 0)) { converged = false; break; }
      converged = !iteration_stalled(state, 0) && solver.converged(state);
    } else if (solver.converged(state)) {
      converged = true;
      break;
    }
  }
  solver.print_final(state, converged);
  return state;
}

} // namespace molresponse_v3::solvers

#endif // MOLRESPONSE_V3_SOLVERS_ITERATE_HPP
