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
//   StopReason Solver::stop_reason(const State&, int iterations, int cap,
//                                  int forced_iterations, bool final_protocol) const;
//   (a solver not yet on rung_exit supplies bool converged(const State&) const
//    instead; none is left in the tree)
//
// Solvers that need optional acceleration (KAIN) or step restriction
// fold those into ::step. The driver itself stays minimal — one loop,
// one exit decision (the shared rung_exit), no policy.
//
// See solvers/convergence_policy.hpp.
// =========================================================================

#include <stdexcept>
#include <functional>
#include <madness/chem/IterationOutput.h>
#include <madness/chem/ConvergenceSettings.h>

namespace molresponse_v3::solvers {

// Adopted solvers own their state and cached equation work. The observer sees
// exactly the evaluated functions, before any update can mutate them.
template <class Solver, class Observer>
madness::iteration::Decision iterate_shared(
    Solver& solver, typename Solver::State& state,
    const madness::iteration::ResolvedSettings& settings, bool final_protocol,
    const madness::iteration::ExpectedTask& expected,
    madness::iteration::IterationInfo info, std::size_t previous_updates,
    const Observer& observe,
    const std::function<void(const madness::iteration::IterationInfo&,
                            const madness::iteration::IterationTiming&)>& observe_timing = {}) {
  namespace it = madness::iteration;
  solver.reset_shared();
  it::Progress progress{0, previous_updates, true};
  for (;;) {
    it::ElapsedTimer<> timer;
    it::IterationTiming timing;
    info.iteration = progress.main_updates;
    const auto metrics = solver.evaluate(state, info);
    timing.steps.emplace_back("evaluate", timer.lap());
    progress.state_valid = !state.diverged;
    const auto decision = it::decide(settings, expected, progress, metrics, final_protocol);
    timing.steps.emplace_back("decide", timer.lap());
    observe(state, metrics, decision, progress);
    timing.steps.emplace_back("checkpoint/report", timer.lap());
    if (decision.action == it::Action::update) {
      solver.advance(state);
      timing.steps.emplace_back("advance", timer.lap());
      ++progress.main_updates;
      ++progress.updates_at_protocol;
    }
    timing.seconds = timer.seconds();
    if (observe_timing) observe_timing(info, timing);
    if (decision.action != it::Action::update) return decision;
  }
}

struct IteratePolicy {
  int max_iters = 25;
  // Explicit diagnostic only: zero preserves the ordinary convergence loop.
  int forced_iterations = 0;
  bool final_protocol = true; // false when a later requested protocol remains
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
    // DESIGN B5: a solver with stop_reason() (the FD and ES solvers) leaves
    // the order of the exits to the shared rung_exit; divergence, its
    // validity check, comes first. The branch below serves a solver without
    // stop_reason(); none is left in the tree.
    if constexpr (requires { solver.stop_reason(state, 0, 0, 0, true); }) {
      if (iteration_diverged(state, 0))
        break;
      state.stop = solver.stop_reason(state, k + 1, policy.max_iters,
                                      policy.forced_iterations, policy.final_protocol);
      if (state.stop == decltype(state.stop)::in_progress)
        continue;
      converged = state.stop == decltype(state.stop)::converged;
      break;
    } else if (policy.forced_iterations || !policy.final_protocol) {
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
