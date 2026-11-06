// Protocol control without numerical kernels: exercise the real drivers and
// shared exit policy with deterministic convergence and divergence evidence.
#include "../solvers/iterate_protocol.hpp"
#include <madness/chem/ConvergenceSettings.h>

#include <iostream>
#include <utility>
#include <vector>

namespace {
namespace cv = madness::convergence;
namespace solvers = molresponse_v3::solvers;
int failures = 0;
void check(bool ok, const char* message) {
  if (!ok) {
    ++failures;
    std::cerr << "FAILED: " << message << '\n';
  }
}

struct Solver {
  struct State {
    int iter = 0; // historical count, never the fresh allowance
    int rung_iter = 0;
    bool diverged = false;
    bool stalled = false;
    cv::StopReason stop = cv::StopReason::in_progress;
  };
  int diverge_at = 0;
  bool checks_pass = true;
  bool transient_convergence = false;
  bool plateau = false;
  bool reported_converged = false;

  [[nodiscard]] State step(State state) const {
    ++state.iter;
    ++state.rung_iter;
    state.diverged = state.rung_iter == diverge_at;
    state.stalled = plateau;
    return state;
  }
  [[nodiscard]] cv::StopReason stop_reason(const State& state, int iterations, int cap,
                            int forced, bool final_protocol = true) const {
    const bool met = checks_pass && state.rung_iter >= 2
                     && (!transient_convergence || state.rung_iter == 2);
    cv::RungSample sample{cv::Decider::engine_check, met, met, state.stalled,
                          iterations, cap, forced};
    sample.final_protocol = final_protocol;
    return cv::rung_exit(sample);
  }
  void refresh_convergence_targets() const {}
  void print_final(const State&, bool converged) { reported_converged = converged; }
};

struct Trace {
  std::vector<std::pair<double, int>> completed;
  [[nodiscard]] static auto prepare() {
    return [](double, Solver&, Solver::State& state) {
      state.rung_iter = 0;
      state.stop = cv::StopReason::in_progress;
    };
  }
  auto save() {
    return [this](double threshold, Solver&, const Solver::State& state) {
      completed.emplace_back(threshold, state.rung_iter);
    };
  }
};
}

int main() {
  {
    Solver solver;
    Trace trace;
    const auto state = solvers::iterate_protocol(
        solver, Solver::State{}, {1e-4, 1e-6}, Trace::prepare(), trace.save(), {5, 0});
    check(trace.completed == std::vector<std::pair<double, int>>{{1e-4, 5}, {1e-6, 2}},
          "coarse stage consumes five iterations and final stage still executes");
    check(state.iter == 7 && solver.reported_converged,
          "final stage retains its convergence stop");
  }
  {
    Solver solver;
    Trace trace;
    solvers::IterateProtocolPolicy policy{5, 0};
    policy.final_protocol = false;
    const auto state = solvers::iterate_protocol(
        solver, Solver::State{}, {1e-4}, Trace::prepare(), trace.save(), policy);
    check(state.iter == 5,
          "one-threshold scheduler call can still be an intermediate stage");
  }
  {
    Solver solver;
    solver.checks_pass = false;
    solver.plateau = true;
    solvers::IteratePolicy policy{5, 0};
    policy.final_protocol = false;
    const auto state = solvers::iterate(solver, Solver::State{}, policy);
    check(state.iter == 5 && state.stop == cv::StopReason::iteration_limit
              && !solver.reported_converged,
          "intermediate plateau cannot shorten the allowance or claim convergence");
  }
  {
    Solver solver;
    solvers::IteratePolicy policy{5, 0};
    policy.final_protocol = false;
    Solver::State saved;
    saved.iter = saved.rung_iter = 37;
    const auto state = solvers::iterate(solver, saved, policy);
    check(state.iter == 42, "restart gets five new iterations despite historical work");
  }
  {
    Solver solver;
    solver.transient_convergence = true;
    solvers::IteratePolicy policy{5, 0};
    policy.final_protocol = false;
    const auto state = solvers::iterate(solver, Solver::State{}, policy);
    check(state.iter == 5 && state.stop == cv::StopReason::iteration_limit,
          "earlier passing checks do not replace the completed-stage verdict");
  }
  {
    Solver solver;
    solver.diverge_at = 3;
    Trace trace;
    const auto state = solvers::iterate_protocol(
        solver, Solver::State{}, {1e-4, 1e-6}, Trace::prepare(), trace.save(), {5, 0});
    check(state.diverged && state.iter == 3 && trace.completed.size() == 1,
          "divergence ends fixed intermediate work without advancing");
    check(!solver.reported_converged, "divergence is never reported as convergence");
  }
  {
    Solver solver;
    solvers::IteratePolicy policy{1, 0};
    policy.final_protocol = false;
    const auto state = solvers::iterate(solver, Solver::State{}, policy);
    check(state.iter == 1 && state.stop == cv::StopReason::iteration_limit,
          "one-iteration allowance is honored");
    policy.max_iters = 5;
    policy.forced_iterations = 3;
    check(solvers::iterate(solver, Solver::State{}, policy).iter == 3,
          "explicit diagnostic count retains precedence");
  }
  std::cout << "iteration protocol failures: " << failures << '\n';
  return failures ? 1 : 0;
}
