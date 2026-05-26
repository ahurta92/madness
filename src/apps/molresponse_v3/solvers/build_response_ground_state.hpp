#ifndef MOLRESPONSE_V3_SOLVERS_BUILD_RESPONSE_GROUND_STATE_HPP
#define MOLRESPONSE_V3_SOLVERS_BUILD_RESPONSE_GROUND_STATE_HPP

// =========================================================================
// Build helpers for response solvers (ES and FD share these):
//
//   ResponseGroundState   = zeroth-order kernel inputs (orbitals, eps,
//                           V_local, fock, Coulomb op, Q-projector).
//                           Rebuilt per protocol step.
//   ESProblem<T, S>       = ResponseGroundState + n_roots.
//   FDProblem<T, S>       = ResponseGroundState + vector<FDPerturbation>.
//                           Lives in solvers/fd_problem.hpp.
//
// Builders:
//   build_response_ground_state_closed_shell(world, gs, c_xc, lo)
//   build_response_ground_state_open_shell  (world, gs, c_xc, lo)
//   build_es_problem_closed_shell           (world, gs, n_roots, c_xc, lo)
//   build_es_problem_open_shell             (world, gs, n_roots, c_xc, lo)
//   build_initial_guess_tda_closed_shell    (world, gs, n_roots) → State
//   build_initial_guess_tda_open_shell      (world, gs, n_roots) → State
//
// Separated from es_solver.hpp because they bring in GroundState and
// ESSolverGuess — keeping those out of the pure architecture header.
// =========================================================================

#include "../ESSolverGuess.hpp"  // make_initial_guess_tda_rhf / _uhf
#include "../GroundState.hpp"
#include "es_solver.hpp"         // ESSolver + ESProblem
#include "response_state.hpp"

#include <madness/chem/SCFOperators.h>

#include <vector>

namespace molresponse_v3 {

/// Build ResponseGroundState for an OpenShell run. Populates both
/// alpha and beta fields. V_local is shared (HF / pure-XC where Vxc
/// is spin-blind — for spin-resolved DFT this would need wiring).
inline ResponseGroundState
build_response_ground_state_open_shell(madness::World &world, GroundState &gs,
                                        double c_xc = 1.0,
                                        double lo   = 1.0e-10) {
  ResponseGroundState t;
  t.amo            = gs.orbitals_alpha();
  t.aeps           = gs.energies_alpha();
  t.V_local_alpha  = gs.V_local();
  t.focka          = gs.focka();
  t.focka_no_diag  = gs.focka_no_diag();
  t.Qa             = gs.Q_alpha();

  t.bmo            = gs.orbitals_beta();
  t.beps           = gs.energies_beta();
  t.V_local_beta   = gs.V_local();         // HF: same potential both spins
  t.fockb          = gs.fockb();
  t.fockb_no_diag  = gs.fockb_no_diag();
  t.Qb             = gs.Q_beta();

  t.coulop = poperatorT(madness::CoulombOperatorPtr(
      world, lo, madness::FunctionDefaults<3>::get_thresh()));

  t.c_xc = c_xc;
  t.lo   = lo;
  return t;
}

/// Build ResponseGroundState for the closed-shell case.
inline ResponseGroundState
build_response_ground_state_closed_shell(madness::World &world, GroundState &gs,
                                          double c_xc = 1.0,
                                          double lo   = 1.0e-10) {
  ResponseGroundState t;
  t.amo            = gs.orbitals_alpha();
  t.aeps           = gs.energies_alpha();
  t.V_local_alpha  = gs.V_local();
  t.focka          = gs.focka();
  t.focka_no_diag  = gs.focka_no_diag();
  t.Qa             = gs.Q_alpha();
  // beta side stays empty — closed-shell.

  t.coulop = poperatorT(madness::CoulombOperatorPtr(
      world, lo, madness::FunctionDefaults<3>::get_thresh()));

  t.c_xc = c_xc;
  t.lo   = lo;
  return t;
}

// ----- ESProblem convenience builders ----------------------------------

template <typename Shell>
inline ESProblem<TDA, Shell>
build_es_problem_tda(madness::World &world, GroundState &gs, int n_roots,
                      double c_xc = 1.0, double lo = 1.0e-10) {
  ESProblem<TDA, Shell> p;
  if constexpr (std::is_same_v<Shell, ClosedShell>)
    p.gs = build_response_ground_state_closed_shell(world, gs, c_xc, lo);
  else
    p.gs = build_response_ground_state_open_shell(world, gs, c_xc, lo);
  p.n_roots = n_roots;
  return p;
}

// ----- Initial-guess adapters ------------------------------------------

/// Adapt the legacy ESSolverGuess output (vector<RealResponseState>)
/// into the new per-root storage (vector<ResponseStateX<ClosedShell>>).
inline ESSolver<TDA, ClosedShell>::State
build_initial_guess_tda_closed_shell(madness::World &world, GroundState &gs,
                                     long n_roots) {
  auto guess = make_initial_guess_tda_rhf(world, gs, n_roots,
                                          /*has_y=*/false);

  ESSolver<TDA, ClosedShell>::State s;
  s.roots.resize(n_roots);
  for (long r = 0; r < n_roots; ++r) {
    s.roots[r].x_alpha = guess[r].x_alpha;
  }
  s.iter = 0;
  return s;
}

/// Adapter for OpenShell TDA — both α and β response components populated.
inline ESSolver<TDA, OpenShell>::State
build_initial_guess_tda_open_shell(madness::World &world, GroundState &gs,
                                   long n_roots) {
  auto guess = make_initial_guess_tda_uhf(world, gs, n_roots);

  ESSolver<TDA, OpenShell>::State s;
  s.roots.resize(n_roots);
  for (long r = 0; r < n_roots; ++r) {
    s.roots[r].x_alpha = guess[r].x_alpha;
    s.roots[r].x_beta  = guess[r].x_beta;
  }
  s.iter = 0;
  return s;
}

// ----- State slicing + oversampled warmup ------------------------------

/// Truncate a TDA-ES State to its lowest n_keep roots. Assumes the
/// input state has already been sorted ascending in omega (call
/// `solver.sort_state_by_omega(...)` first). Returns a fresh State
/// containing only roots[0..n_keep-1]. iter is reset to 0 so the
/// main solver sees a clean run.
template <typename Shell>
inline typename ESSolver<TDA, Shell>::State
slice_state_lowest(const typename ESSolver<TDA, Shell>::State &in,
                   long n_keep) {
  MADNESS_CHECK(n_keep > 0);
  MADNESS_CHECK(static_cast<long>(in.roots.size()) >= n_keep);

  typename ESSolver<TDA, Shell>::State out;
  out.roots.assign(in.roots.begin(), in.roots.begin() + n_keep);

  out.omega = madness::Tensor<double>(n_keep);
  for (long i = 0; i < n_keep; ++i) out.omega(i) = in.omega(i);

  auto trim_double = [&](const std::vector<double> &v) {
    return std::vector<double>(v.begin(),
                               v.begin() + std::min<long>(n_keep, v.size()));
  };
  out.last_bsh_residual     = trim_double(in.last_bsh_residual);
  out.last_density_residual = trim_double(in.last_density_residual);

  if (static_cast<long>(in.rho_alpha_prev.size()) >= n_keep) {
    out.rho_alpha_prev.assign(in.rho_alpha_prev.begin(),
                              in.rho_alpha_prev.begin() + n_keep);
  }

  out.iter     = 0;     // fresh count for the main solver
  out.diverged = false;
  return out;
}

/// Run an over-sampled TDA warm-up and return a state sliced to the
/// lowest `n_roots_final` roots, sorted ascending in omega. Mirrors
/// the legacy iterate_trial → sort → select_functions → main-iter
/// pattern (ExcitedResponse.cpp:83-105):
///
///   * Build problem with n_roots_warmup (~ 2·n_roots_final by
///     convention) so the bundle has room for "extra" trial vectors
///     that converge to higher roots and get filtered out at the
///     downselect step.
///   * Run `warmup_iters` ESSolver step() calls. KAIN is disabled
///     during this pass (`policy.tda_warmup_iters` is overridden
///     to `warmup_iters`); only BSH + step restriction.
///   * Sort by omega, slice to n_roots_final lowest.
///
/// The returned state's KAIN history is empty (it lived in the
/// warmup solver which is now out of scope) — the main solver picks
/// up from a clean bundle. Cheap iters, big payoff for cold guesses.
template <typename Shell>
inline typename ESSolver<TDA, Shell>::State
run_oversampled_tda_warmup(madness::World &world, GroundState &gs,
                            long n_roots_final, long n_roots_warmup,
                            int warmup_iters,
                            ConvergencePolicy base_policy,
                            double c_xc = 1.0, double lo = 1.0e-10,
                            PrintLevel print_level = PrintLevel::Normal) {
  MADNESS_CHECK(n_roots_warmup >= n_roots_final);
  MADNESS_CHECK(warmup_iters >= 0);

  if (n_roots_warmup == n_roots_final && warmup_iters == 0) {
    // No-op fast path — caller wants neither oversample nor warmup.
    if constexpr (std::is_same_v<Shell, ClosedShell>)
      return build_initial_guess_tda_closed_shell(world, gs, n_roots_final);
    else
      return build_initial_guess_tda_open_shell(world, gs, n_roots_final);
  }

  // Force-disable KAIN during warmup — iterate_trial-style pure BSH
  // power iteration. tda_warmup_iters is also zeroed so we don't
  // double-count it (the warmup IS the warmup).
  ConvergencePolicy warmup_policy = base_policy;
  warmup_policy.kain              = false;
  warmup_policy.tda_warmup_iters  = 0;

  if (world.rank() == 0 && print_level >= PrintLevel::Normal) {
    print("\n=== TDA warmup ===  n_roots_warmup =", n_roots_warmup,
          "  n_roots_final =", n_roots_final,
          "  warmup_iters =", warmup_iters);
  }

  // Build the oversampled problem and the warm-up solver. Use a
  // fresh solver instance — its KAIN state lives only here.
  auto problem = build_es_problem_tda<Shell>(world, gs, n_roots_warmup,
                                              c_xc, lo);
  ESSolver<TDA, Shell> warmup_solver(world, std::move(problem),
                                      warmup_policy, print_level);

  typename ESSolver<TDA, Shell>::State state;
  if constexpr (std::is_same_v<Shell, ClosedShell>)
    state = build_initial_guess_tda_closed_shell(world, gs, n_roots_warmup);
  else
    state = build_initial_guess_tda_open_shell(world, gs, n_roots_warmup);

  for (int i = 0; i < warmup_iters; ++i) {
    state = warmup_solver.step(std::move(state));
    if (state.diverged) {
      if (world.rank() == 0)
        print("[WARMUP] aborted at iter", state.iter,
              "— bundle diverged; falling back to last-good state");
      break;
    }
  }

  warmup_solver.sort_state_by_omega(state);

  if (world.rank() == 0 && print_level >= PrintLevel::Normal) {
    print("[WARMUP] final omegas (oversampled, sorted ascending):");
    print(state.omega);
    print("[WARMUP] selecting lowest", n_roots_final, "for main iter");
  }

  return slice_state_lowest<Shell>(state, n_roots_final);
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_BUILD_RESPONSE_GROUND_STATE_HPP
