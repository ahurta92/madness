#ifndef MOLRESPONSE_V3_SOLVERS_BUILD_RESPONSE_GROUND_STATE_HPP
#define MOLRESPONSE_V3_SOLVERS_BUILD_RESPONSE_GROUND_STATE_HPP

// =========================================================================
// Build helpers for response solvers (ES and FD share these):
//   - build_response_ground_state_closed_shell(world, gs, n_roots, c_xc, lo)
//   - build_response_ground_state_open_shell(world, gs, n_roots, c_xc, lo)
//   - build_initial_guess_tda_closed_shell(world, gs, n_roots) -> State
//   - build_initial_guess_tda_open_shell(world, gs, n_roots)   -> State
//
// Separated from es_solver.hpp because they bring in GroundState and
// ESSolverGuess (the legacy guess routine) — keeping those out of the
// pure architecture header.
//
// TODO: drop the `n_roots` argument from the builders — it's a
// solver-level concern, not a ground-state property. The builders
// currently set ResponseGroundState::n_roots but the solver ignores it.
// =========================================================================

#include "../ESSolverGuess.hpp"  // make_initial_guess_tda_rhf
#include "../GroundState.hpp"
#include "es_solver.hpp"
#include "response_state.hpp"

#include <madness/chem/SCFOperators.h>

#include <vector>

namespace molresponse_v3 {

/// Build an ResponseGroundState for an OpenShell run. Populates both alpha and
/// beta fields. V_local is shared (HF / pure-XC where Vxc is spin-blind
/// — for spin-resolved DFT this would need wiring).
inline ResponseGroundState
build_response_ground_state_open_shell(madness::World &world, GroundState &gs,
                               int n_roots, double c_xc = 1.0,
                               double lo = 1.0e-10) {
  ResponseGroundState t;
  t.phi0_alpha          = gs.orbitals_alpha();
  t.eps_alpha           = gs.energies_alpha();
  t.V_local_alpha       = gs.V_local();
  t.fock_full_alpha     = gs.focka();
  t.fock_no_diag_alpha  = gs.focka_no_diag();
  t.Q_alpha             = gs.Q_alpha();

  t.phi0_beta           = gs.orbitals_beta();
  t.eps_beta            = gs.energies_beta();
  t.V_local_beta        = gs.V_local();         // HF: same potential both spins
  t.fock_full_beta      = gs.fockb();
  t.fock_no_diag_beta   = gs.fockb_no_diag();
  t.Q_beta              = gs.Q_beta();

  t.coulop = poperatorT(madness::CoulombOperatorPtr(
      world, lo, madness::FunctionDefaults<3>::get_thresh()));

  t.c_xc    = c_xc;
  t.lo      = lo;
  t.n_roots = n_roots;
  return t;
}

/// Build ResponseGroundState from a prepared GroundState for the closed-shell case.
inline ResponseGroundState
build_response_ground_state_closed_shell(madness::World &world, GroundState &gs,
                                 int n_roots, double c_xc = 1.0,
                                 double lo = 1.0e-10) {
  ResponseGroundState t;
  t.phi0_alpha          = gs.orbitals_alpha();
  t.eps_alpha           = gs.energies_alpha();
  t.V_local_alpha       = gs.V_local();
  t.fock_full_alpha     = gs.focka();
  t.fock_no_diag_alpha  = gs.focka_no_diag();
  t.Q_alpha             = gs.Q_alpha();
  // beta side stays empty — closed-shell.

  // Shared Coulomb (Poisson) operator. The existing solvers build this
  // alongside the kernels; we replicate the convention here.
  t.coulop = poperatorT(madness::CoulombOperatorPtr(
      world, lo, madness::FunctionDefaults<3>::get_thresh()));

  t.c_xc    = c_xc;
  t.lo      = lo;
  t.n_roots = n_roots;
  return t;
}

/// Adapt the legacy ESSolverGuess output (vector<RealResponseState>)
/// into the new per-root storage (vector<ResponseStateX<ClosedShell>>).
/// Returns a fresh State.iter == 0 with omega and residuals empty.
inline ESSolver<TDA, ClosedShell>::State
build_initial_guess_tda_closed_shell(madness::World &world, GroundState &gs,
                                     long n_roots) {
  // Existing make_initial_guess_tda_rhf gives us a vector of
  // RealResponseState. Each has x_alpha populated and the rest empty.
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

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_BUILD_RESPONSE_GROUND_STATE_HPP
