#ifndef MOLRESPONSE_V3_SOLVERS_CONVERGENCE_POLICY_HPP
#define MOLRESPONSE_V3_SOLVERS_CONVERGENCE_POLICY_HPP

// =========================================================================
// ConvergencePolicy — shared by ESSolver<T,S> and (future) FDSolver<T,S>.
//
// Policy carries USER intent (target dconv, behaviour switches). At
// each protocol level the policy resolves to EFFECTIVE TARGETS by
// applying floors against the current FunctionDefaults<3>::thresh —
// you can't ask for a tighter BSH residual than the protocol can
// support, so the policy floors at K * thresh.
//
//   bsh_target     = max(B * dconv_user,  F * thresh)
//   density_target = max(    dconv_user,  D * thresh)
//
// Defaults match the discipline that landed in
// src/apps/molresponse/ExcitedResponse.cpp this session
// (50× thresh BSH floor with WARNING when violated).
// =========================================================================

#include <algorithm>

namespace molresponse_v3 {

struct ConvergencePolicy {
  // User-facing target.
  double dconv_user = 1.0e-4;

  // Floor factors against the current protocol thresh.
  double bsh_thresh_factor   = 50.0;  // bsh_target  >= 50 * thresh
  double bsh_user_factor     =  5.0;  // bsh_target  >= 5  * dconv_user
  double density_thresh_factor = 1.0; // density_target >= 1 * thresh

  // ESSolver memory mode. When true, step() streams Lambda assembly
  // and DROPS V0x/E0x/gamma after the subspace rotation, then
  // RECOMPUTES the rotated pieces to assemble Theta. Costs ~2× the
  // per-root kernel work (V0x/E0x/gamma evaluated twice per iter)
  // but cuts peak memory from ~8M Storage to ~3-4M. Set true when
  // n_occ × n_roots × k³ × leaves is the memory bottleneck.
  // Default false: current "keep pieces, rotate them, assemble Theta
  // from rotated pieces" layout.
  bool stream_theta = false;

  // Diverging-residual bail-out. Triggers only on a runaway residual
  // (BSH residual > guard). The legacy ES guard was 2.0 in normalised
  // units, but FD's iter-1 residual from a "x = perturbation" guess is
  // O(perturbation norm) which can easily exceed 2.0 for closed-shell
  // first-row molecules — so 2.0 false-fires on iter 1 in FD. Bumped
  // to 1e3 (effectively off) until we wire a relative-growth check
  // (compare iter k vs iter k-1, bail if it grew by >2×).
  double explosion_guard = 1.0e3;

  // Minimum iters before we even check convergence (lets KAIN warm up).
  int min_iters_before_conv = 0;

  // ---- KAIN acceleration + step restriction ----
  // KAIN is enabled by default. The solver allocates an
  // XNonlinearSolver over the flat state vector; size of the KAIN
  // subspace history is kain_maxsub. Disable by setting kain=false.
  bool kain = true;
  int  kain_maxsub = 10;
  // Step-restriction: if ||x_new − x_old|| > maxrotn after KAIN,
  // damp the step toward x_old. Set < 0 to disable.
  double maxrotn = 0.5;

  struct Targets {
    double bsh_residual;     // ‖x_old − x_new‖ cap
    double density_residual; // ‖ρ_new − ρ_old‖ cap
  };

  Targets effective_for_thresh(double thresh) const {
    Targets t;
    t.bsh_residual = std::max(bsh_user_factor * dconv_user,
                              bsh_thresh_factor * thresh);
    t.density_residual = std::max(dconv_user,
                                  density_thresh_factor * thresh);
    return t;
  }
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_CONVERGENCE_POLICY_HPP
