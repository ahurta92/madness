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

  // Cluster-unmix threshold factor in rs::diagonalize. The legacy
  // path uses 100·thresh for TDA (loose) and 10·thresh for Full/RPA
  // (tighter — clusters with smaller separation get polar-decomp
  // unmixed). Tighter unmix is more aggressive at freezing eigenvector
  // mixing within a cluster, which stabilises slot identity across
  // iters when omega's are closely spaced (the Li ω=0.04 regime).
  // Default 100 matches the TDA convention; set to 10 for Full or any
  // run with closely-clustered omegas.
  double cluster_unmix_factor = 100.0;

  // Top-of-iter orthonormalization mode for the bundle.
  //
  //   GramSchmidt  — 2-pass Schmidt orthogonalization in increasing slot
  //                  order, then per-slot normalize. Slot 0 is anchored
  //                  (unchanged); higher slots are projected against
  //                  lower ones. Cheap (m² inner products), order-
  //                  dependent — when two states are nearly co-linear
  //                  after BSH+KAIN, the high-index slot loses its
  //                  identity to the low-index one. Default — matches
  //                  legacy iterate_trial discipline.
  //
  //   Lowdin       — symmetric Löwdin: X' = X · S^{-1/2}. All slots
  //                  perturbed equally; no privileged index. Drops
  //                  directions with eigenvalue < 10·thresh (linear-
  //                  dependence detector). Minimises ‖X' − X‖_F among
  //                  orthonormalizations. Recommended when [ROT-SLOTS]
  //                  shows iter-to-iter REORDERED dominance perm.
  enum class OrthonormalizeMode { GramSchmidt, Lowdin };
  OrthonormalizeMode orthonormalize_mode = OrthonormalizeMode::GramSchmidt;

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
  // Blowup safeguard threshold for KAIN coefficient max-abs. When the
  // KAIN solve returns coefficients with |c|max exceeding this value
  // (after rcond escalation has been exhausted), the per-state history
  // is cleared and the iteration falls back to raw BSH for that step.
  //
  // SCF::update_subspace uses 3.0 — appropriate for orbital optimization
  // where large |c| genuinely indicates instability (Brillouin theorem
  // etc.). RESPONSE iterations are linear; large |c| just means strong
  // acceleration toward a polarizable mode (e.g. Li at α≈170). The
  // SCF-style 3.0 trips bailout every iter for such systems and kills
  // KAIN entirely. Bumped to a permissive 100 by default so KAIN can
  // actually accelerate strongly-polarizable response iterations;
  // back-set to 3.0 if you want strict SCF semantics.
  double kain_cmax_cap = 100.0;
  // Warmup oversampling factor — used by the run_oversampled_tda_warmup
  // helper. The warm-up phase runs with ceil(warmup_oversample_factor *
  // n_roots) trial states; after warmup completes the lowest n_roots
  // are kept for the main iteration. Mirrors legacy iterate_trial's
  // "trial bundle is 2× requested states, then select_functions picks
  // the N lowest" pattern (ExcitedResponse.cpp:83-105 and select_
  // functions). Default 1.0 (no oversampling — back-compat). Set to
  // 2.0 for cold CIS-style guesses; the extra trial vectors absorb
  // higher-root contamination so the kept N converge cleanly.
  // Only takes effect if tda_warmup_iters > 0.
  double warmup_oversample_factor = 1.0;

  // TDA-warmup iters: KAIN is disabled for the first `tda_warmup_iters`
  // step() calls, then turns on. Mirrors legacy iterate_trial's
  // "no-KAIN BSH-only filter" pre-pass (ExcitedResponse.cpp:456-666):
  // a few power-iteration steps clean up a cold guess before KAIN's
  // per-slot history starts recording. Without this, KAIN memorizes
  // the unstable iter-0/iter-1 slot ordering and slot identity
  // flickers across iters when the rotation re-permutes. Default 0
  // (KAIN on from iter 1, current behaviour). Recommended 5-10 for
  // CIS-style guesses, 0 for restart from converged lower-protocol
  // state. (Affects ESSolver only; ignored by FDSolver.)
  int  tda_warmup_iters = 0;
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
