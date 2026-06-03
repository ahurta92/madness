#ifndef MOLRESPONSE_V3_SOLVERS_CONVERGENCE_POLICY_HPP
#define MOLRESPONSE_V3_SOLVERS_CONVERGENCE_POLICY_HPP

// =========================================================================
// ConvergencePolicy — shared by ESSolver<T,S> and (future) FDSolver<T,S>.
//
// Policy carries USER intent (target dconv, behaviour switches). At
// each protocol level the policy resolves to EFFECTIVE TARGETS using the
// SCF convention (SCF.cc::solve, lines 2143 + 2382):
//
//   dconv          = max(thresh, dconv_user)          // can't beat the protocol
//   density_target = density_residual_factor * dconv  // SCF uses max(5,natom)
//   bsh_target     = bsh_residual_factor     * dconv  // SCF uses 5.0
//
// `dconv = max(thresh, dconv_user)` is the key idea: at a given protocol
// you cannot resolve a residual below ~thresh (past that you're just
// resolving noise until the next, deeper protocol refines the wavelet
// basis), so dconv floors at thresh; and you never ask for tighter than
// the user wants, so it also floors at dconv_user. BOTH gates then ride at
// ~5× dconv (SCF.cc:2382 uses 5*dconv for bsh and dconv*max(5,natom) for
// density). The 5× headroom is what lets a single coarse protocol reach
// "converged" rather than stalling right at the thresh noise floor — at
// thresh == dconv_user the targets are 5*thresh, comfortably above the
// per-protocol floor.
//
// TODO(study): the factor 5.0 and dconv_user defaults are inherited from
// SCF; response may want different values. Increasing thresh (deeper
// protocol) is what actually buys resolution — within a protocol the
// residual floors at noise after a few iters. A sweep over
// bsh_residual_factor / protocol ladder is planned to pick response-
// appropriate values.
// =========================================================================

#include <algorithm>

namespace molresponse_v3 {

struct ConvergencePolicy {
  // User-facing density-convergence target. The effective target at each
  // protocol is a multiple of max(thresh, dconv_user) — see
  // effective_for_thresh.
  double dconv_user = 1.0e-4;

  // Both gates ride at a multiple of the (protocol-floored) dconv,
  // matching SCF::solve (SCF.cc:2382):
  //   bsh     < 5.0 * dconv
  //   density < dconv * max(5, natom)   -> 5.0 for small molecules
  // Using 5.0 for both gives each gate ~5x headroom over the per-protocol
  // noise floor, so a single coarse protocol can actually reach
  // "converged" instead of stalling right at thresh. (SCF's density gate
  // loosens further with atom count; we keep a flat factor here since the
  // response per-state density residual isn't naturally natom-scaled —
  // revisit in the planned convergence study.) Exposed so the study can
  // sweep them independently.
  double bsh_residual_factor     = 5.0;
  double density_residual_factor = 5.0;

  // Cluster-unmix threshold factor in rs::diagonalize. The legacy
  // path uses 100·thresh for TDA (loose) and 10·thresh for Full/RPA
  // (tighter — clusters with smaller separation get polar-decomp
  // unmixed). Tighter unmix is more aggressive at freezing eigenvector
  // mixing within a cluster, which stabilises slot identity across
  // iters when omega's are closely spaced (the Li ω=0.04 regime).
  // Default 100 matches the TDA convention; set to 10 for Full or any
  // run with closely-clustered omegas.
  double cluster_unmix_factor = 100.0;

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

  // Step-restriction granularity:
  //   PerOrbital — clamp each response function independently (default; the
  //                legacy SCF::do_step_restriction behaviour). A single runaway
  //                orbital is damped without touching the well-behaved ones.
  //   PerState   — measure ||v_new - v_old|| over the WHOLE flattened state (all
  //                functions of a root / channel) and scale the ENTIRE vector by
  //                one factor maxrotn/||diff||. Damps a coherent rotation of the
  //                whole state (useful when a near-degenerate root rotates as a
  //                block rather than one orbital running away).
  enum class StepRestrictMode { PerOrbital, PerState };
  StepRestrictMode step_restrict_mode = StepRestrictMode::PerOrbital;

  struct Targets {
    double bsh_residual;     // ‖x_old − x_new‖ cap
    double density_residual; // ‖ρ_new − ρ_old‖ cap
  };

  Targets effective_for_thresh(double thresh) const {
    // SCF convention (SCF.cc:2143): can't resolve below the protocol's
    // thresh, won't over-tighten past the user's request.
    const double dconv = std::max(thresh, dconv_user);
    Targets t;
    t.density_residual = density_residual_factor * dconv;  // SCF.cc:2382 (da < dconv*max(5,natom))
    t.bsh_residual     = bsh_residual_factor * dconv;      // SCF.cc:2382 (bsh < 5*dconv)
    return t;
  }
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_CONVERGENCE_POLICY_HPP
