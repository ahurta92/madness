#ifndef MOLRESPONSE_V3_KERNELS_TPA_HPP
#define MOLRESPONSE_V3_KERNELS_TPA_HPP

// ===========================================================================
// tpa — two-photon-absorption transition moments (single residue of the
// quadratic response), closed-shell / TDA. See TPA_SCOPING.md (§2).
//
// 2PA is "beta with the C-channel replaced by the excited-state eigenvector
// X_f", evaluated at the two-photon resonance omega_1 = omega_2 = omega_f/2.
// The residue of the C-channel linear response at omega -> omega_f IS X_f, a
// HOMOGENEOUS eigenvector: it has no driving operator, so in the C-channel the
// perturbation operator VC_op -> 0. Everything else reuses the Dalton-validated
// beta contraction core verbatim:
//
//   S_ab = beta_abc( A = mu_a-response@(omega_f/2),
//                    VBC = compute_vbc( B = mu_b-response@(omega_f/2), C = X_f,
//                                       VB_op = mu_b, VC_op = 0 ),
//                    B, C = X_f, VA_op = mu_a )
//
// This is the legacy form, kept for tools/tpa_from_dalton.cpp: the eigenvector
// is fed in as the C leg of the driven source. Production 2PA takes the
// residue in the A channel instead (tpa_pq_spec_sym in calc_executor.hpp;
// tpa_moment_residue below is the compute_vbc-based form of the same
// composition, also used by that tool). Since 2026-09-09
// compute_vbc is the spec build of the equation's source (kernels/vbc.hpp), so
// numbers from this arm differ from those recorded before that date.
//
// The rotationally-invariant strength delta^|| (parallel linear polarization,
// F=G=H=2) is the quantity Dalton reports (two_photon_strengths); S_ab is the
// raw tensor.
// ===========================================================================

#include "beta.hpp"   // beta::beta_abc
#include "vbc.hpp"    // vbc::compute_vbc, ResponseGroundState, ResponseStateXY

#include <madness/chem/DerivedQuantities.h>
#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>

#include <array>
#include <cmath>

namespace molresponse_v3::tpa {

/// The 3x3 two-photon transition-moment tensor S_ab for one excited state f
/// (closed-shell / TDA). Inputs are already-converged states:
///   Xf        — the ES eigenvector for root f, wrapped as XY (TDA: y_alpha=0).
///   mu_resp   — the dipole response mu_a at omega_f/2, per Cartesian axis a.
///   mu_op     — the raw dipole operators (MomentFunctor), per axis.
/// No electronic solve: pure contraction via beta_abc + compute_vbc.
/// The three residue VBC quadratic sources for one excited state f (C = X_f
/// homogeneous, VC_op = 0), one per B-axis. This is the TPA "second-order
/// perturbation vector" — the analogue of beta's VBC source. assemble_tpa saves
/// these to disk (mirroring beta/Raman vbc_states) for reuse + inspection.
inline std::array<ResponseStateXY<ClosedShell>, 3>
tpa_sources(madness::World &world, const ResponseGroundState &g0,
            const ResponseStateXY<ClosedShell> &Xf,
            const std::array<ResponseStateXY<ClosedShell>, 3> &mu_resp,
            const std::array<madness::real_function_3d, 3> &mu_op) {
  madness::real_function_3d zero_op = madness::copy(mu_op[0]);  // VC_op = 0
  zero_op.scale(0.0);
  std::array<ResponseStateXY<ClosedShell>, 3> vbc_b;
  for (int b = 0; b < 3; ++b)
    vbc_b[b] = vbc::compute_vbc<ClosedShell>(world, g0, mu_resp[b], Xf,
                                             mu_op[b], zero_op);
  return vbc_b;
}

/// Contract PRE-BUILT sources into the 3x3 S tensor (symmetrized). Kept separate
/// from tpa_sources so assemble_tpa can build the sources once, save them, then
/// contract — the same build->save->contract shape as the beta path.
inline madness::Tensor<double>
tpa_moment(madness::World &world, const ResponseGroundState &g0,
           const ResponseStateXY<ClosedShell> &Xf,
           const std::array<ResponseStateXY<ClosedShell>, 3> &mu_resp,
           const std::array<madness::real_function_3d, 3> &mu_op,
           const std::array<ResponseStateXY<ClosedShell>, 3> &vbc_b) {
  using namespace madness;
  Tensor<double> S(3L, 3L);
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      S(a, b) = beta::beta_abc<ClosedShell>(world, g0, mu_resp[a], vbc_b[b],
                                            mu_resp[b], Xf, mu_op[a]);

  // Degenerate TPA (two identical photons) => S_ab must be symmetric; the raw
  // beta_abc role split (A vs the B-in-VBC channel) computes ONE ordering, so
  // symmetrize to the physical S_ab = <0|mu_a|p><p|mu_b|f> + (a<->b). Validated
  // against Dalton, whose S is symmetric (upper triangle only).
  Tensor<double> Ssym(3L, 3L);
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      Ssym(a, b) = 0.5 * (S(a, b) + S(b, a));
  return Ssym;
}

/// Convenience: build sources + contract in one call (no save).
inline madness::Tensor<double>
tpa_moment(madness::World &world, const ResponseGroundState &g0,
           const ResponseStateXY<ClosedShell> &Xf,
           const std::array<ResponseStateXY<ClosedShell>, 3> &mu_resp,
           const std::array<madness::real_function_3d, 3> &mu_op) {
  return tpa_moment(world, g0, Xf, mu_resp, mu_op,
                    tpa_sources(world, g0, Xf, mu_resp, mu_op));
}

/// SINGLE-RESIDUE 2PA moment (Parker JCTC 2018 / DALTON QRSMO), the corrected
/// channel assignment (TPA_SCOPING §5l-m). The residue of beta is taken in the
/// A channel (omega_sigma = omega_B+omega_C -> omega_f), so:
///   * V^{bc} is built from the TWO PHOTON RESPONSES N^b, N^c (both at
///     omega_f/2) WITH their dipole operators — the exact quadratic RHS beta
///     uses (Eq. 19 / compute_vbc);
///   * the A-channel response is REPLACED by the eigenvector X_f, so the
///     mu_A property terms (beta's b2/b3) drop — the contraction is the bare
///     paired metric inner (beta's b1 only):
///       S_bc = prefactor * ( <x_f|V^{bc}.x> + <y_f|V^{bc}.y> )
///   * normalization: the formula assumes the Parker/TDHF convention
///     X^t X - Y^t Y = 1 over spin-adapted vectors == 0.5 over our spatial
///     functions. Our ES solver normalizes the SPATIAL metric to 1, so its
///     eigenvectors carry sqrt(2) vs that convention -> fold C_N into
///     `prefactor` (pinned empirically via tools/tpa_from_dalton).
/// S is symmetric by construction (compute_vbc symmetrizes the (B,C)+(C,B)
/// halves), so only 6 V^{bc} builds per root.
inline madness::Tensor<double>
tpa_moment_residue(madness::World &world, const ResponseGroundState &g0,
                   const ResponseStateXY<ClosedShell> &Xf,
                   const std::array<ResponseStateXY<ClosedShell>, 3> &mu_resp,
                   const std::array<madness::real_function_3d, 3> &mu_op,
                   double prefactor = 1.0) {
  using namespace madness;
  Tensor<double> S(3L, 3L);
  for (int b = 0; b < 3; ++b) {
    for (int c = b; c < 3; ++c) {
      const double t0 = wall_time();
      if (world.rank() == 0) {
        printf("  [TPA residue] pair %c%c (%d/6): building V^{bc} "
               "(2 orderings x 3 two-electron builds)...\n",
               "xyz"[b], "xyz"[c], b * (5 - b) / 2 + c + 1);
        fflush(stdout);
      }
      auto V = vbc::compute_vbc<ClosedShell>(
          world, g0, mu_resp[static_cast<size_t>(b)],
          mu_resp[static_cast<size_t>(c)], mu_op[static_cast<size_t>(b)],
          mu_op[static_cast<size_t>(c)]);
      const double mx = inner(Xf.x_alpha, V.x_alpha);
      const double my = inner(Xf.y_alpha, V.y_alpha);
      // Diagnostic: both metric pairings. C_N derivation (TPA_SCOPING §5m)
      // forces + if the ES solver's y-sign convention matches the FD one that
      // validated beta's b1; the M− column decides that empirically (the
      // ~6% √2-vs-measured residual is exactly y-pairing sized).
      if (world.rank() == 0)
        printf("  [TPA residue] pair %c%c (%d/6): M+ = %+.6f   M- = %+.6f   "
               "(x-part %+.6f, y-part %+.6f)   [%.1f s]\n",
               "xyz"[b], "xyz"[c], b * (5 - b) / 2 + c + 1, mx + my, mx - my,
               mx, my, wall_time() - t0);
               fflush(stdout);
      S(b, c) = prefactor * (mx + my);
      S(c, b) = S(b, c);
    }
  }
  return S;
}

/// ONE-ELECTRON part of the single-residue contraction: the mu-operator terms
/// of V^{bc} only (vbc_half_spec's [A] = -Q(v*c) and [M] = c*<phi|v|phi>,
/// both orderings) — no two-electron builds, so this costs inner products.
/// Validated against DALTON's A2B/A2C/X2 families (<=1.3% on every d-aug-QZ
/// ledger element). Production 2PA is S = C_N*(S_1e + S_E3corr), C_N=sqrt(2)
/// (verified 6/6 vs DALTON 2026-07-22, verify_e3_k6.py).
inline madness::Tensor<double>
tpa_moment_residue_1e(madness::World &world, const ResponseGroundState &g0,
                      const ResponseStateXY<ClosedShell> &Xf,
                      const std::array<ResponseStateXY<ClosedShell>, 3> &mu_resp,
                      const std::array<madness::real_function_3d, 3> &mu_op,
                      double prefactor = 1.0) {
  using namespace madness;
  const vecfuncT &phi = g0.amo;
  std::array<Tensor<double>, 3> mv;   // <phi_i | mu_a | phi_j>
  for (int a = 0; a < 3; ++a)
    mv[a] = matrix_inner(world, phi, mul(world, mu_op[a], phi, true));
  // one ordered half: v = mu_b acting on the C-state (mirrors vbc_half_spec)
  auto half = [&](int b, int c) {
    const auto &cx = mu_resp[static_cast<size_t>(c)].x_alpha;
    const auto &cy = mu_resp[static_cast<size_t>(c)].y_alpha;
    const auto &v = mu_op[static_cast<size_t>(b)];
    vecfuncT vx = g0.Qa(mul(world, v, cx, true)); scale(world, vx, -1.0);
    vecfuncT vy = g0.Qa(mul(world, v, cy, true)); scale(world, vy, -1.0);
    gaxpy(world, 1.0, vx, 1.0, transform(world, cx, mv[b], true));
    gaxpy(world, 1.0, vy, 1.0, transform(world, cy, mv[b], true));
    return inner(Xf.x_alpha, vx) + inner(Xf.y_alpha, vy);
  };
  Tensor<double> S(3L, 3L);
  for (int b = 0; b < 3; ++b)
    for (int c = b; c < 3; ++c)
      S(b, c) = S(c, b) = prefactor * (half(b, c) + half(c, b));
  return S;
}

/// Full Monson–McClain two-photon observables from the S tensor, matching
/// Dalton QRSMO (rspvec.F:2803-2821). `omega_ex` = excitation energy (a.u.);
/// both photons are at omega_ex/2. sigma in GM assumes Dalton's 0.1 eV FWHM.
/// The formulas are in madness/chem/DerivedQuantities.h
/// (two_photon_strengths, two_photon_cross_section_gm).
struct Observables {
  double Df, Dg, D_linear, D_circular, R, sigma_linear_gm, sigma_circular_gm;
};

inline Observables
observables(const madness::Tensor<double> &S, double omega_ex) {
  madness::derived::Tensor2 s{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) s[i][j] = S(i, j);
  const auto d = madness::derived::two_photon_strengths(s);
  Observables o;
  o.Df = d.Df;
  o.Dg = d.Dg;
  o.D_linear = d.D_linear;
  o.D_circular = d.D_circular;
  o.R = d.R;
  o.sigma_linear_gm = madness::derived::two_photon_cross_section_gm(o.D_linear, omega_ex);
  o.sigma_circular_gm = madness::derived::two_photon_cross_section_gm(o.D_circular, omega_ex);
  return o;
}

} // namespace molresponse_v3::tpa

#endif // MOLRESPONSE_V3_KERNELS_TPA_HPP
