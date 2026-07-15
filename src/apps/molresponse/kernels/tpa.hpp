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
// This is a CANDIDATE residue (the "VC_op=0 homogeneous C-channel" reduction);
// the exact surviving-term set + the X_f normalization are pinned numerically
// against the Dalton .QUADRA reference (refs/dalton_tpa.json). Because it reuses
// beta_abc/compute_vbc unchanged, a mismatch is a PHYSICS/normalization signal,
// not a contraction-code bug (alpha/beta/ES all already validate).
//
// The rotationally-invariant strength delta^|| (parallel linear polarization,
// F=G=H=2) is the quantity Dalton reports (two_photon_strengths); S_ab is the
// raw tensor.
// ===========================================================================

#include "beta.hpp"   // beta::beta_abc
#include "vbc.hpp"    // vbc::compute_vbc, ResponseGroundState, ResponseStateXY

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

/// Rotationally-invariant two-photon strength for parallel linear polarization:
///   delta^|| = Sum_ab ( F S_aa S_bb + G S_ab S_ab + H S_ab S_ba ), F=G=H=2.
/// (Real S here — TDA, undamped, off-resonance.) This is what Dalton's
/// two_photon_strengths reports, so it is the validation target.
inline double
delta_parallel(const madness::Tensor<double> &S) {
  double d = 0.0;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      d += 2.0 * S(a, a) * S(b, b) + 2.0 * S(a, b) * S(a, b) +
           2.0 * S(a, b) * S(b, a);
  // 1/30 = the isotropic average over molecular orientations (Monson–McClain,
  // JCP 53:29 1970). Without it delta is ~30x too large. Matches Dalton's
  // two_photon_strengths (D = 2*Df + 4*Dg for the symmetric tensor).
  return d / 30.0;
}

/// Full Monson–McClain two-photon observables from the S tensor, matching
/// Dalton QRSMO (rspvec.F:2803-2821). `omega_ex` = excitation energy (a.u.);
/// both photons are at omega_ex/2. sigma in GM assumes Dalton's 0.1 eV FWHM.
struct Observables {
  double Df, Dg, D_linear, D_circular, R, sigma_linear_gm, sigma_circular_gm;
};

inline Observables
observables(const madness::Tensor<double> &S, double omega_ex) {
  double df = 0.0, dg = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      df += S(i, i) * S(j, j);   // Df = sum_ij S_ii S_jj / 30
      dg += S(i, j) * S(i, j);   // Dg = sum_ij S_ij^2   / 30
    }
  Observables o;
  o.Df = df / 30.0;
  o.Dg = dg / 30.0;
  o.D_linear   =  2.0 * o.Df + 4.0 * o.Dg;
  o.D_circular = -2.0 * o.Df + 6.0 * o.Dg;
  const double den = o.Df + 2.0 * o.Dg;
  o.R = (den != 0.0) ? (-o.Df + 3.0 * o.Dg) / den : 0.0;
  // a.u. -> GM (Dalton AU_TO_GM, rspvec.F:1795-6): 8*pi^2 * alpha * a0[pm]^5 /
  // (c[cm/s] * FWHM[au]); FWHM = 0.1 eV baked in. Numerically ~2.170.
  constexpr double PI      = 3.14159265358979324;
  constexpr double ALPHA   = 7.2973525693e-3;   // fine-structure constant
  constexpr double A0_PM   = 52.9177210903;     // Bohr radius in pm
  constexpr double C_CMS   = 2.99792458e10;     // speed of light, cm/s
  constexpr double FWHM_AU = 0.0036749326;      // 0.1 eV in a.u.
  const double AU_TO_GM = 8.0 * PI * PI * ALPHA * std::pow(A0_PM, 5.0) /
                          (C_CMS * FWHM_AU);
  const double ph2 = (0.5 * omega_ex) * (0.5 * omega_ex);
  o.sigma_linear_gm   = o.D_linear   * ph2 * AU_TO_GM;
  o.sigma_circular_gm = o.D_circular * ph2 * AU_TO_GM;
  return o;
}

} // namespace molresponse_v3::tpa

#endif // MOLRESPONSE_V3_KERNELS_TPA_HPP
