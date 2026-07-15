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

namespace molresponse_v3::tpa {

/// The 3x3 two-photon transition-moment tensor S_ab for one excited state f
/// (closed-shell / TDA). Inputs are already-converged states:
///   Xf        — the ES eigenvector for root f, wrapped as XY (TDA: y_alpha=0).
///   mu_resp   — the dipole response mu_a at omega_f/2, per Cartesian axis a.
///   mu_op     — the raw dipole operators (MomentFunctor), per axis.
/// No electronic solve: pure contraction via beta_abc + compute_vbc.
inline madness::Tensor<double>
tpa_moment(madness::World &world, const ResponseGroundState &g0,
           const ResponseStateXY<ClosedShell> &Xf,
           const std::array<ResponseStateXY<ClosedShell>, 3> &mu_resp,
           const std::array<madness::real_function_3d, 3> &mu_op) {
  using namespace madness;

  // Homogeneous C-channel: VC_op = 0 (a valid zero function of the right shape).
  real_function_3d zero_op = copy(mu_op[0]);
  zero_op.scale(0.0);

  // VBC source depends only on the B(=mu_b) / C(=X_f) pair, not on A — build the
  // three (one per B-axis) once and reuse across A.
  std::array<ResponseStateXY<ClosedShell>, 3> vbc_b;
  for (int b = 0; b < 3; ++b)
    vbc_b[b] = vbc::compute_vbc<ClosedShell>(world, g0, mu_resp[b], Xf,
                                             mu_op[b], zero_op);

  Tensor<double> S(3L, 3L);
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      S(a, b) = beta::beta_abc<ClosedShell>(world, g0, mu_resp[a], vbc_b[b],
                                            mu_resp[b], Xf, mu_op[a]);
  return S;
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
  return d;
}

} // namespace molresponse_v3::tpa

#endif // MOLRESPONSE_V3_KERNELS_TPA_HPP
