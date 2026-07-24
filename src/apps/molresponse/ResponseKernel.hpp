#ifndef MOLRESPONSE_V3_RESPONSEKERNEL_HPP
#define MOLRESPONSE_V3_RESPONSEKERNEL_HPP

// =========================================================================
// Legacy ES-guess support surface.
//
// This header once carried the full v2->v3 kernel port (compute_gamma,
// compute_lambda, make_bsh_operators, fd_iteration, the subspace builders,
// ...). Those functions have all been superseded by kernels/ (Kernels<Type,
// Shell>::compute_gamma / compute_V0x / ..., common_ops::make_bsh_operators /
// apply_kinetic, rs:: subspace ops) and had zero call sites, so they were
// removed. What remains is only the small surface still used by the legacy
// excited-state initial guess:
//
//   enum ResponseType        — density/exchange/BSH selector
//   bundle_inner             — metric inner product for a RealResponseState pair
//   orthonormalize_bundle    — Gram-Schmidt across a bundle of RealResponseStates
//   alpha_factor             — polarizability occupancy/y-channel prefactor
//
// Consumers: ESSolverGuess.hpp (orthonormalize_bundle, ResponseType) and
// tests/test_v3_fd_skeleton.cpp (alpha_factor, ResponseType).
// =========================================================================

#include "ResponseFunctions.hpp"  // RealResponseState

#include <madness/mra/mra.h>

#include <cmath>

namespace molresponse_v3 {

using namespace madness;

/// Response type determines density definition, exchange terms, BSH shifts.
/// The storage (ResponseState) is the same for all three — only the
/// algorithm differs.
enum class ResponseType { Static, Full, TDA };

// =========================================================================
// ES bundle helpers — shared by all ES stages (TDA RHF, Full RHF,
// TDA UHF, Full UHF). See docs/11_excited_state_solver_design.md.
// =========================================================================

/// Inner product between two states under the type-appropriate metric:
///   TDA / Static : ⟨a.x|b.x⟩  (positive-definite, identity)
///   Full         : ⟨a.x|b.x⟩ - ⟨a.y|b.y⟩  (symplectic, indefinite)
/// Spin: alpha + beta channels are summed when present.
inline double bundle_inner(World &world, ResponseType type,
                           const RealResponseState &a,
                           const RealResponseState &b) {
  double s = 0.0;
  if (!a.x_alpha.empty() && !b.x_alpha.empty())
    s += madness::inner(a.x_alpha, b.x_alpha);
  if (!a.x_beta.empty() && !b.x_beta.empty())
    s += madness::inner(a.x_beta, b.x_beta);
  if (type == ResponseType::Full) {
    if (!a.y_alpha.empty() && !b.y_alpha.empty())
      s -= madness::inner(a.y_alpha, b.y_alpha);
    if (!a.y_beta.empty() && !b.y_beta.empty())
      s -= madness::inner(a.y_beta, b.y_beta);
  }
  return s;
}

/// Gram-Schmidt orthonormalization of a bundle of response states under
/// the type-appropriate metric.
///
/// For TDA the metric is positive-definite; standard Gram-Schmidt works.
/// For Full the metric is indefinite (S = ⟨x|x⟩ - ⟨y|y⟩) and exact
/// metric-orthonormalization needs hyperbolic Gram-Schmidt. For now we
/// fall back to taking |⟨a|a⟩|^{1/2} as the norm — adequate when the
/// guess is already roughly metric-orthonormal (the legacy code does
/// this same pre-iter pass at TDDFT.cc:1667-1672).
inline void orthonormalize_bundle(World &world, ResponseType type,
                                  std::vector<RealResponseState> &X) {
  long n = static_cast<long>(X.size());
  if (n == 0)
    return;

  auto subtract = [&](RealResponseState &dst, double c,
                      const RealResponseState &src) {
    if (!src.x_alpha.empty())
      gaxpy(world, 1.0, dst.x_alpha, -c, src.x_alpha);
    if (!src.x_beta.empty())
      gaxpy(world, 1.0, dst.x_beta, -c, src.x_beta);
    if (!src.y_alpha.empty())
      gaxpy(world, 1.0, dst.y_alpha, -c, src.y_alpha);
    if (!src.y_beta.empty())
      gaxpy(world, 1.0, dst.y_beta, -c, src.y_beta);
  };
  auto rescale = [&](RealResponseState &s, double f) {
    if (!s.x_alpha.empty())
      scale(world, s.x_alpha, f);
    if (!s.x_beta.empty())
      scale(world, s.x_beta, f);
    if (!s.y_alpha.empty())
      scale(world, s.y_alpha, f);
    if (!s.y_beta.empty())
      scale(world, s.y_beta, f);
  };

  for (long i = 0; i < n; i++) {
    for (long j = 0; j < i; j++) {
      double c = bundle_inner(world, type, X[j], X[i]);
      subtract(X[i], c, X[j]);
    }
    double norm2_val = bundle_inner(world, type, X[i], X[i]);
    double n_abs = std::sqrt(std::abs(norm2_val));
    if (n_abs > 1e-12)
      rescale(X[i], 1.0 / n_abs);
  }
}

// =========================================================================
// Property computation helpers
// =========================================================================

/// Alpha factor for polarizability: alpha = factor * (ip_xa + ip_ya + ip_xb +
/// ip_yb)
///
/// The physical formula is:  α = -Σ_k n_k [⟨x_k|vp_k⟩ + ⟨y_k|vp_k⟩]
/// where n_k = 2 (restricted) or 1 (unrestricted per spin channel).
///
/// For Static/TDA, y is not solved (y=x implied) so ip_ya = ip_yb = 0
/// — the y_factor of 2 in af absorbs the missing y contribution.
/// Empirically validated <1% vs Dalton (d-aug-cc-pVQZ) on H2 / H2O /
/// C2H4 for Static ClosedShell, ~0.6% on Li UHF Static.
///
///   factor = -(spin_factor × y_factor)
///
///   spin_factor: 2 restricted (n_k=2), 1 unrestricted (n_k=1)
///   y_factor:    2 Static/TDA (y=x implied),
///                1 Full (y is solved separately)
///
///   Restricted Static:     -(2×2) = -4.0
///   Restricted Dynamic:    -(2×1) = -2.0
///   Unrestricted Static:   -(1×2) = -2.0
///   Unrestricted Dynamic:  -(1×1) = -1.0
///
/// Practical caveat: the density-scale factor inside `compute_density`
/// (e.g. 4 for Static ClosedShell, 2 for Full ClosedShell) enters the
/// Coulomb feedback inside the response kernel and so is part of the self-
/// consistent equation. The (density-factor, alpha-factor) pair has to
/// be chosen consistently — if you halve af without compensating
/// somewhere in the kernel, α will shift (confirmed empirically:
/// F=4 + af=-2 halved α on H2).
inline double alpha_factor(ResponseType type, bool restricted) {
  double spin_factor = restricted ? 2.0 : 1.0;
  double y_factor = (type == ResponseType::Full) ? 1.0 : 2.0;
  return -(spin_factor * y_factor);
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_RESPONSEKERNEL_HPP
