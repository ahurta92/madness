#ifndef MOLRESPONSE_V3_KERNELS_TWO_ELECTRON_HPP
#define MOLRESPONSE_V3_KERNELS_TWO_ELECTRON_HPP

// ===========================================================================
// Unified two-electron kernel (doc 15 / refactor): ONE response-state-level
// `compute_g(S1, S2, S3)` that both the linear gamma and the VBC quadratic
// source are special cases of.
//
//   gamma involves three response states: the first two (S1, S2) build the
//   density AND supply the exchange bra/kets; the third (S3) is the state we
//   apply to. For linear Full closed-shell, compute_gamma(chi) == compute_g(
//   chi, Phi, Phi) with Phi = the ground state dressed as {phi, phi}:
//     rho   = 2 * ( chi.x*phi + chi.y*phi )            (rho_factor = 2)
//     J     = coulop(rho) applied to phi
//     K.x   = K(phi, chi.x)(phi) + K(chi.y, phi)(phi)
//     K.y   = K(phi, chi.y)(phi) + K(chi.x, phi)(phi)
//   result = Q( J - c_xc*K ).
//
// THIS HEADER (step 1): the Full / XY kernel (ResponseStateXY<Shell>), generic
// over Shell (closed + open). It reproduces Kernels<Full,Shell>::compute_gamma
// term-for-term (so the equivalence test pins it), and returns the density it
// built so the linear iteration reuses it for the Δρ check. TDA/Static (X-only)
// kernels follow in a later step.
//
// Exchange bra/ket order is taken VERBATIM from kernels/full.hpp::compute_gamma
// (Dalton-validated): apply_exchange(bra, ket, apply_to, lo).
// ===========================================================================

#include "tags.hpp"
#include "tda.hpp"   // ResponseGroundState, common_ops::{dot, apply_exchange}
#include "../solvers/response_state.hpp"   // ResponseStateXY<Shell>

#include <madness/mra/mra.h>

#include <type_traits>
#include <utility>

namespace molresponse_v3::two_electron {

/// Unified Full/XY two-electron kernel. S1, S2 form the density + exchange
/// bra/kets; S3 is applied to. Returns (result = Q(J - c_xc*K) applied to S3,
/// rho = the perturbed density built from S1*S2). `rho_factor` is the
/// closed-shell spin factor (2) or 1 per-spin for open shell.
template <typename Shell>
inline std::pair<ResponseStateXY<Shell>, madness::real_function_3d>
compute_g(madness::World &world, const ResponseGroundState &g0,
          const ResponseStateXY<Shell> &S1,
          const ResponseStateXY<Shell> &S2,
          const ResponseStateXY<Shell> &S3,
          double rho_factor) {
  using namespace madness;

  // --- density: rho_factor * sum over channels (S1.c * S2.c) --------------
  real_function_3d rho = common_ops::dot(world, S1.x_alpha, S2.x_alpha);
  rho += common_ops::dot(world, S1.y_alpha, S2.y_alpha);
  if constexpr (std::is_same_v<Shell, OpenShell>) {
    rho += common_ops::dot(world, S1.x_beta, S2.x_beta);
    rho += common_ops::dot(world, S1.y_beta, S2.y_beta);
  }
  rho.scale(rho_factor);
  rho.truncate();

  auto J = apply(*g0.coulop, rho);

  ResponseStateXY<Shell> out;

  // --- alpha block: J*S3 + cross-channel exchange (matches full.hpp) -------
  out.x_alpha = mul(world, J, S3.x_alpha, true);
  out.y_alpha = mul(world, J, S3.y_alpha, true);
  if (g0.c_xc > 0.0) {
    // out.x: K(S2.x, S1.x)(S3.x) + K(S1.y, S2.y)(S3.x)
    auto kx1 = common_ops::apply_exchange(world, S2.x_alpha, S1.x_alpha, S3.x_alpha, g0.lo);
    gaxpy(world, 1.0, out.x_alpha, -g0.c_xc, kx1);
    auto kx2 = common_ops::apply_exchange(world, S1.y_alpha, S2.y_alpha, S3.x_alpha, g0.lo);
    gaxpy(world, 1.0, out.x_alpha, -g0.c_xc, kx2);
    // out.y: K(S2.y, S1.y)(S3.y) + K(S1.x, S2.x)(S3.y)
    auto ky1 = common_ops::apply_exchange(world, S2.y_alpha, S1.y_alpha, S3.y_alpha, g0.lo);
    gaxpy(world, 1.0, out.y_alpha, -g0.c_xc, ky1);
    auto ky2 = common_ops::apply_exchange(world, S1.x_alpha, S2.x_alpha, S3.y_alpha, g0.lo);
    gaxpy(world, 1.0, out.y_alpha, -g0.c_xc, ky2);
  }
  out.x_alpha = g0.Qa(out.x_alpha);
  out.y_alpha = g0.Qa(out.y_alpha);
  truncate(world, out.x_alpha);
  truncate(world, out.y_alpha);

  // --- beta block (open shell only): same structure, bmo channels + Qb -----
  if constexpr (std::is_same_v<Shell, OpenShell>) {
    out.x_beta = mul(world, J, S3.x_beta, true);
    out.y_beta = mul(world, J, S3.y_beta, true);
    if (g0.c_xc > 0.0) {
      auto kxb1 = common_ops::apply_exchange(world, S2.x_beta, S1.x_beta, S3.x_beta, g0.lo);
      gaxpy(world, 1.0, out.x_beta, -g0.c_xc, kxb1);
      auto kxb2 = common_ops::apply_exchange(world, S1.y_beta, S2.y_beta, S3.x_beta, g0.lo);
      gaxpy(world, 1.0, out.x_beta, -g0.c_xc, kxb2);
      auto kyb1 = common_ops::apply_exchange(world, S2.y_beta, S1.y_beta, S3.y_beta, g0.lo);
      gaxpy(world, 1.0, out.y_beta, -g0.c_xc, kyb1);
      auto kyb2 = common_ops::apply_exchange(world, S1.x_beta, S2.x_beta, S3.y_beta, g0.lo);
      gaxpy(world, 1.0, out.y_beta, -g0.c_xc, kyb2);
    }
    out.x_beta = g0.Qb(out.x_beta);
    out.y_beta = g0.Qb(out.y_beta);
    truncate(world, out.x_beta);
    truncate(world, out.y_beta);
  }

  return {std::move(out), std::move(rho)};
}

} // namespace molresponse_v3::two_electron

#endif // MOLRESPONSE_V3_KERNELS_TWO_ELECTRON_HPP
