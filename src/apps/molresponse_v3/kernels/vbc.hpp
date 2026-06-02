#ifndef MOLRESPONSE_V3_KERNELS_VBC_HPP
#define MOLRESPONSE_V3_KERNELS_VBC_HPP

// ===========================================================================
// VBC — the quadratic source for hyperpolarizability beta (2n+1 / VBC
// contraction), closed-shell. Ported from molresponse_v2 SimpleVBCComputer, but
// rebuilt on v3 primitives so the two-electron kernel `compute_g` shares the
// SAME Coulomb / exchange path as the linear kernel
// Kernels<Full, ClosedShell>::compute_gamma (doc 15, beta path). beta is pure
// 2n+1 contraction — no explicit x(2) solve.
//
// THIS HEADER (beta-i, foundation): compute_g + make_zeta only. compute_vbc_i /
// compute_vbc (which build the full quadratic source from two converged
// first-order states) land on top once these conventions are verified.
//
// Convention note (the whole correctness risk): v3 carries the closed-shell
// factor 2 IN the density (compute_density: rho1 = 2*sum phi0*(x+y)), so the
// Coulomb term is J[rho]*phi with no extra factor and exchange is scaled by
// c_xc (= 1 for HF). v2's compute_g instead used 2J - K on an unscaled density
// — identical for HF. We follow the v3 convention so compute_g is consistent
// with the already-Dalton-validated linear kernel.
// ===========================================================================

#include "tags.hpp"
#include "tda.hpp"   // ResponseGroundState, common_ops::{dot, apply_exchange}

#include <madness/mra/mra.h>

#include <utility>
#include <vector>

namespace molresponse_v3::vbc {

using vecfuncT = std::vector<madness::real_function_3d>;

/// (J[rho] - c_xc*K) applied to (phix, phiy), where the perturbed density
///   rho = 2*( sum_p Aleft[p]*Aright[p] + sum_p Bleft[p]*Bright[p] )
/// carries the closed-shell factor 2 (matching compute_density). J = coulop(rho)
/// thus needs no further factor; exchange is scaled by c_xc. The Y / conjugate
/// channel swaps the exchange bra/ket, exactly as v2 compute_g.
///
/// v2 mapping: result_x = 2*Jx - Kx, result_y = 2*Jy - Ky (HF). Here the 2 is
/// folded into rho and the K coefficient is c_xc (= 1 for HF).
inline std::pair<vecfuncT, vecfuncT>
compute_g(madness::World &world, const ResponseGroundState &g0,
          const vecfuncT &Aleft, const vecfuncT &Aright,
          const vecfuncT &Bleft, const vecfuncT &Bright,
          const vecfuncT &phix,  const vecfuncT &phiy) {
  using namespace madness;

  // rho = 2*(Aleft*Aright + Bleft*Bright)   [closed-shell factor 2 in density]
  real_function_3d rho = common_ops::dot(world, Aleft, Aright);
  rho += common_ops::dot(world, Bleft, Bright);
  rho.scale(2.0);
  rho.truncate();

  // Coulomb: J = coulop(rho), applied to each phix / phiy.
  auto J  = apply(*g0.coulop, rho);
  auto Jx = mul(world, J, phix, false);
  auto Jy = mul(world, J, phiy, false);
  world.gop.fence();

  // Exchange (HF/hybrid): K(bra, ket)(apply_to). X channel uses (Aleft,Aright)
  // and (Bleft,Bright); Y channel swaps bra/ket (conjugate), as in v2.
  auto Kx = common_ops::apply_exchange(world, Aleft, Aright, phix, g0.lo);
  {
    auto k2 = common_ops::apply_exchange(world, Bleft, Bright, phix, g0.lo);
    gaxpy(world, 1.0, Kx, 1.0, k2);
  }
  auto Ky = common_ops::apply_exchange(world, Aright, Aleft, phiy, g0.lo);
  {
    auto k2 = common_ops::apply_exchange(world, Bright, Bleft, phiy, g0.lo);
    gaxpy(world, 1.0, Ky, 1.0, k2);
  }
  world.gop.fence();

  // g = J - c_xc*K   (J already carries the factor 2 via rho).
  vecfuncT gx = Jx;
  gaxpy(world, 1.0, gx, -g0.c_xc, Kx);
  vecfuncT gy = Jy;
  gaxpy(world, 1.0, gy, -g0.c_xc, Ky);
  return {std::move(gx), std::move(gy)};
}

/// Occupied-space relaxation term:
///   zeta_BC[p] = - sum_q phi0[q] <y_B[q] | x_C[p]>.
/// v2 make_zeta_bc: -1 * transform(phi0, <y_B|x_C>). Cheap (one matrix_inner +
/// transform); negligible vs. compute_g.
inline vecfuncT
make_zeta(madness::World &world, const vecfuncT &by, const vecfuncT &cx,
          const vecfuncT &phi0) {
  using namespace madness;
  auto mat = matrix_inner(world, by, cx);
  mat.scale(-1.0);
  return transform(world, phi0, mat, true);
}

} // namespace molresponse_v3::vbc

#endif // MOLRESPONSE_V3_KERNELS_VBC_HPP
