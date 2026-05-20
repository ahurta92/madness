#ifndef MOLRESPONSE_V3_KERNELS_FULL_HPP
#define MOLRESPONSE_V3_KERNELS_FULL_HPP

// =========================================================================
// Kernels<Full, ClosedShell>
//
// Frequency-dependent response with X and Y amplitudes (omega != 0).
// Closed-shell, so Storage = ResponseStateXY<ClosedShell> ({x_alpha,
// y_alpha}). The response density couples X and Y:
//
//     rho1 = sum_p phi0[p] * (x_alpha[p] + y_alpha[p])
//
// (closed shell, so factor of 1 — the legacy spin-summed factor of 2
// is absorbed into the source term and into the BSH apply on -2·theta;
// see compute_density below for the exact convention used here, which
// matches molresponse_v3/ResponseKernel.hpp::compute_rho1).
//
// All per-Storage methods return Storage. The BSH apply is paired:
// theta.x_alpha is convolved with BSH(+omega) onto new x_alpha;
// theta.y_alpha is convolved with BSH(-omega) onto new y_alpha. Both
// share the same orbital-energy shift logic.
//
// Open-shell Full will reuse the same algorithmic shape with the
// alpha/beta blocks duplicated — defer until Phase 3 of the roadmap.
// =========================================================================

#include "tags.hpp"
#include "tda.hpp"          // ResponseGroundState, detail_tda::*
#include "../solvers/response_state.hpp"

#include <madness/chem/SCFOperators.h>
#include <madness/chem/projector.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

namespace molresponse_v3 {

template <>
struct Kernels<Full, ClosedShell> {

  using vecfuncT = std::vector<madness::real_function_3d>;
  using State    = ResponseStateXY<ClosedShell>;

  /// Singlet closed-shell Full response density:
  ///   rho1 = 2 * sum_p phi0[p] * (x_alpha[p] + y_alpha[p])
  ///        = spin × (X + Y contributions) × Σφ·(x+y)
  ///
  /// Factor 2 (spin sum) matches the Coulomb coefficient 4(pa|qb) on
  /// the X+Y combination in the singlet (A+B)·(X+Y) sub-block: in
  /// real space, J[2·Σφ(x+y)]·φ = 2·J[Σφ(x+y)]·φ pairs with the
  /// matrix-element 4(pa|qb)(X+Y).
  ///
  /// Self-consistency at Y=X: rho1|Y=X = 2·Σφ·(x+x) = 4·Σφ·x =
  /// `Kernels<Static, ClosedShell>::compute_density`.
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &target,
                  const State    &state) {
    auto sum_xy = madness::copy(world, state.x_alpha);
    gaxpy(world, 1.0, sum_xy, 1.0, state.y_alpha);
    auto rho1 = dot(world, target.phi0_alpha, sum_xy);
    rho1.scale(2.0);
    rho1.truncate();
    return rho1;
  }

  /// gamma_X = Q( J[rho1]*phi0 - c_xc * (K[phi0, x](phi0) + K[y, phi0](phi0)) )
  /// gamma_Y = Q( J[rho1]*phi0 - c_xc * (K[phi0, y](phi0) + K[x, phi0](phi0)) )
  ///
  /// Cross-channel exchange is what couples X and Y in the RPA equations
  /// — without it the result reduces to TDA on each side. Mirrors
  /// legacy ResponseKernel.hpp::compute_gamma Full case.
  /// The Coulomb piece J[rho1]*phi0 is the SAME on both sides.
  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &target,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*target.coulop, rho1);

    // Coulomb piece — same for X and Y.
    auto gamma_x = mul(world, J_rho, target.phi0_alpha, true);
    auto gamma_y = mul(world, J_rho, target.phi0_alpha, true);

    if (target.c_xc > 0.0) {
      // --- X-channel exchange: K[phi0, x](phi0) + K[y, phi0](phi0) -------
      madness::Exchange<double, 3> Kpx(world, target.lo);
      Kpx.set_bra_and_ket(target.phi0_alpha, state.x_alpha);
      Kpx.set_algorithm(madness::Exchange<double, 3>::
                            ExchangeAlgorithm::multiworld_efficient_row);
      auto Kpx_phi = Kpx(target.phi0_alpha);
      gaxpy(world, 1.0, gamma_x, -target.c_xc, Kpx_phi);

      madness::Exchange<double, 3> Kyp(world, target.lo);
      Kyp.set_bra_and_ket(state.y_alpha, target.phi0_alpha);
      Kyp.set_algorithm(madness::Exchange<double, 3>::
                            ExchangeAlgorithm::multiworld_efficient_row);
      auto Kyp_phi = Kyp(target.phi0_alpha);
      gaxpy(world, 1.0, gamma_x, -target.c_xc, Kyp_phi);

      // --- Y-channel exchange: K[phi0, y](phi0) + K[x, phi0](phi0) -------
      madness::Exchange<double, 3> Kpy(world, target.lo);
      Kpy.set_bra_and_ket(target.phi0_alpha, state.y_alpha);
      Kpy.set_algorithm(madness::Exchange<double, 3>::
                            ExchangeAlgorithm::multiworld_efficient_row);
      auto Kpy_phi = Kpy(target.phi0_alpha);
      gaxpy(world, 1.0, gamma_y, -target.c_xc, Kpy_phi);

      madness::Exchange<double, 3> Kxp(world, target.lo);
      Kxp.set_bra_and_ket(state.x_alpha, target.phi0_alpha);
      Kxp.set_algorithm(madness::Exchange<double, 3>::
                            ExchangeAlgorithm::multiworld_efficient_row);
      auto Kxp_phi = Kxp(target.phi0_alpha);
      gaxpy(world, 1.0, gamma_y, -target.c_xc, Kxp_phi);
    }

    gamma_x = target.Q_alpha(gamma_x);
    gamma_y = target.Q_alpha(gamma_y);
    truncate(world, gamma_x);
    truncate(world, gamma_y);
    return State{std::move(gamma_x), std::move(gamma_y)};
  }

  /// V0·x and V0·y. V_local and K0 each act on x and y independently.
  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &target,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Vx = mul_sparse(world, target.V_local_alpha, state.x_alpha, vtol);
    auto Vy = mul_sparse(world, target.V_local_alpha, state.y_alpha, vtol);

    if (target.c_xc > 0.0) {
      madness::Exchange<double, 3> K0(world, target.lo);
      K0.set_bra_and_ket(target.phi0_alpha, target.phi0_alpha);
      K0.set_algorithm(madness::Exchange<double, 3>::
                           ExchangeAlgorithm::multiworld_efficient_row);
      auto k0x = K0(state.x_alpha);
      auto k0y = K0(state.y_alpha);
      gaxpy(world, 1.0, Vx, -target.c_xc, k0x);
      gaxpy(world, 1.0, Vy, -target.c_xc, k0y);
    }
    return State{std::move(Vx), std::move(Vy)};
  }

  /// E0·x and E0·y with the off-diagonal Fock matrix (diag eps absorbed
  /// into BSH).
  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &target,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Ex = transform(world, state.x_alpha,
                        target.fock_no_diag_alpha, vtol, true);
    auto Ey = transform(world, state.y_alpha,
                        target.fock_no_diag_alpha, vtol, true);
    return State{std::move(Ex), std::move(Ey)};
  }

  /// Theta = V0x - E0x + gamma   (both x and y components).
  static State
  compute_theta(madness::World &world,
                const State    &V0x,
                const State    &E0x,
                const State    &gamma) {
    auto Tx = madness::copy(world, V0x.x_alpha);
    gaxpy(world, 1.0, Tx, -1.0, E0x.x_alpha);
    gaxpy(world, 1.0, Tx,  1.0, gamma.x_alpha);
    truncate(world, Tx);

    auto Ty = madness::copy(world, V0x.y_alpha);
    gaxpy(world, 1.0, Ty, -1.0, E0x.y_alpha);
    gaxpy(world, 1.0, Ty,  1.0, gamma.y_alpha);
    truncate(world, Ty);

    return State{std::move(Tx), std::move(Ty)};
  }

  /// Paired BSH apply.
  ///   new_x = Q( BSH(+omega) * (-2 * (theta_x + shift_+ * x_alpha)) )
  ///   new_y = Q( BSH(-omega) * (-2 * (theta_y + shift_- * y_alpha)) )
  /// Each side computes its own shift independently because the LUMO
  /// + omega vs LUMO - omega check produces different shifts.
  static State
  bsh_apply(madness::World &world,
            const ResponseGroundState &target,
            const State    &state,
            const State    &theta,
            double          omega) {
    // --- X side: BSH at +omega ---------------------------------------
    const double shift_p = detail_tda::bsh_shift(target.eps_alpha, omega);
    auto rhs_x = madness::copy(world, theta.x_alpha);
    gaxpy(world, 1.0, rhs_x, shift_p, state.x_alpha);
    scale(world, rhs_x, -2.0);
    truncate(world, rhs_x);
    auto bsh_p = detail_tda::make_bsh_operators(world, target.eps_alpha,
                                                omega, target.lo);
    auto new_x = apply(world, bsh_p, rhs_x);
    new_x      = target.Q_alpha(new_x);
    truncate(world, new_x);

    // --- Y side: BSH at -omega ---------------------------------------
    const double shift_m = detail_tda::bsh_shift(target.eps_alpha, -omega);
    auto rhs_y = madness::copy(world, theta.y_alpha);
    gaxpy(world, 1.0, rhs_y, shift_m, state.y_alpha);
    scale(world, rhs_y, -2.0);
    truncate(world, rhs_y);
    auto bsh_m = detail_tda::make_bsh_operators(world, target.eps_alpha,
                                                -omega, target.lo);
    auto new_y = apply(world, bsh_m, rhs_y);
    new_y      = target.Q_alpha(new_y);
    truncate(world, new_y);

    return State{std::move(new_x), std::move(new_y)};
  }

  /// || (x_old, y_old) - (x_new, y_new) ||, combined norm.
  static double
  compute_residual_norm(madness::World &world,
                        const State    &state_old,
                        const State    &state_new) {
    auto dx = sub(world, state_old.x_alpha, state_new.x_alpha);
    auto dy = sub(world, state_old.y_alpha, state_new.y_alpha);
    const double nx = madness::norm2(world, dx);
    const double ny = madness::norm2(world, dy);
    return std::sqrt(nx * nx + ny * ny);
  }
};

// =========================================================================
// Kernels<Full, OpenShell>
//
// UHF dynamic response. Storage = ResponseStateXY<OpenShell>
// ({x_α, y_α, x_β, y_β}). Shared Coulomb across spins, per-spin
// exchange + BSH. Cross-channel exchange couples X with Y same as
// the closed-shell Full kernel.
//
// Density convention (legacy v3 unrestricted Full):
//     rho = 1 · 1 · (Σ (x_α + y_α)·φ_α + Σ (x_β + y_β)·φ_β)
//         = Σ_spin Σ_p (x_p + y_p) · φ_p
//
// (No spin or X+Y prefactor — for unrestricted, spin_factor=1; and
// for Full y_factor=1.)
// =========================================================================
template <>
struct Kernels<Full, OpenShell> {

  using vecfuncT = std::vector<madness::real_function_3d>;
  using State    = ResponseStateXY<OpenShell>;

  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &target,
                  const State    &state) {
    // alpha (x+y)·phi
    auto sum_xy_a = madness::copy(world, state.x_alpha);
    gaxpy(world, 1.0, sum_xy_a, 1.0, state.y_alpha);
    auto rho_a = dot(world, target.phi0_alpha, sum_xy_a);
    // beta (x+y)·phi
    auto sum_xy_b = madness::copy(world, state.x_beta);
    gaxpy(world, 1.0, sum_xy_b, 1.0, state.y_beta);
    auto rho_b = dot(world, target.phi0_beta, sum_xy_b);
    auto rho = rho_a + rho_b;
    rho.truncate();
    return rho;
  }

  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &target,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*target.coulop, rho1);

    // --- alpha block ---
    auto gamma_ax = mul(world, J_rho, target.phi0_alpha, true);
    auto gamma_ay = mul(world, J_rho, target.phi0_alpha, true);
    if (target.c_xc > 0.0) {
      // K[φα, xα](φα)
      madness::Exchange<double, 3> Kpx_a(world, target.lo);
      Kpx_a.set_bra_and_ket(target.phi0_alpha, state.x_alpha);
      Kpx_a.set_algorithm(madness::Exchange<double, 3>::
                              ExchangeAlgorithm::multiworld_efficient_row);
      auto Kpx_a_phi = Kpx_a(target.phi0_alpha);
      gaxpy(world, 1.0, gamma_ax, -target.c_xc, Kpx_a_phi);

      // K[yα, φα](φα) — cross-channel for X
      madness::Exchange<double, 3> Kyp_a(world, target.lo);
      Kyp_a.set_bra_and_ket(state.y_alpha, target.phi0_alpha);
      Kyp_a.set_algorithm(madness::Exchange<double, 3>::
                              ExchangeAlgorithm::multiworld_efficient_row);
      auto Kyp_a_phi = Kyp_a(target.phi0_alpha);
      gaxpy(world, 1.0, gamma_ax, -target.c_xc, Kyp_a_phi);

      // K[φα, yα](φα)
      madness::Exchange<double, 3> Kpy_a(world, target.lo);
      Kpy_a.set_bra_and_ket(target.phi0_alpha, state.y_alpha);
      Kpy_a.set_algorithm(madness::Exchange<double, 3>::
                              ExchangeAlgorithm::multiworld_efficient_row);
      auto Kpy_a_phi = Kpy_a(target.phi0_alpha);
      gaxpy(world, 1.0, gamma_ay, -target.c_xc, Kpy_a_phi);

      // K[xα, φα](φα) — cross-channel for Y
      madness::Exchange<double, 3> Kxp_a(world, target.lo);
      Kxp_a.set_bra_and_ket(state.x_alpha, target.phi0_alpha);
      Kxp_a.set_algorithm(madness::Exchange<double, 3>::
                              ExchangeAlgorithm::multiworld_efficient_row);
      auto Kxp_a_phi = Kxp_a(target.phi0_alpha);
      gaxpy(world, 1.0, gamma_ay, -target.c_xc, Kxp_a_phi);
    }
    gamma_ax = target.Q_alpha(gamma_ax);
    gamma_ay = target.Q_alpha(gamma_ay);
    truncate(world, gamma_ax);
    truncate(world, gamma_ay);

    // --- beta block ---
    auto gamma_bx = mul(world, J_rho, target.phi0_beta, true);
    auto gamma_by = mul(world, J_rho, target.phi0_beta, true);
    if (target.c_xc > 0.0) {
      madness::Exchange<double, 3> Kpx_b(world, target.lo);
      Kpx_b.set_bra_and_ket(target.phi0_beta, state.x_beta);
      Kpx_b.set_algorithm(madness::Exchange<double, 3>::
                              ExchangeAlgorithm::multiworld_efficient_row);
      auto Kpx_b_phi = Kpx_b(target.phi0_beta);
      gaxpy(world, 1.0, gamma_bx, -target.c_xc, Kpx_b_phi);

      madness::Exchange<double, 3> Kyp_b(world, target.lo);
      Kyp_b.set_bra_and_ket(state.y_beta, target.phi0_beta);
      Kyp_b.set_algorithm(madness::Exchange<double, 3>::
                              ExchangeAlgorithm::multiworld_efficient_row);
      auto Kyp_b_phi = Kyp_b(target.phi0_beta);
      gaxpy(world, 1.0, gamma_bx, -target.c_xc, Kyp_b_phi);

      madness::Exchange<double, 3> Kpy_b(world, target.lo);
      Kpy_b.set_bra_and_ket(target.phi0_beta, state.y_beta);
      Kpy_b.set_algorithm(madness::Exchange<double, 3>::
                              ExchangeAlgorithm::multiworld_efficient_row);
      auto Kpy_b_phi = Kpy_b(target.phi0_beta);
      gaxpy(world, 1.0, gamma_by, -target.c_xc, Kpy_b_phi);

      madness::Exchange<double, 3> Kxp_b(world, target.lo);
      Kxp_b.set_bra_and_ket(state.x_beta, target.phi0_beta);
      Kxp_b.set_algorithm(madness::Exchange<double, 3>::
                              ExchangeAlgorithm::multiworld_efficient_row);
      auto Kxp_b_phi = Kxp_b(target.phi0_beta);
      gaxpy(world, 1.0, gamma_by, -target.c_xc, Kxp_b_phi);
    }
    gamma_bx = target.Q_beta(gamma_bx);
    gamma_by = target.Q_beta(gamma_by);
    truncate(world, gamma_bx);
    truncate(world, gamma_by);

    return State{std::move(gamma_ax), std::move(gamma_ay),
                 std::move(gamma_bx), std::move(gamma_by)};
  }

  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &target,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Vxa = mul_sparse(world, target.V_local_alpha, state.x_alpha, vtol);
    auto Vya = mul_sparse(world, target.V_local_alpha, state.y_alpha, vtol);
    auto Vxb = mul_sparse(world, target.V_local_beta,  state.x_beta,  vtol);
    auto Vyb = mul_sparse(world, target.V_local_beta,  state.y_beta,  vtol);
    if (target.c_xc > 0.0) {
      madness::Exchange<double, 3> K0_a(world, target.lo);
      K0_a.set_bra_and_ket(target.phi0_alpha, target.phi0_alpha);
      K0_a.set_algorithm(madness::Exchange<double, 3>::
                             ExchangeAlgorithm::multiworld_efficient_row);
      auto k0_ax = K0_a(state.x_alpha);
      auto k0_ay = K0_a(state.y_alpha);
      gaxpy(world, 1.0, Vxa, -target.c_xc, k0_ax);
      gaxpy(world, 1.0, Vya, -target.c_xc, k0_ay);

      madness::Exchange<double, 3> K0_b(world, target.lo);
      K0_b.set_bra_and_ket(target.phi0_beta, target.phi0_beta);
      K0_b.set_algorithm(madness::Exchange<double, 3>::
                             ExchangeAlgorithm::multiworld_efficient_row);
      auto k0_bx = K0_b(state.x_beta);
      auto k0_by = K0_b(state.y_beta);
      gaxpy(world, 1.0, Vxb, -target.c_xc, k0_bx);
      gaxpy(world, 1.0, Vyb, -target.c_xc, k0_by);
    }
    return State{std::move(Vxa), std::move(Vya),
                 std::move(Vxb), std::move(Vyb)};
  }

  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &target,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Exa = transform(world, state.x_alpha,
                         target.fock_no_diag_alpha, vtol, true);
    auto Eya = transform(world, state.y_alpha,
                         target.fock_no_diag_alpha, vtol, true);
    auto Exb = transform(world, state.x_beta,
                         target.fock_no_diag_beta,  vtol, true);
    auto Eyb = transform(world, state.y_beta,
                         target.fock_no_diag_beta,  vtol, true);
    return State{std::move(Exa), std::move(Eya),
                 std::move(Exb), std::move(Eyb)};
  }

  static State
  compute_theta(madness::World &world,
                const State    &V0x,
                const State    &E0x,
                const State    &gamma) {
    auto theta_axy = [&](const vecfuncT &Vx, const vecfuncT &Ex,
                         const vecfuncT &g) {
      auto T = madness::copy(world, Vx);
      gaxpy(world, 1.0, T, -1.0, Ex);
      gaxpy(world, 1.0, T,  1.0, g);
      truncate(world, T);
      return T;
    };
    auto Txa = theta_axy(V0x.x_alpha, E0x.x_alpha, gamma.x_alpha);
    auto Tya = theta_axy(V0x.y_alpha, E0x.y_alpha, gamma.y_alpha);
    auto Txb = theta_axy(V0x.x_beta,  E0x.x_beta,  gamma.x_beta);
    auto Tyb = theta_axy(V0x.y_beta,  E0x.y_beta,  gamma.y_beta);
    return State{std::move(Txa), std::move(Tya),
                 std::move(Txb), std::move(Tyb)};
  }

  static State
  bsh_apply(madness::World &world,
            const ResponseGroundState &target,
            const State    &state,
            const State    &theta,
            double          omega) {
    auto apply_side = [&](const madness::Tensor<double> &eps,
                          const madness::QProjector<double, 3> &Q,
                          const vecfuncT &theta_v,
                          const vecfuncT &x_v,
                          double signed_omega) {
      const double shift = detail_tda::bsh_shift(eps, signed_omega);
      auto rhs = madness::copy(world, theta_v);
      gaxpy(world, 1.0, rhs, shift, x_v);
      scale(world, rhs, -2.0);
      truncate(world, rhs);
      auto bsh = detail_tda::make_bsh_operators(world, eps, signed_omega,
                                                target.lo);
      auto out = apply(world, bsh, rhs);
      out = Q(out);
      truncate(world, out);
      return out;
    };
    auto new_xa = apply_side(target.eps_alpha, target.Q_alpha,
                             theta.x_alpha, state.x_alpha,  omega);
    auto new_ya = apply_side(target.eps_alpha, target.Q_alpha,
                             theta.y_alpha, state.y_alpha, -omega);
    auto new_xb = apply_side(target.eps_beta,  target.Q_beta,
                             theta.x_beta,  state.x_beta,   omega);
    auto new_yb = apply_side(target.eps_beta,  target.Q_beta,
                             theta.y_beta,  state.y_beta,  -omega);
    return State{std::move(new_xa), std::move(new_ya),
                 std::move(new_xb), std::move(new_yb)};
  }

  static double
  compute_residual_norm(madness::World &world,
                        const State    &state_old,
                        const State    &state_new) {
    auto dxa = sub(world, state_old.x_alpha, state_new.x_alpha);
    auto dya = sub(world, state_old.y_alpha, state_new.y_alpha);
    auto dxb = sub(world, state_old.x_beta,  state_new.x_beta);
    auto dyb = sub(world, state_old.y_beta,  state_new.y_beta);
    const double na  = madness::norm2(world, dxa);
    const double nya = madness::norm2(world, dya);
    const double nb  = madness::norm2(world, dxb);
    const double nyb = madness::norm2(world, dyb);
    return std::sqrt(na*na + nya*nya + nb*nb + nyb*nyb);
  }
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_KERNELS_FULL_HPP
