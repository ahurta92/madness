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
                  const ResponseGroundState &g0,
                  const State    &state) {
    auto sum_xy = madness::copy(world, state.x_alpha);
    gaxpy(world, 1.0, sum_xy, 1.0, state.y_alpha);
    auto rho1 = dot(world, g0.amo, sum_xy);
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
                const ResponseGroundState &g0,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*g0.coulop, rho1);

    // Coulomb piece — same for X and Y.
    auto gamma_x = mul(world, J_rho, g0.amo, true);
    auto gamma_y = mul(world, J_rho, g0.amo, true);

    if (g0.c_xc > 0.0) {
      // --- X-channel exchange: K[phi0, x](phi0) + K[y, phi0](phi0) -------
      auto Kpx_phi = detail_tda::apply_exchange(world, g0.amo, state.x_alpha, g0.amo, g0.lo);
      gaxpy(world, 1.0, gamma_x, -g0.c_xc, Kpx_phi);

      auto Kyp_phi = detail_tda::apply_exchange(world, state.y_alpha, g0.amo, g0.amo, g0.lo);
      gaxpy(world, 1.0, gamma_x, -g0.c_xc, Kyp_phi);

      // --- Y-channel exchange: K[phi0, y](phi0) + K[x, phi0](phi0) -------
      auto Kpy_phi = detail_tda::apply_exchange(world, g0.amo, state.y_alpha, g0.amo, g0.lo);
      gaxpy(world, 1.0, gamma_y, -g0.c_xc, Kpy_phi);

      auto Kxp_phi = detail_tda::apply_exchange(world, state.x_alpha, g0.amo, g0.amo, g0.lo);
      gaxpy(world, 1.0, gamma_y, -g0.c_xc, Kxp_phi);
    }

    gamma_x = g0.Qa(gamma_x);
    gamma_y = g0.Qa(gamma_y);
    truncate(world, gamma_x);
    truncate(world, gamma_y);
    return State{std::move(gamma_x), std::move(gamma_y)};
  }

  /// V0·x and V0·y. V_local and K0 each act on x and y independently.
  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Vx = mul_sparse(world, g0.V_local_alpha, state.x_alpha, vtol);
    auto Vy = mul_sparse(world, g0.V_local_alpha, state.y_alpha, vtol);

    if (g0.c_xc > 0.0) {
      auto k0x = detail_tda::apply_exchange(world, g0.amo, g0.amo, state.x_alpha, g0.lo);
      auto k0y = detail_tda::apply_exchange(world, g0.amo, g0.amo, state.y_alpha, g0.lo);
      gaxpy(world, 1.0, Vx, -g0.c_xc, k0x);
      gaxpy(world, 1.0, Vy, -g0.c_xc, k0y);
    }
    return State{std::move(Vx), std::move(Vy)};
  }

  /// E0·x and E0·y with the off-diagonal Fock matrix (diag eps absorbed
  /// into BSH).
  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Ex = transform(world, state.x_alpha,
                        g0.focka_no_diag, vtol, true);
    auto Ey = transform(world, state.y_alpha,
                        g0.focka_no_diag, vtol, true);
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
            const ResponseGroundState &g0,
            const State    &state,
            const State    &theta,
            double          omega) {
    // --- X side: BSH at +omega ---------------------------------------
    const double shift_p = detail_tda::bsh_shift(g0.aeps, omega);
    auto rhs_x = madness::copy(world, theta.x_alpha);
    gaxpy(world, 1.0, rhs_x, shift_p, state.x_alpha);
    scale(world, rhs_x, -2.0);
    truncate(world, rhs_x);
    auto bsh_p = detail_tda::make_bsh_operators(world, g0.aeps,
                                                omega, g0.lo);
    auto new_x = apply(world, bsh_p, rhs_x);
    new_x      = g0.Qa(new_x);
    truncate(world, new_x);

    // --- Y side: BSH at -omega ---------------------------------------
    const double shift_m = detail_tda::bsh_shift(g0.aeps, -omega);
    auto rhs_y = madness::copy(world, theta.y_alpha);
    gaxpy(world, 1.0, rhs_y, shift_m, state.y_alpha);
    scale(world, rhs_y, -2.0);
    truncate(world, rhs_y);
    auto bsh_m = detail_tda::make_bsh_operators(world, g0.aeps,
                                                -omega, g0.lo);
    auto new_y = apply(world, bsh_m, rhs_y);
    new_y      = g0.Qa(new_y);
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
                  const ResponseGroundState &g0,
                  const State    &state) {
    // alpha (x+y)·phi
    auto sum_xy_a = madness::copy(world, state.x_alpha);
    gaxpy(world, 1.0, sum_xy_a, 1.0, state.y_alpha);
    auto rho_a = dot(world, g0.amo, sum_xy_a);
    // beta (x+y)·phi
    auto sum_xy_b = madness::copy(world, state.x_beta);
    gaxpy(world, 1.0, sum_xy_b, 1.0, state.y_beta);
    auto rho_b = dot(world, g0.bmo, sum_xy_b);
    auto rho = rho_a + rho_b;
    rho.truncate();
    return rho;
  }

  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &g0,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*g0.coulop, rho1);

    // --- alpha block ---
    auto gamma_ax = mul(world, J_rho, g0.amo, true);
    auto gamma_ay = mul(world, J_rho, g0.amo, true);
    if (g0.c_xc > 0.0) {
      // K[φα, xα](φα)
      auto Kpx_a_phi = detail_tda::apply_exchange(world, g0.amo, state.x_alpha, g0.amo, g0.lo);
      gaxpy(world, 1.0, gamma_ax, -g0.c_xc, Kpx_a_phi);

      // K[yα, φα](φα) — cross-channel for X
      auto Kyp_a_phi = detail_tda::apply_exchange(world, state.y_alpha, g0.amo, g0.amo, g0.lo);
      gaxpy(world, 1.0, gamma_ax, -g0.c_xc, Kyp_a_phi);

      // K[φα, yα](φα)
      auto Kpy_a_phi = detail_tda::apply_exchange(world, g0.amo, state.y_alpha, g0.amo, g0.lo);
      gaxpy(world, 1.0, gamma_ay, -g0.c_xc, Kpy_a_phi);

      // K[xα, φα](φα) — cross-channel for Y
      auto Kxp_a_phi = detail_tda::apply_exchange(world, state.x_alpha, g0.amo, g0.amo, g0.lo);
      gaxpy(world, 1.0, gamma_ay, -g0.c_xc, Kxp_a_phi);
    }
    gamma_ax = g0.Qa(gamma_ax);
    gamma_ay = g0.Qa(gamma_ay);
    truncate(world, gamma_ax);
    truncate(world, gamma_ay);

    // --- beta block ---
    auto gamma_bx = mul(world, J_rho, g0.bmo, true);
    auto gamma_by = mul(world, J_rho, g0.bmo, true);
    if (g0.c_xc > 0.0) {
      auto Kpx_b_phi = detail_tda::apply_exchange(world, g0.bmo, state.x_beta, g0.bmo, g0.lo);
      gaxpy(world, 1.0, gamma_bx, -g0.c_xc, Kpx_b_phi);

      auto Kyp_b_phi = detail_tda::apply_exchange(world, state.y_beta, g0.bmo, g0.bmo, g0.lo);
      gaxpy(world, 1.0, gamma_bx, -g0.c_xc, Kyp_b_phi);

      auto Kpy_b_phi = detail_tda::apply_exchange(world, g0.bmo, state.y_beta, g0.bmo, g0.lo);
      gaxpy(world, 1.0, gamma_by, -g0.c_xc, Kpy_b_phi);

      auto Kxp_b_phi = detail_tda::apply_exchange(world, state.x_beta, g0.bmo, g0.bmo, g0.lo);
      gaxpy(world, 1.0, gamma_by, -g0.c_xc, Kxp_b_phi);
    }
    gamma_bx = g0.Qb(gamma_bx);
    gamma_by = g0.Qb(gamma_by);
    truncate(world, gamma_bx);
    truncate(world, gamma_by);

    return State{std::move(gamma_ax), std::move(gamma_ay),
                 std::move(gamma_bx), std::move(gamma_by)};
  }

  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Vxa = mul_sparse(world, g0.V_local_alpha, state.x_alpha, vtol);
    auto Vya = mul_sparse(world, g0.V_local_alpha, state.y_alpha, vtol);
    auto Vxb = mul_sparse(world, g0.V_local_beta,  state.x_beta,  vtol);
    auto Vyb = mul_sparse(world, g0.V_local_beta,  state.y_beta,  vtol);
    if (g0.c_xc > 0.0) {
      auto k0_ax = detail_tda::apply_exchange(world, g0.amo, g0.amo, state.x_alpha, g0.lo);
      auto k0_ay = detail_tda::apply_exchange(world, g0.amo, g0.amo, state.y_alpha, g0.lo);
      gaxpy(world, 1.0, Vxa, -g0.c_xc, k0_ax);
      gaxpy(world, 1.0, Vya, -g0.c_xc, k0_ay);

      auto k0_bx = detail_tda::apply_exchange(world, g0.bmo, g0.bmo, state.x_beta, g0.lo);
      auto k0_by = detail_tda::apply_exchange(world, g0.bmo, g0.bmo, state.y_beta, g0.lo);
      gaxpy(world, 1.0, Vxb, -g0.c_xc, k0_bx);
      gaxpy(world, 1.0, Vyb, -g0.c_xc, k0_by);
    }
    return State{std::move(Vxa), std::move(Vya),
                 std::move(Vxb), std::move(Vyb)};
  }

  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Exa = transform(world, state.x_alpha,
                         g0.focka_no_diag, vtol, true);
    auto Eya = transform(world, state.y_alpha,
                         g0.focka_no_diag, vtol, true);
    auto Exb = transform(world, state.x_beta,
                         g0.fockb_no_diag,  vtol, true);
    auto Eyb = transform(world, state.y_beta,
                         g0.fockb_no_diag,  vtol, true);
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
            const ResponseGroundState &g0,
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
                                                g0.lo);
      auto out = apply(world, bsh, rhs);
      out = Q(out);
      truncate(world, out);
      return out;
    };
    auto new_xa = apply_side(g0.aeps, g0.Qa,
                             theta.x_alpha, state.x_alpha,  omega);
    auto new_ya = apply_side(g0.aeps, g0.Qa,
                             theta.y_alpha, state.y_alpha, -omega);
    auto new_xb = apply_side(g0.beps,  g0.Qb,
                             theta.x_beta,  state.x_beta,   omega);
    auto new_yb = apply_side(g0.beps,  g0.Qb,
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

// Interface enforcement — currently FD-only (FDSolver<Full, *>). When
// ESSolver<Full, *> lands, promote these to MV3_ASSERT_ES_KERNEL after
// adding compute_T0x / compute_E0x_full / compute_lambda to the Full
// kernels.
#include "kernel_interface.hpp"
MV3_ASSERT_FD_KERNEL(::molresponse_v3::Kernels<::molresponse_v3::Full,
                                                ::molresponse_v3::ClosedShell>);
MV3_ASSERT_FD_KERNEL(::molresponse_v3::Kernels<::molresponse_v3::Full,
                                                ::molresponse_v3::OpenShell>);

#endif // MOLRESPONSE_V3_KERNELS_FULL_HPP
