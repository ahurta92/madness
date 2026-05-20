#ifndef MOLRESPONSE_V3_KERNELS_STATIC_HPP
#define MOLRESPONSE_V3_KERNELS_STATIC_HPP

// =========================================================================
// Kernels<Static, ClosedShell>
//
// Frequency-dependent response at omega = 0 (the "static" case). For
// closed shells, Y == X by construction, so Storage = ResponseStateX
// (X only). The per-root building blocks (compute_density, gamma, V0x,
// E0x, theta, bsh_apply, compute_residual_norm) are byte-identical to
// Kernels<TDA, ClosedShell>'s — what's absent here is everything used
// to build the subspace eigenproblem (T0x, E0x_full, compute_lambda).
// The FD solver doesn't diagonalize, so omitting those keeps the
// signature surface honest: a Static-FD step cannot accidentally call
// a routine it has no business calling.
//
// The bsh_apply signature takes omega like TDA's so FDSolver can pass
// channel.omega uniformly (== 0.0 here).
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
struct Kernels<Static, ClosedShell> {

  using vecfuncT = std::vector<madness::real_function_3d>;
  using State    = ResponseStateX<ClosedShell>;

  /// Singlet closed-shell static response density:
  ///   rho1 = 4 * sum_p amo[p] * x_alpha[p]
  ///        = (spin × 2) · (Y=X × 2) · Σ_p φ_p · x_p
  ///
  /// Factor 4 matches the Coulomb coefficient in the singlet (A+B)
  /// matrix element 4(pa|qb)(X+Y) at Y=X, and is the convention used
  /// by legacy v3 ResponseKernel.hpp::compute_response_density for the
  /// Static case. Self-consistent with Full's density formula at Y=X:
  ///   2·Σφ·(x+y) |_{y=x} = 4·Σφ·x.
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &g0,
                  const State    &state) {
    auto rho1 = dot(world, g0.amo, state.x_alpha);
    rho1.scale(4.0);
    rho1.truncate();
    return rho1;
  }

  /// gamma = Q( J[rho1]*phi0 - c_xc * (K[phi0, x](phi0) + K[x, phi0](phi0)) ).
  ///
  /// Static = the Y=X limit of Full, so BOTH exchange terms remain.
  /// (TDA drops the second one by neglecting Y; that's the only
  /// difference between Static and TDA gammas at the kernel level.)
  /// Mirrors molresponse_v3/ResponseKernel.hpp::compute_gamma Static case.
  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &g0,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*g0.coulop, rho1);
    auto gamma = mul(world, J_rho, g0.amo, true);
    if (g0.c_xc > 0.0) {
      // K[phi0, x_alpha](phi0)
      auto Kpx_phi = detail_tda::apply_exchange(world, g0.amo, state.x_alpha, g0.amo, g0.lo);
      gaxpy(world, 1.0, gamma, -g0.c_xc, Kpx_phi);

      // K[x_alpha, phi0](phi0) — the second term that distinguishes
      // Static from TDA.
      auto Kxp_phi = detail_tda::apply_exchange(world, state.x_alpha, g0.amo, g0.amo, g0.lo);
      gaxpy(world, 1.0, gamma, -g0.c_xc, Kxp_phi);
    }
    gamma = g0.Qa(gamma);
    truncate(world, gamma);
    return State{std::move(gamma)};
  }

  /// V0·x = V_local * x - c_xc * K[phi0, phi0](x).
  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Vx = mul_sparse(world, g0.V_local_alpha, state.x_alpha, vtol);
    if (g0.c_xc > 0.0) {
      auto k0x = detail_tda::apply_exchange(world, g0.amo, g0.amo, state.x_alpha, g0.lo);
      gaxpy(world, 1.0, Vx, -g0.c_xc, k0x);
    }
    return State{std::move(Vx)};
  }

  /// E0·x with the OFF-DIAGONAL Fock only (diag eps absorbed into BSH).
  /// In a canonical basis this is zero; kept for localized-orbital paths.
  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    return State{transform(world, state.x_alpha,
                           g0.focka_no_diag, vtol, true)};
  }

  /// Theta = V0x - E0x + gamma.   (BSH driver input; FD adds the
  /// perturbation source on top of this, see FDSolver::step().)
  static State
  compute_theta(madness::World &world,
                const State    &V0x,
                const State    &E0x,
                const State    &gamma) {
    auto T = madness::copy(world, V0x.x_alpha);
    gaxpy(world, 1.0, T, -1.0, E0x.x_alpha);
    gaxpy(world, 1.0, T,  1.0, gamma.x_alpha);
    truncate(world, T);
    return State{std::move(T)};
  }

  /// new_x = Q( BSH(omega) * (-2 * (theta + shift * x)) ).
  /// `omega` is passed for shape uniformity with TDA / Full kernels;
  /// FDSolver<Static,*> always passes 0.0.
  static State
  bsh_apply(madness::World &world,
            const ResponseGroundState &g0,
            const State    &state,
            const State    &theta,
            double          omega) {
    const double shift = detail_tda::bsh_shift(g0.aeps, omega);
    auto rhs = madness::copy(world, theta.x_alpha);
    gaxpy(world, 1.0, rhs, shift, state.x_alpha);
    scale(world, rhs, -2.0);
    truncate(world, rhs);

    auto bsh_ops = detail_tda::make_bsh_operators(world, g0.aeps,
                                                  omega, g0.lo);
    auto new_x   = apply(world, bsh_ops, rhs);
    new_x        = g0.Qa(new_x);
    truncate(world, new_x);
    return State{std::move(new_x)};
  }

  /// || x_old - x_new ||.
  static double
  compute_residual_norm(madness::World &world,
                        const State    &state_old,
                        const State    &state_new) {
    auto diff = sub(world, state_old.x_alpha, state_new.x_alpha);
    return madness::norm2(world, diff);
  }
};

// =========================================================================
// Kernels<Static, OpenShell>
//
// UHF static response. Storage = ResponseStateX<OpenShell> ({x_alpha,
// x_beta}). The Coulomb response density is spin-summed (single
// real_function_3d carrying total δρ), while exchange and BSH are
// per-spin because exchange is same-spin and eps_α ≠ eps_β.
//
// Density convention (matches legacy v3 ResponseKernel.hpp::
// compute_response_density for unrestricted Static):
//
//     rho = 2 * (Σ x_α · φ0_α + Σ x_β · φ0_β)
//         = (spin_factor=1 for unrestricted)
//           × (y_factor=2 for Static)
//           × Σ_spin Σ_p x_p · φ_p
//
// V_local is shared between spins (it's the spin-blind ground-state
// local potential V_nuc + V_J + V_xc^total). For HF (c_xc=1, no XC
// kernel) this is exact; for hybrid DFT with spin-resolved Vxc the
// per-spin V_local fields would need wiring (open TODO).
// =========================================================================
template <>
struct Kernels<Static, OpenShell> {

  using vecfuncT = std::vector<madness::real_function_3d>;
  using State    = ResponseStateX<OpenShell>;

  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &g0,
                  const State    &state) {
    auto rho_a = dot(world, g0.amo, state.x_alpha);
    auto rho_b = dot(world, g0.bmo,  state.x_beta);
    auto rho   = rho_a + rho_b;
    rho.scale(2.0);
    rho.truncate();
    return rho;
  }

  /// gamma_α = Q_α( J[rho]*φ_α - c_xc(K[φ_α,x_α](φ_α) + K[x_α,φ_α](φ_α)) )
  /// gamma_β = Q_β( J[rho]*φ_β - c_xc(K[φ_β,x_β](φ_β) + K[x_β,φ_β](φ_β)) )
  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &g0,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*g0.coulop, rho1);

    // --- alpha block ---
    auto gamma_a = mul(world, J_rho, g0.amo, true);
    if (g0.c_xc > 0.0) {
      auto Kpx_a_phi = detail_tda::apply_exchange(world, g0.amo, state.x_alpha, g0.amo, g0.lo);
      gaxpy(world, 1.0, gamma_a, -g0.c_xc, Kpx_a_phi);

      auto Kxp_a_phi = detail_tda::apply_exchange(world, state.x_alpha, g0.amo, g0.amo, g0.lo);
      gaxpy(world, 1.0, gamma_a, -g0.c_xc, Kxp_a_phi);
    }
    gamma_a = g0.Qa(gamma_a);
    truncate(world, gamma_a);

    // --- beta block (same shape, beta orbitals) ---
    auto gamma_b = mul(world, J_rho, g0.bmo, true);
    if (g0.c_xc > 0.0) {
      auto Kpx_b_phi = detail_tda::apply_exchange(world, g0.bmo, state.x_beta, g0.bmo, g0.lo);
      gaxpy(world, 1.0, gamma_b, -g0.c_xc, Kpx_b_phi);

      auto Kxp_b_phi = detail_tda::apply_exchange(world, state.x_beta, g0.bmo, g0.bmo, g0.lo);
      gaxpy(world, 1.0, gamma_b, -g0.c_xc, Kxp_b_phi);
    }
    gamma_b = g0.Qb(gamma_b);
    truncate(world, gamma_b);

    return State{std::move(gamma_a), std::move(gamma_b)};
  }

  /// V0·x for each spin. V_local is shared; K0 is same-spin.
  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    // Alpha
    auto Vx_a = mul_sparse(world, g0.V_local_alpha, state.x_alpha, vtol);
    if (g0.c_xc > 0.0) {
      auto k0_ax = detail_tda::apply_exchange(world, g0.amo, g0.amo, state.x_alpha, g0.lo);
      gaxpy(world, 1.0, Vx_a, -g0.c_xc, k0_ax);
    }
    // Beta
    auto Vx_b = mul_sparse(world, g0.V_local_beta, state.x_beta, vtol);
    if (g0.c_xc > 0.0) {
      auto k0_bx = detail_tda::apply_exchange(world, g0.bmo, g0.bmo, state.x_beta, g0.lo);
      gaxpy(world, 1.0, Vx_b, -g0.c_xc, k0_bx);
    }
    return State{std::move(Vx_a), std::move(Vx_b)};
  }

  /// E0·x using the per-spin off-diagonal Fock.
  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Ea = transform(world, state.x_alpha,
                        g0.focka_no_diag, vtol, true);
    auto Eb = transform(world, state.x_beta,
                        g0.fockb_no_diag,  vtol, true);
    return State{std::move(Ea), std::move(Eb)};
  }

  /// Theta = V0x - E0x + gamma  for each spin.
  static State
  compute_theta(madness::World &world,
                const State    &V0x,
                const State    &E0x,
                const State    &gamma) {
    auto Ta = madness::copy(world, V0x.x_alpha);
    gaxpy(world, 1.0, Ta, -1.0, E0x.x_alpha);
    gaxpy(world, 1.0, Ta,  1.0, gamma.x_alpha);
    truncate(world, Ta);

    auto Tb = madness::copy(world, V0x.x_beta);
    gaxpy(world, 1.0, Tb, -1.0, E0x.x_beta);
    gaxpy(world, 1.0, Tb,  1.0, gamma.x_beta);
    truncate(world, Tb);

    return State{std::move(Ta), std::move(Tb)};
  }

  /// Per-spin BSH apply. Each spin gets its own shift, BSH ops, and Q.
  static State
  bsh_apply(madness::World &world,
            const ResponseGroundState &g0,
            const State    &state,
            const State    &theta,
            double          omega) {
    // alpha
    const double shift_a = detail_tda::bsh_shift(g0.aeps, omega);
    auto rhs_a = madness::copy(world, theta.x_alpha);
    gaxpy(world, 1.0, rhs_a, shift_a, state.x_alpha);
    scale(world, rhs_a, -2.0);
    truncate(world, rhs_a);
    auto bsh_a = detail_tda::make_bsh_operators(world, g0.aeps,
                                                omega, g0.lo);
    auto new_xa = apply(world, bsh_a, rhs_a);
    new_xa      = g0.Qa(new_xa);
    truncate(world, new_xa);

    // beta
    const double shift_b = detail_tda::bsh_shift(g0.beps, omega);
    auto rhs_b = madness::copy(world, theta.x_beta);
    gaxpy(world, 1.0, rhs_b, shift_b, state.x_beta);
    scale(world, rhs_b, -2.0);
    truncate(world, rhs_b);
    auto bsh_b = detail_tda::make_bsh_operators(world, g0.beps,
                                                omega, g0.lo);
    auto new_xb = apply(world, bsh_b, rhs_b);
    new_xb      = g0.Qb(new_xb);
    truncate(world, new_xb);

    return State{std::move(new_xa), std::move(new_xb)};
  }

  /// Combined-norm residual: || (x_α, x_β)_old − (x_α, x_β)_new ||.
  static double
  compute_residual_norm(madness::World &world,
                        const State    &state_old,
                        const State    &state_new) {
    auto da = sub(world, state_old.x_alpha, state_new.x_alpha);
    auto db = sub(world, state_old.x_beta,  state_new.x_beta);
    const double na = madness::norm2(world, da);
    const double nb = madness::norm2(world, db);
    return std::sqrt(na * na + nb * nb);
  }
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_KERNELS_STATIC_HPP
