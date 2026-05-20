#ifndef MOLRESPONSE_V3_KERNELS_TDA_HPP
#define MOLRESPONSE_V3_KERNELS_TDA_HPP

// =========================================================================
// Kernels<TDA, ClosedShell>  and  Kernels<TDA, OpenShell>
//
// Per-state building blocks for Tamm-Dancoff response. Naming follows
// thesis / molresponse convention:
//
//   Lambda = T·x + V0·x - E0_full·x + gamma     (has kinetic, used to
//                                                 build subspace A)
//   Theta  = V0·x - E0_nodi·x + gamma            (no kinetic, used as
//                                                 BSH input; in FD an
//                                                 extra V_perturbation
//                                                 source term is added)
//
// where E0_full = fock_full and E0_nodi = fock_no_diag. In a canonical
// orbital basis fock_full is diagonal (= diag(eps)) and fock_no_diag is
// zero. The choice is dictated by what BSH absorbs: BSH(omega) inverts
// (T - eps - omega), so theta carries the OFF-diagonal Fock only; the
// subspace matrix A_ij = <x_i|Lambda x_j> uses the full Fock so the
// eigenvalues of (A, S) are the omegas.
//
// Per-state operations (V0x, T0x, E0x_full, E0x_nodi, gamma) are
// individually accessible so:
//   - ES iteration:    build V0x/T0x/E0x_full/gamma, assemble Lambda,
//                      diagonalize, rotate {V0x, E0x_nodi, gamma} by U,
//                      assemble Theta from rotated pieces, BSH apply.
//   - FD iteration:    build V0x/E0x_nodi/gamma, assemble Theta + V_p,
//                      BSH apply.
// Both share the same V0x / E0x_nodi / gamma kernel bodies.
// =========================================================================

#include "../solvers/response_state.hpp"
#include "tags.hpp"

#include <madness/chem/SCFOperators.h>
#include <madness/chem/projector.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/tensor/tensor.h>

namespace molresponse_v3 {

using poperatorT = std::shared_ptr<madness::real_convolution_3d>;

/// Read-only problem definition. Closed-shell runs leave *_beta empty.
struct ResponseGroundState {
  std::vector<madness::real_function_3d> phi0_alpha;
  madness::Tensor<double>                eps_alpha;
  std::vector<madness::real_function_3d> phi0_beta;
  madness::Tensor<double>                eps_beta;

  madness::real_function_3d V_local_alpha;
  madness::real_function_3d V_local_beta;

  madness::Tensor<double>   fock_full_alpha;       // for Lambda
  madness::Tensor<double>   fock_full_beta;
  madness::Tensor<double>   fock_no_diag_alpha;    // for Theta
  madness::Tensor<double>   fock_no_diag_beta;

  poperatorT                              coulop;
  madness::QProjector<double, 3>          Q_alpha;
  madness::QProjector<double, 3>          Q_beta;
  double                                  c_xc = 1.0;
  double                                  lo   = 1.0e-10;

  int                                     n_roots = 0;
};

namespace detail_tda {

inline std::vector<madness::real_function_3d>
apply_kinetic(madness::World &world,
              const std::vector<madness::real_function_3d> &v) {
  if (v.empty()) return {};
  std::vector<madness::real_function_3d> result;
  for (int d = 0; d < 3; ++d) {
    madness::real_derivative_3d D(world, d);
    auto dv = apply(world, D, v);
    auto dv2 = apply(world, D, dv);
    if (result.empty()) result = std::move(dv2);
    else gaxpy(world, 1.0, result, 1.0, dv2);
  }
  scale(world, result, -0.5);
  truncate(world, result);
  return result;
}

inline double bsh_shift(const madness::Tensor<double> &eps, double omega) {
  constexpr double guard = 0.05;
  double lumo_shifted = eps(eps.size() - 1) + omega;
  return (lumo_shifted >= 0.0) ? -guard - lumo_shifted : 0.0;
}

inline std::vector<poperatorT>
make_bsh_operators(madness::World &world,
                   const madness::Tensor<double> &eps, double omega,
                   double lo) {
  const double tol = madness::FunctionDefaults<3>::get_thresh();
  const double shift = bsh_shift(eps, omega);
  std::vector<poperatorT> ops(eps.size());
  for (long p = 0; p < eps.size(); ++p) {
    const double mu = std::sqrt(-2.0 * (eps(p) + omega + shift));
    ops[p] = poperatorT(madness::BSHOperatorPtr3D(world, mu, lo, tol));
  }
  return ops;
}

} // namespace detail_tda

// =========================================================================
// Kernels<TDA, ClosedShell>
// =========================================================================
template <>
struct Kernels<TDA, ClosedShell> {

  using vecfuncT = std::vector<madness::real_function_3d>;
  using State    = ResponseStateX<ClosedShell>;

  // ---- per-state operator applications (the shared building blocks) ----

  /// rho1 = 2 * sum_p phi0_alpha[p] * x_alpha[p].
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &target,
                  const State    &state) {
    auto rho1 = dot(world, target.phi0_alpha, state.x_alpha);
    rho1.scale(2.0);
    rho1.truncate();
    return rho1;
  }

  /// gamma = Q( J[rho1]*phi0_alpha - c_xc * K[phi0_alpha, x_alpha](phi0_alpha) ).
  ///
  /// Returns Storage so the (Type, Shell) uniform shape holds: ES Full
  /// will return gamma with both x_alpha and y_alpha components; OpenShell
  /// will also fill the *_beta components. TDA-ClosedShell carries
  /// x_alpha only.
  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &target,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*target.coulop, rho1);
    auto gamma = mul(world, J_rho, target.phi0_alpha, true);

    if (target.c_xc > 0.0) {
      madness::Exchange<double, 3> Kpx(world, target.lo);
      Kpx.set_bra_and_ket(target.phi0_alpha, state.x_alpha);
      Kpx.set_algorithm(madness::Exchange<double, 3>::
                            ExchangeAlgorithm::multiworld_efficient_row);
      auto Kpx_phi = Kpx(target.phi0_alpha);
      gaxpy(world, 1.0, gamma, -target.c_xc, Kpx_phi);
    }
    gamma = target.Q_alpha(gamma);
    truncate(world, gamma);
    return State{std::move(gamma)};
  }

  /// V0·x = V_local * x - c_xc * K[phi0, phi0](x).
  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &target,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Vx = mul_sparse(world, target.V_local_alpha, state.x_alpha, vtol);
    if (target.c_xc > 0.0) {
      madness::Exchange<double, 3> K0(world, target.lo);
      K0.set_bra_and_ket(target.phi0_alpha, target.phi0_alpha);
      K0.set_algorithm(madness::Exchange<double, 3>::
                           ExchangeAlgorithm::multiworld_efficient_row);
      auto k0x = K0(state.x_alpha);
      gaxpy(world, 1.0, Vx, -target.c_xc, k0x);
    }
    return State{std::move(Vx)};
  }

  /// T·x — kinetic energy applied to the response orbitals.
  static State
  compute_T0x(madness::World &world,
              const ResponseGroundState &/*target*/,
              const State    &state) {
    return State{detail_tda::apply_kinetic(world, state.x_alpha)};
  }

  /// E0·x with the FULL Fock matrix (canonical: diag eps; localized:
  /// diag + off-diag). Used to assemble Lambda for the subspace matrix.
  static State
  compute_E0x_full(madness::World &world,
                   const ResponseGroundState &target,
                   const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    return State{transform(world, state.x_alpha,
                           target.fock_full_alpha, vtol, true)};
  }

  /// E0·x with the OFF-DIAGONAL Fock only. Used to assemble Theta for
  /// the BSH driver, where the diagonal eps is absorbed into the BSH
  /// Green's function. In a canonical basis this is the zero vector.
  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &target,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    return State{transform(world, state.x_alpha,
                           target.fock_no_diag_alpha, vtol, true)};
  }

  // ---- assemblies (pure linear combinations of the pieces above) ----

  /// Lambda = T0x + V0x - E0x_full + gamma.   (subspace matrix input)
  static State
  compute_lambda(madness::World &world,
                 const State    &T0x,
                 const State    &V0x,
                 const State    &E0x_full,
                 const State    &gamma) {
    auto L = madness::copy(world, T0x.x_alpha);
    gaxpy(world, 1.0, L,  1.0, V0x.x_alpha);
    gaxpy(world, 1.0, L, -1.0, E0x_full.x_alpha);
    gaxpy(world, 1.0, L,  1.0, gamma.x_alpha);
    truncate(world, L);
    return State{std::move(L)};
  }

  /// Theta = V0x - E0x + gamma.               (BSH driver input)
  /// FD response adds a perturbation source term: Theta_FD = Theta + V_p.
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

  // ---- BSH and residual ----

  /// new_x = Q( BSH(omega) * (-2 * (theta + shift * x)) ).
  /// The shift matches the one inside make_bsh_operators(omega).
  static State
  bsh_apply(madness::World &world,
            const ResponseGroundState &target,
            const State    &state,
            const State    &theta,
            double          omega) {
    const double shift = detail_tda::bsh_shift(target.eps_alpha, omega);

    auto rhs = madness::copy(world, theta.x_alpha);
    gaxpy(world, 1.0, rhs, shift, state.x_alpha);
    scale(world, rhs, -2.0);
    truncate(world, rhs);

    auto bsh_ops = detail_tda::make_bsh_operators(world, target.eps_alpha,
                                                  omega, target.lo);
    auto new_x   = apply(world, bsh_ops, rhs);
    new_x        = target.Q_alpha(new_x);
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
// Kernels<TDA, OpenShell>
//
// UHF TDA — same shape as Static OpenShell but with ONE exchange term
// per spin block (no Y coupling), plus T0x / E0x_full / Lambda used by
// the ES subspace eigenvalue problem.
//
// Density convention: rho = Σ x_α·φ_α + Σ x_β·φ_β  (NO extra factor —
// TDA drops Y, so we don't get the Y=X doubling that Static carries).
// At the restricted limit (α=β orbitals/X equal), this reduces to
// 2·Σ x·φ = Kernels<TDA, ClosedShell>'s density. ✓
//
// V_local is spin-blind (matches HF; spin-resolved Vxc DFT is a TODO).
// =========================================================================
template <>
struct Kernels<TDA, OpenShell> {

  using vecfuncT = std::vector<madness::real_function_3d>;
  using State    = ResponseStateX<OpenShell>;

  /// rho = Σ_p φ_α_p · x_α_p + Σ_p φ_β_p · x_β_p.
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &target,
                  const State    &state) {
    auto rho_a = dot(world, target.phi0_alpha, state.x_alpha);
    auto rho_b = dot(world, target.phi0_beta,  state.x_beta);
    auto rho1  = rho_a + rho_b;
    rho1.truncate();
    return rho1;
  }

  /// gamma_α = Q_α( J[rho]*φ_α - c_xc K[φ_α, x_α](φ_α) )
  /// gamma_β = Q_β( J[rho]*φ_β - c_xc K[φ_β, x_β](φ_β) )
  /// (Single exchange term per spin — that's what distinguishes TDA
  /// from Static, which has K[φ,x](φ) + K[x,φ](φ) per spin.)
  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &target,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*target.coulop, rho1);

    // alpha
    auto gamma_a = mul(world, J_rho, target.phi0_alpha, true);
    if (target.c_xc > 0.0) {
      madness::Exchange<double, 3> Kpx_a(world, target.lo);
      Kpx_a.set_bra_and_ket(target.phi0_alpha, state.x_alpha);
      Kpx_a.set_algorithm(madness::Exchange<double, 3>::
                              ExchangeAlgorithm::multiworld_efficient_row);
      auto Kpx_a_phi = Kpx_a(target.phi0_alpha);
      gaxpy(world, 1.0, gamma_a, -target.c_xc, Kpx_a_phi);
    }
    gamma_a = target.Q_alpha(gamma_a);
    truncate(world, gamma_a);

    // beta
    auto gamma_b = mul(world, J_rho, target.phi0_beta, true);
    if (target.c_xc > 0.0) {
      madness::Exchange<double, 3> Kpx_b(world, target.lo);
      Kpx_b.set_bra_and_ket(target.phi0_beta, state.x_beta);
      Kpx_b.set_algorithm(madness::Exchange<double, 3>::
                              ExchangeAlgorithm::multiworld_efficient_row);
      auto Kpx_b_phi = Kpx_b(target.phi0_beta);
      gaxpy(world, 1.0, gamma_b, -target.c_xc, Kpx_b_phi);
    }
    gamma_b = target.Q_beta(gamma_b);
    truncate(world, gamma_b);

    return State{std::move(gamma_a), std::move(gamma_b)};
  }

  /// V0·x = V_local · x − c_xc · K[φ, φ](x), per spin.
  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &target,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;

    auto Vx_a = mul_sparse(world, target.V_local_alpha, state.x_alpha, vtol);
    if (target.c_xc > 0.0) {
      madness::Exchange<double, 3> K0_a(world, target.lo);
      K0_a.set_bra_and_ket(target.phi0_alpha, target.phi0_alpha);
      K0_a.set_algorithm(madness::Exchange<double, 3>::
                             ExchangeAlgorithm::multiworld_efficient_row);
      auto k0_ax = K0_a(state.x_alpha);
      gaxpy(world, 1.0, Vx_a, -target.c_xc, k0_ax);
    }
    auto Vx_b = mul_sparse(world, target.V_local_beta, state.x_beta, vtol);
    if (target.c_xc > 0.0) {
      madness::Exchange<double, 3> K0_b(world, target.lo);
      K0_b.set_bra_and_ket(target.phi0_beta, target.phi0_beta);
      K0_b.set_algorithm(madness::Exchange<double, 3>::
                             ExchangeAlgorithm::multiworld_efficient_row);
      auto k0_bx = K0_b(state.x_beta);
      gaxpy(world, 1.0, Vx_b, -target.c_xc, k0_bx);
    }
    return State{std::move(Vx_a), std::move(Vx_b)};
  }

  /// T·x per spin (kinetic energy on response orbitals — diagonal in spin).
  static State
  compute_T0x(madness::World &world,
              const ResponseGroundState &/*target*/,
              const State    &state) {
    auto Ta = detail_tda::apply_kinetic(world, state.x_alpha);
    auto Tb = detail_tda::apply_kinetic(world, state.x_beta);
    return State{std::move(Ta), std::move(Tb)};
  }

  /// E0·x with the FULL Fock per spin (Lambda assembly: T0 + V0 − E0_full + γ).
  static State
  compute_E0x_full(madness::World &world,
                   const ResponseGroundState &target,
                   const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Ea = transform(world, state.x_alpha,
                        target.fock_full_alpha, vtol, true);
    auto Eb = transform(world, state.x_beta,
                        target.fock_full_beta,  vtol, true);
    return State{std::move(Ea), std::move(Eb)};
  }

  /// E0·x with the OFF-DIAGONAL Fock per spin (Theta assembly).
  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &target,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Ea = transform(world, state.x_alpha,
                        target.fock_no_diag_alpha, vtol, true);
    auto Eb = transform(world, state.x_beta,
                        target.fock_no_diag_beta,  vtol, true);
    return State{std::move(Ea), std::move(Eb)};
  }

  /// Lambda = T0x + V0x - E0x_full + gamma   per spin.
  static State
  compute_lambda(madness::World &world,
                 const State    &T0x,
                 const State    &V0x,
                 const State    &E0x_full,
                 const State    &gamma) {
    auto La = madness::copy(world, T0x.x_alpha);
    gaxpy(world, 1.0, La,  1.0, V0x.x_alpha);
    gaxpy(world, 1.0, La, -1.0, E0x_full.x_alpha);
    gaxpy(world, 1.0, La,  1.0, gamma.x_alpha);
    truncate(world, La);

    auto Lb = madness::copy(world, T0x.x_beta);
    gaxpy(world, 1.0, Lb,  1.0, V0x.x_beta);
    gaxpy(world, 1.0, Lb, -1.0, E0x_full.x_beta);
    gaxpy(world, 1.0, Lb,  1.0, gamma.x_beta);
    truncate(world, Lb);

    return State{std::move(La), std::move(Lb)};
  }

  /// Theta = V0x - E0x + gamma   per spin.
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

  /// Per-spin BSH apply, each side with its own ε and shift.
  static State
  bsh_apply(madness::World &world,
            const ResponseGroundState &target,
            const State    &state,
            const State    &theta,
            double          omega) {
    // alpha
    const double shift_a = detail_tda::bsh_shift(target.eps_alpha, omega);
    auto rhs_a = madness::copy(world, theta.x_alpha);
    gaxpy(world, 1.0, rhs_a, shift_a, state.x_alpha);
    scale(world, rhs_a, -2.0);
    truncate(world, rhs_a);
    auto bsh_a = detail_tda::make_bsh_operators(world, target.eps_alpha,
                                                omega, target.lo);
    auto new_xa = apply(world, bsh_a, rhs_a);
    new_xa      = target.Q_alpha(new_xa);
    truncate(world, new_xa);

    // beta
    const double shift_b = detail_tda::bsh_shift(target.eps_beta, omega);
    auto rhs_b = madness::copy(world, theta.x_beta);
    gaxpy(world, 1.0, rhs_b, shift_b, state.x_beta);
    scale(world, rhs_b, -2.0);
    truncate(world, rhs_b);
    auto bsh_b = detail_tda::make_bsh_operators(world, target.eps_beta,
                                                omega, target.lo);
    auto new_xb = apply(world, bsh_b, rhs_b);
    new_xb      = target.Q_beta(new_xb);
    truncate(world, new_xb);

    return State{std::move(new_xa), std::move(new_xb)};
  }

  /// Combined-norm residual: sqrt(||Δx_α||² + ||Δx_β||²).
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

#endif // MOLRESPONSE_V3_KERNELS_TDA_HPP
