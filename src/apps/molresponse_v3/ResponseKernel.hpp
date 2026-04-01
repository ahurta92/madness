#ifndef MOLRESPONSE_V3_RESPONSEKERNEL_HPP
#define MOLRESPONSE_V3_RESPONSEKERNEL_HPP

#include "GroundState.hpp"
#include "ResponseFunctions.hpp"

#include <madness/chem/SCFOperators.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

namespace molresponse_v3 {

using namespace madness;

/// Response type determines density definition, exchange terms, BSH shifts.
/// The storage (ResponseState) is the same for all three — only the
/// algorithm differs.
enum class ResponseType { Static, Full, TDA };

// =========================================================================
// Response density
// =========================================================================

/// Compute response density from one spin channel (no spin factor).
///
/// Returns the bare spin-channel density:
///   Static: sum(x_i * phi_i)        [y implied = x]
///   Full:   sum((x_i + y_i) * phi_i) [x,y independent]
///   TDA:    sum(x_i * phi_i)          [y = 0]
///
/// The spin factor (2x for restricted Static/TDA) is applied by
/// compute_response_density, NOT here.
inline real_function_3d compute_spin_density(
    World& world,
    ResponseType type,
    const vector_real_function_3d& x,
    const vector_real_function_3d& y,
    const vector_real_function_3d& phi) {

    auto xphi = mul(world, x, phi, true);

    switch (type) {
    case ResponseType::Static:
        return sum(world, xphi, true);

    case ResponseType::Full: {
        auto yphi = mul(world, y, phi, true);
        return sum(world, xphi, true) + sum(world, yphi, true);
    }

    case ResponseType::TDA:
        return sum(world, xphi, true);
    }
    MADNESS_EXCEPTION("Unknown ResponseType", 0);
}

/// Compute total response density handling both spins.
///
/// Restricted Static: rho = 2 * sum(x_alpha * phi_alpha)  [factor 2 from y=x]
/// Unrestricted Static: rho = sum(x_alpha * phi_alpha) + sum(x_beta * phi_beta)
/// Restricted Full: rho = sum((x+y)_alpha * phi_alpha)    [no extra factor]
/// Unrestricted Full: rho = sum above for alpha + beta
inline real_function_3d compute_response_density(
    World& world,
    ResponseType type,
    const RealResponseState& state,
    const GroundState& gs) {

    auto rho_alpha = compute_spin_density(
        world, type, state.x_alpha, state.y_alpha, gs.orbitals_alpha());

    if (gs.is_spin_restricted()) {
        // For restricted Static/TDA: factor 2 from y=x (B-matrix dropped)
        // For restricted Full: no extra factor (x+y already in spin_density)
        if (type == ResponseType::Static || type == ResponseType::TDA) {
            rho_alpha = 2.0 * rho_alpha;
        }
        return rho_alpha;
    }

    // Unrestricted: sum alpha + beta (no doubling — each spin explicit)
    auto rho_beta = compute_spin_density(
        world, type, state.x_beta, state.y_beta, gs.orbitals_beta());
    return rho_alpha + rho_beta;
}

// =========================================================================
// Response coupling (gamma)
// =========================================================================

/// Compute response coupling potential for one spin channel.
///
/// gamma = Q * [2*J[rho]*phi - c_xc * exchange_terms]
///
/// Exchange terms depend on ResponseType:
///   Static: K[phi,x](phi) + K[x,phi](phi)
///   Full x: K[phi,x](phi) + K[y,phi](phi)
///   TDA:    K[phi,x](phi)  [single term]
///
/// For Full y-channel, call with x and y swapped.
///
/// Spin handling: phi, x, y must all be from the SAME spin channel.
/// Coulomb uses rho_total (sum over both spins). Exchange is same-spin
/// only — this matches SCF::apply_potential(ispin) which constructs
/// Exchange(world, calc, ispin) using only the ispin orbital set.
///
/// ispin parameter (0=alpha, 1=beta) is reserved for future XC kernel
/// (fxc) which needs spin-dependent second derivatives.
inline vector_real_function_3d compute_gamma(
    World& world,
    ResponseType type,
    const vector_real_function_3d& x,
    const vector_real_function_3d& y,
    const vector_real_function_3d& phi,
    const real_function_3d& rho_total,
    const poperatorT& coulop,
    const QProjector<double, 3>& Q,
    double c_xc,
    double lo,
    int ispin = 0) {

    long n = phi.size();

    // Coulomb: J[rho] * phi
    auto J_rho = apply(*coulop, rho_total);
    auto Jphi = mul(world, J_rho, phi, true);
    J_rho.clear();

    // Start with 2 * J * phi
    scale(world, Jphi, 2.0);
    auto gamma = std::move(Jphi);

    // Exchange terms (scaled by c_xc)
    if (c_xc > 0.0) {
        Exchange<double, 3> Kpx(world, lo);
        Kpx.set_bra_and_ket(phi, x);
        Kpx.set_algorithm(Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
        auto Kpx_phi = Kpx(phi);  // K[phi,x](phi)

        switch (type) {
        case ResponseType::Static: {
            // K[phi,x](phi) + K[x,phi](phi)
            Exchange<double, 3> Kxp(world, lo);
            Kxp.set_bra_and_ket(x, phi);
            Kxp.set_algorithm(Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
            auto Kxp_phi = Kxp(phi);  // K[x,phi](phi)
            gaxpy(world, 1.0, gamma, -c_xc, Kpx_phi);
            gaxpy(world, 1.0, gamma, -c_xc, Kxp_phi);
            break;
        }
        case ResponseType::Full: {
            // For x-channel: K[phi,x](phi) + K[y,phi](phi)
            Exchange<double, 3> Kyp(world, lo);
            Kyp.set_bra_and_ket(y, phi);
            Kyp.set_algorithm(Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
            auto Kyp_phi = Kyp(phi);  // K[y,phi](phi)
            gaxpy(world, 1.0, gamma, -c_xc, Kpx_phi);
            gaxpy(world, 1.0, gamma, -c_xc, Kyp_phi);
            break;
        }
        case ResponseType::TDA: {
            // K[phi,x](phi) only
            gaxpy(world, 1.0, gamma, -c_xc, Kpx_phi);
            break;
        }
        }
    }

    // TODO: XC kernel fxc for DFT (hybrid and pure DFT functionals)
    // The remaining correlation beyond c_xc*K is captured by the XC kernel
    // second derivative fxc. For unrestricted, ispin selects spin-dependent
    // fxc (alpha or beta channel). See XCOperator(world, scf, ispin).
    // if (has_fxc) {
    //     XCOperator<double,3> xc_kernel(world, xc_string, spin_restricted,
    //                                     arho, brho, ispin);
    //     auto vxc_response = xc_kernel.apply_kernel(rho_response);
    //     gaxpy(world, 1.0, gamma, 1.0, mul(world, vxc_response, phi));
    // }

    // Q-project
    gamma = Q(gamma);
    truncate(world, gamma);
    return gamma;
}

// =========================================================================
// BSH operators
// =========================================================================

/// Compute BSH level-shift to keep mu^2 positive.
inline double compute_bsh_shift(const Tensor<double>& eps, double omega) {
    constexpr double guard = 0.05;
    double lumo_shifted = eps(eps.size() - 1) + omega;
    return (lumo_shifted >= 0.0) ? -guard - lumo_shifted : 0.0;
}

/// Create BSH operators for one spin channel.
///
/// mu_p = sqrt(-2 * (eps_p + omega + shift))
inline std::vector<poperatorT> make_bsh_operators(
    World& world,
    const Tensor<double>& orbital_energies,
    double omega,
    double lo) {

    double tol = FunctionDefaults<3>::get_thresh();
    double shift = compute_bsh_shift(orbital_energies, omega);
    long n = orbital_energies.size();

    std::vector<poperatorT> ops(n);
    for (long p = 0; p < n; p++) {
        double mu = std::sqrt(-2.0 * (orbital_energies(p) + omega + shift));
        ops[p] = poperatorT(BSHOperatorPtr3D(world, mu, lo, tol));
    }
    return ops;
}

/// Create BSH operators for y-channel (Full TDDFT only).
/// Uses -omega: mu_p = sqrt(-2 * (eps_p - omega))
inline std::vector<poperatorT> make_bsh_operators_y(
    World& world,
    const Tensor<double>& orbital_energies,
    double omega,
    double lo) {

    // y-channel: eps_p - omega < 0 guaranteed for bound states, no shift needed
    return make_bsh_operators(world, orbital_energies, -omega, lo);
}

// =========================================================================
// FD iteration step (one spin channel)
// =========================================================================

/// One BSH update step for frequency-dependent response (one spin channel).
///
/// Computes:
///   k0x = K[phi,phi](x)
///   v0x = V_local * x - c_xc * k0x
///   eps_x = transform(x, fock_no_diag)
///   theta = -2 * (v0x - eps_x + gamma + v_p)
///   x_new = Q * BSH(theta)
///
/// gamma must be precomputed via compute_gamma().
inline vector_real_function_3d fd_iteration_step(
    World& world,
    const vector_real_function_3d& x,
    const vector_real_function_3d& phi,
    const vector_real_function_3d& v_perturbation,
    const real_function_3d& V_local,
    const Tensor<double>& fock_no_diag,
    const vector_real_function_3d& gamma,
    const QProjector<double, 3>& Q,
    const std::vector<poperatorT>& bsh_ops,
    double c_xc,
    double lo) {

    double vtol = FunctionDefaults<3>::get_thresh() * 0.1;

    // Ground-state exchange: K[phi,phi](x)
    vector_real_function_3d k0x;
    if (c_xc > 0.0) {
        Exchange<double, 3> K0(world, lo);
        K0.set_bra_and_ket(phi, phi);
        K0.set_algorithm(Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
        k0x = K0(x);
    } else {
        k0x = zero_functions<double, 3>(world, x.size());
    }

    // V_local * x
    auto Vx = mul_sparse(world, V_local, x, vtol);

    // V0*x = V_local*x - c_xc * K0*x
    gaxpy(world, 1.0, Vx, -c_xc, k0x);
    k0x.clear();

    // Off-diagonal Fock coupling: sum_j F_ij * x_j
    auto eps_x = transform(world, x, fock_no_diag, vtol, true);

    // Theta = -2 * (V0*x - eps*x + gamma + v_p)
    // Assemble: theta = V0x - eps_x + gamma + v_p
    gaxpy(world, 1.0, Vx, -1.0, eps_x);  // Vx = V0x - eps_x
    gaxpy(world, 1.0, Vx, 1.0, gamma);   // Vx += gamma
    gaxpy(world, 1.0, Vx, 1.0, v_perturbation);  // Vx += v_p
    scale(world, Vx, -2.0);              // theta = -2 * (...)
    truncate(world, Vx);

    // Apply BSH
    auto x_new = apply(world, bsh_ops, Vx);

    // Q-project
    x_new = Q(x_new);
    truncate(world, x_new);

    return x_new;
}

/// Full FD iteration step handling both spins and x/y channels.
///
/// For Static/TDA: updates x only.
/// For Full: updates both x and y (with swapped exchange coupling for y).
/// For unrestricted: applies alpha and beta channels separately with
///   shared total response density but separate Fock/exchange/BSH.
inline RealResponseState fd_iteration(
    World& world,
    ResponseType type,
    const RealResponseState& current,
    const RealResponseState& perturbation,
    const GroundState& gs,
    const poperatorT& coulop,
    const std::vector<poperatorT>& bsh_alpha_x,
    const std::vector<poperatorT>& bsh_alpha_y,   // empty for Static/TDA
    const std::vector<poperatorT>& bsh_beta_x,    // empty if restricted
    const std::vector<poperatorT>& bsh_beta_y,    // empty if restricted or Static/TDA
    double omega) {

    // Total response density (shared across spins for Coulomb)
    auto rho = compute_response_density(world, type, current, gs);

    double c_xc = gs.hf_exchange_coefficient();
    double lo = gs.params().lo();
    RealResponseState result;

    // --- Alpha x-channel (ispin=0) ---
    auto gamma_alpha_x = compute_gamma(world, type,
        current.x_alpha, current.y_alpha,
        gs.orbitals_alpha(), rho, coulop, gs.Q_alpha(), c_xc, lo, /*ispin=*/0);

    result.x_alpha = fd_iteration_step(world,
        current.x_alpha, gs.orbitals_alpha(),
        perturbation.x_alpha,
        gs.V_local(), gs.focka_no_diag(),
        gamma_alpha_x, gs.Q_alpha(), bsh_alpha_x, c_xc, lo);

    // --- Alpha y-channel (Full only, ispin=0) ---
    if (type == ResponseType::Full) {
        // For y-channel: exchange terms are K[phi,y]+K[x,phi] (swapped from x)
        auto gamma_alpha_y = compute_gamma(world, type,
            current.y_alpha, current.x_alpha,  // swap x and y
            gs.orbitals_alpha(), rho, coulop, gs.Q_alpha(), c_xc, lo, /*ispin=*/0);

        result.y_alpha = fd_iteration_step(world,
            current.y_alpha, gs.orbitals_alpha(),
            perturbation.y_alpha.empty() ? perturbation.x_alpha : perturbation.y_alpha,
            gs.V_local(), gs.focka_no_diag(),
            gamma_alpha_y, gs.Q_alpha(), bsh_alpha_y, c_xc, lo);
    }

    // --- Beta channels (unrestricted only, ispin=1) ---
    if (!gs.is_spin_restricted()) {
        auto gamma_beta_x = compute_gamma(world, type,
            current.x_beta, current.y_beta,
            gs.orbitals_beta(), rho, coulop, gs.Q_beta(), c_xc, lo, /*ispin=*/1);

        result.x_beta = fd_iteration_step(world,
            current.x_beta, gs.orbitals_beta(),
            perturbation.x_beta,
            gs.V_local(), gs.fockb_no_diag(),
            gamma_beta_x, gs.Q_beta(), bsh_beta_x, c_xc, lo);

        if (type == ResponseType::Full) {
            auto gamma_beta_y = compute_gamma(world, type,
                current.y_beta, current.x_beta,
                gs.orbitals_beta(), rho, coulop, gs.Q_beta(), c_xc, lo, /*ispin=*/1);

            result.y_beta = fd_iteration_step(world,
                current.y_beta, gs.orbitals_beta(),
                perturbation.y_beta.empty() ? perturbation.x_beta : perturbation.y_beta,
                gs.V_local(), gs.fockb_no_diag(),
                gamma_beta_y, gs.Q_beta(), bsh_beta_y, c_xc, lo);
        }
    }

    return result;
}

// =========================================================================
// ES building blocks (excited-state specific)
// =========================================================================

/// Compute Lambda: potential operator matrix on response space.
///
/// Lambda = (V0*x - eps_diag*x + gamma)
///
/// Unlike FD theta which uses V_local (no kinetic), Lambda uses the full
/// Fock matrix for the diagonal term. This is used for subspace rotation.
inline vector_real_function_3d compute_lambda(
    World& world,
    ResponseType type,
    const vector_real_function_3d& x,
    const vector_real_function_3d& y,
    const vector_real_function_3d& phi,
    const real_function_3d& V_local,
    const Tensor<double>& fock_no_diag,
    const real_function_3d& rho_response,
    const poperatorT& coulop,
    const QProjector<double, 3>& Q,
    double c_xc,
    double lo) {

    double vtol = FunctionDefaults<3>::get_thresh() * 0.1;

    // Ground-state potential
    vector_real_function_3d k0x;
    if (c_xc > 0.0) {
        Exchange<double, 3> K0(world, lo);
        K0.set_bra_and_ket(phi, phi);
        K0.set_algorithm(Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
        k0x = K0(x);
    } else {
        k0x = zero_functions<double, 3>(world, x.size());
    }

    auto Vx = mul_sparse(world, V_local, x, vtol);
    gaxpy(world, 1.0, Vx, -c_xc, k0x);  // V0*x

    // Off-diagonal Fock
    auto eps_x = transform(world, x, fock_no_diag, vtol, true);

    // Response coupling
    auto gamma = compute_gamma(world, type, x, y, phi,
                                rho_response, coulop, Q, c_xc, lo);

    // Lambda = V0*x - eps*x + gamma
    gaxpy(world, 1.0, Vx, -1.0, eps_x);
    gaxpy(world, 1.0, Vx, 1.0, gamma);
    truncate(world, Vx);

    return Vx;
}

/// Compute subspace potential matrix: A_ij = <x_i | lambda_j>
inline Tensor<double> compute_subspace_matrix(
    World& world,
    const std::vector<vector_real_function_3d>& states,
    const std::vector<vector_real_function_3d>& lambda) {

    long m = states.size();
    Tensor<double> A(m, m);
    for (long i = 0; i < m; i++) {
        for (long j = 0; j < m; j++) {
            A(i, j) = madness::inner(states[i], lambda[j]);
        }
    }
    return A;
}

/// Compute overlap matrix: S_ij = <x_i | x_j>
/// For Full TDDFT with symplectic metric: S_ij = <x_i|x_j> - <y_i|y_j>
inline Tensor<double> compute_overlap_matrix(
    World& world,
    ResponseType type,
    const std::vector<vector_real_function_3d>& states_x,
    const std::vector<vector_real_function_3d>& states_y) {

    long m = states_x.size();
    Tensor<double> S(m, m);
    for (long i = 0; i < m; i++) {
        for (long j = 0; j < m; j++) {
            S(i, j) = madness::inner(states_x[i], states_x[j]);
            if (type == ResponseType::Full && !states_y.empty()) {
                S(i, j) -= madness::inner(states_y[i], states_y[j]);
            }
        }
    }
    return S;
}

/// Diagonalize generalized eigenvalue problem A*U = S*U*diag(omega).
/// Returns {eigenvalues, transformation_matrix}.
inline std::pair<Tensor<double>, Tensor<double>> diagonalize_subspace(
    const Tensor<double>& A,
    const Tensor<double>& S) {

    Tensor<double> evals;
    Tensor<double> U;
    sygvp(World::get_default(), A, S, 1, U, evals);
    return {evals, U};
}

// =========================================================================
// Property computation helpers
// =========================================================================

/// Alpha factor for polarizability: alpha_AB = factor * <x_A | v_B>
inline double alpha_factor(ResponseType type, bool restricted) {
    if (restricted) {
        return (type == ResponseType::Full) ? -2.0 : -4.0;
    } else {
        return (type == ResponseType::Full) ? -1.0 : -2.0;
    }
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_RESPONSEKERNEL_HPP
