#pragma once

// ResponseKernels.hpp — shared physics building blocks for excited-state and
// frequency-dependent response calculations.
//
// All functions here are free functions that depend only on their explicit
// arguments. No solver state, no restart logic, no metadata, no MPI topology
// decisions.
//
// The two problem types use these kernels differently:
//
//   Frequency-dependent response (single state, fixed omega):
//     potentials = compute_response_potentials(world, state, gs)
//     theta      = -2 * (potentials.v0 - epsilon + potentials.gamma + vp)
//     x_new      = BSH(omega) * theta;  Qhat * x_new
//
//   Excited-state response (M-state bundle, variational omega):
//     for each i: potentials[i] = compute_response_potentials(world, states[i], gs)
//     build_rotation_matrices(states, lambdas, S, A)
//     diag = diagonalize_bundle(world, S, A)       // -> U, omega
//     states   = rotate_bundle(world, states, diag.U)
//     lambdas  = rotate_bundle(world, lambdas, diag.U)
//     theta_i  = -2 * lambdas[i]    // vp = 0 for excited states
//     x_new_i  = BSH(omega_i) * theta_i;  Qhat * x_new_i
//
// NOTE: The numerical code extracted here from ExcitedStateBundleSolver.cpp
// has not yet been fully validated against the legacy molresponse results.
// Treat as a faithful port requiring cross-validation.

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/tensor/tensor.h>
#include <madness/world/world.h>

#include "GroundStateData.hpp"
#include "ResponseState.hpp"   // must precede ResponseVector.hpp to break circular include
#include "ResponseVector.hpp"
#include "ResponseVectorKernels.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

using namespace madness;

// ============================================================================
// Inner products and norms
// ============================================================================

// Plain (Euclidean) inner product over all active channels.
//
//   ⟨Φ|Φ'⟩_plain = ⟨flat|flat'⟩ = Σ_k Σ_i ⟨channel_k_i | channel'_k_i⟩
//
// flat already collects all channels in contiguous order, so this reduces to
// a single vmra inner() call.  Used to build the A matrix in the excited-state
// subspace rotation (see build_rotation_matrices).
// Unrestricted support is not yet implemented.
template <typename ResponseType>
[[nodiscard]] double response_plain_inner(madness::World &world,
                                          const ResponseType &lhs,
                                          const ResponseType &rhs) {
    static_assert(!response_is_unrestricted_v<ResponseType>,
                  "response_plain_inner: unrestricted not yet implemented");
    return inner(world, lhs.flat, rhs.flat).sum();
}

// Symplectic metric inner product for the TDDFT generalized eigenvalue problem.
//
// The metric arises from the Casida equation structure:
//
//   ⟨Φ|Φ'⟩_metric = ⟨x_α|x'_α⟩ − ⟨y_α|y'_α⟩
//
// The MINUS sign on the y block is NOT arbitrary — it reflects the symplectic
// structure of the Casida (TDDFT) response matrix:
//
//   ⎡ A   B  ⎤ ⎡X⎤       ⎡ 1   0 ⎤ ⎡X⎤
//   ⎣ B*  A* ⎦ ⎣Y⎦  = ω  ⎣−1   0 ⎦ ⎣Y⎦
//
// For TDA (StaticRestricted) B = 0 and there is no y block, so the metric
// reduces to a plain overlap: ⟨x|x'⟩.
// Used in build_rotation_matrices to form the overlap matrix S.
// Unrestricted support is not yet implemented.
template <typename ResponseType>
[[nodiscard]] double response_metric_inner(madness::World &world,
                                           const ResponseType &lhs,
                                           const ResponseType &rhs) {
    static_assert(!response_is_unrestricted_v<ResponseType>,
                  "response_metric_inner: unrestricted not yet implemented");
    double metric = inner(world, lhs.x_alpha, rhs.x_alpha).sum();
    if constexpr (response_has_y_channel_v<ResponseType>)
        metric -= inner(world, lhs.y_alpha, rhs.y_alpha).sum();
    return metric;
}

template <typename ResponseType>
[[nodiscard]] double response_metric_norm2(madness::World &world,
                                           const ResponseType &response) {
    return response_metric_inner(world, response, response);
}

// ============================================================================
// Flat-vector utilities
// ============================================================================

// Norm of a flat function vector: sqrt(<flat|flat>).
[[nodiscard]] inline double
state_norm(madness::World &world, const vector_real_function_3d &flat) {
    return std::sqrt(std::max(0.0, inner(world, flat, flat).sum()));
}

// Replace any uninitialised entries in flat with zero functions.
inline void ensure_initialized_flat(madness::World &world,
                                    vector_real_function_3d &flat) {
    if (flat.empty()) return;
    auto zeros = zero_functions_compressed<double, 3>(
        world, static_cast<int>(flat.size()));
    for (size_t i = 0; i < flat.size(); ++i)
        if (!flat[i].is_initialized()) flat[i] = zeros[i];
}

// Scale all active channels by 1/sqrt(metric_norm2).
// The metric is <x|x> - <y|y>, so this normalises to "metric-unit" length.
// Unrestricted not yet implemented.
template <typename ResponseType>
void normalize_response_metric(madness::World &world, ResponseType &state) {
    static_assert(!response_is_unrestricted_v<ResponseType>,
                  "normalize_response_metric: unrestricted not yet implemented");
    const double metric = response_metric_norm2(world, state);
    if (!std::isfinite(metric) || metric <= 1.0e-14) return;
    scale(world, state.flat, 1.0 / std::sqrt(metric), false);
    state.sync();
    world.gop.fence();
}

// Project all active channels onto the complement of the ground-state orbitals.
// Applies gs.Qhat to the flat vector (Qhat already operates correctly on
// concatenated channel vectors).
// Unrestricted not yet implemented.
template <typename ResponseType>
void project_response_channels(madness::World &world, ResponseType &state,
                                const GroundStateData &gs) {
    static_assert(!response_is_unrestricted_v<ResponseType>,
                  "project_response_channels: unrestricted not yet implemented");
    state.flat = gs.Qhat(state.flat);
    state.sync();
}

// Compute the one-particle response density from the first-order density matrix.
//
// From the first-order density matrix:
//   γ^C(r,r') = Σ_i [ x_i^C(r) φ_i†(r')  +  φ_i(r) y_i^C(r') ]
//
// the response density ρ¹(r) is obtained by taking r → r':
//   ρ¹(r) = Σ_i [ x_i^C(r) φ_i*(r)  +  y_i^{C*}(r) φ_i(r) ]
//
// For TDA (StaticRestricted): y ≡ 0, only the first term contributes.
// For full TDDFT (DynamicRestricted): both x and y contribute.
//
// This density is used for convergence tracking (change in ρ¹ per iteration).
// Unrestricted not yet implemented.
template <typename ResponseType>
[[nodiscard]] madness::real_function_3d
compute_response_density(madness::World &world, const ResponseType &state,
                         const GroundStateData &gs) {
    static_assert(!response_is_unrestricted_v<ResponseType>,
                  "compute_response_density: unrestricted not yet implemented");
    const size_t n = std::min(state.x_alpha.size(), gs.orbitals.size());
    madness::real_function_3d rho{madness::real_factory_3d(world)};
    for (size_t i = 0; i < n; ++i) {
        rho += state.x_alpha[i] * gs.orbitals[i];
        if constexpr (response_has_y_channel_v<ResponseType>)
            rho += state.y_alpha[i] * gs.orbitals[i];
    }
    rho.truncate();
    return rho;
}

// ============================================================================
// Kinetic energy applied to a flat response vector: T*flat = -0.5 * nabla^2
//
// Used in iterate_excited to build the full (T + V - E + Gamma)*x lambda that
// gives a positive-definite A matrix for the subspace rotation, matching the
// legacy Lambda_X = T0X + V0X - E0X + gamma formula.
//
// NOTE: The kinetic contribution is only needed for building the rotation
// matrix A.  The BSH theta = (V - E + Gamma)*x should NOT include T (T is
// handled implicitly by the BSH operator itself).
// ============================================================================
[[nodiscard]] inline vector_real_function_3d
apply_kinetic_flat(madness::World &world, const vector_real_function_3d &flat) {
    if (flat.empty()) return {};
    madness::real_derivative_3d Dx(world, 0);
    madness::real_derivative_3d Dy(world, 1);
    madness::real_derivative_3d Dz(world, 2);

    // apply(world, D, f) reconstructs internally; no explicit copy needed
    auto dvx = apply(world, Dx, flat);
    auto dvy = apply(world, Dy, flat);
    auto dvz = apply(world, Dz, flat);
    auto dvx2 = apply(world, Dx, dvx);
    auto dvy2 = apply(world, Dy, dvy);
    auto dvz2 = apply(world, Dz, dvz);

    // T = -0.5 * (d²/dx² + d²/dy² + d²/dz²)
    auto result = dvx2;
    gaxpy(world, 1.0, result, 1.0, dvy2);
    gaxpy(world, 1.0, result, 1.0, dvz2);
    scale(world, result, -0.5);
    truncate(world, result);
    return result;
}

// ============================================================================
// BSH operators for excited-state response
//
// In the BSH (Bound-State Helmholtz) integral-equation form of the response
// equations, each orbital p contributes one Green's function operator per
// flat slot.  The operator parameter mu is derived from:
//
//   x-channel:  k_p^x = sqrt(-2(ε_p + ω))   [positive-frequency]
//   y-channel:  k_p^y = sqrt(-2(ε_p - ω))   [negative-frequency; TDA: not used]
//
// where ε_p is the ground-state orbital energy and ω is the excitation energy.
//
// Numerical stabilization: if ε_p + freq > 0 (unbound state), a shift is
// applied to keep the BSH denominator negative (well-conditioned):
//   shift = -(ε_p + freq + 0.05)
//   mu = sqrt(-2 * max(ε_p + freq + shift, -1e-8))
//
// Returns one BSH operator per flat slot; orbital_shifts_out carries the
// per-slot shifts so the caller can apply them to the source term.
//
// Channel-frequency mapping:
//   StaticRestricted:   [x_α: +omega]                      (N operators)
//   DynamicRestricted:  [x_α: +omega, y_α: -omega]         (2N operators)
// ============================================================================
template <typename ResponseType>
[[nodiscard]] std::vector<poperatorT>
make_excited_bsh_operators(madness::World &world, const ResponseType &state,
                           double omega, const GroundStateData &gs,
                           std::vector<double> &orbital_shifts_out) {
    static_assert(!response_is_unrestricted_v<ResponseType>,
                  "make_excited_bsh_operators: unrestricted not yet implemented");
    constexpr double lo           = 1.0e-8;
    constexpr double shift_factor = 0.05;
    const double thresh = madness::FunctionDefaults<3>::get_thresh();
    const size_t n = gs.orbitals.size();

    std::vector<poperatorT> ops;
    orbital_shifts_out.clear();

    auto append_block = [&](double freq, bool apply_shift) {
        for (size_t p = 0; p < n; ++p) {
            const double e = gs.getEnergies()(long(p));
            double shift = 0.0;
            if (apply_shift && (e + freq) > 0.0)
                shift = -(e + freq + shift_factor);
            const double denom      = e + freq + shift;
            const double stabilized = std::min(denom, -1.0e-8);
            const double mu         = std::sqrt(std::max(1.0e-16, -2.0 * stabilized));
            ops.emplace_back(poperatorT(BSHOperatorPtr3D(world, mu, lo, thresh)));
            orbital_shifts_out.push_back(shift);
        }
    };

    if constexpr (!response_has_y_channel_v<ResponseType>) {
        // x-only channel (StaticRestricted, TDARestricted): one block at +omega.
        append_block(omega, true);
    } else {
        // x and y channels (DynamicRestricted): x block at +omega, y block at -omega.
        append_block( omega, true);
        append_block(-omega, false);
    }

    // Guard against size mismatch (e.g. partially-initialised state)
    if (ops.size() != state.flat.size()) {
        orbital_shifts_out.resize(state.flat.size(), 0.0);
        if (ops.size() > state.flat.size()) ops.resize(state.flat.size());
    }
    return ops;
}

// ============================================================================
// Epsilon: Hamiltonian_no_diag applied to response channels
//
// Returns a flat vector in the same channel order as ResponseType::flat:
//   StaticRestricted:    [x_alpha...]
//   DynamicRestricted:   [x_alpha..., y_alpha...]
//   StaticUnrestricted:  [x_alpha..., x_beta...]
//   DynamicUnrestricted: [x_alpha..., y_alpha..., x_beta..., y_beta...]
//
// Used inside compute_response_potentials to form lambda = v0 - epsilon + gamma.
// ============================================================================
// Apply Hamiltonian_no_diag to each spin-orbital channel separately.
// H_no_diag is Norb x Norb, so it must be applied to x and y independently —
// never to the concatenated flat vector.
// Unrestricted (alpha+beta) Hamiltonian coupling is not yet implemented.
template <typename ResponseType>
[[nodiscard]] vector_real_function_3d
apply_hamiltonian_no_diag(madness::World &world,
                          const ResponseType &response,
                          const GroundStateData &gs) {
    if constexpr (response_is_unrestricted_v<ResponseType>) {
        MADNESS_EXCEPTION("apply_hamiltonian_no_diag: unrestricted not yet implemented", 1);
    } else if constexpr (!response_has_y_channel_v<ResponseType>) {
        // x-only channel (StaticRestricted, TDARestricted, any future 1-channel restricted type).
        return transform(world, response.x_alpha, gs.Hamiltonian_no_diag, true);
    } else {
        // x and y channels (DynamicRestricted): transform each with the same Norb x Norb H.
        auto eps   = transform(world, response.x_alpha, gs.Hamiltonian_no_diag, true);
        auto eps_y = transform(world, response.y_alpha, gs.Hamiltonian_no_diag, true);
        eps.insert(eps.end(), eps_y.begin(), eps_y.end());
        return eps;
    }
}

// ============================================================================
// Response potentials
//
// Encapsulates the right-hand side terms of the BSH integral-equation form
// of the coupled response equations.  For orbital p:
//
//   RHS_p = V^0 x_p  −  Σ_{i≠p} ε_{ip} x_i  +  g_p'[γ^C] φ_p  +  V_p^C
//         = v0_p      −  epsilon_p             +  gamma_p         +  vp_p
//
// Returned fields:
//
//   v0     = V^0 x_p = (V_local − c_xc K) x_p
//            Local potential (Coulomb + XC local) minus ground-state exchange.
//            No kinetic energy term — T is handled implicitly by the BSH operator.
//
//   epsilon = Σ_{i≠p} ε_{ip} x_i,   ε_{ip} = ∫ φ_i† F^0 φ_p dr
//            Off-diagonal Fock matrix coupling.  Arises from orbital localization;
//            essential for efficiency in MRA (avoids canonical orbitals).
//            Computed by apply_hamiltonian_no_diag().
//
//   gamma  = g_p'[γ^C] φ_p
//            Response exchange-correlation kernel (Γ_X in the paper).
//            Linear coupling of the response density to the XC functional.
//
//   lambda = v0 − epsilon + gamma
//            Full effective Hamiltonian action on x_p.  Used to build the
//            A matrix in the excited-state subspace rotation:
//              A_{ij} = ½ (⟨Φ_i|λ_j⟩_plain + ⟨Φ_j|λ_i⟩_plain)
//
// BSH source term:
//   Frequency-dependent: theta = −2 * (lambda + vp)   [vp = perturbation]
//   Excited-state:       theta = −2 * lambda_rotated   [vp = 0; vp absorbed into rotation]
//
// T is NOT included here.  The BSH operator Ĝ(k_p) accounts for kinetic energy
// implicitly through its Green's function definition.
// ============================================================================
template <typename ResponseType>
struct ResponseStatePotentials {
    ResponseType lambda;  // v0 - epsilon + gamma (for building S, A)
    ResponseType v0;      // (V_local - c_xc*K)*x, no kinetic term
    ResponseType gamma;   // response exchange-correlation coupling
};

// Compute {lambda, v0, gamma} for a single response state.
// Uses compute_ground_exchange and compute_gamma_response (direct per-vector kernels).
template <typename ResponseType>
[[nodiscard]] ResponseStatePotentials<ResponseType>
compute_response_potentials(madness::World &world,
                            const ResponseType &state,
                            const GroundStateData &gs) {
    const auto &flat = state.flat;
    const size_t n   = flat.size();

    if (n == 0) {
        return {state, state, state};
    }

    auto k0 = compute_ground_exchange(world, state, gs.orbitals);
    const double c_xc  = gs.xcf_.hf_exchange_coefficient();
    auto gx = compute_gamma_response(world, state, gs.orbitals, gs.Qhat, c_xc);
    auto v0_flat       = gs.V_local * flat - c_xc * k0.flat;
    auto epsilon_flat  = apply_hamiltonian_no_diag(world, state, gs);
    auto lambda_flat   = v0_flat - epsilon_flat + gx;

    // Reconstruct typed response from a flat vector.
    // flat -> channels via sync(), which splits by channel size.
    auto from_flat = [&](vector_real_function_3d f) -> ResponseType {
        auto r  = state;
        r.flat  = std::move(f);
        r.sync();
        return r;
    };

    ResponseStatePotentials<ResponseType> result;
    result.v0     = from_flat(std::move(v0_flat));
    result.gamma  = from_flat(std::move(gx));
    result.lambda = from_flat(std::move(lambda_flat));
    return result;
}

// ============================================================================
// Excited-state subspace rotation
// ============================================================================

// Build overlap (S) and energy (A) matrices for the M-state excited-state bundle.
//
// These are the S and A matrices of the generalized eigenvalue problem
// solved in each excited-state iteration:
//
//   A c = ω S c
//
// where c are the subspace mixing coefficients and ω are the excitation energies.
//
//   S_ij = ⟨Φ_i|Φ_j⟩_metric
//        = ⟨x_i|x_j⟩ − ⟨y_i|y_j⟩   (symplectic metric; see response_metric_inner)
//
//   A_ij = ½ (⟨Φ_i|λ_j⟩_plain + ⟨Φ_j|λ_i⟩_plain)   (symmetrised)
//
// where λ_j = v0_j − ε_j + γ_j is the effective Hamiltonian action on state j
// (see compute_response_potentials).
//
// Both S and A are symmetrised before being passed to diagonalize_bundle.
// S and A are inputs to diagonalize_bundle.
template <typename ResponseType>
void build_rotation_matrices(madness::World &world,
                             const std::vector<ResponseType> &states,
                             const std::vector<ResponseType> &lambdas,
                             madness::Tensor<double> &S,
                             madness::Tensor<double> &A) {
    const size_t n = states.size();
    S = madness::Tensor<double>(n, n);
    A = madness::Tensor<double>(n, n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            const double sij = response_metric_inner(world, states[i], states[j]);
            S(i, j) = sij;
            S(j, i) = sij;
            const double aij = response_plain_inner(world, states[i], lambdas[j]);
            const double aji = response_plain_inner(world, states[j], lambdas[i]);
            A(i, j) = aij;
            A(j, i) = aji;
        }
    }
    S = 0.5 * (S + transpose(S));
    A = 0.5 * (A + transpose(A));
}

// SVD conditioning of S when it is near-singular (near-linear dependence
// between states). Floors small singular values to avoid sygvp failure.
struct OverlapConditioningResult {
    madness::Tensor<double> overlap;
    size_t floored_singular_values = 0;
    double floor                   = 0.0;
    double min_singular_value      = 0.0;
    bool applied                   = false;
    bool success                   = false;
};

[[nodiscard]] inline OverlapConditioningResult
condition_overlap_matrix(const madness::Tensor<double> &overlap,
                         double thresh) {
    OverlapConditioningResult result;
    const size_t n = static_cast<size_t>(overlap.dim(0));
    if (n == 0) {
        result.success  = true;
        result.overlap  = overlap;
        return result;
    }

    madness::Tensor<double> L, sigma, R;
    auto S_copy = copy(overlap);
    try {
        svd(S_copy, L, sigma, R);
    } catch (...) {
        return result;
    }

    const double floor = std::max(1.0e-12, 10.0 * std::max(thresh, 1.0e-14));
    result.floor = floor;
    result.min_singular_value = (sigma.dim(0) > 0) ? sigma(long(0)) : 0.0;
    for (int64_t i = 1; i < sigma.dim(0); ++i)
        result.min_singular_value = std::min(result.min_singular_value, sigma(i));

    std::vector<double> clipped(n, floor);
    for (size_t i = 0; i < n && i < static_cast<size_t>(sigma.dim(0)); ++i) {
        const double s = sigma(long(i));
        if (!std::isfinite(s) || s < floor)
            ++result.floored_singular_values;
        else
            clipped[i] = s;
    }

    madness::Tensor<double> conditioned(n, n);
    for (size_t r = 0; r < n; ++r)
        for (size_t c = 0; c < n; ++c) {
            double v = 0.0;
            for (size_t k = 0; k < n; ++k)
                v += L(long(r), long(k)) * clipped[k] * L(long(c), long(k));
            conditioned(long(r), long(c)) = v;
        }

    result.overlap = 0.5 * (conditioned + transpose(conditioned));
    result.applied = result.floored_singular_values > 0;
    result.success = true;
    return result;
}

// Rotate a response bundle by U (state-index rotation).
// Mirrors the legacy response_space transform: for each output state i,
//   rotated[i].flat = sum_j U(j,i) * states[j].flat
// flat encodes all active channels for a state, so this single gaxpy loop
// handles x, y (and eventually alpha/beta) uniformly.
// Unrestricted support is not yet implemented.
template <typename ResponseType>
[[nodiscard]] std::vector<ResponseType>
rotate_bundle(madness::World &world,
              const std::vector<ResponseType> &states,
              const madness::Tensor<double>   &U) {
    static_assert(!response_is_unrestricted_v<ResponseType>,
                  "rotate_bundle: unrestricted not yet implemented");
    const size_t M       = states.size();
    if (M == 0) return states;
    const size_t n_slots = states[0].flat.size();

    std::vector<ResponseType> rotated = states;
    for (size_t i = 0; i < M; ++i) {
        auto result_flat = zero_functions_compressed<double, 3>(world, n_slots);
        for (size_t j = 0; j < M; ++j)
            gaxpy(world, 1.0, result_flat, U(long(j), long(i)), states[j].flat);
        rotated[i].flat = std::move(result_flat);
    }
    for (auto &s : rotated)
        s.sync();
    world.gop.fence();
    return rotated;
}

// Result of diagonalize_bundle.
struct DiagonalizeResult {
    bool                     success    = false;
    madness::Tensor<double>  U;
    std::vector<double>      omega;
    // Mapping from sorted output slot -> original slot before energy ordering.
    // Needed by ExcitedResponse driver to track root permutations.
    std::vector<size_t>      slot_order;
};

// Solve the generalized eigenvalue problem S*U = A*U*diag(omega) for the
// excited-state subspace. S is the metric matrix, A is the energy matrix.
//
// Post-processing applied to the eigenvector matrix U:
//   1. Swap columns to maximise diagonal (keeps roots tracking the same states).
//   2. Fix phases so diagonal elements of U are positive.
//   3. Within near-degenerate clusters, apply SVD to remove arbitrary mixing.
//   4. Sort columns by ascending eigenvalue.
//
// Falls back to SVD conditioning of S and retries if sygvp throws.
// Returns DiagonalizeResult::success = false if both attempts fail.
[[nodiscard]] inline DiagonalizeResult
diagonalize_bundle(madness::World        &world,
                   madness::Tensor<double> S,
                   madness::Tensor<double> A,
                   int                    print_level = 0) {
    DiagonalizeResult result;
    const size_t n = static_cast<size_t>(S.dim(0));
    if (n == 0) return result;

    // Floor near-zero diagonal elements so S is positive definite.
    const double diag_floor =
        std::max(1.0e-12, 10.0 * FunctionDefaults<3>::get_thresh());
    for (size_t i = 0; i < n; ++i)
        if (!std::isfinite(S(i, i)) || std::abs(S(i, i)) < diag_floor)
            S(i, i) = diag_floor;
    S = 0.5 * (S + transpose(S));
    A = 0.5 * (A + transpose(A));
    const auto A_copy = copy(A);

    madness::Tensor<double> eigenvectors(n, n);
    madness::Tensor<double> eigenvalues(n);
    bool conditioned = false;

    try {
        sygvp(world, A, S, 1, eigenvectors, eigenvalues);
    } catch (...) {
        const auto cond =
            condition_overlap_matrix(S, FunctionDefaults<3>::get_thresh());
        if (!cond.success) return result;
        if (print_level > 0 && world.rank() == 0)
            madness::print("DIAGONALIZE_RETRY reason=condition_overlap_svd"
                           " floored=", cond.floored_singular_values);
        try {
            auto S2 = cond.overlap;
            auto A2 = copy(A_copy);
            sygvp(world, A2, S2, 1, eigenvectors, eigenvalues);
            conditioned = cond.applied;
        } catch (...) {
            return result;
        }
    }

    auto U = eigenvectors;
    result.omega.resize(n);
    for (size_t i = 0; i < n; ++i)
        result.omega[i] = eigenvalues(long(i));

    // Step 1: swap columns to keep large coefficients near the diagonal.
    bool switched = true;
    while (switched) {
        switched = false;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                const double sold = U(i, i) * U(i, i) + U(j, j) * U(j, j);
                const double snew = U(i, j) * U(i, j) + U(j, i) * U(j, i);
                if (snew <= sold) continue;
                for (size_t row = 0; row < n; ++row)
                    std::swap(U(long(row), long(i)), U(long(row), long(j)));
                std::swap(result.omega[i], result.omega[j]);
                switched = true;
            }
        }
    }

    // Step 2: fix phases.
    for (size_t i = 0; i < n; ++i)
        if (U(long(i), long(i)) < 0.0)
            for (size_t row = 0; row < n; ++row)
                U(long(row), long(i)) *= -1.0;

    // Step 3: SVD rotation within near-degenerate clusters.
    const double thresh_degen = FunctionDefaults<3>::get_thresh();
    size_t ilo = 0;
    while (ilo + 1 < n) {
        size_t ihi = ilo;
        while ((ihi + 1) < n) {
            const double tol = thresh_degen * 10.0 *
                               std::max(std::abs(result.omega[ilo]), 1.0);
            if (std::abs(result.omega[ilo] - result.omega[ihi + 1]) >= tol)
                break;
            ++ihi;
        }
        const size_t nclus = ihi - ilo + 1;
        if (nclus > 1) {
            madness::Tensor<double> q(nclus, nclus);
            for (size_t r = 0; r < nclus; ++r)
                for (size_t c = 0; c < nclus; ++c)
                    q(long(r), long(c)) = U(long(ilo + r), long(ilo + c));
            madness::Tensor<double> W, sigma, VH;
            try {
                svd(q, W, sigma, VH);
                madness::Tensor<double> q_rot(nclus, nclus);
                for (size_t r = 0; r < nclus; ++r)
                    for (size_t c = 0; c < nclus; ++c) {
                        double v = 0.0;
                        for (size_t k = 0; k < nclus; ++k)
                            v += W(long(r), long(k)) * VH(long(k), long(c));
                        q_rot(long(c), long(r)) = v;  // transpose(W * VH)
                    }
                madness::Tensor<double> block(n, nclus);
                for (size_t row = 0; row < n; ++row)
                    for (size_t c = 0; c < nclus; ++c)
                        block(long(row), long(c)) = U(long(row), long(ilo + c));
                for (size_t row = 0; row < n; ++row)
                    for (size_t c = 0; c < nclus; ++c) {
                        double v = 0.0;
                        for (size_t k = 0; k < nclus; ++k)
                            v += block(long(row), long(k)) *
                                 q_rot(long(k), long(c));
                        U(long(row), long(ilo + c)) = v;
                    }
            } catch (...) {
                if (print_level > 0 && world.rank() == 0)
                    madness::print("DIAGONALIZE_WARN"
                                   " reason=degenerate_cluster_svd_exception"
                                   " ilo=", ilo, " ihi=", ihi);
            }
        }
        ilo = ihi + 1;
    }

    // Step 4: sort columns by ascending eigenvalue.
    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return result.omega[a] < result.omega[b]; });

    madness::Tensor<double> sorted_U(n, n);
    std::vector<double> sorted_omega(n);
    for (size_t col = 0; col < n; ++col) {
        const size_t src     = order[col];
        sorted_omega[col]    = result.omega[src];
        for (size_t row = 0; row < n; ++row)
            sorted_U(long(row), long(col)) = U(long(row), long(src));
        if (sorted_U(long(col), long(col)) < 0.0)
            for (size_t row = 0; row < n; ++row)
                sorted_U(long(row), long(col)) *= -1.0;
    }

    if (conditioned && print_level > 1 && world.rank() == 0)
        madness::print("DIAGONALIZE_RECOVERED reason=condition_overlap_svd n=", n);

    result.U          = std::move(sorted_U);
    result.omega      = std::move(sorted_omega);
    result.slot_order = std::move(order);
    result.success    = true;
    return result;
}
