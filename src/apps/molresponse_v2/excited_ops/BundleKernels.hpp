#pragma once
// excited_ops/BundleKernels.hpp — Generic building blocks for
// excited-state bundle iteration.
//
// All functions here are generic templates dispatched via if constexpr inside
// the kernels in ResponseKernels.hpp.  No per-type specialisations are needed:
//   TDARestrictedResponse   — x-only channel, TDA density (Σ xᵢφᵢ)
//   DynamicRestrictedResponse — x+y channels, full TDDFT density (Σ(xᵢ+yᵢ)φᵢ)
//
// The main iteration loop in ExcitedResponse.cpp (iterate_excited<R>) calls
// these helpers in the order listed below, which maps directly onto the
// algorithm steps visible in the loop:
//
//   1. compute_excited_potentials   — V₀, γ, λ per root
//   2. diagonalize_and_rotate_bundle — S + A → ω, U; rotate states + pots
//   3. snapshot_bundle              — deep copy of bundle for KAIN / residual
//   4. compute_response_densities   — pre-step densities for Δρ metric
//   5. ExcitedStateStep             — θ = V₀−ε·x+γ; x_new = G(μ)·θ  (per root)
//   6. apply_step_project_normalize — step restriction, Q-project, normalize
//
// Bundle inner products (transpose_bundle, transpose_x, transpose_y,
// inner, inner_x, inner_y) and build_rotation_matrices live in namespace
// excited_ops here.  metric_inner is NOT defined here — it is a type-specific
// overload in global namespace (found by ADL at instantiation time):
//   TDARestrictedResponse     → RestrictedTDABundleOps.hpp    S = inner_x(⋯)
//   DynamicRestrictedResponse → RestrictedFullBundleOps.hpp   S = inner_x(⋯) − inner_y(⋯)
//
// Keeping inner/transpose_bundle namespaced avoids conflict with
// madness::inner(World&, vector<Function<T,N>>&, ...) from vmra.h —
// response_all() is a catch-all template so SFINAE cannot distinguish them.
//
// Helper types:
//   BundlePotentials<R>             — {lambda, v0, gamma} per root

#include "../ExcitedResponse.hpp"
#include "../ResponseBundle.hpp"
#include "../ResponseDebugLoggerMacros.hpp"
#include "../ResponseKernels.hpp"
#include "../ResponseSolverConstants.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace excited_ops {

// ============================================================================
// Terminology used in this namespace (short forms):
//
//   Chi       std::vector<R>   — M response vectors; the X_space analogue
//   flat      response_all(s)  — flat function vector for one state s:
//                                 TDA  → [x₀…xₙ]       (N slots)
//                                 Full → [x₀…xₙ, y₀…yₙ] (2N slots)
//   Chi_T     [slot × state]   — Chi transposed: Chi_T[p][i] = flat(Chi[i])[p]
//
// ── transpose_bundle ─────────────────────────────────────────────────────────
// Reorganise Chi from [state × slot] to [slot × state] layout.
// Works for any R: TDA uses N flat slots, Full uses 2N.
// ============================================================================

template <typename R>
[[nodiscard]] std::vector<madness::vector_real_function_3d>
transpose_bundle(const std::vector<R> &Chi)
{
    if (Chi.empty()) return {};
    const size_t M = Chi.size();
    const size_t P = response_all(Chi[0]).size();
    std::vector<madness::vector_real_function_3d> Chi_T(
        P, madness::vector_real_function_3d(M));
    for (size_t i = 0; i < M; ++i) {
        const auto &flat = response_all(Chi[i]);
        for (size_t p = 0; p < P; ++p)
            Chi_T[p][i] = flat[p];
    }
    return Chi_T;
}

// ── inner(World, Chi, Gamma) → Tensor<double> ────────────────────────────────
// Bundle-level plain inner product matrix: result(i,j) = ⟨Chi[i] | Gamma[j]⟩
//
// Mirrors inner(X_space, X_space) in legacy molresponse/x_space.cc:
//   1. Transpose Chi and Gamma to [slot × state] layout.
//   2. Accumulate matrix_inner(Chi_T[p], Gam_T[p]) over all flat slots p.
//
// Lives in namespace excited_ops (not global) to avoid ambiguity with
// madness::inner(World&, vector<Function<T,N>>&, ...) from vmra.h
// (response_all is a catch-all template; SFINAE cannot distinguish the types).
template <typename R>
[[nodiscard]] madness::Tensor<double>
inner(madness::World       &world,
      const std::vector<R> &Chi,
      const std::vector<R> &Gamma)
{
    const long M = static_cast<long>(Chi.size());
    const long N = static_cast<long>(Gamma.size());
    if (M == 0 || N == 0) return madness::Tensor<double>(M, N);
    const auto Chi_T = transpose_bundle(Chi);
    const auto Gam_T = transpose_bundle(Gamma);
    madness::Tensor<double> result(M, N);
    for (size_t p = 0; p < Chi_T.size(); ++p) {
        result += matrix_inner(world, Chi_T[p], Gam_T[p]);
        world.gop.fence();
    }
    return result;
}

// ── transpose_x / transpose_y ────────────────────────────────────────────────
// Channel-specific [state × slot] → [slot × state] transposes.
// transpose_x uses the x-channel; transpose_y uses the y-channel.
// transpose_y requires a y-channel (full TDDFT only); enforced via static_assert.

template <typename R>
[[nodiscard]] std::vector<madness::vector_real_function_3d>
transpose_x(const std::vector<R> &Chi)
{
    if (Chi.empty()) return {};
    const size_t M  = Chi.size();
    const size_t Nx = response_x(Chi[0]).size();
    std::vector<madness::vector_real_function_3d> Chi_T(
        Nx, madness::vector_real_function_3d(M));
    for (size_t i = 0; i < M; ++i) {
        const auto &x = response_x(Chi[i]);
        for (size_t p = 0; p < Nx; ++p)
            Chi_T[p][i] = x[p];
    }
    return Chi_T;
}

template <typename R>
[[nodiscard]] std::vector<madness::vector_real_function_3d>
transpose_y(const std::vector<R> &Chi)
{
    static_assert(response_has_y_channel_v<R>,
                  "transpose_y: response type has no y-channel");
    if (Chi.empty()) return {};
    const size_t M  = Chi.size();
    const size_t Ny = response_y(Chi[0]).size();
    std::vector<madness::vector_real_function_3d> Chi_T(
        Ny, madness::vector_real_function_3d(M));
    for (size_t i = 0; i < M; ++i) {
        const auto &y = response_y(Chi[i]);
        for (size_t p = 0; p < Ny; ++p)
            Chi_T[p][i] = y[p];
    }
    return Chi_T;
}

// ── inner_x / inner_y ────────────────────────────────────────────────────────
// Channel-specific M×N inner product matrices.
//   inner_x: result(i,j) = Σ_p ⟨Chi[i].x[p] | Gamma[j].x[p]⟩
//   inner_y: result(i,j) = Σ_p ⟨Chi[i].y[p] | Gamma[j].y[p]⟩  (y-channel only)
//
// Used by type-specific metric_inner overloads:
//   TDA  →  return inner_x(world, Chi, Gamma);
//   Full →  return inner_x(world, Chi, Gamma) - inner_y(world, Chi, Gamma);

template <typename R>
[[nodiscard]] madness::Tensor<double>
inner_x(madness::World       &world,
        const std::vector<R> &Chi,
        const std::vector<R> &Gamma)
{
    const long M = static_cast<long>(Chi.size());
    const long N = static_cast<long>(Gamma.size());
    if (M == 0 || N == 0) return madness::Tensor<double>(M, N);
    const auto Chi_T = transpose_x(Chi);
    const auto Gam_T = transpose_x(Gamma);
    madness::Tensor<double> result(M, N);
    for (size_t p = 0; p < Chi_T.size(); ++p) {
        result += matrix_inner(world, Chi_T[p], Gam_T[p]);
        world.gop.fence();
    }
    return result;
}

template <typename R>
[[nodiscard]] madness::Tensor<double>
inner_y(madness::World       &world,
        const std::vector<R> &Chi,
        const std::vector<R> &Gamma)
{
    static_assert(response_has_y_channel_v<R>,
                  "inner_y: response type has no y-channel");
    const long M = static_cast<long>(Chi.size());
    const long N = static_cast<long>(Gamma.size());
    if (M == 0 || N == 0) return madness::Tensor<double>(M, N);
    const auto Chi_T = transpose_y(Chi);
    const auto Gam_T = transpose_y(Gamma);
    madness::Tensor<double> result(M, N);
    for (size_t p = 0; p < Chi_T.size(); ++p) {
        result += matrix_inner(world, Chi_T[p], Gam_T[p]);
        world.gop.fence();
    }
    return result;
}

// ── metric_gram_schmidt ───────────────────────────────────────────────────────
// Metric-aware modified Gram-Schmidt orthonormalization of a response bundle.
//
// Mirrors gram_schmidt(Chi.X, Chi.Y) + normalize(Chi) in the reference
// update_x_space_excited, called at the START of every iteration.
//
// Inner product used:
//   TDA  (y ≡ 0): ⟨s_i|s_j⟩  = ⟨x_i|x_j⟩
//   Full (x+y):   ⟨s_i|s_j⟩  = ⟨x_i|x_j⟩ − ⟨y_i|y_j⟩   (symplectic metric)
//
// Steps:
//   1. Q-project each state onto virtual (unoccupied) subspace
//   2. Modified Gram-Schmidt: subtract projections onto earlier states
//   3. Final metric normalisation pass
//
// States whose metric norm falls below tiny_norm are left untouched (near-zero
// states are handled by the subsequent S-matrix conditioning in sygvp).
template <typename R>
void metric_gram_schmidt(madness::World       &world,
                         std::vector<R>       &states,
                         const GroundStateData &gs)
{
    if (states.empty()) return;
    constexpr double tiny_norm = 1.0e-11;

    // 1. Q-project all states
    for (auto &s : states)
        project_response_channels(world, s, gs);
    world.gop.fence();

    // 2. Modified Gram-Schmidt with metric inner product
    for (size_t i = 0; i < states.size(); ++i) {
        const double ni2 = response_metric_norm2(world, states[i]);
        if (ni2 <= tiny_norm * tiny_norm) continue;
        scale(world, response_all(states[i]), 1.0 / std::sqrt(ni2), false);
        states[i].sync();
        for (size_t j = i + 1; j < states.size(); ++j) {
            const double proj = response_metric_inner(world, states[i], states[j]);
            gaxpy(world, 1.0, response_all(states[j]), -proj,
                  response_all(states[i]), false);
            states[j].sync();
        }
    }

    // 3. Final normalization pass
    for (auto &s : states)
        normalize_response_metric(world, s);
    world.gop.fence();
}

// ── build_rotation_matrices ───────────────────────────────────────────────────
// Build the M×M overlap (S) and energy (A) matrices for the generalized
// eigenvalue problem  A c = ω S c  in each excited-state iteration.
//
//   S_ij = ⟨Φ_i|Φ_j⟩_metric   — dispatched to a type-specific metric_inner
//                                  overload (global namespace, type-specific file):
//                                  TDA  → ⟨x|x'⟩        (ops/TDARestrictedOps.hpp)
//                                  Full → ⟨x|x'⟩−⟨y|y'⟩  (RestrictedFullBundleOps.hpp)
//   A_ij = ½(⟨Φ_i|λ_j⟩ + ⟨Φ_j|λ_i⟩)  (symmetrised; λ includes kinetic term)
//
// λ must already include the kinetic term (see diagonalize_and_rotate_bundle).
template <typename R>
void build_rotation_matrices(madness::World          &world,
                             const std::vector<R>    &states,
                             const std::vector<R>    &lambdas,
                             madness::Tensor<double> &S,
                             madness::Tensor<double> &A)
{
    S = metric_inner(world, states, states);
    S = 0.5 * (S + madness::transpose(S));

    A = inner(world, states, lambdas);
    A = 0.5 * (A + madness::transpose(A));
}

// ============================================================================
// BundlePotentials<R> — per-root potential vectors
// ============================================================================

template <typename ResponseType>
struct BundlePotentials {
    std::vector<ResponseType> lambda;  // V₀ − ε·x + γ  (for building S, A)
    std::vector<ResponseType> v0;      // (V_local − c_xc·K₀)·x
    std::vector<ResponseType> gamma;   // Q̂·[response XC coupling]
};

// ============================================================================
// 1. compute_excited_potentials<R>
//    Compute {V₀, γ, λ} for every root in the bundle.
//    Wraps compute_response_potentials<R> from ResponseKernels.hpp.
// ============================================================================

template <typename R>
[[nodiscard]] BundlePotentials<R>
compute_excited_potentials(madness::World                  &world,
                           const std::vector<R>            &states,
                           const GroundStateData           &gs)
{
    BundlePotentials<R> pots;
    pots.lambda.reserve(states.size());
    pots.v0.reserve(states.size());
    pots.gamma.reserve(states.size());
    for (const auto &state : states) {
        auto p = compute_response_potentials(world, state, gs);
        pots.lambda.push_back(std::move(p.lambda));
        pots.v0.push_back(std::move(p.v0));
        pots.gamma.push_back(std::move(p.gamma));
    }
    return pots;
}

// ============================================================================
// 2. diagonalize_and_rotate_bundle<R>
//    Build the M×M overlap (S) and energy (A) matrices, solve A·c = S·c·Ω,
//    then rotate states and potentials by the eigenvector matrix U.
//    Updates omega in place; returns false if diagonalisation fails.
//
//    Matrix construction (delegated to build_rotation_matrices):
//      S_ij = ⟨Φ_i|Φ_j⟩_metric           (symplectic: ⟨x|x'⟩ − ⟨y|y'⟩)
//      A_ij = ½(⟨Φ_i|λ_j⟩ + ⟨Φ_j|λ_i⟩)   (symmetrised; λ includes −T term)
// ============================================================================

template <typename R>
[[nodiscard]] bool
diagonalize_and_rotate_bundle(madness::World        &world,
                              std::vector<R>        &states,
                              BundlePotentials<R>   &potentials,
                              std::vector<double>   &omega,
                              int                    print_level,
                              ResponseDebugLogger   *logger)
{
    // Add kinetic term to lambda before building A:
    //   Λ_i = T·x_i + λ_i   (T acts on the STATE, not on lambda)
    //
    // λ_i = v0_i − ε_i + γ_i  (no kinetic; T is implicit in BSH).
    // For the subspace A matrix we need the full Lambda_X from the reference:
    //   Λ_i = T·x_i + v0_i − ε_i + γ_i = T·x_i + λ_i
    // so kinetic_flat must be applied to states[i], NOT to potentials.lambda[i].
    std::vector<R> lambdas_with_kinetic;
    lambdas_with_kinetic.reserve(potentials.lambda.size());
    for (size_t i = 0; i < potentials.lambda.size(); ++i) {
        const auto &lam = potentials.lambda[i];
        auto kinetic_flat = apply_kinetic_flat(world, response_all(states[i]));  // T * x_i
        auto lam_all = response_all(lam);
        madness::gaxpy(world, 1.0, lam_all, 1.0, kinetic_flat);   // λ_i + T·x_i
        lambdas_with_kinetic.push_back(
            clone_with_all(world, lam, std::move(lam_all)));
    }
    world.gop.fence();

    // Build S and A matrices
    madness::Tensor<double> S, A;
    DEBUG_TIMED_BLOCK(world, logger, "rotate.build_matrices", {
        build_rotation_matrices(world, states, lambdas_with_kinetic, S, A);
    });

    if (print_level >= 2 && world.rank() == 0) {
        madness::print("  S matrix:");  madness::print(S);
        madness::print("  A matrix:");  madness::print(A);
    }

    // Diagonalise A·c = S·c·Ω  →  omega, U
    auto diag = diagonalize_bundle(world, S, A, print_level);
    if (!diag.success) {
        if (print_level > 0 && world.rank() == 0)
            madness::print("iterate_excited: diagonalisation failed");
        return false;
    }

    omega = diag.omega;
    if (print_level >= 2 && world.rank() == 0) {
        madness::print("  omega after diag:");
        for (size_t i = 0; i < omega.size(); ++i)
            madness::print("    root=", i,
                           " omega=", omega[i],
                           " eV=",    omega[i] * 27.2114);
    }
    DEBUG_LOG_VALUE(world, logger, "omega.post_diag", omega);

    // Rotate states and all three potential vectors by U
    DEBUG_TIMED_BLOCK(world, logger, "rotate.states",  { states = rotate_bundle(world, states, diag.U); });
    DEBUG_TIMED_BLOCK(world, logger, "rotate.lambda",  { potentials.lambda = rotate_bundle(world, potentials.lambda, diag.U); });
    DEBUG_TIMED_BLOCK(world, logger, "rotate.v0",      { potentials.v0     = rotate_bundle(world, potentials.v0,     diag.U); });
    DEBUG_TIMED_BLOCK(world, logger, "rotate.gamma",   { potentials.gamma  = rotate_bundle(world, potentials.gamma,  diag.U); });

    // Normalize each state after rotation (mirrors normalize(Chi) inside
    // diagonalizeFockMatrix / diagonalizeFullResponseMatrix in the reference).
    // sygvp eigenvectors are orthonormal in coefficient space; renormalize in
    // function space with the metric norm sqrt(<x|x> - <y|y>).
    for (auto &s : states)
        normalize_response_metric(world, s);
    world.gop.fence();

    return true;
}

// ============================================================================
// 3. snapshot_bundle<R>
//    Deep-copy the current bundle states for use as the KAIN "previous" bundle
//    and for computing residuals (||updated − prev||).
// ============================================================================

template <typename R>
[[nodiscard]] ResponseBundle<R>
snapshot_bundle(madness::World &world, const std::vector<R> &states)
{
    std::vector<R> snapped;
    snapped.reserve(states.size());
    for (const auto &s : states)
        snapped.push_back(clone_with_all(world, s,
                              madness::copy(world, response_all(s), false)));
    world.gop.fence();
    return ResponseBundle<R>(std::move(snapped));
}

// ============================================================================
// 4. compute_response_densities<R>
//    Compute the first-order response density ρ¹ for each root.
//    Delegates to compute_response_density<R> in ResponseKernels.hpp.
// ============================================================================

template <typename R>
[[nodiscard]] std::vector<madness::real_function_3d>
compute_response_densities(madness::World       &world,
                           const std::vector<R> &states,
                           const GroundStateData &gs)
{
    std::vector<madness::real_function_3d> densities;
    densities.reserve(states.size());
    for (const auto &s : states)
        densities.push_back(compute_response_density(world, s, gs));
    return densities;
}

// ============================================================================
// 5. ExcitedStateStep<R>
//    One BSH iteration for a single excited-state root.
//
//    Residual function:
//      θ_p = V₀_p − ε_{ip}·x_i + γ_p
//
//    BSH update (per orbital slot p, absorbing level-shift correction):
//      x_p^new = G(μ_p) · [−2·(θ_p + shift_p · x_p)]
//      μ_p = √(−2·(ε_p + ω)),  shift applied when ε_p + ω ≥ 0
//
//    Channel layout:
//      TDA  (x-only):   N operators at +ω
//      Full (x + y):    N operators at +ω  followed by N at −ω
// ============================================================================

template <typename R>
[[nodiscard]] R
ExcitedStateStep(madness::World    &world,
                 const R           &state,
                 const R           &v0,
                 const R           &gamma,
                 double             omega,
                 const GroundStateData &gs)
{
    // θ_p = V₀_p − ε_{ip}·x_i + γ_p
    auto epsilon_flat = apply_hamiltonian_no_diag(world, state, gs);
    auto theta_all    = response_all(v0) - epsilon_flat + response_all(gamma);

    // Build BSH operators; collect per-slot level-shift corrections
    std::vector<double> shifts;
    auto bsh_ops = make_excited_bsh_operators(world, state, omega, gs, shifts);
    for (size_t p = 0; p < theta_all.size() && p < shifts.size(); ++p)
        theta_all[p] = ResponseSolverConstants::k_bsh_residual_prefactor
                       * (theta_all[p] + shifts[p] * response_all(state)[p]);

    auto updated_all = madness::apply(world, bsh_ops, theta_all);
    madness::truncate(world, updated_all);
    world.gop.fence();
    return clone_with_all(world, state, std::move(updated_all));
}

// ============================================================================
// 6. apply_step_project_normalize<R>
//    Post-KAIN update: per-root step restriction, Q-projection, normalisation.
//    Updates `states` in place from `candidate` (KAIN bundle or raw BSH step).
//    Returns per-root density changes Δρ for diagnostics.
//
//    Step restriction: if ||candidate[i] − prev[i]|| > max_rotation,
//      candidate is linearly mixed back toward prev to cap the step size.
// ============================================================================

template <typename R>
[[nodiscard]] std::vector<double>
apply_step_project_normalize(
        madness::World              &world,
        std::vector<R>              &states,           // updated in place
        const ResponseBundle<R>     &candidate,        // KAIN or raw BSH result
        const ResponseBundle<R>     &prev,             // previous states
        const std::vector<madness::real_function_3d> &densities_before,
        const GroundStateData       &gs,
        double                       max_rotation)
{
    const size_t M = states.size();
    std::vector<double> density_changes(M, 0.0);

    for (size_t i = 0; i < M; ++i) {
        auto candidate_flat =
            madness::copy(world, response_all(candidate[i]), false);
        world.gop.fence();

        // Step restriction: clamp step to max_rotation
        if (max_rotation > 0.0) {
            auto step_flat = sub(world, candidate_flat, response_all(prev[i]));
            const double step_norm = state_norm(world, step_flat);
            if (step_norm > max_rotation) {
                const double s = max_rotation / step_norm;
                madness::gaxpy(world, s, candidate_flat,
                               1.0 - s, response_all(prev[i]));
            }
        }

        // Assign, project onto unoccupied subspace, normalize
        assign_all_and_sync(states[i], std::move(candidate_flat));
        project_response_channels(world, states[i], gs);
        normalize_response_metric(world, states[i]);

        // Δρ for diagnostics
        const auto rho_after  = compute_response_density(world, states[i], gs);
        const auto drho       = rho_after - densities_before[i];
        density_changes[i]    = std::sqrt(std::max(0.0, drho.inner(drho)));
    }
    world.gop.fence();
    return density_changes;
}

// ============================================================================
// Diagnostics helper
// ============================================================================

template <typename R>
inline void
print_potential_diagnostics(madness::World               &world,
                             const std::vector<R>         &states,
                             const BundlePotentials<R>    &potentials,
                             int                           print_level)
{
    if (print_level < 2) return;

    // inner() is a collective: ALL ranks must participate regardless of rank.
    // Only move the rank-0 guard to the print statements below.
    const auto S_gamma  = inner(world, states, potentials.gamma);
    const auto S_v0     = inner(world, states, potentials.v0);
    const auto S_lambda = inner(world, states, potentials.lambda);

    if (world.rank() != 0) return;

    for (size_t i = 0; i < states.size(); ++i)
        madness::print("  POT_DIAG root=", i,
                       " <x|gamma>=",  S_gamma (long(i), long(i)),
                       " <x|v0>=",     S_v0    (long(i), long(i)),
                       " <x|lambda>=", S_lambda(long(i), long(i)));

    // Full off-diagonal coupling at highest verbosity
    if (print_level >= 3) {
        madness::print("  <Chi|Gamma> matrix:");  madness::print(S_gamma);
        madness::print("  <Chi|V0> matrix:");     madness::print(S_v0);
        madness::print("  <Chi|Lambda> matrix:"); madness::print(S_lambda);
    }
}

} // namespace excited_ops
