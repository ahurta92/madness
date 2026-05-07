#ifndef MOLRESPONSE_V3_ESSOLVER_HPP
#define MOLRESPONSE_V3_ESSOLVER_HPP

#include "ESSolverGuess.hpp"
#include "FDSolver.hpp"        // for PrintLevel
#include "GroundState.hpp"
#include "ResponseFunctions.hpp"
#include "ResponseKernel.hpp"

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace molresponse_v3 {

using namespace madness;

/// Result of one excited-state bundle solve.
struct ESSolveResult {
    std::vector<RealResponseState> roots;  // converged response states (orthonormal)
    Tensor<double> omegas;                 // excitation energies (sorted ascending)
    std::vector<double> residuals;         // per-root final BSH residual norm
    int iterations = 0;
    bool converged = false;
};

/// Excited-state response solver — Stage 1: closed-shell TDA only.
///
/// Iteration outline (mirrors molresponse_legacy/iterate_excited.cc with v3
/// kernel helpers):
///   1. orthonormalize bundle (TDA metric is identity → standard GS)
///   2. for each root s: rho_s = density(X_s),  Λ_s = compute_lambda(...)
///   3. A = ⟨X|Λ⟩,  S = ⟨X|X⟩,  (ω, U) = sygvp(A, S)
///   4. rotate X and Λ by U
///   5. per-root: theta = -2(Λ - ω·X);  x_new = Q · BSH(ω) · theta;
///                step-restrict; residual
///   6. converged if max_res < residual_target  AND  max_dω < dconv
///
/// Future stages extend this to Full (y channel + symplectic metric) and
/// to UHF (alpha + beta).  Per-root KAIN is deferred (no acceleration in
/// stage 1 — we want to observe vanilla convergence first).
inline ESSolveResult es_solve(
    World& world,
    ResponseType type,
    long num_roots,
    GroundState& gs,
    long maxiter = 25,
    double dconv = 1e-4,
    double maxrotn = 0.5,
    PrintLevel print_level = PrintLevel::Normal,
    const std::vector<RealResponseState>* initial_guess = nullptr) {

    // -- Stage 1 invariants --
    MADNESS_CHECK(type == ResponseType::TDA);
    MADNESS_CHECK(gs.is_spin_restricted());
    MADNESS_CHECK(num_roots >= 1);

    const double thresh = FunctionDefaults<3>::get_thresh();
    const double lo = gs.params().lo();
    const double c_xc = gs.hf_exchange_coefficient();
    auto coulop = poperatorT(CoulombOperatorPtr(world, lo, 0.001 * thresh));

    // Convergence targets — same shape as fd_solve (FDSolver.hpp:75-77).
    const double density_target = std::max(thresh * 10.0, dconv)
        * std::max(5.0, static_cast<double>(gs.molecule().natom()));
    const double residual_target = density_target * 10.0;
    const double omega_target = std::max(dconv, thresh * 10.0);

    // ---- Initial guess and ω estimates ----
    std::vector<RealResponseState> X;
    if (initial_guess && static_cast<long>(initial_guess->size()) == num_roots) {
        X = *initial_guess;
        // Re-project in case the protocol changed since the guess was made.
        for (auto& s : X) project_response_state(world, s);
        if (print_level >= PrintLevel::Normal && world.rank() == 0)
            print("ES: starting from supplied guess (", X.size(), " roots)");
    } else {
        // Guess-refinement phase (mirrors molresponse_legacy/calc_runner.cc:609-694
        // and TDDFT::iterate_guess at TDDFT.cc:3539).
        //
        // Random Gaussian-envelope noise has huge kinetic energy, so its initial
        // ⟨X|Λ|X⟩ is dominated by ⟨X|T|X⟩ — the eigenvalues come out at tens
        // of hartree, far from the 0.5-hartree physical excitation range.  A
        // few BSH iterations on a *larger* trial bundle (2N), sorted by ω
        // afterwards and truncated to the lowest N, brings the guess into the
        // physically meaningful range before the main iteration starts.
        const long n_trial = 2 * num_roots;
        auto X_trial = make_initial_guess_tda_rhf(world, gs, n_trial);
        if (print_level >= PrintLevel::Normal && world.rank() == 0)
            print("ES: refining ", n_trial, " trial states (5 iterations) "
                  "before main solve");
        auto guess_result = es_solve(
            world, type, n_trial, gs,
            /*maxiter=*/5, dconv, maxrotn,
            (print_level >= PrintLevel::Verbose) ? PrintLevel::Normal
                                                  : PrintLevel::Silent,
            &X_trial);
        // Lowest ω first (sygvp returns ascending eigenvalues).
        X.assign(guess_result.roots.begin(),
                 guess_result.roots.begin() + num_roots);
        if (print_level >= PrintLevel::Normal && world.rank() == 0) {
            print("ES: guess-refinement done. Selected ω:");
            for (long s = 0; s < num_roots; s++)
                print("    seed", s, "  ω=", guess_result.omegas(s));
        }
    }

    Tensor<double> omegas(num_roots);
    for (long s = 0; s < num_roots; s++) omegas(s) = 0.3 + 0.05 * s;

    ESSolveResult result;
    result.iterations = 0;
    result.residuals.assign(num_roots, 0.0);

    for (long iter = 0; iter < maxiter; iter++) {
        result.iterations = iter + 1;

        // 1. Pre-iter Gram-Schmidt (no-op for converged bundle, important
        //    for stages with the symplectic metric — kept here so the
        //    code path stays uniform across types).
        orthonormalize_bundle(world, type, X);

        // 2. Build full Λ_s for the SUBSPACE matrix (needs T·X and full
        //    Fock; this is what gives the right excitation energy as the
        //    sygvp eigenvalue).
        std::vector<vector_real_function_3d> Lambda_sub(num_roots);
        for (long s = 0; s < num_roots; s++) {
            auto rho = compute_response_density(world, type, X[s], gs);
            Lambda_sub[s] = compute_lambda_subspace(
                world, type,
                X[s].x_alpha, X[s].y_alpha, gs.orbitals_alpha(),
                gs.V_local(), gs.focka(), rho, coulop,
                gs.Q_alpha(), c_xc, lo);
        }

        // 3. Build A and S in the current bundle, then diagonalize.
        std::vector<vector_real_function_3d> X_flat(num_roots);
        std::vector<vector_real_function_3d> Y_flat;  // empty for TDA
        for (long s = 0; s < num_roots; s++) X_flat[s] = X[s].x_alpha;
        Tensor<double> A = compute_subspace_matrix(world, X_flat, Lambda_sub);
        Tensor<double> S = compute_overlap_matrix(world, type, X_flat, Y_flat);
        auto [omega_new, U] = diagonalize_subspace(A, S);

        if (world.rank() == 0) {
            print("  A (subspace lambda):");  print(A);
            print("  S (overlap):");          print(S);
            print("  ω(new):");                print(omega_new);
            print("  U:");                     print(U);
        }

        // 4. Rotate bundle by U so each root is now an eigenvector.
        rotate_bundle(world, X, U);

        // 5. Recompute Θ_s on the *rotated* X for the BSH driver. (Θ is
        //    NOT linear in X across the bundle in the same simple way as
        //    Λ once the subspace has been collapsed to eigenvectors —
        //    just recompute it. This mirrors legacy iterate_excited.cc
        //    where Compute_Theta_X is called after compute_new_omegas_transform.)
        std::vector<double> per_root_res(num_roots, 0.0);
        std::vector<RealResponseState> X_next(num_roots);
        for (long s = 0; s < num_roots; s++) {
            auto rho_s = compute_response_density(world, type, X[s], gs);
            auto Theta_s = compute_lambda(
                world, type,
                X[s].x_alpha, X[s].y_alpha, gs.orbitals_alpha(),
                gs.V_local(), gs.focka_no_diag(), rho_s, coulop,
                gs.Q_alpha(), c_xc, lo);

            // Compute the same shift used inside make_bsh_operators so we
            // can apply it to theta as well — without that, the BSH update
            // is wrong whenever ε+ω > 0 (which is most early iterations
            // when ω from a noisy random guess is huge).
            const double shift_s =
                compute_bsh_shift(gs.energies_alpha(), omega_new(s));
            auto bsh_s = make_bsh_operators(
                world, gs.energies_alpha(), omega_new(s), lo);

            auto theta = build_es_theta(world, Theta_s, X[s].x_alpha, shift_s);
            auto x_bsh = apply(world, bsh_s, theta);
            x_bsh = gs.Q_alpha()(x_bsh);
            truncate(world, x_bsh);

            // Per-root residual = ||x_new - X[s].x_alpha||
            auto res_vec = sub(world, x_bsh, X[s].x_alpha);
            double res_norm = norm2(world, res_vec);
            per_root_res[s] = res_norm;

            // Step restriction.
            if (res_norm > maxrotn) {
                double f = maxrotn / res_norm;
                gaxpy(world, f, x_bsh, 1.0 - f, X[s].x_alpha);
            }

            X_next[s].x_alpha = std::move(x_bsh);
        }

        // 6. Compute max residual and ω-change for convergence check.
        double max_res = 0.0, max_domega = 0.0;
        for (long s = 0; s < num_roots; s++) {
            max_res = std::max(max_res, per_root_res[s]);
            max_domega = std::max(
                max_domega, std::abs(omega_new(s) - omegas(s)));
        }

        // Report.
        if (print_level >= PrintLevel::Normal && world.rank() == 0) {
            print("\nES iter ", iter, "  max_res=", max_res,
                  "  max_dω=", max_domega);
            for (long s = 0; s < num_roots; s++) {
                print("    root", s, "  ω=", omega_new(s),
                      "  res=", per_root_res[s]);
            }
        }

        X = std::move(X_next);
        omegas = copy(omega_new);
        result.residuals = per_root_res;
        result.omegas = copy(omegas);

        if (iter > 0 && max_res < residual_target && max_domega < omega_target) {
            result.converged = true;
            if (print_level >= PrintLevel::Normal && world.rank() == 0)
                print("ES CONVERGED at iter ", iter);
            break;
        }
    }

    result.roots = std::move(X);
    return result;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_ESSOLVER_HPP
