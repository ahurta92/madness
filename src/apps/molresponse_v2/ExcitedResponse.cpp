// ExcitedResponse.cpp — TDDFT/HF excited-state solver implementation.
//
// Two-layer design (mirrors the header):
//
//   Layer 1  iterate_excited<R>  (top of this file)
//     Pure numerics.  Takes (world, states, omega, gs, params), runs the
//     diagonalise-rotate-BSH-KAIN loop, returns ExcitedIterDiagnostics.
//     No files, no restart, no MPI topology decisions.
//
//   Layer 2  ExcitedResponse::Impl
//     Protocol driver: loads ground state, builds a localized-Gaussian
//     fresh guess, calls iterate_excited<R>(), and prints per-iteration
//     diagnostics.  Restart/checkpoint is not yet implemented here.

#include "ExcitedResponse.hpp"

#include <madness/chem/SCF.h>
#include <madness/mra/mra.h>
#include <madness/mra/nonlinsol.h>
#include <madness/mra/vmra.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <random>

using namespace madness;

// ============================================================================
// Anonymous namespace — helpers used only in this translation unit
// ============================================================================
namespace {

// ── LocalizedGaussianGuess ───────────────────────────────────────────────────
//
// Gaussian function centered at `origin` with polynomial prefactor:
//   f(r) = x^p0 * y^p1 * z^p2 * exp(-exponent * |r - origin|^2)

template <std::size_t NDIM>
class LocalizedGaussianGuess
    : public FunctionFunctorInterface<double, NDIM> {
    using coordT = Vector<double, NDIM>;
    coordT           origin_;
    double           exponent_;
    std::vector<int> powers_;
public:
    LocalizedGaussianGuess(const coordT      &origin,
                           double             exponent,
                           std::vector<int>   powers)
        : origin_(origin), exponent_(exponent), powers_(std::move(powers)) {}

    double operator()(const coordT &xyz) const override {
        double radial2  = 0.0;
        double prefactor = 1.0;
        for (std::size_t i = 0; i < NDIM; ++i) {
            const double x = xyz[i] - origin_[i];
            radial2 += x * x;
            if (powers_[i] > 0)
                prefactor *= std::pow(x, powers_[i]);
        }
        return prefactor * std::exp(-exponent_ * radial2);
    }
};

// ── Gram-Schmidt project and orthonormalize ──────────────────────────────────
//
// Projects each state onto the complement of the ground-state occupied space
// (via QProjector), then applies modified Gram-Schmidt orthonormalization.
// States whose norm falls below `tiny_norm` after projection are removed.

void project_and_orthonormalize(
    World                                     &world,
    std::vector<vector_real_function_3d>      &states,
    const QProjector<double, 3>               &Qhat)
{
    if (states.empty()) return;

    // Project.
    for (auto &s : states)
        s = Qhat(s);

    // Remove near-zero states.
    constexpr double tiny_norm = 1.0e-11;
    std::vector<vector_real_function_3d> kept;
    kept.reserve(states.size());
    for (auto &s : states)
        if (state_norm(world, s) > tiny_norm)
            kept.push_back(std::move(s));
    states = std::move(kept);
    if (states.empty()) return;

    // Modified Gram-Schmidt.
    for (size_t i = 0; i < states.size(); ++i) {
        const double ni = state_norm(world, states[i]);
        if (ni <= tiny_norm) continue;
        scale(world, states[i], 1.0 / ni, false);
        for (size_t j = i + 1; j < states.size(); ++j) {
            const double proj = inner(world, states[i], states[j]).sum();
            gaxpy(world, 1.0, states[j], -proj, states[i], false);
        }
    }
    // Final renormalization.
    for (auto &s : states) {
        const double ni = state_norm(world, s);
        if (ni > tiny_norm)
            scale(world, s, 1.0 / ni, false);
        truncate(world, s, FunctionDefaults<3>::get_thresh());
    }
    world.gop.fence();
}

} // namespace

// ============================================================================
// Helpers — exposed via ExcitedResponse.hpp for use in tests and drivers
// ============================================================================

std::vector<double>
estimate_initial_omega(World                                        &world,
                       const std::vector<vector_real_function_3d>   &x_states,
                       const Tensor<double>                         &orbital_energies)
{
    (void)world;  // not needed but kept for API symmetry
    std::vector<double> omega(x_states.size(), 0.0);
    for (size_t i = 0; i < x_states.size(); ++i) {
        double weighted = 0.0, denom = 0.0;
        for (size_t orb = 0; orb < x_states[i].size(); ++orb) {
            const double n2    = std::max(0.0, x_states[i][orb].inner(x_states[i][orb]));
            const double abs_e = std::abs(orbital_energies(long(orb)));
            weighted += abs_e * n2;
            denom    += n2;
        }
        omega[i] = (denom > 1.0e-14)
                       ? weighted / denom + 0.02 * double(i + 1)
                       : 0.10 + 0.05 * double(i + 1);
    }
    return omega;
}

std::vector<vector_real_function_3d>
build_fresh_guess_x_states(World              &world,
                            size_t              n_gen,
                            size_t              n_orbitals,
                            const Molecule     &mol,
                            const QProjector<double, 3> &Qhat,
                            unsigned            seed)
{
    const auto &atoms = mol.get_atoms();
    if (atoms.empty() || n_orbitals == 0) return {};

    std::mt19937                           rng(seed);
    std::uniform_real_distribution<double> exp_dist(0.04, 1.20);
    std::uniform_real_distribution<double> shift_dist(-0.35, 0.35);
    std::uniform_int_distribution<int>     axis_dist(0, 2);
    std::uniform_int_distribution<int>     power_dist(0, 2);

    std::vector<vector_real_function_3d> states;
    states.reserve(n_gen);

    for (size_t i = 0; i < n_gen; ++i) {
        const auto &atom = atoms[i % atoms.size()];
        coord_3d center  = atom.get_coords();
        center[0] += shift_dist(rng);
        center[1] += shift_dist(rng);
        center[2] += shift_dist(rng);

        std::vector<int> powers(3, 0);
        powers[axis_dist(rng)] = power_dist(rng);
        const double exponent = exp_dist(rng);

        auto gauss = real_factory_3d(world).functor(real_functor_3d(
            new LocalizedGaussianGuess<3>(center, exponent, powers)));

        auto state = zero_functions_compressed<double, 3>(
            world, static_cast<int>(n_orbitals));
        state[i % n_orbitals] = gauss;
        states.push_back(std::move(state));
    }

    project_and_orthonormalize(world, states, Qhat);
    return states;
}

// ============================================================================
// iterate_excited<R>(ResponseBundle) — bundle-level KAIN kernel
//
// Algorithm per iteration:
//   1. Compute {lambda, v0, gamma} for every state in the bundle.
//   2. Build overlap (S) and energy (A) matrices; diagonalise A·c = ω·S·c.
//   3. Rotate the bundle and potentials by the eigenvector matrix U.
//   4. Snapshot states after rotation (prev_bundle for KAIN residuals).
//   5. For each state i: compute BSH update → updated_bundle[i].
//   6. Form residual_bundle = updated_bundle − prev_bundle.
//   7. Bundle-level KAIN: kain_bundle = solver.update(prev_bundle, residual_bundle).
//   8. Choose candidate: kain (iter >= 2) or BSH-updated (iter 0–1).
//   9. For each state i: step-restrict candidate, assign, project, normalize.
//  10. Convergence: max_i ||updated[i] − prev[i]|| <= dconv.
// ============================================================================

template <typename R>
ExcitedIterDiagnostics
iterate_excited(World                     &world,
                ResponseBundle<R>         &bundle,
                std::vector<double>       &omega,
                const GroundStateData     &gs,
                const ExcitedSolverParams &params)
{
    static_assert(!response_is_unrestricted_v<R>,
                  "iterate_excited: unrestricted not yet implemented");

    auto &states = bundle.states();
    const size_t M = states.size();
    ExcitedIterDiagnostics result;
    if (M == 0) return result;
    if (omega.size() < M) omega.resize(M, 0.0);

    const size_t N = bundle.num_orbitals();

    // Bundle-level KAIN — one solver for all M states jointly.
    bundle_solver<R> solver(ResponseBundleAllocator<R>{world, M, N}, false);
    solver.set_maxsub(static_cast<int>(params.maxsub));

    // ── Main iteration loop ──────────────────────────────────────────────────
    for (size_t iter = 0; iter < params.maxiter; ++iter) {

        // 1. Potentials for all states in the response space.
        std::vector<R> v0s, gammas, lambdas;
        v0s.reserve(M); gammas.reserve(M); lambdas.reserve(M);
        for (size_t i = 0; i < M; ++i) {
            auto pots = compute_response_potentials(world, states[i], gs);
            lambdas.push_back(std::move(pots.lambda));
            v0s.push_back(std::move(pots.v0));
            gammas.push_back(std::move(pots.gamma));
        }

        // Diagnostic inner products (print_level >= 2).
        if (params.print_level >= 2 && world.rank() == 0) {
            for (size_t i = 0; i < M; ++i) {
                auto xgx  = inner(world, states[i].flat, gammas[i].flat);
                auto xv0x = inner(world, states[i].flat, v0s[i].flat);
                auto xlx  = inner(world, states[i].flat, lambdas[i].flat);
                madness::print("  DIAG state=", i,
                               " <x|gamma>=", xgx.sum(),
                               " <x|v0>=", xv0x.sum(),
                               " <x|lambda>=", xlx.sum());
            }
        }

        // 2. Subspace rotation matrices.
        //    A_ij includes kinetic T so eigenvalues are physical excitation
        //    energies.  BSH theta does NOT include T (implicit in BSH kernel).
        std::vector<R> lambdas_rot;
        lambdas_rot.reserve(M);
        for (size_t i = 0; i < M; ++i) {
            auto T_flat = apply_kinetic_flat(world, states[i].flat);
            R lr = lambdas[i];
            assign_flat_and_sync(lr, add(world, lr.flat, T_flat));
            lambdas_rot.push_back(std::move(lr));
        }

        Tensor<double> S, A;
        build_rotation_matrices(world, states, lambdas_rot, S, A);

        if (params.print_level >= 2 && world.rank() == 0) {
            madness::print("  S matrix:"); madness::print(S);
            madness::print("  A matrix:"); madness::print(A);
        }

        // 3. Diagonalise.
        auto diag = diagonalize_bundle(world, S, A, params.print_level);
        if (!diag.success) {
            if (params.print_level > 0 && world.rank() == 0)
                madness::print("iterate_excited: diagonalisation failed at iter=", iter);
            break;
        }
        omega = diag.omega;

        if (params.print_level >= 2 && world.rank() == 0) {
            madness::print("  omega after diag:");
            for (size_t i = 0; i < M; ++i)
                madness::print("    state=", i, " omega=", omega[i],
                               " eV=", omega[i] * 27.2114);
        }

        // 4. Rotate response space and potentials.
        states = rotate_bundle(world, states, diag.U);
        v0s    = rotate_bundle(world, v0s,    diag.U);
        gammas = rotate_bundle(world, gammas, diag.U);

        // 5a. Snapshot states after rotation (for KAIN residuals and density change).
        //     prev_flats[i] = deep copy of states[i].flat.
        std::vector<vector_real_function_3d> prev_flats;
        prev_flats.reserve(M);
        for (const auto &s : states)
            prev_flats.push_back(copy(world, s.flat, false));
        world.gop.fence();

        // Pre-compute response densities before BSH (for drho tracking).
        std::vector<real_function_3d> rho_prevs;
        rho_prevs.reserve(M);
        for (size_t i = 0; i < M; ++i)
            rho_prevs.push_back(compute_response_density(world, states[i], gs));

        // 5b. BSH loop — compute updated flat for each state.
        std::vector<vector_real_function_3d> updated_flats;
        updated_flats.reserve(M);
        std::vector<double> residuals(M, 0.0);

        for (size_t i = 0; i < M; ++i) {
            auto eps_flat   = apply_hamiltonian_no_diag(world, states[i], gs);
            auto theta_flat = v0s[i].flat - eps_flat + gammas[i].flat;

            std::vector<double> shifts;
            auto bsh_ops = make_excited_bsh_operators(
                world, states[i], omega[i], gs, shifts);

            for (size_t p = 0; p < theta_flat.size() && p < shifts.size(); ++p)
                theta_flat[p] = -2.0 * (theta_flat[p] + shifts[p] * states[i].flat[p]);

            auto updated_flat = apply(world, bsh_ops, theta_flat);
            truncate(world, updated_flat);
            world.gop.fence();

            auto res_flat = sub(world, updated_flat, prev_flats[i]);
            residuals[i]  = state_norm(world, res_flat);

            updated_flats.push_back(std::move(updated_flat));
        }
        world.gop.fence();

        // 6. Build prev_bundle and updated_bundle for KAIN.
        {
            std::vector<R> prev_states, upd_states;
            prev_states.reserve(M);
            upd_states.reserve(M);
            for (size_t i = 0; i < M; ++i) {
                R ps(N);
                ps.flat = std::move(prev_flats[i]);
                ps.sync();
                prev_states.push_back(std::move(ps));

                R us(N);
                us.flat = std::move(updated_flats[i]);
                us.sync();
                upd_states.push_back(std::move(us));
            }
            ResponseBundle<R> prev_bundle(std::move(prev_states));
            ResponseBundle<R> updated_bundle(std::move(upd_states));

            // 7. Bundle KAIN update.
            auto residual_bundle = updated_bundle - prev_bundle;
            auto kain_bundle     = solver.update(prev_bundle, residual_bundle);

            // 8. Choose candidate.
            const ResponseBundle<R> &candidate_bundle =
                (iter >= 2) ? kain_bundle : updated_bundle;

            // 9. Step-restrict, assign, project, normalize.
            std::vector<double> drhos(M, 0.0);
            for (size_t i = 0; i < M; ++i) {
                auto candidate_flat = copy(world, candidate_bundle[i].flat, false);
                world.gop.fence();

                if (params.max_rotation > 0.0) {
                    auto step_flat      = sub(world, candidate_flat, prev_bundle[i].flat);
                    const double step_norm = state_norm(world, step_flat);
                    if (step_norm > params.max_rotation) {
                        const double s = params.max_rotation / step_norm;
                        gaxpy(world, s, candidate_flat, (1.0 - s), prev_bundle[i].flat);
                    }
                }

                assign_flat_and_sync(states[i], std::move(candidate_flat));
                project_response_channels(world, states[i], gs);
                normalize_response_metric(world, states[i]);

                const auto rho_new = compute_response_density(world, states[i], gs);
                const auto drho    = rho_new - rho_prevs[i];
                drhos[i] = std::sqrt(std::max(0.0, drho.inner(drho)));
            }
            world.gop.fence();

            // 10. Diagnostics.
            const double max_res  = *std::max_element(residuals.begin(), residuals.end());
            const double max_drho = *std::max_element(drhos.begin(), drhos.end());
            result.iteration_max_residuals.push_back(max_res);
            result.iteration_max_density_changes.push_back(max_drho);
            result.iterations_used = iter + 1;

            if (params.print_level > 0 && world.rank() == 0) {
                madness::print("iterate_excited"
                               " iter=", iter,
                               " max_res=", max_res,
                               " max_drho=", max_drho);
                for (size_t i = 0; i < M; ++i)
                    madness::print("  state=", i,
                                   " omega=", omega[i],
                                   " eV=",    omega[i] * 27.2114,
                                   " res=",   residuals[i]);
            }

            if (max_res <= params.dconv) {
                result.converged = true;
                break;
            }
        }
    }

    // Propagate last-iteration max to per-state vectors.
    result.residual_norms.assign(M, 0.0);
    result.density_change_norms.assign(M, 0.0);
    if (!result.iteration_max_residuals.empty()) {
        const double last_res  = result.iteration_max_residuals.back();
        const double last_drho = result.iteration_max_density_changes.empty()
                                     ? 0.0
                                     : result.iteration_max_density_changes.back();
        std::fill(result.residual_norms.begin(), result.residual_norms.end(), last_res);
        std::fill(result.density_change_norms.begin(),
                  result.density_change_norms.end(), last_drho);
    }
    return result;
}

// Explicit instantiations (restricted variants only; unrestricted throws via
// static_assert at compile time).
template ExcitedIterDiagnostics
iterate_excited<StaticRestrictedResponse>(
    World &, ResponseBundle<StaticRestrictedResponse> &,
    std::vector<double> &, const GroundStateData &, const ExcitedSolverParams &);

template ExcitedIterDiagnostics
iterate_excited<TDARestrictedResponse>(
    World &, ResponseBundle<TDARestrictedResponse> &,
    std::vector<double> &, const GroundStateData &, const ExcitedSolverParams &);

template ExcitedIterDiagnostics
iterate_excited<DynamicRestrictedResponse>(
    World &, ResponseBundle<DynamicRestrictedResponse> &,
    std::vector<double> &, const GroundStateData &, const ExcitedSolverParams &);

// ============================================================================
// ExcitedResponse::Impl — protocol driver
//
// Loads ground-state data lazily (cached per k/thresh), builds a
// localized-Gaussian fresh guess, dispatches to iterate_excited<R>().
// No restart/checkpoint yet — that is the next migration step.
// ============================================================================

class ExcitedResponse::Impl {
public:
    explicit Impl(ExcitedSolverConfig config) : config_(std::move(config)) {}

    [[nodiscard]] ExcitedProtocolResult
    solve_protocol(World &world, const ExcitedProtocolInput &input) {
        ExcitedProtocolResult result;
        result.attempted = true;

        ensure_ground_state(world);

        const size_t n_orb    = static_cast<size_t>(gs_->getNumOrbitals());
        const size_t n_states = input.num_states;
        const size_t n_gen    = std::max<size_t>(2 * n_states, 4);

        // Build fresh Gaussian guess (x_alpha channel only).
        auto x_states = build_fresh_guess_x_states(
            world, n_gen, n_orb,
            gs_->getMolecule(), gs_->Qhat,
            42u + static_cast<unsigned>(input.protocol_index));

        if (x_states.size() < n_states) {
            if (world.rank() == 0)
                madness::print("ExcitedResponse: not enough linearly independent guess"
                               " states (got", x_states.size(),
                               "need", n_states, ") — returning uncoverged");
            result.stage_status = "insufficient_guess";
            return result;
        }
        // Truncate to exactly num_states.
        if (x_states.size() > n_states)
            x_states.resize(n_states);

        ExcitedSolverParams params;
        params.maxiter      = input.maxiter;
        params.maxsub       = input.maxsub;
        params.dconv        = input.dconv;
        params.print_level  = config_.print_level;
        params.tda          = input.tda;
        params.max_rotation = 0.25;

        std::vector<double> omega =
            estimate_initial_omega(world, x_states, gs_->getEnergies());

        ExcitedIterDiagnostics diag;
        if (input.tda) {
            auto states =
                pack_guess_states<TDARestrictedResponse>(world, x_states, n_orb);
            diag = iterate_excited(world, states, omega, *gs_, params);
        } else {
            auto states =
                pack_guess_states<DynamicRestrictedResponse>(world, x_states, n_orb);
            diag = iterate_excited(world, states, omega, *gs_, params);
        }

        if (world.rank() == 0) {
            const char *variant = input.tda ? "TDA (TDARestricted)"
                                            : "Full (DynamicRestricted)";
            madness::print("\nExcitedResponse::solve_protocol summary"
                           " protocol=", input.protocol_index,
                           " variant=", variant);
            madness::print("  converged=",   diag.converged,
                           " iterations=",   diag.iterations_used,
                           " dconv=",        input.dconv);
            madness::print("  Excitation energies:");
            for (size_t i = 0; i < omega.size(); ++i)
                madness::print("    state", i,
                               "  omega =", omega[i], " au",
                               "  =", omega[i] * 27.2114, " eV");
        }

        result.converged      = diag.converged;
        result.iterations     = diag.iterations_used;
        result.energies       = omega;
        result.residual_norms = diag.residual_norms;
        result.stage_status   = diag.converged ? "converged" : "not_converged";
        result.response_variant = input.tda ? "TDARestricted" : "DynamicRestricted";
        return result;
    }

private:
    ExcitedSolverConfig           config_;
    std::unique_ptr<GroundStateData> gs_;
    int                           gs_k_      = -1;
    double                        gs_thresh_ = std::numeric_limits<double>::infinity();

    void ensure_ground_state(World &world) {
        const int    k      = FunctionDefaults<3>::get_k();
        const double thresh = FunctionDefaults<3>::get_thresh();

        if (gs_ != nullptr && gs_k_ == k &&
            std::abs(gs_thresh_ - thresh) <= 1.0e-14)
            return;

        if (gs_ == nullptr)
            gs_ = std::make_unique<GroundStateData>(
                world, config_.archive_file, Molecule());

        gs_->prepareOrbitals(world, k, thresh);

        const double vtol  = std::max(1.0e-12, 0.1 * thresh);
        const auto   coulop = CoulombOperator(world, 1.0e-8, thresh);
        gs_->computePreliminaries(world, coulop, vtol, "moldft.fock.json");

        gs_k_      = k;
        gs_thresh_ = thresh;
        world.gop.fence();
    }
};

// ============================================================================
// ExcitedResponse public interface
// ============================================================================

ExcitedResponse::ExcitedResponse(ExcitedSolverConfig config)
    : impl_(std::make_unique<Impl>(std::move(config))) {}

ExcitedResponse::~ExcitedResponse() = default;

ExcitedProtocolResult
ExcitedResponse::solve_protocol(World &world, const ExcitedProtocolInput &input) {
    return impl_->solve_protocol(world, input);
}
