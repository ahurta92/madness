#pragma once

// ExcitedResponse.hpp — TDDFT/HF excited-state solver.
//
// Two-layer design:
//
//   Layer 1 — Algorithm kernel (top of ExcitedResponse.cpp):
//     iterate_excited<R>(world, states, omega, gs, params)
//     Pure numerics: potentials → S,A → diagonalise → rotate → BSH → KAIN.
//     No files, no restart, no MPI topology decisions.
//
//   Layer 2 — Protocol driver (ExcitedResponse class):
//     Handles restart/checkpoint, resolves the response variant via a single
//     std::visit, calls iterate_excited<R>, records results via ResponseRecord2.
//
// I/O types (ExcitedProtocolInput, ExcitedProtocolResult, ExcitedRootDescriptor)
// are defined in ResponseRecord.hpp so that ResponseRecord2 can use them
// directly without a circular dependency.

#include "GroundStateData.hpp"
#include "ResponseKernels.hpp"
#include "ResponseRecord.hpp"
#include "ResponseVector.hpp"
// ResponseBundle.hpp must come after ResponseVector.hpp (ResponseState.hpp circular dependency)
#include "ResponseBundle.hpp"

#include <madness/world/world.h>

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

class ResponseDebugLogger;

// ============================================================================
// Parameters for the iterate_excited kernel
// ============================================================================

struct ExcitedSolverParams {
    size_t maxiter      = 20;
    size_t maxsub       = 8;      // KAIN subspace size
    double dconv        = 1.0e-4; // convergence target (residual norm)
    double max_rotation = 0.25;   // step-restriction threshold (0 = disabled)
    int    print_level  = 0;
    bool   tda          = false;
    ResponseDebugLogger *debug_logger = nullptr;
};

// Per-iteration diagnostics returned by iterate_excited.
struct ExcitedIterDiagnostics {
    bool   converged             = false;
    size_t iterations_used       = 0;
    std::vector<double> residual_norms;
    std::vector<double> density_change_norms;
    std::vector<double> iteration_max_residuals;
    std::vector<double> iteration_max_density_changes;
};

// ============================================================================
// Helper functions exposed for testing and direct use
// ============================================================================

/// Build localized Gaussian guess states for the x-channel.
///
/// Generates n_gen trial functions using atoms from mol.  After generation the
/// states are projected onto Q (unoccupied subspace) and orthonormalized.
std::vector<madness::vector_real_function_3d>
build_fresh_guess_x_states(madness::World                 &world,
                            size_t                          n_gen,
                            size_t                          n_orbitals,
                            const madness::Molecule        &mol,
                            const madness::QProjector<double,3> &Qhat,
                            unsigned                        seed = 42u);

/// Estimate initial excitation energies from orbital energies and guess norms.
std::vector<double>
estimate_initial_omega(madness::World                                        &world,
                       const std::vector<madness::vector_real_function_3d>   &x_states,
                       const madness::Tensor<double>                         &orbital_energies);

/// Pack raw x_states into typed response vectors.
///
/// For full response (DynamicRestricted), the y-channel is initialized to zero.
/// States are padded or trimmed to n_orbitals.
template <typename R>
inline std::vector<R>
pack_guess_states(madness::World                                          &world,
                  const std::vector<madness::vector_real_function_3d>    &x_states,
                  size_t                                                   n_orbitals)
{
    std::vector<R> states;
    states.reserve(x_states.size());
    for (const auto &xs : x_states) {
        R rv;
        auto &x = response_x(rv);
        x = copy(world, xs, true);
        if (x.size() < n_orbitals)
            x.resize(n_orbitals,
                zero_functions_compressed<double, 3>(world, 1)[0]);
        x.resize(n_orbitals);

        if constexpr (response_has_y_channel_v<R>)
            response_y(rv) = zero_functions_compressed<double, 3>(
                world, static_cast<int>(n_orbitals));

        rv.flatten();
        states.push_back(std::move(rv));
    }
    return states;
}

// ============================================================================
// Algorithm kernel — pure numerics, no I/O, no restart
//
// Bundle overload (primary): uses a single bundle-level KAIN solver to
// accelerate all M states jointly, capturing inter-state couplings from the
// subspace rotation step.
//
// Vector overload (backward-compat wrapper): wraps the bundle overload.
//
// Both are restricted to closed-shell restricted variants for now.
// ============================================================================

/// Primary overload: operates on ResponseBundle<R> with bundle-level KAIN.
template <typename R>
ExcitedIterDiagnostics
iterate_excited(madness::World            &world,
                ResponseBundle<R>         &bundle,
                std::vector<double>       &omega,
                const GroundStateData     &gs,
                const ExcitedSolverParams &params);

/// Backward-compatible wrapper: delegates to the bundle overload.
template <typename R>
inline ExcitedIterDiagnostics
iterate_excited(madness::World            &world,
                std::vector<R>            &states,
                std::vector<double>       &omega,
                const GroundStateData     &gs,
                const ExcitedSolverParams &params)
{
    ResponseBundle<R> bundle(std::move(states));
    auto diag = iterate_excited(world, bundle, omega, gs, params);
    states = std::move(bundle.states());
    return diag;
}

// ============================================================================
// Config for the ExcitedResponse driver
// ============================================================================

struct ExcitedSolverConfig {
    std::string         archive_file;
    std::string         output_prefix = "excited";
    std::vector<double> protocols;
    int                 print_level   = 0;
};

// ============================================================================
// ExcitedResponse — protocol driver
//
// Owns all restart/checkpoint state. solve_protocol() is called once per
// protocol level by the outer response manager.
// ============================================================================

class ExcitedResponse {
public:
    explicit ExcitedResponse(ExcitedSolverConfig config);
    ~ExcitedResponse();

    [[nodiscard]] ExcitedProtocolResult
    solve_protocol(madness::World              &world,
                   const ExcitedProtocolInput  &input);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};
