#pragma once
// ResponseSolver.hpp — facade for the response solver subsystem.
//
// Shared infrastructure (KAIN allocator + solver typedef) lives here.
// All per-type ops (alpha_factor, compute_density, make_bsh_operators,
// CoupledResponseEquations) are defined in the ops/ sub-headers below and
// are available to any translation unit that includes this file.
//
// To understand the math for a specific response type, open its ops header:
//   ops/StaticRestrictedOps.hpp   — static (ω=0) restricted
//   ops/TDARestrictedOps.hpp      — TDA excited-state restricted
//   ops/DynamicRestrictedOps.hpp  — full TDDFT restricted
//   ops/StaticUnrestrictedOps.hpp — static unrestricted (stubs, phase 5+)
//   ops/DynamicUnrestrictedOps.hpp — dynamic unrestricted (stubs, phase 5+)

#include "GroundStateData.hpp"
#include "ResponseDebugLogger.hpp"
#include "ResponseManager.hpp"
#include "ResponseVector.hpp"
#include <madness/mra/nonlinsol.h>

// ---------------------------------------------------------------------------
// Shared KAIN infrastructure
// ---------------------------------------------------------------------------

struct response_vector_allocator {
    World& world;
    const size_t n_orbitals;
    response_vector_allocator(World& world, size_t n_orbitals)
        : world(world), n_orbitals(n_orbitals) {}
    vector_real_function_3d operator()() {
        return zero_functions<double, 3>(world, static_cast<int>(n_orbitals));
    }
};

using response_solver = XNonlinearSolver<vector_real_function_3d, double,
                                         response_vector_allocator>;

// ---------------------------------------------------------------------------
// Per-type ops modules
// Each header defines: alpha_factor, compute_density, make_bsh_operators,
// and CoupledResponseEquations for its concrete response type.
// ---------------------------------------------------------------------------
#include "ops/StaticRestrictedOps.hpp"
#include "ops/TDARestrictedOps.hpp"
#include "ops/DynamicRestrictedOps.hpp"
#include "ops/StaticUnrestrictedOps.hpp"
#include "ops/DynamicUnrestrictedOps.hpp"
