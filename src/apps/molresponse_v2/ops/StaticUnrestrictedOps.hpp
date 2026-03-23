#pragma once
// ops/StaticUnrestrictedOps.hpp
//
// Ops module for static (ω=0) unrestricted (open-shell) response.
//
// Status: NOT YET IMPLEMENTED (phase 5+).
//   All functions throw std::runtime_error at runtime.
//   Compile-time dispatch via overload resolution is fully wired up.
//
// When implementing:
//   - alpha/beta spin channels are stored separately
//   - alpha_factor remains −2.0
//   - compute_density: ρ = Σ(xα·φα + xβ·φβ)  (no spin doubling — each channel is explicit)
//   - make_bsh_operators: 2N operators [bsh_α | bsh_β], exponents μ = √(−2·ε_p)
//   - CoupledResponseEquations: separate α and β equations, shared Coulomb J[ρ]
//   See docs/response_refactor_plan.md, Step 4 / phase 5 notes.

#include "ResponseVector.hpp"
#include "ResponseManager.hpp"
#include "GroundStateData.hpp"
#include "ResponseDebugLogger.hpp"
#include "functypedefs.h"

constexpr double alpha_factor(const StaticUnrestrictedResponse&) { return -2.0; }

inline real_function_3d
compute_density(World& /*world*/, const StaticUnrestrictedResponse& /*rvec*/,
                const vector_real_function_3d& /*phi0*/) {
    MADNESS_EXCEPTION("compute_density: StaticUnrestrictedResponse not implemented", 1);
}

inline std::vector<poperatorT>
make_bsh_operators(World& /*world*/, const ResponseManager& /*rm*/,
                   double /*freq*/, const Tensor<double>& /*eps*/, int /*n*/,
                   ResponseDebugLogger& /*lgr*/,
                   const StaticUnrestrictedResponse& /* tag */) {
    MADNESS_EXCEPTION("make_bsh_operators: StaticUnrestrictedResponse not implemented", 1);
}

inline vector_real_function_3d
CoupledResponseEquations(World& /*world*/, const GroundStateData& /*g_s*/,
                          const StaticUnrestrictedResponse& /*vecs*/,
                          const vector_real_function_3d& /*v_p*/,
                          const std::vector<poperatorT>& /*bsh*/,
                          const ResponseManager& /*rm*/,
                          ResponseDebugLogger& /*lgr*/) {
    MADNESS_EXCEPTION("CoupledResponseEquations: StaticUnrestrictedResponse not implemented", 1);
}
