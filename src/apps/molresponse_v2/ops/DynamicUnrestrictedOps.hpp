#pragma once
// ops/DynamicUnrestrictedOps.hpp
//
// Ops module for dynamic (ω≠0) unrestricted (open-shell) response — full TDDFT.
//
// Status: NOT YET IMPLEMENTED (phase 5+).
//   All functions throw std::runtime_error at runtime.
//   Compile-time dispatch via overload resolution is fully wired up.
//
// When implementing:
//   - alpha/beta spin channels: flat layout [xα | yα | xβ | yβ] → 4N slots
//   - alpha_factor: −2.0
//   - make_bsh_operators: 4N operators [bsh_xα | bsh_yα | bsh_xβ | bsh_yβ]
//       μ^x = √(−2·(ε+ω)), μ^y = √(−2·(ε−ω)) per spin channel
//   - CoupledResponseEquations: coupled α and β equations with shared Coulomb J[ρ¹]
//       ρ¹ = Σ[(xα+yα)·φα + (xβ+yβ)·φβ]
//   See docs/response_refactor_plan.md, Step 4 / phase 5 notes.

#include "ResponseVector.hpp"
#include "ResponseManager.hpp"
#include "GroundStateData.hpp"
#include "ResponseDebugLogger.hpp"
#include "functypedefs.h"

constexpr double alpha_factor(const DynamicUnrestrictedResponse&) { return -2.0; }

inline real_function_3d
compute_density(World& /*world*/, const DynamicUnrestrictedResponse& /*rvec*/,
                const vector_real_function_3d& /*phi0*/) {
    MADNESS_EXCEPTION("compute_density: DynamicUnrestrictedResponse not implemented", 1);
}

inline std::vector<poperatorT>
make_bsh_operators(World& /*world*/, const ResponseManager& /*rm*/,
                   double /*freq*/, const Tensor<double>& /*eps*/, int /*n*/,
                   ResponseDebugLogger& /*lgr*/,
                   const DynamicUnrestrictedResponse& /* tag */) {
    MADNESS_EXCEPTION("make_bsh_operators: DynamicUnrestrictedResponse not implemented", 1);
}

inline vector_real_function_3d
CoupledResponseEquations(World& /*world*/, const GroundStateData& /*g_s*/,
                          const DynamicUnrestrictedResponse& /*vecs*/,
                          const vector_real_function_3d& /*v_p*/,
                          const std::vector<poperatorT>& /*bsh*/,
                          const ResponseManager& /*rm*/,
                          ResponseDebugLogger& /*lgr*/) {
    MADNESS_EXCEPTION("CoupledResponseEquations: DynamicUnrestrictedResponse not implemented", 1);
}
