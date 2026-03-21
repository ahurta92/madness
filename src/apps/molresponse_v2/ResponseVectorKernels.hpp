#pragma once

// ResponseVectorKernels.hpp — Direct (non-MacroTask) per-vector kernel
// implementations.
//
// Provides compute_ground_exchange<R> and compute_gamma_response<R> as direct
// free-function replacements for the ResponseComputeGroundExchange and
// ResponseComputeGammaX MacroTask wrappers.  Physics is identical; the
// MacroTask scaffolding (per-orbital batching, serialisation, subworld
// dispatch) is removed because molresponse_v2 processes one response vector
// at a time so the overhead has no architectural benefit.
//
// Exchange operator convention (matches response_macrotask.hpp K lambda):
//   Exchange::set_bra_and_ket(bra, ket) applied to f:
//     K(f)[r] = Σ_i bra_i(r) ∫ ket_i(r') f(r') / |r-r'| dr'
//
// No physics changes — results are numerically identical to the macrotask path.

#include <madness/chem/SCF.h>
#include <madness/chem/SCFOperators.h>
#include <madness/chem/projector.h>
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>

#include "ResponseVector.hpp"

using namespace madness;

// ============================================================================
// K(bra, ket) — Exchange operator factory
//
// Returns a ready-to-apply Exchange<double,3> with bra/ket set and the
// algorithm defaulting to multiworld_efficient_row.
//
// Usage:  K(world, phi0, x)(phi0)   reads as   K[φ₀,x](φ₀)
// ============================================================================
[[nodiscard]] inline Exchange<double, 3>
K(madness::World& world,
  const vector_real_function_3d& bra,
  const vector_real_function_3d& ket,
  Exchange<double, 3>::ExchangeAlgorithm alg =
      Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row) {
    constexpr double lo = 1.e-10;
    Exchange<double, 3> ex(world, lo);
    ex.set_bra_and_ket(bra, ket);
    ex.set_algorithm(alg);
    return ex;
}

// ============================================================================
// compute_ground_exchange<R>
//
// Applies K[φ₀,φ₀] to each response channel independently.
//
// No Q-projection (mirrors the commented-out Q(g) in
// ResponseComputeGroundExchange).
//
// The Exchange operator's bra/ket (φ₀, N functions) must match the size of
// the vector it acts on, so each channel (x_alpha, y_alpha) is processed
// separately — not as a single flat application.
//
// Returns R with:
//   TDA/Static (x-only):  result.x_alpha = K[φ₀,φ₀](x_alpha)
//   TDDFT (x+y):          result.x_alpha = K[φ₀,φ₀](x_alpha)
//                          result.y_alpha = K[φ₀,φ₀](y_alpha)
// ============================================================================
template <typename R>
[[nodiscard]] R
compute_ground_exchange(madness::World& world,
                        const R& response,
                        const vector_real_function_3d& phi0) {
    if constexpr (response_is_unrestricted_v<R>) {
        MADNESS_EXCEPTION(
            "compute_ground_exchange: unrestricted not yet implemented", 1);
    }

    R result(response.num_orbitals());
    if (phi0.empty()) return result;

    auto k0 = K(world, phi0, phi0);
    result.x_alpha = k0(response.x_alpha);
    if constexpr (response_has_y_channel_v<R>) {
        result.y_alpha = k0(response.y_alpha);
    }
    result.flatten();
    return result;
}

// ============================================================================
// compute_gamma_response<R>
//
// Computes the full response XC coupling vector, Q-projected.
//
// Per orbital p, three branches by response type:
//
//   TDA (excited-state, y≡0) — legacy compute_gamma_tda:
//     ρ = Σ x_i·φ₀_i   (no spin factor; 2× on J)
//     gx[p] = 2·J[ρ]·φ₀_p − c_xc·K[φ₀,x](φ₀_p)
//
//   Static (freq-dependent, y≡x) — legacy compute_gamma_static:
//     ρ = 2·Σ x_i·φ₀_i  (factor-of-2 from y≡x)
//     gx[p] = 2·J[ρ]·φ₀_p − c_xc·(K[φ₀,x](φ₀_p) + K[x,φ₀](φ₀_p))
//
//   TDDFT x-channel:
//     ρ_channel = Σ(x_i+y_i)·φ₀_i  (shared with y-channel — computed once)
//     K_a = K[φ₀,x](φ₀_p)   set_bra_and_ket(φ₀,  x)
//     K_b = K[y,φ₀](φ₀_p)   set_bra_and_ket(y,   φ₀)
//
//   TDDFT y-channel:
//     (same ρ_channel and J)
//     K_a = K[φ₀,y](φ₀_p)   set_bra_and_ket(φ₀,  y)
//     K_b = K[x,φ₀](φ₀_p)   set_bra_and_ket(x,   φ₀)
//
// For TDDFT, J[ρ_channel] is computed once and shared — this avoids the 2N
// independent Coulomb evaluations that the per-orbital macrotask performed.
//
// Returns:
//   TDA:           N functions  (Q-projected)
//   Static:        N functions  (Q-projected)
//   TDDFT:        2N functions  in flat [gx|gy] layout (Q-projected)
// ============================================================================
template <typename R>
[[nodiscard]] vector_real_function_3d
compute_gamma_response(madness::World& world,
                       const R& response,
                       const vector_real_function_3d& phi0,
                       const madness::QProjector<double, 3>& Q,
                       double c_xc = 1.0) {
    if constexpr (response_is_unrestricted_v<R>) {
        MADNESS_EXCEPTION(
            "compute_gamma_response: unrestricted not yet implemented", 1);
    }

    if (phi0.empty()) return {};

    constexpr double lo    = 1.e-10;
    const double     thresh = FunctionDefaults<3>::get_thresh();
    const auto&      x     = response.x_alpha;

    if constexpr (!response_has_y_channel_v<R>) {
        auto xphi = mul(world, x, phi0, true);

        if constexpr (std::is_same_v<R, TDARestrictedResponse>) {
            // ---- TDA (excited-state, y≡0) ----
            // ρ = Σ x_i·φ_i, single K
            auto rho   = sum(world, xphi, true);
            auto J_rho = apply(CoulombOperator(world, lo, thresh), rho);

            auto gx = 2.0 * (J_rho * phi0) - c_xc * K(world, phi0, x)(phi0);
            return Q(gx);

        } else {
            // ---- Static (freq-dependent, y≡x) ----
            // ρ = 2·Σ x_i·φ_i, both K terms
            auto rho   = 2.0 * sum(world, xphi, true);
            auto J_rho = apply(CoulombOperator(world, lo, thresh), rho);

            //  K[φ₀,x](φ₀) + K[x,φ₀](φ₀)
            auto gx = 2.0 * (J_rho * phi0)
                     - c_xc * (K(world, phi0, x)(phi0) + K(world, x, phi0)(phi0));
            return Q(gx);
        }

    } else {
        // ---- TDDFT (x and y channels) ----
        const auto& y = response.y_alpha;

        // ρ¹ = Σ(x_i+y_i)·φ₀_i — computed once, shared across both channels
        auto xphi  = mul(world, x, phi0, true);
        auto yphi  = mul(world, y, phi0, true);
        auto rho1  = sum(world, xphi, true) + sum(world, yphi, true);
        auto J_rho = apply(CoulombOperator(world, lo, thresh), rho1);
        auto Jphi  = J_rho * phi0;  // N functions (shared)

        // x-channel:  K[φ₀,x](φ₀) + K[y,φ₀](φ₀)
        auto gx = 2.0 * Jphi
                 - c_xc * (K(world, phi0, x)(phi0) + K(world, y, phi0)(phi0));

        // y-channel:  K[φ₀,y](φ₀) + K[x,φ₀](φ₀)
        auto gy = 2.0 * Jphi
                 - c_xc * (K(world, phi0, y)(phi0) + K(world, x, phi0)(phi0));

        // Assemble flat [gx | gy], Q-project
        gx.insert(gx.end(), gy.begin(), gy.end());
        return Q(gx);
    }
}
