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
// compute_ground_exchange<R>
//
// Applies K[φ₀,φ₀] to all active response channels simultaneously.
//
// No Q-projection (mirrors the commented-out Q(g) in
// ResponseComputeGroundExchange).
//
// Returns:
//   TDA/Static (x-only):    N functions  — K[φ₀,φ₀](x_p) for p in [0..N-1]
//   TDDFT (x+y channels):  2N functions  — K[φ₀,φ₀](flat_p) for p in [0..2N-1]
// ============================================================================
template <typename R>
[[nodiscard]] vector_real_function_3d
compute_ground_exchange(madness::World& world,
                        const R& response,
                        const vector_real_function_3d& phi0) {
    if constexpr (response_is_unrestricted_v<R>) {
        MADNESS_EXCEPTION(
            "compute_ground_exchange: unrestricted not yet implemented", 1);
    }

    if (phi0.empty()) return {};

    constexpr double lo = 1.e-10;
    Exchange<double, 3> k0(world, lo);
    k0.set_bra_and_ket(phi0, phi0);
    k0.set_algorithm(
        Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);

    if constexpr (!response_has_y_channel_v<R>) {
        return k0(response.x_alpha);  // N functions
    } else {
        return k0(response.flat);     // 2N functions
    }
}

// ============================================================================
// compute_gamma_response<R>
//
// Computes the full response XC coupling vector, Q-projected.
//
// Per orbital p:
//   gx[p] = 2·J[ρ_channel]·φ₀_p  −  K_a(φ₀_p)  −  K_b(φ₀_p)
//
// Channel definitions and set_bra_and_ket arguments:
//
//   TDA/Static (x-only):
//     ρ_channel = 2·Σ x_i·φ₀_i
//     K_a = K[φ₀,x](φ₀_p)   set_bra_and_ket(φ₀,  x)
//     K_b = K[x,φ₀](φ₀_p)   set_bra_and_ket(x,   φ₀)
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
//   TDA/Static:    N functions  (Q-projected)
//   TDDFT:        2N functions  in flat [gx|gy] layout (Q-projected)
// ============================================================================
template <typename R>
[[nodiscard]] vector_real_function_3d
compute_gamma_response(madness::World& world,
                       const R& response,
                       const vector_real_function_3d& phi0,
                       const madness::QProjector<double, 3>& Q) {
    if constexpr (response_is_unrestricted_v<R>) {
        MADNESS_EXCEPTION(
            "compute_gamma_response: unrestricted not yet implemented", 1);
    }

    if (phi0.empty()) return {};

    constexpr double lo = 1.e-10;
    const double thresh = FunctionDefaults<3>::get_thresh();

    if constexpr (!response_has_y_channel_v<R>) {
        // ---- TDA / Static (x-only) ----

        // ρ = 2·Σ x_i·φ₀_i
        auto xphi  = mul(world, response.x_alpha, phi0, true);
        auto rho   = 2.0 * sum(world, xphi, true);
        auto J_rho = apply(CoulombOperator(world, lo, thresh), rho);
        auto Jphix = J_rho * phi0;  // N functions

        // K_a = K[φ₀,x](φ₀_p)  →  set_bra_and_ket(φ₀, x)
        Exchange<double, 3> ka_op(world, lo);
        ka_op.set_bra_and_ket(phi0, response.x_alpha);
        ka_op.set_algorithm(
            Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);

        // K_b = K[x,φ₀](φ₀_p)  →  set_bra_and_ket(x, φ₀)
        Exchange<double, 3> kb_op(world, lo);
        kb_op.set_bra_and_ket(response.x_alpha, phi0);
        kb_op.set_algorithm(
            Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);

        auto gx = 2.0 * Jphix - ka_op(phi0) - kb_op(phi0);
        return Q(gx);

    } else {
        // ---- TDDFT (x and y channels) ----

        // ρ¹ = Σ(x_i+y_i)·φ₀_i — computed once, shared across both channels
        auto xphi  = mul(world, response.x_alpha, phi0, true);
        auto yphi  = mul(world, response.y_alpha, phi0, true);
        auto rho1  = sum(world, xphi, true) + sum(world, yphi, true);
        auto J_rho = apply(CoulombOperator(world, lo, thresh), rho1);
        auto Jphix = J_rho * phi0;  // N functions (shared)

        // x-channel:
        //   K_a = K[φ₀,x](φ₀_p)  →  set_bra_and_ket(φ₀, x)
        //   K_b = K[y,φ₀](φ₀_p)  →  set_bra_and_ket(y,  φ₀)
        Exchange<double, 3> kax(world, lo);
        kax.set_bra_and_ket(phi0, response.x_alpha);
        kax.set_algorithm(
            Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
        Exchange<double, 3> kbx(world, lo);
        kbx.set_bra_and_ket(response.y_alpha, phi0);
        kbx.set_algorithm(
            Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
        auto gx = 2.0 * Jphix - kax(phi0) - kbx(phi0);

        // y-channel:
        //   K_a = K[φ₀,y](φ₀_p)  →  set_bra_and_ket(φ₀, y)
        //   K_b = K[x,φ₀](φ₀_p)  →  set_bra_and_ket(x,  φ₀)
        Exchange<double, 3> kay(world, lo);
        kay.set_bra_and_ket(phi0, response.y_alpha);
        kay.set_algorithm(
            Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
        Exchange<double, 3> kby(world, lo);
        kby.set_bra_and_ket(response.x_alpha, phi0);
        kby.set_algorithm(
            Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
        auto gy = 2.0 * Jphix - kay(phi0) - kby(phi0);

        // Assemble flat [gx | gy], Q-project
        gx.insert(gx.end(), gy.begin(), gy.end());
        return Q(gx);
    }
}
