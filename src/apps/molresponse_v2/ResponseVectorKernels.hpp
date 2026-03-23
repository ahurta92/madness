#pragma once

// ResponseVectorKernels.hpp — exchange operator factory and response XC kernels.
//
// The K() factory is used by both the per-type ops headers (ops/*.hpp) and the
// excited-state path.
//
// compute_ground_exchange<R> and compute_gamma_response<R> are retained here
// because compute_response_potentials<R> in ResponseKernels.hpp uses them for
// the excited-state bundle iteration.  The linear-response ops headers
// (ops/*.hpp) inline the same math directly instead of calling these helpers,
// which keeps each response type self-contained for the BSH iteration loop.

#include <madness/chem/SCF.h>
#include <madness/chem/SCFOperators.h>
#include <madness/chem/projector.h>
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>

#include "ResponseSolverConstants.hpp"
#include "ResponseVector.hpp"

using namespace madness;

// ============================================================================
// K(bra, ket) — Exchange operator factory
//
// Returns a ready-to-apply Exchange<double,3> with bra/ket set.
// Usage:  K(world, phi0, x)(phi0)   reads as   K[φ₀,x](φ₀)
// ============================================================================
[[nodiscard]] inline Exchange<double, 3>
K(madness::World& world,
  const vector_real_function_3d& bra,
  const vector_real_function_3d& ket,
  Exchange<double, 3>::ExchangeAlgorithm alg =
      Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row) {
    Exchange<double, 3> ex(world, ResponseSolverConstants::k_coulomb_lo);
    ex.set_bra_and_ket(bra, ket);
    ex.set_algorithm(alg);
    return ex;
}

// ============================================================================
// compute_ground_exchange<R>
//
// Applies K[φ₀,φ₀] to each response channel independently.
// Used by compute_response_potentials<R> in ResponseKernels.hpp (excited-state path).
// The linear-response ops headers inline this logic directly.
//
// Returns R with:
//   x-only (Static/TDA):  result.x_alpha = K[φ₀,φ₀](x_alpha)
//   x+y    (TDDFT):       result.x_alpha = K[φ₀,φ₀](x_alpha)
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
// Computes the response XC coupling vector, Q-projected.
// Used by compute_response_potentials<R> in ResponseKernels.hpp (excited-state path).
// The linear-response ops headers inline this logic directly.
//
// Branches:
//   TDA  (TDARestrictedResponse):      ρ = Σ xᵢφᵢ,           single K
//   Static (StaticRestrictedResponse): ρ = 2·Σ xᵢφᵢ,         both K terms
//   TDDFT (DynamicRestrictedResponse): ρ = Σ(xᵢ+yᵢ)φᵢ, shared J, per-channel K
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

    const double thresh = FunctionDefaults<3>::get_thresh();
    const auto&      x     = response.x_alpha;

    if constexpr (!response_has_y_channel_v<R>) {
        auto xphi = mul(world, x, phi0, true);

        if constexpr (std::is_same_v<R, TDARestrictedResponse>) {
            // TDA: ρ = Σ xᵢφᵢ, single K term
            auto rho   = sum(world, xphi, true);
            auto J_rho = apply(CoulombOperator(world, ResponseSolverConstants::k_coulomb_lo, thresh), rho);
            auto gx = 2.0 * (J_rho * phi0) - c_xc * K(world, phi0, x)(phi0);
            return Q(gx);
        } else {
            // Static: ρ = 2·Σ xᵢφᵢ, both K terms
            auto rho   = 2.0 * sum(world, xphi, true);
            auto J_rho = apply(CoulombOperator(world, ResponseSolverConstants::k_coulomb_lo, thresh), rho);
            auto gx = 2.0 * (J_rho * phi0)
                     - c_xc * (K(world, phi0, x)(phi0) + K(world, x, phi0)(phi0));
            return Q(gx);
        }
    } else {
        // TDDFT: ρ = Σ(xᵢ+yᵢ)φᵢ, shared J, per-channel asymmetric K
        const auto& y = response.y_alpha;
        auto xphi  = mul(world, x, phi0, true);
        auto yphi  = mul(world, y, phi0, true);
        auto rho1  = sum(world, xphi, true) + sum(world, yphi, true);
        auto J_rho = apply(CoulombOperator(world, ResponseSolverConstants::k_coulomb_lo, thresh), rho1);
        auto Jphi  = J_rho * phi0;

        auto gx = 2.0 * Jphi
                 - c_xc * (K(world, phi0, x)(phi0) + K(world, y, phi0)(phi0));
        auto gy = 2.0 * Jphi
                 - c_xc * (K(world, phi0, y)(phi0) + K(world, x, phi0)(phi0));
        gx.insert(gx.end(), gy.begin(), gy.end());
        return Q(gx);
    }
}
