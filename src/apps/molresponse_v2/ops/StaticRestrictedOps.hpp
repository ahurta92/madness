#pragma once
// ops/StaticRestrictedOps.hpp
//
// Self-contained ops module for static (ω=0) restricted (closed-shell) response.
//
// ---- Equation solved --------------------------------------------------------
//   (A)·x = v_p       (B-matrix dropped: y ≡ 0 in the ω→0 limit)
//
//   θ_p = −2 [ V₀ x_p  −  ε_{ip} x_i  +  Γˢᵗᵃᵗⁱᶜ[x]_p  +  v_p^C ]
//
//   Ground-state potential acting on x:
//     V₀ x_p  =  V_local · x_p  −  c_xc · K[φ₀,φ₀](x_p)
//
//   Off-diagonal Fock coupling:
//     ε_{ip} x_i  =  transform(x, H_no_diag)
//
//   Response coupling (static kernel — both exchange terms, factor-2 density):
//     ρ  =  2 · Σᵢ xᵢ · φᵢ               (y ≡ x in ω→0 limit)
//     Γˢᵗᵃᵗⁱᶜ[x]  =  Q̂ · [ 2·J[ρ]·φ  −  c_xc·(K[φ₀,x]·φ + K[x,φ₀]·φ) ]
//
//   BSH update:   x_p^new  =  G(μ_p) · θ_p,   μ_p = √(−2·ε_p)
//   Final:        result   =  Q̂ · x^new
//
// ---- alpha_factor -----------------------------------------------------------
//   α_{AB}(0) = −4 · ⟨x_A | V_B⟩
//   Factor −4 = −2 (spin degeneracy) × 2 (y≡x, B-matrix dropped)
// -----------------------------------------------------------------------------

#include "ResponseVectorKernels.hpp"   // K() factory, CoulombOperator, using namespace madness
#include "ResponseDebugLoggerMacros.hpp"
#include "ResponseSolverConstants.hpp"
#include "ResponseSolverUtils.hpp"
#include "GroundStateData.hpp"
#include "ResponseManager.hpp"
#include "functypedefs.h"

// ---- alpha_factor -----------------------------------------------------------

constexpr double alpha_factor(const StaticRestrictedResponse&) { return -4.0; }

// ---- compute_density --------------------------------------------------------
// ρ¹(r) = 2 · Σᵢ xᵢ(r) · φᵢ*(r)
// Factor 2: spin degeneracy (restricted, α=β).
// Uses x-only density (B-matrix / y channel dropped at ω=0).

inline real_function_3d
compute_density(World& world, const StaticRestrictedResponse& rvec,
                const vector_real_function_3d& phi0) {
    auto xphi = mul(world, rvec.x_alpha, phi0, true);
    return ResponseSolverConstants::k_restricted_spin_factor * sum(world, xphi, true);
}

// ---- make_bsh_operators -----------------------------------------------------
// Returns N operators, one per occupied orbital.
// BSH exponent: μ_p = √(−2·ε_p)   (ω = 0, no frequency shift)
// A level-shift is applied when ε_p + ω ≥ 0 to keep operators bounded.

inline std::vector<poperatorT>
make_bsh_operators(World& world, const ResponseManager& response_manager,
                   const double freq, const Tensor<double>& orbital_energies,
                   const int n, ResponseDebugLogger& /*logger*/,
                   const StaticRestrictedResponse& /* tag */) {
    const double x_shifts = ResponseSolverUtils::compute_bsh_x_shift(orbital_energies, n, freq);
    return ResponseSolverUtils::make_bsh_operators_response(
        world, x_shifts, freq, orbital_energies, response_manager.params().lo());
}

// ---- CoupledResponseEquations -----------------------------------------------
// One BSH iteration step.  All response math for this type lives here.
//
// Ground exchange (no Q-projection):
//   k0x  =  K[φ₀,φ₀](x)
//
// Response coupling (static branch — inlined from compute_gamma_response):
//   ρ    =  2·Σ xᵢ·φᵢ
//   g_x  =  Q̂·[ 2·J[ρ]·φ  −  c_xc·(K[φ₀,x]·φ + K[x,φ₀]·φ) ]

inline vector_real_function_3d
CoupledResponseEquations(World& world, const GroundStateData& g_s,
                          const StaticRestrictedResponse& vecs,
                          const vector_real_function_3d& v_p,
                          const std::vector<poperatorT>& bsh_x,
                          const ResponseManager& /*response_manager*/,
                          ResponseDebugLogger& logger) {
    const auto   c_xc  = g_s.xcf_.hf_exchange_coefficient();
    const auto&  x     = vecs.x_alpha;
    const auto&  phi0  = g_s.orbitals;

    // --- Ground exchange: K₀[x] = K[φ₀,φ₀](x) ---
    vector_real_function_3d k0x;
    DEBUG_TIMED_BLOCK(world, &logger, "g0_task",
                      { k0x = K(world, phi0, phi0)(x); });

    // --- Response coupling (static kernel) ---
    // ρ = 2·Σ xᵢ·φᵢ   (y ≡ x at ω=0, factor 2 from spin degeneracy)
    // Γˢᵗᵃᵗⁱᶜ[x] = Q̂·[ 2·J[ρ]·φ  −  c_xc·(K[φ₀,x]·φ + K[x,φ₀]·φ) ]
    vector_real_function_3d g_x;
    DEBUG_TIMED_BLOCK(world, &logger, "gx_task", {
        const double thresh = FunctionDefaults<3>::get_thresh();
        auto xphi  = mul(world, x, phi0, true);
        auto rho   = ResponseSolverConstants::k_restricted_spin_factor * sum(world, xphi, true);
        auto J_rho = apply(CoulombOperator(world, ResponseSolverConstants::k_coulomb_lo, thresh), rho);
        g_x = g_s.Qhat(2.0 * (J_rho * phi0)
                       - c_xc * (K(world, phi0, x)(phi0) + K(world, x, phi0)(phi0)));
    });

    // --- Residual and BSH update ---
    auto v_local  = g_s.V_local * x;
    auto v0x      = v_local - c_xc * k0x;
    auto epsilonx = transform(world, x, g_s.Hamiltonian_no_diag, true);

    auto thetax = ResponseSolverConstants::k_bsh_residual_prefactor * (v0x - epsilonx + g_x + v_p);
    truncate(world, thetax);
    auto rsh = apply(world, bsh_x, thetax);
    rsh = g_s.Qhat(rsh);  // project out occupied space
    return rsh;
}
