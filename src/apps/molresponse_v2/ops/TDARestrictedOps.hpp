#pragma once
// ops/TDARestrictedOps.hpp
//
// Self-contained ops module for TDA (Tamm-Dancoff Approximation) restricted
// (closed-shell) response — used for excited-state calculations.
//
// ---- Equation solved --------------------------------------------------------
//   (A)·x = ω·x    (B-matrix dropped: y ≡ 0 as an approximation)
//
//   θ_p = −2 [ V₀ x_p  −  ε_{ip} x_i  +  Γᵀᴰᴬ[x]_p  +  v_p^C ]
//
//   Ground-state potential acting on x:
//     V₀ x_p  =  V_local · x_p  −  c_xc · K[φ₀,φ₀](x_p)
//
//   Off-diagonal Fock coupling:
//     ε_{ip} x_i  =  transform(x, H_no_diag)
//
//   Response coupling (TDA kernel — half-density, single exchange term):
//     ρ  =  Σᵢ xᵢ · φᵢ                       (no ×2 spin factor — TDA convention)
//     Γᵀᴰᴬ[x]  =  Q̂ · [ 2·J[ρ]·φ  −  c_xc · K[φ₀,x]·φ ]
//
//   Contrast with static kernel: ρ = 2·Σ xᵢφᵢ and two K terms (K[φ₀,x] + K[x,φ₀]).
//   TDA drops both the factor-2 density and the K[x,φ₀] cross term.
//
//   BSH update:   x_p^new  =  G(μ_p) · θ_p,   μ_p = √(−2·(ε_p + ω))
//   Final:        result   =  Q̂ · x^new
//
// ---- alpha_factor -----------------------------------------------------------
//   α_{AB}(ω) = −4 · ⟨x_A | V_B⟩
//   Same prefactor as StaticRestricted: spin factor −2 × y≡0 factor 2.
//
// ---- compute_density note ---------------------------------------------------
//   The convergence-tracking density uses the full spin factor:
//     ρ¹ = 2 · Σᵢ xᵢ · φᵢ
//   This differs from the XC coupling density above (no ×2) by convention:
//   the TDA kernel uses a half-density for the response coupling, but the
//   observable density (for polarizability and convergence) includes the spin
//   factor as normal.
// -----------------------------------------------------------------------------

#include "ResponseVectorKernels.hpp"   // K() factory, CoulombOperator, using namespace madness
#include "ResponseDebugLoggerMacros.hpp"
#include "ResponseSolverConstants.hpp"
#include "ResponseSolverUtils.hpp"
#include "GroundStateData.hpp"
#include "ResponseManager.hpp"
#include "functypedefs.h"

// ---- alpha_factor -----------------------------------------------------------

constexpr double alpha_factor(const TDARestrictedResponse&) { return -4.0; }

// ---- compute_density --------------------------------------------------------
// ρ¹(r) = 2 · Σᵢ xᵢ(r) · φᵢ*(r)
// Factor 2: spin degeneracy (restricted).  x-only channel (y ≡ 0).

inline real_function_3d
compute_density(World& world, const TDARestrictedResponse& rvec,
                const vector_real_function_3d& phi0) {
    auto xphi = mul(world, rvec.x_alpha, phi0, true);
    return ResponseSolverConstants::k_restricted_spin_factor * sum(world, xphi, true);
}

// ---- make_bsh_operators -----------------------------------------------------
// Returns N operators, one per occupied orbital.
// BSH exponent: μ_p = √(−2·(ε_p + ω))   (ω = excitation energy estimate)
// Level-shift applied when ε_p + ω ≥ 0.

inline std::vector<poperatorT>
make_bsh_operators(World& world, const ResponseManager& response_manager,
                   const double freq, const Tensor<double>& orbital_energies,
                   const int n, ResponseDebugLogger& /*logger*/,
                   const TDARestrictedResponse& /* tag */) {
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
// Response coupling (TDA kernel — inlined from compute_gamma_response TDA branch):
//   ρ    =  Σ xᵢ·φᵢ   (no ×2 — TDA convention; cf. static: ρ = 2·Σ xᵢφᵢ)
//   g_x  =  Q̂·[ 2·J[ρ]·φ  −  c_xc · K[φ₀,x]·φ ]
//   (single exchange term; cf. static: K[φ₀,x] + K[x,φ₀])

inline vector_real_function_3d
CoupledResponseEquations(World& world, const GroundStateData& g_s,
                          const TDARestrictedResponse& vecs,
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

    // --- Response coupling (TDA kernel) ---
    // ρ = Σ xᵢ·φᵢ   (no ×2 — TDA uses half-density for XC coupling)
    // Γᵀᴰᴬ[x] = Q̂·[ 2·J[ρ]·φ  −  c_xc · K[φ₀,x]·φ ]
    vector_real_function_3d g_x;
    DEBUG_TIMED_BLOCK(world, &logger, "gx_task", {
        const double thresh = FunctionDefaults<3>::get_thresh();
        auto xphi  = mul(world, x, phi0, true);
        auto rho   = sum(world, xphi, true);               // no ×2 — TDA convention
        auto J_rho = apply(CoulombOperator(world, ResponseSolverConstants::k_coulomb_lo, thresh), rho);
        g_x = g_s.Qhat(2.0 * (J_rho * phi0) - c_xc * K(world, phi0, x)(phi0));
    });                                                     // single K term

    // --- Residual and BSH update ---
    auto v_local  = g_s.V_local * x;
    auto v0x      = v_local - c_xc * k0x;
    auto epsilonx = transform(world, x, g_s.Hamiltonian_no_diag, true);

    auto thetax = ResponseSolverConstants::k_bsh_residual_prefactor * (v0x - epsilonx + g_x + v_p);
    truncate(world, thetax);
    auto rsh = apply(world, bsh_x, thetax);
    rsh = g_s.Qhat(rsh);
    return rsh;
}
