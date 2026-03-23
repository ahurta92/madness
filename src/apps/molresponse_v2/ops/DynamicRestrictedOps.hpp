#pragma once
// ops/DynamicRestrictedOps.hpp
//
// Self-contained ops module for dynamic (ω≠0) restricted (closed-shell) response
// — full Casida / TDDFT coupled equations.
//
// ---- Equations solved -------------------------------------------------------
//   [ A   B  ] [ x ]       [ v_p   ]
//   [ B*  A* ] [ y ]  = ω  [ v_p*  ]
//
//   Residual per orbital p in flat [x|y] layout (2N slots):
//     θ_p^x = −2 [ V₀ x_p  −  ε_{ip} x_i  +  Γˣ[x,y]_p  +  v_p^C   ]
//     θ_p^y = −2 [ V₀ y_p  −  ε_{ip} y_i  +  Γʸ[x,y]_p  +  v_p^{C*} ]
//
//   Ground-state potential (same for both channels):
//     V₀ f_p  =  V_local · f_p  −  c_xc · K[φ₀,φ₀](f_p)
//
//   Response coupling — TDDFT kernel (shared density, asymmetric exchange):
//     ρ¹  =  Σᵢ (xᵢ + yᵢ) · φᵢ                      (computed once, shared)
//     J   =  J[ρ¹] · φ                                (N functions, shared)
//     Γˣ[x,y]  =  Q̂ · [ 2·J  −  c_xc·(K[φ₀,x]·φ + K[y,φ₀]·φ) ]
//     Γʸ[x,y]  =  Q̂ · [ 2·J  −  c_xc·(K[φ₀,y]·φ + K[x,φ₀]·φ) ]
//   Assembled flat: g = [Γˣ | Γʸ]
//
//   BSH updates:
//     x_p^new  =  G(μ_p^x) · θ_p^x,   μ_p^x = √(−2·(ε_p + ω))
//     y_p^new  =  G(μ_p^y) · θ_p^y,   μ_p^y = √(−2·(ε_p − ω))
//   Returns 2N operators: [bsh_x | bsh_y]
//
// ---- alpha_factor -----------------------------------------------------------
//   α_{AB}(ω) = −2 · ⟨x_A | V_B⟩
//   Factor −2: spin degeneracy (restricted); both x and y channels contribute,
//   halving the prefactor relative to the x-only (TDA/static) case.
// -----------------------------------------------------------------------------

#include "ResponseVectorKernels.hpp"   // K() factory, CoulombOperator, using namespace madness
#include "ResponseDebugLoggerMacros.hpp"
#include "ResponseSolverConstants.hpp"
#include "ResponseSolverUtils.hpp"
#include "GroundStateData.hpp"
#include "ResponseManager.hpp"
#include "functypedefs.h"

// ---- alpha_factor -----------------------------------------------------------

constexpr double alpha_factor(const DynamicRestrictedResponse&) { return -2.0; }

// ---- compute_density --------------------------------------------------------
// ρ¹(r) = 2 · Σᵢ [xᵢ(r) + yᵢ(r)] · φᵢ*(r)
// Flat layout [x|y] is multiplied against doubled φ₀ to cover both channels.
// Factor 2: spin degeneracy (restricted).

inline real_function_3d
compute_density(World& world, const DynamicRestrictedResponse& rvec,
                const vector_real_function_3d& phi0) {
    auto phi_phi = phi0;
    phi_phi.insert(phi_phi.end(), phi0.begin(), phi0.end());
    auto xphi = mul(world, rvec.flat, phi_phi, true);
    return ResponseSolverConstants::k_restricted_spin_factor * sum(world, xphi, true);
}

// ---- make_bsh_operators -----------------------------------------------------
// Returns 2N operators: [bsh_x[0..N-1] | bsh_y[0..N-1]]
// x-channel BSH exponent: μ_p^x = √(−2·(ε_p + ω))
// y-channel BSH exponent: μ_p^y = √(−2·(ε_p − ω))
// Level-shift applied to x-channel when ε_p + ω ≥ 0.
// No shift on y-channel: ε_p − ω < 0 is guaranteed for bound states below ω.

inline std::vector<poperatorT>
make_bsh_operators(World& world, const ResponseManager& response_manager,
                   const double freq, const Tensor<double>& orbital_energies,
                   const int n, ResponseDebugLogger& /*logger*/,
                   const DynamicRestrictedResponse& /* tag */) {
    const double x_shifts = ResponseSolverUtils::compute_bsh_x_shift(orbital_energies, n, freq);
    auto bsh_x = ResponseSolverUtils::make_bsh_operators_response(
        world, x_shifts, freq, orbital_energies, response_manager.params().lo());
    auto bsh_y = ResponseSolverUtils::make_bsh_operators_response(
        world, 0.0, -freq, orbital_energies, response_manager.params().lo());
    bsh_x.insert(bsh_x.end(), bsh_y.begin(), bsh_y.end());
    return bsh_x;
}

// ---- CoupledResponseEquations -----------------------------------------------
// One BSH iteration step.  All response math for this type lives here.
//
// Ground exchange (applied independently per channel, no Q-projection):
//   k0x  =  K[φ₀,φ₀](x)
//   k0y  =  K[φ₀,φ₀](y)
//   k0_flat = [k0x | k0y]
//
// Response coupling (TDDFT kernel — inlined from compute_gamma_response TDDFT branch):
//   ρ¹   =  Σᵢ(xᵢ+yᵢ)·φᵢ   (shared density, one Coulomb evaluation)
//   J    =  J[ρ¹]·φ           (N functions, reused for both channels)
//   Γˣ   =  Q̂·[2J − c_xc·(K[φ₀,x]·φ + K[y,φ₀]·φ)]
//   Γʸ   =  Q̂·[2J − c_xc·(K[φ₀,y]·φ + K[x,φ₀]·φ)]
//   g_x  =  [Γˣ | Γʸ]   (flat 2N)

inline vector_real_function_3d
CoupledResponseEquations(World& world, const GroundStateData& g_s,
                          const DynamicRestrictedResponse& vecs,
                          const vector_real_function_3d& v_p,
                          const std::vector<poperatorT>& bsh_x,
                          const ResponseManager& /*response_manager*/,
                          ResponseDebugLogger& logger) {
    const auto  c_xc  = g_s.xcf_.hf_exchange_coefficient();
    const auto& xvec  = vecs.x_alpha;
    const auto& yvec  = vecs.y_alpha;
    const auto& all_x = vecs.flat;
    const auto& phi0  = g_s.orbitals;

    // --- Ground exchange: K[φ₀,φ₀] applied per channel ---
    vector_real_function_3d k0_flat;
    DEBUG_TIMED_BLOCK(world, &logger, "g0_task", {
        auto k0  = K(world, phi0, phi0);
        auto k0x = k0(xvec);
        auto k0y = k0(yvec);
        k0_flat  = k0x;
        k0_flat.insert(k0_flat.end(), k0y.begin(), k0y.end());
    });

    // --- Response coupling (TDDFT kernel) ---
    // ρ¹ = Σ(xᵢ+yᵢ)·φᵢ — one Coulomb evaluation, shared by both channels
    // x-channel: Γˣ = Q̂·[2J − c_xc·(K[φ₀,x]·φ + K[y,φ₀]·φ)]
    // y-channel: Γʸ = Q̂·[2J − c_xc·(K[φ₀,y]·φ + K[x,φ₀]·φ)]
    vector_real_function_3d g_x;
    DEBUG_TIMED_BLOCK(world, &logger, "gx_task", {
        const double thresh = FunctionDefaults<3>::get_thresh();
        auto xphi  = mul(world, xvec, phi0, true);
        auto yphi  = mul(world, yvec, phi0, true);
        auto rho1  = sum(world, xphi, true) + sum(world, yphi, true);
        auto J_rho = apply(CoulombOperator(world, ResponseSolverConstants::k_coulomb_lo, thresh), rho1);
        auto Jphi  = J_rho * phi0;                                      // N functions, shared

        auto gx = 2.0 * Jphi
                 - c_xc * (K(world, phi0, xvec)(phi0) + K(world, yvec, phi0)(phi0));
        auto gy = 2.0 * Jphi
                 - c_xc * (K(world, phi0, yvec)(phi0) + K(world, xvec, phi0)(phi0));
        gx.insert(gx.end(), gy.begin(), gy.end());
        g_x = g_s.Qhat(gx);
    });

    // --- Residual and BSH update ---
    auto v_local  = g_s.V_local * all_x;
    auto v0x      = v_local - c_xc * k0_flat;
    auto epsilonx = transform(world, xvec, g_s.Hamiltonian_no_diag, true);
    auto epsilony = transform(world, yvec, g_s.Hamiltonian_no_diag, true);
    epsilonx.insert(epsilonx.end(), epsilony.begin(), epsilony.end());

    auto thetax = ResponseSolverConstants::k_bsh_residual_prefactor * (v0x - epsilonx + g_x + v_p);
    truncate(world, thetax);
    auto rsh = apply(world, bsh_x, thetax);
    rsh = g_s.Qhat(rsh);
    return rsh;
}
