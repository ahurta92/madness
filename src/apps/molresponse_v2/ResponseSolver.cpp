#include "ResponseSolver.hpp"

#include "ResponseDebugLogger.hpp"
#include "ResponseVectorKernels.hpp"
#include "ResponseDebugLoggerMacros.hpp"
#include "ResponseManager.hpp"
#include "ResponseSolverUtils.hpp"
#include "ResponseVector.hpp"
#include "functypedefs.h"
#include "projector.h"

const double restricted_density_factor = 2.0;

// ============================================================================
// compute_density — response density from x (and y) channels
// ============================================================================

// StaticRestrictedResponse (ω=0 frequency response, x-only channel):
//   ρ¹(r) = 2 Σ_i x_i(r) φ_i*(r)
// Factor 2: spin degeneracy (restricted, α=β), so total ρ = 2·ρ_α.
// Used for convergence tracking and polarizability assembly.
// NOTE: Uses TDA-like x-only density (B-matrix dropped). For the full ω=0
//       static polarizability including B, use DynamicRestrictedResponse at ω=0.
real_function_3d compute_density(World& world, const StaticRestrictedResponse& rvec,
                                 const vector_real_function_3d& phi0) {
    auto xphi = mul(world, rvec.x_alpha, phi0, true);
    return restricted_density_factor * sum(world, xphi, true);
}

// TDARestrictedResponse (TDA excited-state, x-only channel):
//   ρ¹(r) = 2 Σ_i x_i(r) φ_i*(r)
// Same x-only kernel as StaticRestrictedResponse. Kept as a separate overload
// so the type system enforces semantic distinction between ω=0 and TDA.
real_function_3d compute_density(World& world, const TDARestrictedResponse& rvec,
                                 const vector_real_function_3d& phi0) {
    auto xphi = mul(world, rvec.x_alpha, phi0, true);
    return restricted_density_factor * sum(world, xphi, true);
}

// DynamicRestrictedResponse (full TDDFT, x and y channels):
//   ρ¹(r) = 2 Σ_i [x_i(r) + y_i(r)] φ_i*(r)
// Stored in flat layout [x|y], so the doubled phi0 matches the 2N flat slots.
// Factor 2: spin degeneracy (restricted).
real_function_3d compute_density(World& world, const DynamicRestrictedResponse& rvec,
                                 const vector_real_function_3d& phi0) {
    auto phi_phi = phi0;
    phi_phi.insert(phi_phi.end(), phi0.begin(), phi0.end());

    auto xphi = mul(world, rvec.flat, phi_phi, true);
    return restricted_density_factor * sum(world, xphi, true);
}

real_function_3d compute_density(World& world, const StaticUnrestrictedResponse& rvec,
                                 const vector_real_function_3d& phi0) {
    throw std::runtime_error("compute_density for StaticUnrestrictedResponse not implemented");
    // auto phi_phi = phi0;
    // phi_phi.insert(phi_phi.end(), phi0.begin(), phi0.end());
    // auto xphi = mul(world, rvec.flat, phi_phi, true);
    // return 2.0 * sum(world, xphi, true);
}

real_function_3d compute_density(World& world, const DynamicUnrestrictedResponse& rvec,
                                 const vector_real_function_3d& phi0) {
    throw std::runtime_error("compute_density for DynamicUnrestrictedResponse not implemented");
    // auto phi_phi = phi0;
    // phi_phi.insert(phi_phi.end(), phi0.begin(), phi0.end());
    // auto &x = rvec.flat;
    //
    // auto xphi = mul(world, x, phi_phi, true);
}

// ============================================================================
// make_bsh_operators — build per-orbital BSH Green's functions
// ============================================================================
//
// StaticRestrictedResponse (ω=0, x-channel only):
//   BSH exponent per orbital p:  μ_p = √(-2·ε_p)   (ω=0, no frequency shift)
//   Returns N operators (one per occupied orbital).
//   A level-shift is applied when any ε_p + ω ≥ 0 to keep operators bounded.
std::vector<poperatorT> make_bsh_operators(World& world, const ResponseManager& response_manager, const double freq,
                                           const Tensor<double>& orbital_energies, const int n,
                                           ResponseDebugLogger& logger, const StaticRestrictedResponse& /* vecs */) {
    auto bsh_x = std::vector<poperatorT>(n);
    double x_shifts = 0.0;
    const double shift_factor = 0.05;
    if ((orbital_energies[static_cast<long>(n) - 1] + freq) >= 0.0) {
        x_shifts = -shift_factor - (freq + orbital_energies[static_cast<long>(n) - 1]);
    }
    bsh_x = ResponseSolverUtils::make_bsh_operators_response(world, x_shifts, freq, orbital_energies,
                                                             response_manager.params().lo());
    return bsh_x;
}

// TDARestrictedResponse (TDA excited-state, x-channel only):
//   Same BSH exponent form as StaticRestrictedResponse but called for ω > 0
//   (the excitation energy estimate omega[i] is passed as freq).
//   Returns N operators (one per occupied orbital).
std::vector<poperatorT> make_bsh_operators(World& world, const ResponseManager& response_manager, const double freq,
                                           const Tensor<double>& orbital_energies, const int n,
                                           ResponseDebugLogger& logger, const TDARestrictedResponse& /* vecs */) {
    auto bsh_x = std::vector<poperatorT>(n);
    double x_shifts = 0.0;
    const double shift_factor = 0.05;
    if ((orbital_energies[static_cast<long>(n) - 1] + freq) >= 0.0) {
        x_shifts = -shift_factor - (freq + orbital_energies[static_cast<long>(n) - 1]);
    }
    bsh_x = ResponseSolverUtils::make_bsh_operators_response(world, x_shifts, freq, orbital_energies,
                                                             response_manager.params().lo());
    return bsh_x;
}

// DynamicRestrictedResponse (full TDDFT, x and y channels):
//   x-channel BSH exponent:  μ_p^x = √(-2·(ε_p + ω))
//   y-channel BSH exponent:  μ_p^y = √(-2·(ε_p - ω))
//   Returns 2N operators: [bsh_x[0..N-1] | bsh_y[0..N-1]].
std::vector<poperatorT> make_bsh_operators(World& world, const ResponseManager& response_manager, const double freq,
                                           const Tensor<double>& orbital_energies, const int n,
                                           ResponseDebugLogger& logger, const DynamicRestrictedResponse& /* vecs */) {
    auto bsh_x = std::vector<poperatorT>(2 * n);
    double x_shifts = 0.0;
    const double shift_factor = 0.05;
    if ((orbital_energies[static_cast<long>(n) - 1] + freq) >= 0.0) {
        x_shifts = -shift_factor - (freq + orbital_energies[static_cast<long>(n) - 1]);
    }
    bsh_x = ResponseSolverUtils::make_bsh_operators_response(world, x_shifts, freq, orbital_energies,
                                                             response_manager.params().lo());
    auto bsh_y = ResponseSolverUtils::make_bsh_operators_response(world, 0.0, -freq, orbital_energies,
                                                                  response_manager.params().lo());

    bsh_x.insert(bsh_x.end(), bsh_y.begin(), bsh_y.end());
    return bsh_x;
}

// ============================================================================
// CoupledResponseEquations — one BSH iteration step
// ============================================================================
//
// StaticRestrictedResponse (ω=0, x-channel only):
//   Solves: (A)·x = v   (TDA-like: B-matrix dropped, y≡0)
//
//   Residual RHS per orbital p:
//     θ_p = −2 [ V₀ x_p − ε_{ip} x_i  +  Γ_x[x]  +  v_p^C ]
//   where:
//     V₀ x_p  = (V_local − c_xc·K) x_p     local ground-state potential on x
//     ε_{ip}  = Σ_{i≠p} H_no_diag_{ip} x_i  off-diagonal Fock coupling
//     Γ_x[x]  = response XC kernel (x-only density, tda=true macrotask flag)
//     v_p^C   = Q̂ v^C φ_p               projected perturbation
//
//   BSH update: x_p^new = G(μ_p)·θ_p,   μ_p = √(−2·ε_p)
//   NOTE: uses x-only density (TDA/A-only kernel), not the full A+B static kernel.
vector_real_function_3d CoupledResponseEquations(World& world, const GroundStateData& g_s,
                                                 const StaticRestrictedResponse& vecs,
                                                 const vector_real_function_3d& v_p,
                                                 const std::vector<poperatorT>& bsh_x,
                                                 const ResponseManager& response_manager, ResponseDebugLogger& logger) {

    auto c_xc = g_s.xcf_.hf_exchange_coefficient();
    vector_real_function_3d k_0;
    vector_real_function_3d g_x;

    DEBUG_TIMED_BLOCK(world, &logger, "g0_task",
                      { k_0 = compute_ground_exchange(world, vecs, g_s.orbitals); });
    DEBUG_TIMED_BLOCK(world, &logger, "gx_task",
                      { g_x = compute_gamma_response(world, vecs, g_s.orbitals, g_s.Qhat); });

    auto v_local = g_s.V_local * vecs.x_alpha;
    auto v0x = v_local - c_xc * k_0;
    auto epsilonx = transform(world, vecs.x_alpha, g_s.Hamiltonian_no_diag, true);

    auto thetax = -2.0 * (v0x - epsilonx + g_x + v_p);
    truncate(world, thetax);
    auto rsh = apply(world, bsh_x, thetax);
    rsh = g_s.Qhat(rsh); // project out the ground state density from the response
    return rsh;
}

// TDARestrictedResponse (TDA excited-state, x-channel only):
//   Identical kernel to StaticRestrictedResponse — same x-only density, same
//   macrotask tda=true flag.  Kept as a separate overload for semantic clarity:
//   the caller's type makes it explicit that this is an excited-state TDA solve,
//   not a frequency-domain ω=0 solve.
//
//   BSH exponent per orbital p: μ_p = √(−2·(ε_p + ω))  (ω = excitation energy estimate)
vector_real_function_3d CoupledResponseEquations(World& world, const GroundStateData& g_s,
                                                 const TDARestrictedResponse& vecs,
                                                 const vector_real_function_3d& v_p,
                                                 const std::vector<poperatorT>& bsh_x,
                                                 const ResponseManager& response_manager, ResponseDebugLogger& logger) {

    auto c_xc = g_s.xcf_.hf_exchange_coefficient();
    vector_real_function_3d k_0;
    vector_real_function_3d g_x;

    DEBUG_TIMED_BLOCK(world, &logger, "g0_task",
                      { k_0 = compute_ground_exchange(world, vecs, g_s.orbitals); });
    DEBUG_TIMED_BLOCK(world, &logger, "gx_task",
                      { g_x = compute_gamma_response(world, vecs, g_s.orbitals, g_s.Qhat); });

    auto v_local = g_s.V_local * vecs.x_alpha;
    auto v0x = v_local - c_xc * k_0;
    auto epsilonx = transform(world, vecs.x_alpha, g_s.Hamiltonian_no_diag, true);

    auto thetax = -2.0 * (v0x - epsilonx + g_x + v_p);
    truncate(world, thetax);
    auto rsh = apply(world, bsh_x, thetax);
    rsh = g_s.Qhat(rsh);
    return rsh;
}

// DynamicRestrictedResponse (full TDDFT, ω≠0, x and y channels):
//   Solves the Casida coupled equations:
//     [ A   B  ][X]       [V ]
//     [ B*  A* ][Y]  = ω  [V*]
//
//   Residual RHS in flat layout [x|y] (2N slots):
//     θ_p^x = −2 [ V₀ x_p − ε_{ip} x_i  +  Γ_xy[x,y]_x  +  v_p^C   ]
//     θ_p^y = −2 [ V₀ y_p − ε_{ip} y_i  +  Γ_xy[x,y]_y  +  v_p^{C†} ]
//   where Γ_xy uses the full response density: ρ¹ = Σ(x_i+y_i)φ_i  (tda=false)
//
//   BSH updates:
//     x_p^new = G(μ^x)·θ^x,   μ^x = √(−2·(ε_p+ω))
//     y_p^new = G(μ^y)·θ^y,   μ^y = √(−2·(ε_p−ω))
vector_real_function_3d CoupledResponseEquations(World& world, const GroundStateData& g_s,
                                                 const DynamicRestrictedResponse& vecs,
                                                 const vector_real_function_3d& v_p,
                                                 const std::vector<poperatorT>& bsh_x,
                                                 const ResponseManager& response_manager, ResponseDebugLogger& logger) {
    const auto& xvec = vecs.x_alpha;
    const auto& yvec = vecs.y_alpha;
    const auto& all_x = vecs.flat;

    vector_real_function_3d k_0;
    vector_real_function_3d g_x;

    DEBUG_TIMED_BLOCK(world, &logger, "g0_task",
                      { k_0 = compute_ground_exchange(world, vecs, g_s.orbitals); });
    DEBUG_TIMED_BLOCK(world, &logger, "gx_task",
                      { g_x = compute_gamma_response(world, vecs, g_s.orbitals, g_s.Qhat); });

    auto c_xc = g_s.xcf_.hf_exchange_coefficient();
    auto v_local = g_s.V_local * all_x;
    auto v0x = v_local - c_xc * k_0;
    auto epsilonx = transform(world, xvec, g_s.Hamiltonian_no_diag, true);
    auto epsilony = transform(world, yvec, g_s.Hamiltonian_no_diag, true);
    epsilonx.insert(epsilonx.end(), epsilony.begin(), epsilony.end());
    auto thetax = -2.0 * (v0x - epsilonx + g_x + v_p);
    truncate(world, thetax);
    auto rsh = apply(world, bsh_x, thetax);
    truncate(world, thetax);
    rsh = g_s.Qhat(rsh);
    return rsh;
}

// vector_real_function_3d
// ResponseSolverPolicy<StaticRestrictedResponse>::CoupledResponseEquations(
//     World &world, const GroundStateData &gs, const StaticRestrictedResponse
//     &vecs, const vector_real_function_3d &vp, const std::vector<poperatorT>
//     &bsh_x, const ResponseManager &rm, ResponseDebugLogger &logger) {
//
//   auto &x = vecs.x_alpha;
//   auto &all_x = vecs.flat;
//
//   auto num_orbitals = gs.orbitals.size();
//
//   std::vector<int> state_index;
//   std::vector<int> ii;
//
//   int i = 0;
//   int si = 0;
//   for (int j = 0; j < num_orbitals; j++) {
//     state_index.push_back(si);
//     ii.push_back(i++);
//   }
//
//   ResponseComputeGroundExchange t0;
//   MacroTask g0_task(world, t0);
//   ResponseComputeGammaX tresponse;
//   MacroTask gx_task(world, tresponse);
//
//   auto c_xc = gs.xcf_.hf_exchange_coefficient();
//   vector_real_function_3d k0;
//   vector_real_function_3d gx;
//
//   DEBUG_TIMED_BLOCK(world, &logger, "g0_task", { k0 = g0_task(ii,
//   state_index, all_x, gs.orbitals, true); }); DEBUG_TIMED_BLOCK(world,
//   &logger, "gx_task", { gx = gx_task(ii, state_index, all_x, gs.orbitals,
//   true); });
//
//   auto v_local = gs.V_local * x;
//   auto v0x = v_local - c_xc * k0;
//   auto epsilonx = transform(world, x, gs.Hamiltonian_no_diag, true);
//
//   auto thetax = -2.0 * (v0x - epsilonx + gx + vp);
//   truncate(world, thetax);
//   auto rsh = apply(world, bsh_x, thetax);
//   rsh = gs.Qhat(rsh); // project out the ground state density from the
//   response return rsh;
// }
//
// vector_real_function_3d
// ResponseSolverPolicy<DynamicRestrictedResponse>::CoupledResponseEquations(
//     World &world, const GroundStateData &gs, const DynamicRestrictedResponse
//     &vecs, const vector_real_function_3d &vp, const std::vector<poperatorT>
//     &bsh_x, const ResponseManager &rm, ResponseDebugLogger &logger) {
//
//   auto &x = vecs.x_alpha;
//   auto &y = vecs.y_alpha;
//   auto &all_x = vecs.flat;
//
//   auto num_orbitals = gs.orbitals.size();
//
//   std::vector<int> state_index;
//   std::vector<int> ii;
//
//   int i = 0;
//   for (int j = 0; j < 2 * num_orbitals; j++) {
//     state_index.push_back(0);
//     ii.push_back(i++);
//   }
//
//   ResponseComputeGroundExchange t0;
//   MacroTask g0_task(world, t0);
//   ResponseComputeGammaX tresponse;
//   MacroTask gx_task(world, tresponse);
//
//   vector_real_function_3d k0;
//   vector_real_function_3d gx;
//
//   DEBUG_TIMED_BLOCK(world, &logger, "g0_task", { k0 = g0_task(ii,
//   state_index, all_x, gs.orbitals, false); }); DEBUG_TIMED_BLOCK(world,
//   &logger, "gx_task", { gx = gx_task(ii, state_index, all_x, gs.orbitals,
//   false); });
//
//   auto c_xc = gs.xcf_.hf_exchange_coefficient();
//   auto v_local = gs.V_local * all_x;
//   auto v0x = v_local - c_xc * k0;
//   auto epsilonx = transform(world, x, gs.Hamiltonian_no_diag, true);
//   auto epsilony = transform(world, y, gs.Hamiltonian_no_diag, true);
//   epsilonx.insert(epsilonx.end(), epsilony.begin(), epsilony.end());
//   auto thetax = -2.0 * (v0x - epsilonx + gx + vp);
//   truncate(world, thetax);
//   auto rsh = apply(world, bsh_x, thetax);
//   truncate(world, thetax);
//   rsh = gs.Qhat(rsh);
//   return rsh;
// }
