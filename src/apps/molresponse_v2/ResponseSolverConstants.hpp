#pragma once
// ResponseSolverConstants.hpp — named physics and numerical constants used
// across the response ops modules.  All values carry documented justification.

namespace ResponseSolverConstants {

// Short-range cutoff for Coulomb and BSH operators (atomic units).
// 1e-10 is the standard MADNESS choice; smaller values increase cost with
// negligible accuracy gain for response-property calculations.
constexpr double k_coulomb_lo = 1.e-10;

// Prefactor in the BSH residual vector:
//   θ_p = k_bsh_residual_prefactor * (V₀·x_p − ε·x + Γ[x] + v_p)
// Comes from the BSH Green's function identity: G(μ) = (−∇² + μ²)⁻¹,
// combined with the response equation written as (−∇² + μ²)·x = θ.
constexpr double k_bsh_residual_prefactor = -2.0;

// Spin degeneracy factor for restricted (closed-shell) calculations.
// Restricted: α-spin = β-spin, so each occupied spatial orbital contributes
// twice to the response density and to the polarizability prefactor.
constexpr double k_restricted_spin_factor = 2.0;

} // namespace ResponseSolverConstants
