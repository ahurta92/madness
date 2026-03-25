#pragma once
// excited_ops/RestrictedFullBundleOps.hpp
//
// Full TDDFT (Casida) excited-state bundle ops — restricted closed-shell.
//
// All iteration logic for DynamicRestrictedResponse is handled by the generic
// templates in BundleKernels.hpp.  The per-type dispatch happens
// implicitly via if constexpr branches in ResponseKernels.hpp:
//
//   compute_gamma_response<DynamicRestrictedResponse>:
//     ρ¹ = Σᵢ (xᵢ+yᵢ)·φᵢ     (x+y density; full TDDFT)
//     γ_x = Q̂·[2·J[ρ¹]·φ − c_xc·(K[φ₀,x]·φ + K[y,φ₀]·φ)]
//     γ_y = Q̂·[2·J[ρ¹]·φ − c_xc·(K[φ₀,y]·φ + K[x,φ₀]·φ)]
//     flat = [γ_x | γ_y]
//
//   make_excited_bsh_operators<DynamicRestrictedResponse>:
//     2N operators: N at +ω (x channel), N at −ω (y channel)
//     μ_p^x = √(−2·(ε_p + ω)),  μ_p^y = √(−2·(ε_p − ω))

#include "BundleKernels.hpp"

// ---- metric_inner -----------------------------------------------------------
// Symplectic metric inner product matrix for full TDDFT:
//   S_ij = ⟨x_i|x_j⟩ − ⟨y_i|y_j⟩
//
// The minus sign on the y block reflects the Casida equation structure:
//   ⎡ A   B  ⎤ ⎡X⎤       ⎡ 1   0 ⎤ ⎡X⎤
//   ⎣ B*  A* ⎦ ⎣Y⎦  = ω  ⎣−1   0 ⎦ ⎣Y⎦
//
// Used by excited_ops::build_rotation_matrices (BundleKernels.hpp) to form
// the S matrix of the generalized eigenvalue problem A·c = ω·S·c.

[[nodiscard]] inline madness::Tensor<double>
metric_inner(madness::World&                                       world,
             const std::vector<DynamicRestrictedResponse>&         Chi,
             const std::vector<DynamicRestrictedResponse>&         Gamma)
{
    return excited_ops::inner_x(world, Chi, Gamma) - excited_ops::inner_y(world, Chi, Gamma);
}
