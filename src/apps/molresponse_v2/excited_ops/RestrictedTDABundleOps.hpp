#pragma once
// excited_ops/RestrictedTDABundleOps.hpp
//
// TDA (Tamm-Dancoff Approximation) excited-state bundle ops.
//
// All iteration logic for TDARestrictedResponse is handled by the generic
// templates in BundleKernels.hpp.  The per-type dispatch happens
// implicitly via if constexpr branches in ResponseKernels.hpp:
//
//   compute_gamma_response<TDARestrictedResponse>:
//     ρ¹ = Σᵢ xᵢ·φᵢ          (x-only; no y channel in TDA)
//     γ  = Q̂·[2·J[ρ¹]·φ − c_xc·K[φ₀,x]·φ]
//
//   make_excited_bsh_operators<TDARestrictedResponse>:
//     N operators at +ω  (x channel only)
//     μ_p = √(−2·(ε_p + ω))
//
// metric_inner for TDARestrictedResponse is defined here (not in
// ops/TDARestrictedOps.hpp) so it can delegate to excited_ops::inner_x.
// It is found by ADL when build_rotation_matrices<TDARestrictedResponse>
// is instantiated in ExcitedResponse.cpp.

#include "BundleKernels.hpp"
#include "../ops/TDARestrictedOps.hpp"

// ---- metric_inner -----------------------------------------------------------
// Symplectic metric inner product matrix for TDA: S_ij = ⟨x_i|x_j⟩
//
// TDA has no y-channel (y ≡ 0), so the Casida symplectic metric
//   ⟨Φ_i|Φ_j⟩_metric = ⟨x_i|x_j⟩ − ⟨y_i|y_j⟩
// reduces to a plain x-overlap.

[[nodiscard]] inline madness::Tensor<double>
metric_inner(madness::World&                                    world,
             const std::vector<TDARestrictedResponse>&         Chi,
             const std::vector<TDARestrictedResponse>&         Gamma)
{
    return excited_ops::inner_x(world, Chi, Gamma);
}
