#ifndef MOLRESPONSE_V3_KERNELS_TDA_BATCH_HPP
#define MOLRESPONSE_V3_KERNELS_TDA_BATCH_HPP

// =========================================================================
// Bundle-batched TDA building blocks (closed shell).  Stage-1 of the ES
// root-parallel work (doc 19): turn the M sequential per-root collective
// passes in ESSolver::step_rotate_pieces into ONE pass over the flattened
// M·n_occ vecfunc, for the operations that are pure per-function maps.
//
// Why ONLY these three (V0x, T0x, BSH):
//   They are elementwise over the orbital list — out[i] depends on in[i]
//   alone, with NO cross-function reduction.  Concatenating M roots into a
//   length-(M·n_occ) vector therefore leaves every out[i] bit-for-bit
//   identical to the per-root call; only the fence/task-fill structure
//   changes (M passes → 1).  This is the same reasoning that lets the
//   per-orbital BSH apply at tda.hpp:225-227 batch over n_occ — we extend
//   it across roots.
//
// What is deliberately NOT batched here:
//   - compute_density: `dot` is mul-THEN-sum; batching the mul reorders the
//     internal truncate in the reduction and drifts the result (a batched
//     density variant was tried and rolled back — see common_ops.hpp:48).
//   - compute_gamma exchange: K[φ, x_s](φ) is a DISTINCT operator per root
//     (ket = x_s), so it cannot share one operator across the bundle — this
//     is the embarrassingly-parallel-over-roots piece Stage-2 targets.
//   - E0x / E0x_full: the focka transform mixes orbitals WITHIN a root;
//     batching needs a dense (M·n)×(M·n) block-diagonal matrix = MORE work.
//
// Bodies mirror Kernels<TDA, ClosedShell>::{compute_V0x, compute_T0x,
// bsh_apply} exactly; the --es-batch A/B run in test_v3_es_skeleton pins
// them numerically equal to the per-root reference path.
// =========================================================================

#include "common_ops.hpp"   // apply_kinetic, make_bsh_operators,
                            // apply_ground_exchange, bsh_shift, poperatorT
#include "tda.hpp"          // ResponseGroundState
#include "../solvers/response_state.hpp"  // ResponseStateX<ClosedShell>

#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>

#include <cstddef>
#include <vector>

namespace molresponse_v3::tda_batch {

using State    = ResponseStateX<ClosedShell>;
using vecfuncT = std::vector<madness::real_function_3d>;

/// Concatenate every root's x_alpha into one flat vecfunc of length M·n_occ
/// (root-major: [r0_orb0 … r0_orbN-1, r1_orb0 …]).
inline vecfuncT flatten(const std::vector<State> &roots) {
  vecfuncT f;
  for (const auto &r : roots)
    f.insert(f.end(), r.x_alpha.begin(), r.x_alpha.end());
  return f;
}

/// Scatter a flat [M·n_occ] vecfunc back into out[s].x_alpha (n_occ each).
/// `out` already sized to M; only x_alpha is overwritten.
inline void unflatten_into(const vecfuncT &f, std::size_t n,
                           std::vector<State> &out) {
  for (std::size_t s = 0; s < out.size(); ++s)
    out[s].x_alpha.assign(f.begin() + s * n, f.begin() + (s + 1) * n);
}

/// V0·x over the whole bundle.  Mirror of compute_V0x:
///   V_local·x − c_xc·K0(x),  K0 = K[φ,φ] (one fixed operator → batches).
inline vecfuncT compute_V0x_flat(madness::World &world,
                                 const ResponseGroundState &g0,
                                 const vecfuncT &Xf) {
  const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
  auto Vx = mul_sparse(world, g0.V_local_alpha, Xf, vtol);
  if (g0.c_xc > 0.0) {
    auto k0x = common_ops::apply_ground_exchange(world, g0.K0_alpha, g0.amo,
                                                 Xf, g0.lo);
    gaxpy(world, 1.0, Vx, -g0.c_xc, k0x);
  }
  return Vx;
}

/// T·x over the whole bundle.  Mirror of compute_T0x.
inline vecfuncT compute_T0x_flat(madness::World &world, const vecfuncT &Xf) {
  return common_ops::apply_kinetic(world, Xf);
}

/// BSH apply over the whole bundle.  Mirror of bsh_apply, but each root has
/// its own ω ⇒ its own shift and its own n_occ-long operator set; we build
/// the per-root rhs slice and op set, concatenate to M·n_occ, and run ONE
/// batched apply + Q + truncate.
///   x_new = Q( BSH(ω_s) · (−2·(θ_s + shift_s·x_s)) )   per root s.
inline vecfuncT bsh_apply_flat(madness::World &world,
                               const ResponseGroundState &g0,
                               const vecfuncT &Xf,        // pre-BSH, flat
                               const vecfuncT &ThetaF,    // flat
                               const madness::Tensor<double> &omega,
                               std::size_t n) {
  const std::size_t M = static_cast<std::size_t>(omega.size());
  auto rhs = madness::copy(world, ThetaF);
  std::vector<poperatorT> ops;
  ops.reserve(M * n);
  for (std::size_t s = 0; s < M; ++s) {
    const double shift = common_ops::bsh_shift(g0.aeps, omega(static_cast<long>(s)));
    for (std::size_t p = 0; p < n; ++p)
      rhs[s * n + p].gaxpy(1.0, Xf[s * n + p], shift, false);  // rhs += shift·x
    auto bsh_s = common_ops::make_bsh_operators(world, g0.aeps,
                                                omega(static_cast<long>(s)), g0.lo);
    ops.insert(ops.end(), bsh_s.begin(), bsh_s.end());
  }
  world.gop.fence();              // close the deferred per-slice gaxpys
  scale(world, rhs, -2.0);
  truncate(world, rhs);
  auto nx = apply(world, ops, rhs);   // single batched Green's-function apply
  nx = g0.Qa(nx);
  truncate(world, nx);
  return nx;
}

} // namespace molresponse_v3::tda_batch

#endif // MOLRESPONSE_V3_KERNELS_TDA_BATCH_HPP
