#ifndef MOLRESPONSE_V3_KERNELS_EXCHANGE_CTX_HPP
#define MOLRESPONSE_V3_KERNELS_EXCHANGE_CTX_HPP

// =========================================================================
// Tensor-layer exchange for linear response (doc 28; framework doc 27).
//
// Every exchange is K[a,b](c)_k = Σ_i a_i · g(b,c)_{ik},  g(b,c)_{ik} = Poisson(b_i·c_k).
// We build the shared convolution tensors ONCE, then assemble θ as cheap
// contractions — instead of the op-level compute_V0x/compute_gamma which each
// rebuild their own exchange (so Tx was built twice; no home for the g0 cache).
//
//   g0[i*n+k] = Poisson(φ_i·φ_k)   (Class 1: φ-only, cache-eligible per protocol)
//   Tx[i*n+k] = Poisson(φ_i·x_k)   (Class 2: per response, transient)
//
// Static ClosedShell θ = V0x − E0x + γ with:
//   V0x   = V_local·x − c_xc·groundK(x),  groundK(x)_k = Σ_i φ_i·Tx[i,k]   (= K[φ,φ](x))
//   γ     = Q( J·φ − c_xc·( direct + cross ) )                            (Q on γ ONLY)
//             direct_k = Σ_i x_i·g0[i,k]   (= K[φ,x](φ))
//             cross_k  = Σ_i φ_i·Tx[k,i]   (= K[x,φ](φ) — the TRANSPOSE contraction)
//
// The same Tx feeds BOTH groundK (V0x) and cross (γ) — one build, two contractions.
// This header is the FD ClosedShell slice (Static now; Full = +Ty + the Y block).
// Numerics: the explicit Tx-build + dot-contraction is the SAME math as the
// Exchange operator but a different truncation/accumulation path → A/B-to-thresh,
// not bit-identical. Reference (compute_V0x/compute_gamma) stays the gate.
// =========================================================================

#include "tda.hpp"     // ResponseGroundState, poperatorT, Kernels<> primary template
#include "static.hpp"  // Kernels<Static, ClosedShell> specialization (compute_E0x reuse)

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>

#include <vector>

namespace molresponse_v3 {
namespace exch {

using vecfuncT = std::vector<madness::real_function_3d>;

/// Pair-density convolution tensor T[i*nc+k] = Poisson(b_i · c_k), |b|=nb, |c|=nc.
inline vecfuncT
build_pair_tensor(madness::World &world, const poperatorT &coulop,
                  const vecfuncT &b, const vecfuncT &c, double vtol) {
  vecfuncT prod;
  prod.reserve(b.size() * c.size());
  for (const auto &bi : b) {
    auto pi = mul_sparse(world, bi, c, vtol);     // bi · c_k for all k
    prod.insert(prod.end(), pi.begin(), pi.end());
  }
  madness::truncate(world, prod, vtol);
  auto T = madness::apply(world, *coulop, prod);  // Poisson over the whole bundle, ONE wave
  madness::truncate(world, T, vtol);
  return T;   // length nb*nc, row-major [i*nc+k]
}

/// out_k = Σ_i a_i · T[i*stride + k]   (column contraction; |a|=na, |out|=stride).
inline vecfuncT
contract_col(madness::World &world, const vecfuncT &a, const vecfuncT &T,
             std::size_t stride) {
  const std::size_t na = a.size();
  vecfuncT out(stride);
  for (std::size_t k = 0; k < stride; ++k) {
    vecfuncT col(na);
    for (std::size_t i = 0; i < na; ++i) col[i] = T[i * stride + k];
    out[k] = madness::dot(world, a, col);         // Σ_i a_i · col_i
  }
  return out;
}

/// out_k = Σ_i a_i · T[k*n + i]   (row contraction = transpose; square n×n tensor).
inline vecfuncT
contract_row(madness::World &world, const vecfuncT &a, const vecfuncT &T,
             std::size_t n) {
  vecfuncT out(n);
  for (std::size_t k = 0; k < n; ++k) {
    vecfuncT row(T.begin() + k * n, T.begin() + (k + 1) * n);
    out[k] = madness::dot(world, a, row);
  }
  return out;
}

// ---- Class-1 cache (φ-only) ---------------------------------------------
/// g0[i*n+k] = Poisson(φ_i·φ_k). Built once per protocol (φ fixed); Inc 2 caches
/// it on ResponseGroundState. Standalone callers build it here.
inline vecfuncT
build_g0(madness::World &world, const ResponseGroundState &gs, double vtol) {
  return build_pair_tensor(world, gs.coulop, gs.amo, gs.amo, vtol);
}

// ---- per-response Class-2 context ---------------------------------------
struct ResponseExchangeCtx {
  vecfuncT                  Tx;   // Poisson(φ_i·x_k), n*n
  vecfuncT                  Ty;   // Full only
  madness::real_function_3d J;    // Poisson(ρ)  (Coulomb)
};

/// build_ctx for Static ClosedShell: J + Tx (no Ty).
inline ResponseExchangeCtx
build_ctx_static_cs(madness::World &world, const ResponseGroundState &gs,
                    const ResponseStateX<ClosedShell> &state,
                    const madness::real_function_3d &rho, double vtol) {
  ResponseExchangeCtx ctx;
  ctx.J  = madness::apply(*gs.coulop, rho);
  ctx.Tx = build_pair_tensor(world, gs.coulop, gs.amo, state.x_alpha, vtol);
  return ctx;
}

/// θ = V0x − E0x + γ  for Static ClosedShell, assembled from {g0, ctx}.
/// Mirrors fd_solver's θ EXACTLY: Q applies to the γ block ONLY.
inline ResponseStateX<ClosedShell>
assemble_theta_static_cs(madness::World &world, const ResponseGroundState &gs,
                         const ResponseStateX<ClosedShell> &state,
                         const vecfuncT &g0, const ResponseExchangeCtx &ctx) {
  const double thr  = madness::FunctionDefaults<3>::get_thresh();
  const double vtol = thr * 0.1;
  const std::size_t n = gs.amo.size();

  // --- V0x (NO Q):  V_local·x − c_xc·groundK ---
  auto theta = mul_sparse(world, gs.V_local_alpha, state.x_alpha, vtol);
  if (gs.c_xc > 0.0) {
    auto groundK = contract_col(world, gs.amo, ctx.Tx, n);   // Σ_i φ_i Tx[i,k]
    madness::gaxpy(world, 1.0, theta, -gs.c_xc, groundK);
  }
  // --- − E0x (NO Q): off-diagonal Fock transform (reuse the kernel) ---
  {
    auto E0x = Kernels<Static, ClosedShell>::compute_E0x(world, gs, state);
    madness::gaxpy(world, 1.0, theta, -1.0, E0x.x_alpha);
  }
  // --- + γ (WITH Q):  Q( J·φ − c_xc·(direct + cross) ) ---
  auto g = mul(world, ctx.J, gs.amo, true);                  // J·φ
  if (gs.c_xc > 0.0) {
    auto direct = contract_col(world, state.x_alpha, g0, n); // Σ_i x_i g0[i,k]
    auto cross  = contract_row(world, gs.amo, ctx.Tx, n);    // Σ_i φ_i Tx[k,i]
    madness::gaxpy(world, 1.0, g, -gs.c_xc, direct);
    madness::gaxpy(world, 1.0, g, -gs.c_xc, cross);
  }
  g = gs.Qa(g);
  madness::truncate(world, g, vtol);
  madness::gaxpy(world, 1.0, theta, 1.0, g);                 // θ = V0x − E0x + γ
  madness::truncate(world, theta, thr);
  return ResponseStateX<ClosedShell>{std::move(theta)};
}

} // namespace exch
} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_KERNELS_EXCHANGE_CTX_HPP
