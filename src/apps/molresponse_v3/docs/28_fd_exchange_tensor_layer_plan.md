# 28 — Linear-response (FD) exchange tensor layer: implementation plan

Status: PLAN (propose-diff-first; implements docs 26/27 for the FD path). Scope:
**linear response (FD) ClosedShell — Static + Full**. ES (`assemble_lambda`) and
VBC/β reuse the same layer but are NOT wired here. The op-level kernels
(`compute_V0x`/`compute_gamma`) stay as the gated REFERENCE.

## 0. Why a new layer (the problem)
Today's FD step (`fd_solver.hpp:236-266`) builds θ per response as
`compute_V0x − compute_E0x + compute_gamma`. The kernels are decomposed **by
operation**, but the exchange reuse is **by tensor** (doc 27): `T_x` is needed by
*both* the V0x ground exchange and the γ cross term; `g₀` by *both* γ-direct X and
Y. As black-box op calls, `T_x` is built twice and there's no home for the `g₀`
cache. Fix: a tensor layer below the ops — build the shared convolution tensors
ONCE, then *assemble* θ as cheap contractions.

## 1. The layer

```cpp
// kernels/exchange_ctx.hpp  (new)

// g(b,c)_{ik} = Poisson(b_i · c_k);  K[a,b](c)_k = Σ_i a_i · g(b,c)_{ik}  (doc 27 §1)

// --- Class 1: φ-only, CACHED per protocol on the ground state -------------
struct GroundExchangeCache {            // lives on ResponseGroundState
  vecfuncT g0;                          // g0[i*n+k] = Poisson(φ_i·φ_k), n² (symmetric)
  bool built = false;
};
// built in build_response_ground_state_* (like K0); rebuilt on re-prep.

// --- Class 2: response-dependent, per response, TRANSIENT ------------------
struct ResponseExchangeCtx {
  vecfuncT Tx;                          // Tx[i*n+k] = Poisson(φ_i·x_k), n²
  vecfuncT Ty;                          // Full only (Poisson(φ_i·y_k))
  madness::real_function_3d J;          // Poisson(ρ)  (Coulomb / Hartree)
};

template <typename Type, typename Shell>
ResponseExchangeCtx build_ctx(World&, const ResponseGroundState&,
                              const State&, const real_function_3d& rho);
//   builds J = coulop(rho); Tx = Poisson(φ_i·x_k) ONCE (and Ty for Full).

template <typename Type, typename Shell>
State assemble_theta(World&, const ResponseGroundState&, const State&,
                     const GroundExchangeCache&, const ResponseExchangeCtx&);
//   θ = [V_local·x − c_xc·groundK(Tx)] − E0x + Q( J·φ − c_xc·γ(g0,Tx[,Ty]) )
//   (assemble_lambda<…> — adds T0x for the ES subspace — is the ES extension.)
```

Contraction helpers (cheap mul+reduce, NO convolutions):
- `groundK_k = Σ_i φ_i · Tx[i,k]`         (V0x ground exchange, from Tx)
- `gamma_direct_k = Σ_i x_i · g0[i,k]`    (γ direct, from cached g0)
- `gamma_cross_k  = Σ_i φ_i · Tx[k,i]`    (γ cross, from Tx — the TRANSPOSE contraction)

## 2. Term → tensor → contraction map (FD ClosedShell)

| path | term | tensor | contraction | Q? | in |
|---|---|---|---|---|---|
| both | V0x ground-K (X) | **Tx** | Σ_i φ_i Tx[i,k] | no | θ via V0x |
| Full | V0x ground-K (Y) | **Ty** | Σ_i φ_i Ty[i,k] | no | θ via V0x |
| both | γ direct (X) | **g₀** | Σ_i x_i g₀[i,k] | yes | γ |
| Full | γ direct (Y) | **g₀** | Σ_i y_i g₀[i,k] | yes | γ |
| Static | γ cross | **Tx** | Σ_i φ_i Tx[k,i] | yes | γ |
| Full | γ cross (γ_X) | **Ty** | Σ_i φ_i Ty[k,i] | yes | γ |
| Full | γ cross (γ_Y) | **Tx** | Σ_i φ_i Tx[k,i] | yes | γ |

Shared builds: **Tx serves V0x(X) + γ-cross(Static / Full-γ_Y); Ty serves V0x(Y)
+ γ-cross(Full-γ_X); g₀ serves γ-direct X + Y (cached).** Convolution builds per
response-iter: Static `Tx` (was 2×n²); Full `Tx`+`Ty` (was 4×n²). g₀ once/protocol.

**Q discipline (must match the reference exactly):** Q is applied to the γ block
ONLY (`apply_gamma` does Q+truncate); V0x and E0x are NOT Q-projected. `assemble_theta`
must keep V_local·x − c_xc·groundK and −E0x un-projected, and Q only the
`J·φ − c_xc·(direct+cross)` block. (`compute_V0x` no Q; `compute_gamma` Q.)

## 3. Timings (instrument each type — extends `--es-time`)
New `ES_TIMING`/`FD_TIMING` fields (wall, fenced under the timing flag only):
- `coulomb`     — J = Poisson(ρ) + J·φ (Hartree)
- `g0_build`    — g₀ convolutions (report once/protocol; ~0 per-iter when cached)
- `g0_contract` — γ-direct reductions (Class 1)
- `Tx_build` / `Ty_build` — Class-2 response convolutions (the dominant cost)
- `tx_contract` — groundK + γ-cross reductions
- keep `E0x`, `bsh`, `kain`, `total`.
Baseline = reference path (flag 0). Each increment reports the delta on these fields.

## 4. Implementation increments (each: gated, A/B bit-identity, timed)
- **Inc 1 — the layer + FD rewire (delivers the Tx fusion).** Add
  `kernels/exchange_ctx.hpp` (`build_ctx`, `assemble_theta`, contraction helpers,
  parameterized by `<Static|Full, ClosedShell>`). Rewire `fd_solver.hpp::step` to
  `build_ctx` then `assemble_theta` under a gate (`policy_.exchange_tensor` /
  `--fd-tensor`, 0 = current per-op path = REFERENCE). g₀ built per-iter here (not
  cached yet) so Inc 1 is pure restructure + fusion. Add the §3 timers. A/B
  bit-identity (machine-eps) vs the reference on h2o-α. → measures baseline + the
  Tx-fusion delta.
- **Inc 2 — cache g₀ per protocol.** Build `GroundExchangeCache.g0` in
  `build_response_ground_state_*` (like `K0`); `build_ctx`/`assemble` read it.
  Direct-exchange convolutions → once/protocol. A/B + g0_build/g0_contract delta.
- **Inc 3 (later, with the G strategy + VBC) — batch Tx over states; tile;
  parallelize.** Out of scope here.

## 5. Files / wiring
- NEW `kernels/exchange_ctx.hpp` — the layer (templated on Type/Shell).
- `kernels/tda.hpp` — `ResponseGroundState` gains `GroundExchangeCache g0_cache`.
- `solvers/build_response_ground_state.hpp` — build g₀ (Inc 2).
- `solvers/fd_solver.hpp::step` — gated branch: ctx+assemble vs current.
- `solvers/convergence_policy.hpp` — `exchange_tensor` knob; `main.cpp`/skeleton
  `--fd-tensor`; thread through.
- `tests` — extend the FD skeleton / a `verify_fd_tensor.sh` A/B (flag 0 vs 1).
- Reference path (`compute_V0x`/`compute_gamma`) UNTOUCHED.

## 6. A/B verification
Bit-identity (target machine-eps; `apply_exchange` may tile/truncate slightly
differently than the explicit Tx build → A/B to ~1e-10, and confirm converged α
matches the recorded reference within tolerance). Single node, h2o: Static (ω=0) +
one dynamic ω, Full path too. Gate 0 = reference; reuse the `--es-batch` A/B
discipline (`verify_es_batch.sh` analog for FD).

## 7. Risks / notes
- **Q discipline** (§2) — the single easiest way to break bit-identity; mirror it exactly.
- **Truncation order** — build Tx with the same `vtol` as `apply_exchange`; truncate
  contractions to match. Expect ~1e-10 A/B, not exact.
- **Memory** — Tx (n²) [+Ty] transient per response; tile in blocks if n² is heavy
  (mirror MADNESS exchange tiling). g₀ cache resident = n² distributed (Inc 2).
- **Designed for extension** — `build_ctx`/`assemble_*` templated on `<Type,Shell>`;
  `assemble_lambda` (ES, adds T0x) and the VBC source reuse the SAME ctx + g₀ cache.
  Linear-response FD is the first consumer; ES/VBC wire in after the G strategy.
- **Supersession** — this is the concrete first slice of doc 26 (γ batching) +
  doc 27 (the K[a,b](c) catalog); the G>1 scheduler (doc 25) runs on top.
