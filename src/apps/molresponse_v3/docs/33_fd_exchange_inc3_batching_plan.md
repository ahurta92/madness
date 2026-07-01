# 33 — FD exchange tensor Inc-3: batch / tile / parallelize + perf instrumentation

Status: PLAN (propose-diff-first). Scope: **linear response (FD) ClosedShell —
Static + Full**, building on doc 28 Inc-1/Inc-2 (the `exchange_ctx` tensor layer +
the per-protocol `g0` cache). Inc-3 is doc 28 §4's deferred increment + the board's
exchange "NEXT". The reference kernels (`compute_V0x`/`compute_gamma`/`compute_E0x`)
stay the gate-0 oracle. The `exchange_ctx` shape is pinned in
`operator_contracts.md` ("Tensor-layer exchange").

## 0. Starting point (what Inc-1/2 already deliver)
- `exchange_ctx.hpp`: `build_ctx_*` builds `{Tx[,Ty], J}` once per response-iter;
  `assemble_theta_*` assembles θ as contractions of `{g0, Tx[,Ty]}` (doc 27 §2 catalog).
- **Tx fusion (Inc-1):** one `Tx` build serves both V0x ground-exch (`col(φ,Tx)`)
  and γ-cross (`row(φ,Tx)`); `Ty` serves the Full Y terms. Static Class-2 convs
  2n²→n², Full 4n²→2n².
- **g0 cache (Inc-2):** `g0 = Poisson(φ_i·φ_k)` built once per protocol on
  `ResponseGroundState::g0_alpha`; γ-direct (`col(x/y, g0)`) is contractions-only
  after iter 1. Correct across the protocol climb because `solve_fd_protocol`
  rebuilds `tgt.gs` (fresh-empty `g0_alpha`) on every (re)prepare.
- A/B floor vs the reference: ~1e-3 RELATIVE on a converged channel
  (`verify_fd_tensor.sh`).

What's left (the cost): the **Class-2 Poisson convolutions** (`Tx`,`Ty` — the n²
internode-traffic builds) are still issued as separate fenced waves, and
`build_pair_tensor` materializes all n² products before the apply (a memory peak
that scales with n_occ — the c6h6/naphthalene OOM regime).

## 1. Where batching actually pays for FD (the M=1 reality)
The production path `solve_fd_protocol` solves **one (pert, freq) state per call**
(`calc_executor.hpp:276`, `tgt.responses.resize(1)`) → `M=1` inside `fd_solver::step`.
So "batch Tx over states" has **no FD win in the per-state scheduler** as wired today.
The real Inc-3 wins, in order of FD payoff:

1. **Within-state wave fusion (Full):** `Tx` and `Ty` are two `build_pair_tensor`
   calls, each fences. Fuse into ONE Poisson wave (`[φ_i·x_k, φ_i·y_k]` products →
   one truncate → one apply → split). 2 waves → 1. (Static is already 1 wave.)
2. **Tiling `build_pair_tensor` (memory):** bound the n² product/Poisson peak by
   processing φ-rows (or (i,k) blocks) in tiles, mirroring MADNESS's exchange
   tiling. Directly serves the project scaling goal (>20 occ; the OOM blocker).
3. **Batch over states (M>1):** the doc-26 `compute_gamma_flat` generalization —
   pays for **ES bundles** (M roots), **multi-channel FD** (AXES=xyz / multi-ω),
   and **VBC/β**, and is where the `parallel-runtime` multi-world (docs 31/32)
   localizes the Class-2 comm. Cross-thread; lands after 3a/3b prove the kernel.

## 2. Increments (each: gated, A/B-to-floor, measured as a delta)

### Inc-3a — fuse Tx+Ty into one Poisson wave (Full)
- New `build_pair_tensor_multi(world, coulop, b, {c1,c2,…}, vtol)` (or a 2-ket
  overload) that concatenates the product bundles, truncates once, applies the
  Poisson over the whole bundle in one wave, then splits back to `Tx`,`Ty`.
- `build_ctx_full_cs` calls it instead of two `build_pair_tensor` calls.
- Static unchanged. A/B: Full h2o-α (ω=0.057). Meter: `Tx_build`+`Ty_build` fence
  count 2→1.

### Inc-3b — tile `build_pair_tensor` (memory peak)
- Add a tile size knob (`policy_.exchange_tile` / env, default = no tiling =
  current behavior). When set, build the products + Poisson in row/(i,k) blocks,
  truncating each block before the next, so peak ≈ `tile·n·|φ|` not `n²·|φ|`.
- Pure memory/throughput restructure — **bit-identical to untiled** (same convs,
  same truncation order within a block; choose the block boundary to not change
  per-product truncation). A/B + a `MEMORY_HWM` delta at a large-n fixture
  (c6h6 at a coarse probe protocol — the 15d "measure low-k" trick).

### Inc-3c — batch Tx/cross over states (M>1) [cross-thread; after 3a/3b]
- Generalize `assemble_theta_*` to a flat M-state form (`compute_gamma_flat`,
  doc 26 §3): one Coulomb wave for M `J[ρ_s]`, one Poisson wave for the M·n²
  cross products, state-major layout.
- **No-mixing contract (doc 26 §4):** `s` is NEVER a reduction index. Direct uses
  φ-only g0; cross products carry `(s,i,k)` and the `Σ_i` is bounded to state `s`'s
  slice. FORBIDDEN: handing flattened `[x_0..x_{M-1}]` as one Exchange ket.
- Consumers: ES bundle (`es_solver` step), multi-channel FD, VBC source. Coordinate
  the batched-wave shape with `parallel-runtime` (the interface it distributes).

## 3. Perf instrumentation (reconciles doc 28 §3 with perf-model doc 29)
doc 28 §3 proposed a v3-local `FD_TIMING` line. **perf-model (doc 29) rejected a
v3-local wall-timer** in favor of core `WorldProfile` (compile flag
`WORLD_PROFILE_ENABLE`, runtime env `MADQC_PROFILE_JSON`, zero-effect when off).
Inc-3 adopts that: exchange's contribution is **named `PROFILE_BLOCK` meters** around
its phases, so their CPU/comm/**call-count** attribution flows into perf-model's
profile JSON. The meter names are exchange's slice of perf-model's pinned schema:

| meter | wraps | what the cost model reads |
|---|---|---|
| `coulomb` | `J = Poisson(ρ)` + `J·φ` | per-iter Coulomb cost |
| `g0_build` | `build_g0` convs | once/protocol; ≈0 per-iter when cached |
| `g0_contract` | γ-direct `col(x/y,g0)` reductions | Class-1 contraction cost |
| `Tx_build` / `Ty_build` | the Class-2 Poisson wave(s) | **Tx/wave count** + bytes = the n² conv coefficient |
| `tx_contract` | `col(φ,Tx)` + `row(φ,Tx/Ty)` reductions | Class-2 contraction cost |

- **Tx/tile counts** = the call counts of the `Tx_build`/`Ty_build` blocks (and a
  small explicit tile counter for 3b) — the board's "report Tx/tile counts" obligation.
- Zero-effect when off (the `PROFILE_*` macros compile to nothing without the flag);
  no numerics touched.
- Wall stays the coarse `StateMetrics.wall_s` / `PROTOCOL_START/DONE` layer
  (perf-model PM-3 joins the fine CPU/comm profile to it).

## 4. A/B verification (per increment)
- `verify_fd_tensor.sh` (gate 0 = reference, gate 1 = `--fd-tensor`), convergence-
  gated, RELATIVE `TOL=2e-3`. Static (ω=0) + Full (ω=0.057), single node, h2o.
- 3a/3c: A/B-to-floor (different accumulation order). 3b: bit-identical to untiled
  (assert exact, not just within floor).
- Multi-protocol climb (`PROTOCOL=1e-4,1e-6`) at least once per increment to keep
  the g0-cache-rebuild path covered.
- Do not regress H2O/CH3OH/C2H4 α; aim 3b at bringing c6h6 under the per-task budget.

## 5. Files / wiring
- `kernels/exchange_ctx.hpp` — `build_pair_tensor_multi` (3a), tiled
  `build_pair_tensor` (3b), `PROFILE_BLOCK` meters (3d), the flat M-state form (3c).
- `solvers/convergence_policy.hpp` — `exchange_tile` knob (3b); thread `--fd-tensor-tile`.
- `solvers/fd_solver.hpp::step` — unchanged control flow; the gate-1 branch calls
  the same `assemble_theta_tensor` entry points.
- Reference path (`compute_V0x`/`compute_gamma`/`compute_E0x`) UNTOUCHED.
- perf-model owns the `WorldProfile` JSON emitter + schema; exchange only adds the
  named blocks.

## 6. Risks / notes
- **Q discipline + truncation order** (doc 28 §2/§7) — re-verify after every reshape.
- **3b boundary:** tiling must not change per-product truncation, or it stops being
  bit-identical to untiled (it becomes A/B-to-floor like 3a). Pick the block edge
  on the φ-row (i) index so each `Poisson(φ_i·c_k)` is truncated identically.
- **3c memory:** M·n² cross products transient — tile by state (reuses 3b's tiler).
- **Cross-thread:** 3c's batched-wave shape is the `exchange ↔ parallel-runtime`
  contract; 3d's meter names are the `exchange ↔ perf-model` contract — both pinned
  in `operator_contracts.md`.

## 7. Sequencing
3a (Full wave fusion) → 3b (tiling; the scaling payoff) → 3d (meters; can land with
3a so each step is measured) → 3c (state batching; cross-thread, after the FD kernel
is proven). Land each behind the existing `--fd-tensor` gate with an A/B before the next.
