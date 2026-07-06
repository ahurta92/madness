# 24 — S2: persistent node-subworld ES iteration (design + diff proposal)

Status: PROPOSAL (propose-diff-first — no solver code cut yet). Closed-shell TDA.
Builds on doc 23 (S1 green) + doc 22 §8.

## 0. Premise — every primitive is proven

- **Inc 1** node-aligned subworld (`make_node_aligned_subworld`).
- **Inc 2** φ = one distributed copy/node.
- **keystone** distributed A/S = per-node partial columns + universe allreduce (== direct, 8.9e-16).
- **S1** the GS-dependent build (V0x direct-exch, γ, E0x, T0x → Λ) is bit-identical in a
  subworld (2.5e-12 @ n_occ=5, 3.9e-12 @ c6h6), under the
  `FunctionDefaults<3>::set_default_pmap(subworld)` discipline.

S2 composes them into one distributed ES iteration.

## 1. Architecture decision — a NEW driver, not an `es_solver.hpp` step edit

The per-root build is local, but the **subspace matrix `A_ij=⟨X_i|Λ_j⟩` and the dense
rotation `X_i←Σ_j U_ji X_j` couple ALL roots** — they span subworlds. `ESSolver::step*`
is single-world (every kernel call is `K::foo(world_, …)`). So S2 is **not** a `step()`
tweak; it is a new iteration driver that orchestrates G subworld-local builds + a
universe allreduce. `es_solver.hpp`'s `step_rotate_pieces` / `step_recompute_pieces`
stay **untouched** as the `G=0` reference. The driver REUSES the existing kernels
(`K::compute_*`, `K::bsh_apply`), `rs::inner/metric/diagonalize/transform`,
`ResponseSubspaceKain`, and `node_subworlds.hpp` — it only changes *where* each runs
and adds the M×M allreduce + X re-broadcast.

`step_recompute_pieces` (es_solver.hpp:575) is the phase template: build streamed Λ →
subspace → rotate X → recompute θ on rotated X → BSH → KAIN. The driver mirrors these
phases with the subspace/rotation made distributed.

## 2. The distributed iteration (one `step`)

Persistent per protocol: node-subworld pool (1/node), GS resident in each subworld
(built there under the subworld pmap, S1-style), **X replicated** in every subworld,
roots partitioned round-robin by node index `g` → owned set `S_g` (keystone's
`my_node_index`). `FunctionDefaults<3>::set_default_pmap(subworld)` is in force during
the local phases; restored to universe only for the scalar allreduce + teardown.

```
local  (subworld g, owned j∈S_g):   build Λ_j   (K::compute_T0x/V0x/E0x_full/gamma → assemble_lambda)
local  (subworld g):  Apart[:,S_g] = rs::inner(X_all, Λ_owned)   (X replicated ⇒ all-i local)
                      Spart[:,S_g] = rs::metric(X_all, X_owned)
                      zero non-owned cols; only subworld-rank-0 contributes (keystone)
univ.:  universe.gop.sum(Apart), sum(Spart)            # M×M SCALARS; Λ never crosses
local  (replicated): rs::diagonalize(A,S) → ω,U        # identical inputs ⇒ identical ω,U
local  (subworld g): rs::transform(X_all,U)            # replicated ⇒ local; all subworlds agree
local  (subworld g, owned j):  recompute θ_j on rotated X_j → K::bsh_apply → residual
local  (subworld g):  KAIN on owned roots (per-owned history)
sync:   re-broadcast each owned, updated root → all subworld replicas (full X for next iter)
```

Per-iter cross-world traffic = the M×M `A,S` allreduce + the X re-broadcast
(M·n_occ functions). The 78% build and Λ never leave the subworld.

## 3. Components & files

1. **Knob.** `ConvergencePolicy::es_subworlds` (int, default 0 = current single-World
   path). Plumb: `main.cpp` `--es-subworlds=G` → `ExecutorContext.es_subworlds` →
   `solve_es_tda_closed_shell`. Mirrors how `--es-batch`/`--es-time` thread through.
2. **New: `solvers/es_subworld_iterate.hpp`** — `iterate_subworld_tda_cs(universe, sub,
   info, gs_sub, state, policy, max_iters, …)`: the distributed analog of
   `solvers::iterate`. Owns root partition, replicated-X bookkeeping, the per-iter flow
   in §2, convergence test (reuse `ESSolver::converged`), and re-broadcast. Closed-shell
   TDA only (the scope).
3. **Wiring (one branch) in `solve_es_tda_closed_shell`** (calc_executor.hpp:436): if
   `ctx.es_subworlds > 0`, build the subworld pool + subworld GS + replicate the
   (warmup-seeded) X, run `iterate_subworld_tda_cs` in place of
   `iterate_protocol(solver, …)`; else the exact current path. Warmup/guess
   (`run_oversampled_tda_warmup`) stays single-World; only the main solve distributes.
   Restart load/save (`try_load_es_bundle`/`save_es_roots`) is universe-side, unchanged.

## 4. Staging (each verifiable, low-risk — the Inc1/2/keystone/S1 discipline)

- **S2a — one-iteration bit-identity (standalone test, NO solver change).** New
  `tests/test_es_iter_subworld.cpp`: from a fixed seed state, run ONE distributed step
  (§2) and ONE single-World `ESSolver::step`, compare resulting `ω` and the rotated/
  BSH'd roots (via `A=⟨X_i|Λ_j⟩` or per-root norms) → expect machine-eps. Composes
  S1 (local Λ) + keystone (A/S allreduce); the gate before any solver edit.
- **S2b — driver + knob + wiring.** Implement `es_subworld_iterate.hpp` + the knob +
  the calc_executor branch. Full protocol solve. A/B `verify_es_batch.sh`-style:
  `--es-subworlds=2` vs `0`, **converged roots only**, 1 and 2 nodes → ω match.
- **S2c — robustness.** KAIN-across-rotation correctness; restart parity; (defer
  full-deflation locking — `policy_.lock_converged` stays single-World for now).

## 5. Open questions / risks

- **KAIN history** is per-state; rotation mixes all roots but X is replicated and U is
  identical, so every subworld holds all rotated roots — the owner of slot j runs KAIN
  on j with its node-local history. History never migrates. (Verify in S2b A/B.)
- **Determinism:** A,S are bit-identical on all ranks post-allreduce ⇒ `sygvp` + the
  identity-preserving fixups must be rank-deterministic for ω,U (hence rotation) to
  agree across subworlds. Confirm no rank-dependent tie-break (S2a catches it).
- **Re-broadcast:** gather owned updated roots → scatter to replicas via
  `copy(target_world, f)` across the world boundary (φ/keystone pattern), under the
  set_default_pmap discipline + teardown fences.
- **GS in subworld:** build via `build_response_ground_state_closed_shell(*sub, …)`
  under the subworld pmap (S1 showed this is bit-identical) — no Cloud needed.
- **Scope:** TDA/ClosedShell only; warmup single-World; locking deferred.

## 6. Verification

S2a: `test_es_iter_subworld` (sbatch 2 nodes × 8) → distributed step == single-World
step at machine-eps. S2b: `verify_es_batch.sh --es-subworlds=2` vs `0` on h2o/c2h4,
converged-roots ω parity at 1 and 2 nodes; `--es-time` for the per-iter coupling cost
(M×M allreduce + X re-broadcast) vs the local build — the number that confirms the
restructure pays off at scale.

## 7. Next action

Implement **S2a** (`test_es_iter_subworld`) — standalone, no solver touch — to gate the
driver, then bring the **S2b** diff (`es_subworld_iterate.hpp` + knob + calc_executor
branch) for approval.
