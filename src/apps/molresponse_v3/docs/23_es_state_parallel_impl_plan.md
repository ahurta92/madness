# 23 — ES state-parallel: implementation plan (post-feasibility)

Status: DESIGN (propose-diff-first; no solver code yet). Closed-shell TDA.
Supersedes the *staging* of doc 22 §1–7; doc 22 §8 (persistent design) stands.

## 0. Decision input — feasibility verdict GREEN (2026-06-14)

The larger-system feasibility study (scope B; `notes/2026-06-14_es_node_phi_feasibility.md`,
`node_phi_feasibility.json`) ran 17 sbatch jobs at 2–4 nodes × 8 ranks:

- node-aligned-subworld φ gives `per_rank_node = if_replicated / R` **exactly** (R =
  ranks/node), at every N (5/8/21/34) and both k=6 and the real k=10/thresh=1e-7,
  ratio = 1.0, node-count-independent. → an exact **8× per-rank reduction** of the
  ground-state-φ term at the production 8-ranks/node geometry.
- Projected real per-rank φ (v2 measured MaxRSS ÷ R): c6h6 ~67/8 ≈ 8.4 GB, naphthalene
  < 8.4 GB → both fit the 72 GB budget where v2 OOMs.
- Distributed subspace A/S allreduce stayed machine-eps (8.9e-16) to M=21 with G=4.

→ The persistent-subworld restructure removes the OOM blocker and is worth building.

## 1. The current iteration (es_solver.hpp) — what we are parallelizing

`ESSolver<Type,Shell>::step` dispatches (es_solver.hpp:310) to `step_rotate_pieces`
(default) or `step_recompute_pieces` (`policy_.stream_theta`). One outer iteration:

| Phase | es_solver.hpp | per-root? | couples all M roots? |
|------|----------------|-----------|----------------------|
| 0 project + orthonormalize | 367–373 | yes | no |
| 1 build: density, γ(exchange), V0x, T0x, E0x, E0x_full, Λ assembly | 375–432 | **yes (~78%)** | no |
| 2 subspace `A=⟨X_i\|Λ_j⟩`, `S=⟨X_i\|X_j⟩`, diagonalize→ω,U | 435–461 | no | **yes (<1%)** |
| 3 rotate X,V0x,E0x,γ by U; θ=V0x−E0x+γ | 463–474 | (rotate couples) | **yes** |
| 4 BSH apply + residual | 476–527 | yes (or batched) | no |
| 5 KAIN + explosion guard | 529–544 | yes (per-state) | no |

Timing (h2o, 3 roots, 4 ranks, `--es-time`): build 42% + γ 36% = **78%** per-root;
subspace coupling **0.5%**. The per-root build is the target; the M×M coupling is cheap.
`--es-batch` (tda_batch.hpp) already bundles V0x/T0x/BSH; γ/E0x stay per-root.

## 2. Proven building blocks — reuse, do not rebuild

- `make_node_aligned_subworld` (node_subworlds.hpp:52) — one subworld/node. **Inc 1.**
- φ shipped into a node-subworld = one distributed copy/node (`test_node_phi`). **Inc 2.**
- distributed A/S = per-node partial columns + universe allreduce, functions never cross
  (`test_subspace_allreduce`, == direct to 8.9e-16). **keystone.**
- `step_recompute_pieces` (es_solver.hpp:575) — build Λ, drop, rotate X, **recompute θ**
  on rotated X. This *is* the recompute mechanism the persistent design needs.
- `rs::inner / rs::metric / rs::diagonalize / rs::transform` (response_space_ops.hpp).

## 3. Recommendation — skip heavyweight Inc 3a; go S1 → S2

**Drop doc 22 §1–7 (MacroTask build fan-out) as a perf step.** It gathers
`[Λ|V0x|E0x|γ] = 4·M·n_occ` functions to the universe *every iteration* — strictly
MORE cross-world traffic than the persistent design ever moves, so its cost number is an
unrepresentative upper bound. Its only durable value is **bit-identity of the build
kernels in subworld context** (esp. direct exchange with `K0=nullptr`). Capture that
value far more cheaply with a focused proof (S1), then build the real design (S2).

### S1 — build-in-subworld bit-identity proof  (standalone test, NO solver change) — DONE ✓

**RESULT (2026-06-14): PASS.** `tests/test_es_build_subworld.cpp`, 2 nodes × 8:
n_occ=5/M=3 → `max|A_sub−A_uni| = 2.5e-12`; n_occ=21(c6h6)/M=6 → `3.9e-12`; both
rc=0. cached-vs-direct exchange diagnostic 1.2–1.6e-12. → the GS-dependent build
kernels (V0x direct-exchange, γ, E0x focka-transform, T0x → Λ) are bit-identical in a
node-subworld at machine-eps, and S2's subworld build will match the single-World
(cached-K0) reference. The one practical gotcha (now a hard contract above): the
subworld build needs `FunctionDefaults<3>::set_default_pmap(subworld)` around it —
exactly what MacroTaskQ does (macrotaskq.h:710/853).



`tests/test_es_build_subworld.cpp` (mirrors `test_subspace_allreduce`): ship the full GS
(φ/amo, V_local_alpha, focka, focka_no_diag, aeps, c_xc, lo) into a node-subworld (Inc 2
machinery), build `Λ_j` (and V0x/E0x/γ) for K roots **locally** with direct exchange,
compare to the universe build → expect machine-eps. This is the **one open question the
φ-count and A/S-scalar proofs did not cover**: do `compute_V0x/compute_gamma/compute_E0x`
(GS-dependent kernels, direct exchange) reproduce bit-for-bit in a subworld? Cheap,
gated, throwaway-safe. Register in CMakeLists + `cm_build`; add `verify_es_build_subworld.sh`.

### S2 — persistent-subworld ES iteration  (solver change; propose-diff before coding)

Gated mode in `es_solver.hpp` (`policy_.es_subworlds` + `--es-subworlds=G`, default 0 =
exact current path = reference). Per protocol: build the node-subworld pool once, ship GS
once (Inc 2), **replicate X across subworlds**, partition the M roots round-robin by node.
Per iteration (doc 22 §8.1), with the 78% build staying local:

```
local  (subworld g, owned roots j):  build Λ_j        (reuse step_recompute_pieces ctx)
local  (subworld g):                 A[:,owned] = ⟨X_i|Λ_j⟩ ;  S[:,owned] = ⟨X_i|X_j⟩
universe:                            allreduce M×M A,S  (keystone — Λ NEVER crosses)
local  (replicated A,S ⇒ identical U): diagonalize → ω,U
local  (subworld g):                 rotate X by U ; recompute θ ; BSH ; KAIN  (owned)
sync:                                re-broadcast updated owned roots to all replicas
```

Per-iter cross-world traffic = M×M scalars + X re-broadcast (M·n_occ functions). A/B
verify vs the single-World reference on **converged roots only** (`verify_es_batch.sh`
pattern, `--es-subworlds=G` vs `0`, 1 and 2 nodes).

## 4. Risks / open points (resolve in S2 design)

- **KAIN history** is per-state. Rotation mixes all roots, but X is replicated and U is
  identical everywhere, so every node holds all rotated roots; the **owner of slot j**
  applies KAIN to slot j with its local history. History stays node-local; no migration.
- **Diagonalization determinism** — A,S are bit-identical on all ranks after the
  allreduce, so `sygvp` + the identity-preserving fixups must be rank-deterministic for U
  (hence the rotation) to agree across subworlds. Confirm no rank-dependent tie-breaking.
- **Re-broadcast** = gather owned roots → scatter to subworld replicas via
  `copy(target_world, f)` (as in the φ/keystone tests). Needs the teardown/fence
  discipline (clear all Functions before `finalize()`; subworld before universe).
- **Direct exchange** (`K0=nullptr`) in the subworld build — S1 proves bit-identity.
  (S1 single-rank: cached-vs-direct exchange agree to 1.2e-12, so S2's subworld
  direct-exchange build matches the single-World cached-K0 reference at machine-eps,
  not truncation-level — A/B will be effectively bit-identical.)
- **Subworld pmap (HARD CONTRACT, found in S1):** Functions used in a subworld must
  ARRIVE there (via `copy(subworld, f)` / Cloud), NOT be built fresh with
  `real_factory_3d(subworld)` — a fresh build inherits the GLOBAL (universe) pmap, so
  operator Isends target ranks outside the subworld → `MPI_ERR_RANK: invalid rank`
  (aborts a fence). The real S2 design ships GS via Cloud (handled), so it is safe;
  but any helper that builds intermediates in the subworld must set/inherit the
  subworld pmap. S1's first 2-node run hit exactly this; fixed by shipping the GS in.
- **Comparability** — only converged roots are bit-stable; unconverged are path-chaotic.

## 5. Next action

S1 DONE (PASS, §3). Next: bring the **S2** `es_solver.hpp` diff for explicit approval
before coding it. The proven building blocks now cover every primitive S2 needs —
subworld creation (Inc 1), φ one-copy-per-node (Inc 2), A/S allreduce (keystone),
**subworld build bit-identity + the set_default_pmap discipline (S1)**. S2 wiring:
persistent node-subworld pool + GS shipped once (Cloud or copy) + X replicated +
roots partitioned; per-iter build-local → M×M allreduce → diagonalize → rotate/θ/
BSH/KAIN local → re-broadcast X; gated `--es-subworlds=G` (0 = exact current path).
```
build / run via cm.sh:  cm_build → verify_es_build_subworld.sh (sbatch, 2 nodes × 8)  [S1, green]
```
