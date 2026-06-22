# 32 — FD state-parallel: implementation plan (F1 → F2 → F3)

Status: **PLAN** (propose-diff-first; no solver code yet). Closed-shell FD.
Implements the decision in doc 31. Reuses the doc 19–24 proven primitives.

## 0. Recap of the decision (doc 31)

One FD state per **node-aligned subworld**; FD states are coupling-free, so the
distributed iteration is **pure fan-out → solve → gather → fence** — none of the ES
M×M allreduce / X re-broadcast. Gated `--fd-subworlds=G` (0 = the current
single-World path = the bit-identical reference). Staged F1 (standalone proof, no
solver change) → F2 (gated `run()` + knob + pre-flight abort) → F3 (weak-scale +
auto-selector). Each stage is verifiable and low-risk, the Inc1/keystone/S1 discipline.

## 1. The exact seam (where the fan-out goes)

- **`CalcManager::run()`** (`calc_executor.hpp:903`): per pass — `world.gop.fence()` →
  reload metadata → `expand_converged_es` → `schedule()` → take `waves.front()`.
  The dispatch is **one line**:
  ```cpp
  for (const auto &item : wave) exec.run_protocol(item);   // :963-964
  world.gop.fence();                                        // :965
  ```
  Today every rank solves each item together, one at a time (single-group). **F2
  replaces this loop with a node-subworld partition.**
- **`solve_fd_protocol<Type,Shell>(ExecutorContext &ctx, pert, freq, thresh, action)`**
  (`:254`): uses `ctx.world` / `ctx.gs` / `ctx.calc_dir` / `ctx.policy`; builds the
  source, runs `solvers::iterate_protocol`, and **`save_fd_state(world, …, ctx.calc_dir,
  …)`** (`:397`) writes the converged state + metrics to disk. → **the filesystem is
  the gather**: each subworld writes its states' archives; the universe re-reads them
  next pass. `FdResponseExecutor::run_protocol` (`:763`) dispatches Type×Shell.
- **Determinism contract** (`calc_executor.hpp:856-872`): the schedule must be
  byte-identical on every rank (same `response_metadata.json` in after a fence,
  `std::map` perturbation order out). With node-subworlds this **still holds**: every
  subworld reads the same metadata after the universe fence, computes the same
  schedule, and the **partition is deterministic** (round-robin by node index, from
  `NodeSubworldInfo`) — so each subworld knows its owned items with no cross-talk and
  no mismatched-collective deadlock. All ranks *within* a subworld agree by
  construction (same World).

## 2. Proven primitives + the one hard discipline (reuse, do not rebuild)

- `make_node_aligned_subworld(universe, &info)` (`node_subworlds.hpp:52`, Inc1) — one
  subworld per physical node via `Split_type(SHARED)`; `NodeSubworldInfo` carries
  `n_nodes` + the within-node rank. `verify_one_subworld_per_node` is the invariant check.
- **The `set_default_pmap(subworld)` discipline (HARD CONTRACT, S1).** Functions used
  in a subworld must *arrive* there via `copy(*sub, f)` (or Cloud), NOT be built fresh
  with `real_factory_3d(*sub)` — a fresh build inherits the **universe** pmap → operator
  Isends target ranks outside the subworld → `MPI_ERR_RANK`. Point
  `FunctionDefaults<3>::set_default_pmap(*sub)` around the subworld phase; restore the
  universe pmap **before** `sub.reset()`.
- **Teardown order** (S1 test, lines 217-221): subworld Functions destruct → `fence` →
  `set_default_pmap(universe)` → `sub.reset()` → `universe.gop.fence()` → universe
  Functions clear → `finalize()`.

`tests/test_es_build_subworld.cpp` is the line-by-line template for §4.

## 3. GS provisioning into the subworld (the key wiring decision)

Each subworld needs the ground state in its own World. Two options:

- **(a) Ship-in via `copy(*sub, f)` from a universe-loaded GS** — S1-proven, safe,
  pmap-correct. **RECOMMENDED for F1 and F2.** Load the GS once in the universe
  (`GroundState::from_archive`, as the run_response/madqc path already does), then
  `copy` `amo` / `V_local` / etc. into each node-subworld under the pmap discipline.
  Per the feasibility study this leaves φ as **one distributed copy per node**
  (per-rank φ = |φ|/R) — the whole memory win.
- **(b) Per-subworld `from_archive` load** — simpler (no copy) but hits the
  ParallelInputArchive **nio/nproc np-mismatch** hazard (the archive was written by
  moldft with a fixed writer nproc; reading into an R-rank subworld is the
  `e544ae063` failure class). **Deferred** until archive np-robustness is confirmed.

## 4. F1 — standalone fan-out/gather proof (NO solver change)

New `tests/test_fd_subworld_fanout.cpp` (mirrors `test_es_build_subworld.cpp`):

1. Universe: load a real GS (h2o fixture) once; `set_response_protocol`.
2. Build the node-subworld pool (`make_node_aligned_subworld`).
3. Choose **K independent FD states** (dipole x/y/z at ω=0 → K=3).
4. **Reference:** solve all K in the universe via `solve_fd_protocol` (Fresh) → record
   each converged α component (e.g. `α_AA = −⟨source_A | x_A⟩`, the existing assembly).
5. **Subworld:** under the pmap discipline, ship the GS into each node-subworld;
   partition the K states round-robin by node index; each subworld calls
   `solve_fd_protocol` on a **subworld-bound `ExecutorContext`** (Fresh) for its owned
   states; gather the α scalars to the universe (`universe.gop` — scalars only,
   Functions never cross).
6. **A/B gate:** `max|α_sub − α_ref| < 1e-9`. Converged FD states are a *unique* fixed
   point (linear response), so they are path-stable across decomposition — the
   converged-only A/B discipline (keystone/S1; `verify_es_batch.sh`).
7. **Geometry:** needs **≥2 nodes** (≥2 subworlds) for a real partition; run 2 nodes × R
   (the S1/keystone geometry). Register in `CMakeLists.txt` + the `cm_build` ninja line;
   add `verify_fd_subworld.sh` (sbatch, 2 nodes × 8). Gate before any `run()` edit.

This exercises the exact function F2 fans out, with the gather via the existing
persistence — the only new code is the pool + partition + subworld-bound ctx/GS.

## 5. F2 — gated `run()` variant + knob + pre-flight abort (solver change)

Bring the diff for approval before coding. Pieces:

1. **Knob:** `ExecutorSettings::fd_subworlds` (int, default 0) ← `main.cpp`
   `--fd-subworlds=G` → `ExecutorContext`. Mirrors `--es-batch` / `--es-subworlds`.
2. **`run()` branch:** if `fd_subworlds > 0`, build the node-subworld pool **once**
   (per protocol or per run), ship the GS in (§3a). Replace the `:963-964` loop: split
   `wave`'s **FD/NuclearFD** items round-robin by node index; each node-subworld runs
   `run_protocol(item)` for its owned items on a subworld-bound executor; `universe.gop.fence()`.
   Non-FD items in a mixed wave (ES / DerivedFD) stay single-World for now (FD-only
   fan-out). `fd_subworlds == 0` is the byte-for-byte current path.
3. **Metadata merge (the one genuinely new persistence wiring):** per-state archives
   are distinct files (keyed by pert/freq) → no collision. `response_metadata.json` is
   the shared file → write **per-group shards** `response_metadata.group<gid>.json`,
   merged by **universe rank 0** after the fence (the v2 subgroup pattern; goes through
   the metadata layer — never write the merged file directly).
4. **Pre-flight memory estimate + clean abort:** before allocating the pool, estimate
   per-rank φ + response from the mem model (`rss ∝ n_occ`, ×k-growth) for the chosen
   `(S, R_state)`; abort with a clear message if it overflows the rank/node budget.
   Also gates the "large" regime until multi-node-per-state lands.
5. **A/B:** `--fd-subworlds=2` vs `0`, **converged states only**, 1 and 2 nodes → α match
   (the `verify_es_batch.sh` discipline). `--es-time`-style per-state wall for the win.

## 6. F3 — weak-scale sweep + auto-selector

- **`weak_scale_fd.sh`** (sibling of `strong_scale_es.sh`): hold states/node fixed
  (≈1 state/node, fixed occ/node), grow `(molecule, nodes)` together; capture per-state
  `wall_s` + `MEMORY_HWM`; `aggregate_weak_scale.py` → flat-line / efficiency plot. Also
  calibrates the doc 31 §4 GB thresholds (where per-state mem crosses rank/node budgets
  at the real k=10 protocol).
- **Auto-selector:** `(n_occ, k, available nodes, budget) → (S, R_state)` via the doc 31
  §4 regime model; manual `--fd-subworlds` override stays.

## 7. Risks / open points

- **Metadata write race** → per-group shards merged on rank 0 (§5.3) — the only new
  persistence mechanism; everything else is existing per-state archives.
- **GS np-mismatch** if per-subworld archive load → use ship-in (§3a).
- **Single-node degenerate case:** `n_nodes == 1` → one subworld == universe → the
  partition is a no-op and **must** equal the `G=0` path bit-identically. A free
  regression check (assert it in F1/F2).
- **Determinism:** all subworlds must read identical metadata post-fence (guaranteed by
  the `std::map` order contract, `:868-872`); the round-robin partition is pure.
- **Teardown discipline:** subworld before universe; restore `set_default_pmap`
  before `reset()`; clear all subworld Functions first (S1 lines 217-221).
- **Mixed waves:** a wave with FD + ES + DerivedFD items — F2 fans out only the FD
  subset, ES stays single-World. Confirm the partition + the non-FD remainder both fence
  cleanly on the universe.

## 8. Next action

Bring the **F1 diff** for approval: `tests/test_fd_subworld_fanout.cpp` +
`CMakeLists.txt` registration + the `cm_build` ninja line + `verify_fd_subworld.sh`
(sbatch 2 nodes × 8). **No solver touch** — the gate before the F2 `run()` change.

## Do-not-touch / contracts (carried from doc 31)

Reference kernels (`compute_V0x`/`compute_gamma`/`compute_E0x` — exchange's gate-0
oracle); the single-World `run()` path stays the `G=0` bit-identical reference;
perf-model's meter/profile schema; viz's `dump_mra_trees`/legacy export. Collectives
on all ranks (print-gate only); **never write `response_metadata.json` directly**
(metadata layer only); **USER runs builds/solves on the alloc**; propose-diff-first on
solver/runtime edits; commit with `git -c core.hooksPath=/dev/null`.
