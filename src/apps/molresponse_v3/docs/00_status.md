# molresponse_v3 — Active Workstream Status

**Purpose:** the single page an agent reads first so it does *not* re-derive
context. It says what is in flight, which files are hot (= conflict risk), how to
test each thing, and the standing contracts. Keep it short and current — update
the affected row in the same change that lands the work.

Last updated: 2026-06-08 (branch `molresponse-feature-next`) — WS2 stabilized
(`ca2c21dfc` + `3a7a1dfee`, doc-15 #5). **Ground-up rebuild started:** master
architecture is `docs/16_architecture.md` (L0→L5 layer cake, R0..R5 sequence,
state-parallel LAST). **R0a done:** `run_response` seam
(`orchestrator/response_workflow.hpp`) + `ExecutorContext`/`ExecutorSettings`
split + new `test_run_response` driver — h2o α bit-identical to baseline.
**R1a done** (`145383592`): stage timing → `Output.timing`. **R1b done**
(`ceed4716d`): uniform per-state `wall_s` + `MEMORY_HWM` across FD/ES/VBC
(`StateMetrics.wall_s`; worst-task `rss_gb_max`); validated PASS on a shared-dir
h2o run (all states `wall_s>0`, 84 `MEMORY_HWM` lines). Next: R1c (scheduler
trace + `PROTOCOL_START/DONE` → `Output.diagnostics`), then R0b, then R3 madqc.
Known follow-ups surfaced by the R1b run (NOT R1b bugs): ES stalls unconverged at
1e-4 (doesn't climb); β incomplete when dynamic VBC pairs don't all reach the top
protocol.

Workflow + build/run/validate harness: `cm.sh` in
`/gpfs/scratch/ahurtado/madness_es_bench/` (its `README.md` is the command
catalog). Run on a compute node via the `run-on-allocation` skill.

---

## Active workstreams

### WS1 — ES property printing  *(in progress)*
Post-convergence excited-state transition properties: transition dipole,
oscillator strength, transition quadrupole, dominant-occupied weights, AO
(sto-3g) Mulliken population. Port of legacy `TDDFT::analysis` /
`analyze_vectors`. Closed-shell only.
- **Hot files:** `solvers/es_analysis.hpp` (new), `calc/calc_executor.hpp`
  (wired into the TDA/Full ES solve, ~494/645), `tests/test_calc_manager_run.cpp`
  (`--es-analyze-only`).
- **Test:** `cm_es h2 3` (analysis runs after convergence → writes
  `es_analysis__<key>.json`); analyze-only: `test_calc_manager_run
  --es-analyze-only [--es-full]` on a calc dir with a converged bundle.
- **State:** computes + prints + JSON for TDA and Full; wired into solve path and
  standalone load path. **Not yet** validated against legacy/Dalton numbers, and
  **no `cm_` shorthand** for analyze-only.
- **Gotchas:** every inner/norm/`Function::size()` is collective — keep it on all
  ranks, gate only printing on rank 0 (see `analyze_response_orbitals`). Doc
  comments for `es_analysis_to_json`/`report_es_analysis` are currently tangled.

### WS2 — restart correctness + convergence control  *(stabilizing)*
Restart / reconcile precedence, seeding, divergence handling, convergence knobs.
- **Hot files:** `solvers/convergence_policy.hpp`, `calc/calc_executor.hpp`,
  `calc/calc_manager.hpp`, `solvers/fd_save_load.hpp`, `solvers/iterate*.hpp`.
- **Test:** `cm_idem <mol>` (Skip / restart-safety), `cm_resume <mol>`
  (cross-protocol restart), `cm_smoke`.
- **State:** ES→derived-FD expansion is metadata-driven + restart-safe; expansion
  node id now keys on the true ω (= ωₙ/2), guarded by a DAG id/freq invariant in
  `test_calc_manager_run`. Fresh-FD `x0=0` seed + nuclear-FD source fixes landed.
  Policy defaults `kain_maxsub 5` / `step_restrict PerState` **committed** (`ca2c21dfc`).
- **Reconcile/load unified (`3a7a1dfee`, doc-15 refinement #5):** one shared pure
  helper `best_usable_fd_source_key` (`response_metadata.hpp`) now backs BOTH the
  reconcile verdict (`reconcile_protocol`) and the archive load (`try_load_fd_state`),
  so they can never disagree. Fixes: (a) a coarser-or-equal **partial** (not
  converged, not diverged) is now a usable Restart seed instead of being discarded
  to Fresh; (b) `try_load` now **excludes diverged** snapshots (never seeds a
  blown-up state). Removed dead `has_coarser_converged_fd`. Two new `cm_unit` rows
  cover it. FD-only — ES/VBC coarser helpers unchanged (parallel follow-up).
- **Convergence-control knobs (new):** `--conv-factor=F` / `--bsh-factor` /
  `--density-factor` set the gate target = `F·max(thresh,dconv)` (default 5),
  exposing the previously CLI-unreachable `ConvergencePolicy` factor fields.
  `--accept-at-maxiter` records a **non-diverged** FD that exhausts `--maxiter`
  without hitting target as `converged=true` + an `accepted=true` marker (real
  residual still saved). cm.sh: `RAMAN_CONV`, `RAMAN_ACCEPT`; RUN CONFIG echoes
  `accept_at_maxiter`. **Scheduler contract learned:** a not-converged FD stays
  `Resume` → the wave signature repeats → `run()` halts on "no progress" and the
  node **never climbs to the next protocol** (so a stiff channel that misses
  target at a coarse rung blocks the whole ladder *and* the downstream VBC
  prerequisite gate). `--accept-at-maxiter` is what lets it advance — acceptance
  must fire at **every** rung, not just the last, or the climb never starts.

### WS3 — Raman / VBC save-reload  *(in progress)*
Single-component vibrational Raman = β(dipole; dipole, nuclear); VBC quadratic
source persistence; in-flight kernel-naming refactor.
- **Hot files:** `kernels/two_electron.hpp`, `kernels/vbc.hpp`,
  `solvers/vbc_save_load.hpp`, `kernels/{full,static,tda}.hpp`,
  `calc/calc_executor.hpp` (VBC save/load ~695/1020).
- **Test:** `cm_check_raman <mol>`, `cm_beta <mol>`, and `cm_equiv <mol>`
  (kernel-equivalence gate — run it after any `two_electron`/kernel edit).
- **State:** single-component Raman now runs **end-to-end** on h2o: dipole + nuclear
  FD → VBC quadratic source → β(dipole; dipole, nuclear) contraction → `[RAMAN]`
  tensor values (e.g. A=z → 4.68), recorded PASS via `cm_record`. VBC save/load
  wired. Refactor renames `ExPair→ExchangePair`, `apply_channel_raw→apply_gamma_raw`,
  "channel"→"gamma component" across the kernels — keep the rename complete +
  `cm_equiv` green.
- **Convergence learning:** the nuclear-displacement FD (`MolecularDerivativeFunctor`,
  a cusped ∂V_nuc/∂R source) **floors at ~4e-3 at k6/1e-4** — ~8× the 5e-4 gate —
  while the smooth dipole FD converges to ~3e-4. It is **resolution-limited**, not a
  solver-tuning issue: it only improves by climbing to k8/1e-6 (Restart reprojects
  the accepted k6 state to k8 and keeps iterating). Practical recipe:
  `PROTOCOL=1e-4,1e-6 MAXITER=… RAMAN_ACCEPT=1 cm_check_raman h2o` (see WS2 knobs).
- **Validation gap:** no Raman reference yet — `cm_record` reports PASS as "recorded,
  no ref". **Next:** get a molresponse_v2 Raman value for h2o to gate correctness.

### Cross-cutting — core-lib debug-logging tweak
`src/madness/chem/exchangeoperator.h` + `src/madness/mra/macrotaskq.h`: moved a
`MacroTaskInfo` parser print into a verbosity-gated `set_macro_task_info`. These
are **core MADNESS libraries** (affect all of madness, not just v3) — touch with
care; a full `ninja` (not just the cm targets) is the real check.

---

## Hot-file conflict map (read before parallelizing)

| File | WS1 | WS2 | WS3 | Notes |
|------|----|----|----|-------|
| `calc/calc_executor.hpp` | ✓ | ✓ | ✓ | **Highest conflict risk** — all three touch it. Serialize edits here or isolate in worktrees. |
| `calc/calc_manager.hpp` |   | ✓ |   | scheduling/reconcile |
| `solvers/convergence_policy.hpp` |   | ✓ |   | |
| `solvers/es_analysis.hpp` | ✓ |   |   | isolated — low risk |
| `kernels/*`, `solvers/vbc_save_load.hpp` |   |   | ✓ | run `cm_equiv` after |

---

## Parallel-agent protocol

1. **One workstream per agent.** Read this file + the `cm.sh` README first; do not
   re-explore the tree to rediscover what's above.
2. **Isolate writers of shared files.** If an agent will edit `calc_executor.hpp`
   (or any ✓-in-multiple-columns file), spawn it with `isolation: "worktree"` so
   parallel edits don't collide; reconcile on merge.
3. **Build/test only the cm targets** (`cm_build` → `cm_unit` → the workstream's
   `cm_*` command on the allocation), except the cross-cutting core-lib change,
   which needs a full `ninja`.
4. **Validate via the ledger:** `cm_record <calc_dir> <mol>` → PASS/FAIL vs
   `madness_studies/refs/madness_results.json` (commit-stamped history).

## Standing contracts (do not regress)
- Non-trivial solver/runtime edits: propose a concrete diff + get approval before
  writing.
- Never write `response_metadata.json` directly — go through the metadata layer.
- ES is **closed-shell only**; open-shell ES + Full-RPA open-shell are out of
  scope. ES Gaussian references use **d-aug-cc-pVQZ** (single-aug manufactures a
  phantom ~3% error on diffuse roots) — `madness_studies/refs/dalton_tdhf.json`.
- Commit with `git -c core.hooksPath=/dev/null` (a repo hook corrupts
  `.git/index`).
