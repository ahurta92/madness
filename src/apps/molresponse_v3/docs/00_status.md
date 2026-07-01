<!-- THREAD BANNER (perf-model) — drop on merge to trunk -->
> **THREAD: perf-model** — branch `perf-model`, off `molresponse-feature-next`. Instrumentation + cost model (the measurement arm).
> Mandate, test recipe, do-not-touch, contracts: see the **perf-model** brief in `madness_studies/RELEASE_STATUS.md`.
> Exposes the meter API + profile schema that `exchange` reports into and `parallel-runtime` uses to settle the doc-24-vs-25 fork.
> Log progress under a `## perf-model log` section in this file; the status below is inherited from trunk — append, don't rewrite.

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
h2o run (all states `wall_s>0`, 84 `MEMORY_HWM` lines). **R1c done**
(`ba509aa53`): scheduler trace (per-wave reconcile actions, stop_reason, passes) →
`Output.diagnostics` + `PROTOCOL_START/DONE` lines. **R1 (observability) COMPLETE**
(R1a+R1b+R1c; R1d operation-timing optional). **R0b done** (`db32d7038`): main.cpp
rewritten as a thin run_response app (installed binary drives the seam; `--archive`
CLI); retired legacy test_v3_solver/test_es_solver + dead top-level FDSolver.hpp/
ESSolver.hpp (kept the live ESSolverGuess/ResponseFunctions/ResponseKernel chain).
**L1 (contract+orchestrator) COMPLETE.** **R3a done** (`6463f21d9`): madqc v3
response engine behind `response.engine=v3` (alpha only) — `madqc --wf=response`
drives v3 via `ResponseApplication<molresponse_v3_lib>` (GroundState from the
in-memory SCF; run_response split into run_response_with_ground). Validated by
single-node madqc smoke: h2o α_zz=8.5346 in calc_info. Fixed 2 bugs (GroundState
reloading a nonexistent archive on the in-memory path → `from_memory_` flag;
doubled calc_dir). **R3b done** (`ff475c5bd`): multi-property mapping
(requested_properties → polarizability + hyperpolarizability + single-component
raman + resonant/excited via merge_plans) + assembly now does alpha AND beta (not
XOR). Validated: h2o through madqc engine=v3 yields α_zz=8.5346 AND β_zzz=7.760.
R2 export/viz is a parallel agent (dump_mra_trees etc.). **cm routed through madqc**
(`cm_mq`/`cm_mq_beta`/`cm_mq_es`). **R4 diagnostic-study harness set up** (scratch
`cm_study` + repo `madness_studies/refs/study_analyze.py`): sweeps the production
path over a molecule set at a COARSE probe protocol (single static α, 1 FD state —
the 15d "measure low-k" trick keeps c6h6/naphthalene safe), appends per-state
metrics (coeffs/rss_gb/wall_s/iters) to `refs/r4_study_runs.tsv`, and prints the
empirical `mem_per_task(n_occ,k)` model (rss/coeffs ~ n_occ fit + k6→k8 growth
factor) feeding L2 (15d pre-flight abort / 15c subworld sizing). Fixtures added:
`c6h6` (12 atoms, 21 occ), `naphthalene` (18 atoms, 34 occ) — geometry-only;
cm_mq re-runs SCF. **R4 study ran (SLURM job 2005092, 4 nodes):** Sweep 1 (n_occ
scaling at k6) PASS — clean linear `mem_per_task`: `rss_GB ~ 0.93 + 0.031*n_occ`,
`coeffs/state ~ 3.51e5*n_occ` (2.81 MB/occ-orbital) across h2o(5)/c2h4(8)/c6h6(21)/
naphthalene(34). **Sweep 2 (k-growth k6->k8) surfaced a real bug:** the multi-protocol
climb on the **madqc / in-memory GroundState path crashes** ("tensors do not conform").
Root cause: `GroundState::prepare` reproject is guarded by `original_k_ != target_k`,
valid only because the archive path reloads pristine MOs each call; the from_memory
path skips that reload but mutates `scf_->amo` in place, so climbing back to
`original_k_` skips reprojection and leaves stale coarse-k orbitals. **R3 smokes were
all single-rung — first time the madqc climb ran.** A SECOND bug then surfaced on the
restart-in-place path (rerun madqc in an existing dir): **segfault in
`GroundState::build_fock_matrices`** because madqc validates a valid SCF archive as
`Ok` and `lib_.calc()` constructs the SCF WITHOUT loading MOs (`Applications.hpp`
~179-184/634) → the in-memory `scf_calc->amo` is **empty** → the v3 adapter (which
built `GroundState` from the live SCF) dereferences nothing. **Both bugs share one root:
R3a's in-memory `GroundState` shortcut.** v2 is immune — it loads the ground state from
the moldft archive (`scf_calc->work_dir`+`prefix+.restartdata`, `MolresponseLib.hpp`
~1149). **Real fix (applied, NOT committed):** (1) madqc adapter now builds via
`GroundState::from_archive` exactly like v2 (resolve archive from `scf_calc->work_dir`);
(2) `from_archive` resets `from_memory_=false` so `prepare()` reloads pristine MOs from
disk on each climb (restores pre-R3a behavior); (3) reverted the in-memory pristine-
snapshot band-aid. One mechanism (archive load) fixes climb + restart segfault; compiles
clean. **Validation pending** (SLURM job 2015778: FRESH h2o climb via from_archive +
restart-in-place c2h4/c6h6 maxiter=80). Caveat from Sweep 1: c2h4/c6h6/naphthalene
hit the 25-iter cap at k6 (static α unconverged — memory robust, wall is to-cap). Then
R5 state-parallel. (raman maps single-component only — full tensor
deferred; ES uses default SolidHarmonics guess — VirtualAO/es-guess madqc knob =
follow-up.) Open follow-ups: ES stalls unconverged at 1e-4 (blocks ES/2PA/
resonant-Raman + R2 ES-density export); β incomplete when dynamic VBC don't all climb.
**ES-guess work (active, doc 17):** **A) `ESGuessMode::VirtualAO` DONE** (`18f853182`):
virtual-orbital "NWChem" CIS-diagonal guess — h2o recovers all four roots in order
(incl. the 0.378 root SolidHarmonics missed). B) Dalton restart via
`Dalton_Interface : ES_Interface` / Molden adapter; C) seed from Dalton excitation
vectors — both still on the roadmap.
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

---

## perf-model log
*(measurement-arm thread; design anchor: `docs/29_perf_model_design.md`. Append
newest-first; the status above is inherited from trunk — do not rewrite it.)*

- **2026-07-01 — fit sweep ran (job 2051068); the compute proxy was wrong, fixed.**
  5/7 profiles (k axis clean: h2/h2o at k6 *and* k8, context confirms k=6/8). The
  **rank axis (h2o np2/np4) failed on a SLURM slot config** — `mpirun -np>1` needs
  `ntasks-per-node >= max NP`; was 1. **Fixed** in the sbatch (→ 4). No comm term
  yet (still all NP=1 → bytes=0).
  - **The honest fit exposed a model-form bug.** `coeffs·k` linear: in-sample
    R²=0.9955 but **LOO mean |err| = 122%** (negative intercept, negative predicted
    walls). LOO immediately flagged what R² hid.
  - **Root cause + fix (data-backed):** `apply` is 85–92% of work and its cost
    tracks **apply-call-count × k³** (per-box 3D convolution), NOT coefficient
    count. Across the h2..c2h4 / k6..k8 sweep `coeffs·k` scatters **9.9×** while
    `apply_calls·k³` collapses to **1.9×** (C ≈ 16 µs/(call·k³)). Swapped the proxy
    in `refs/perf_model_fit.py`. New model `wall ≈ c0 + c1·(apply·k³/P) + c3·fences`
    → **LOO 122% → 13.3%** (h2o/c2h4 to 1.6–14%; tiny h2 ~25% = fixed overhead).
  - **Parallel efficiency:** effective threads = Σcpu/wall = **4.3–5.9 of 10**
    (43–59%), rising with problem size — a size-dependent `threads_eff` is the next
    refinement (and the parallel-runtime φ input).
  - **NEXT:** (1) resubmit the slot-fixed job → the **comm/fence term** (rank axis)
    for multi-node prediction; (2) the a-priori predictor for a NEW molecule still
    needs **apply_calls(n_occ, k)** — a separate "how much work" model (tree size ×
    iters) — that's the path to predicting C6H6/naphthalene before launch.

- **2026-06-30 (pm) — broadened fit sweep prepped (ready to submit).**
  `es_bench/perf_model_fit_sweep.sbatch`: 7 shapes across two new axes — **k**
  (h2/h2o at k6 *and* k8, single fixed rung each, no climb) + **rank** (h2o at
  NP=1/2/4, the only rows with non-zero bytes/fences). c2h4 stays k6 (n_occ
  anchor). Reuses the WORLD_PROFILE build (`cm_build` incremental). Upgraded
  `refs/perf_model_fit.py` to make the fit honest: compute proxy is now
  **`coeffs·k / P`** (§9 `N·c/R`, so multi-rank rows are coherent), **drops
  all-zero columns** (kills the R²=1 artifact), reports **in-sample R² with an
  (over)determined warning**, and adds **leave-one-out mean |err|** (the real
  generalization number). Manifest gained an optional 4th `geometry` column so
  labels can differ from the molecule name (h2o_k8). Verified on job 2049760's
  data: 3pts/3params now correctly reports "overdetermined, LOO skipped".
  Submit: `sbatch es_bench/perf_model_fit_sweep.sbatch`.

- **2026-06-30 (pm) — first SLURM sweep ran (job 2049760); context+build validated.**
  `es_bench/perf_model_sweep.sbatch` (exclusive xeonmax node, NP=1, fixed
  PROTOCOL=1e-4/AXES=z/MAXITER=10, 2 reps × {h2,h2o,c2h4}). **All green on the
  infrastructure:** the context-block edits **compile** (first build incl. the
  worldprofile core-lib change); **`context` emits non-null** `{n_threads:10,
  k:6, thresh:1e-4}` on every run; **stability PASS** (call counts bit-identical
  rep1↔rep2, phase_set_diff=0; cpu noise 7–13%); **`apply` = 85–87%** of work
  across all three (the phase to model); no crashes.
  - **The "fit" is NOT real yet — do not trust its coefficients.** 3 shapes, and
    at NP=1 `bytes=0` so only 3 columns are active `[1, coeffs·k, fences]` → the
    design matrix is square → residual 0 → **R²=1.0000 is an artifact
    (interpolation, not prediction)**; the intercept even comes out negative
    (−7.41). Need #shapes > #params **and a held-out point** before the fit means
    anything.
  - **Real finding: thread parallel efficiency is only 30–45%.** cpu_total/wall ⇒
    effective threads ≈ **3.0 / 4.5 / 4.1** for h2/h2o/c2h4 (of 10 nominal) — so
    the `cpu_max/n_threads` wall estimate overpredicts by ~2.2–3.3×. The cost
    model must carry an **effective-parallelism (threads_eff)** term, not nominal
    threads. This directly feeds parallel-runtime's φ/imbalance term (doc 32 §6).
  - **NEXT (broaden so the fit is real):** (1) more shapes — molecule × k (run each
    at a fixed single rung k6 *and* k8, not a climb, to avoid the iter-count
    confound) → variation in both n_occ and k; (2) a **rank axis** (NP=1,2,4 on one
    shape) → strong-scaling + populates the currently-zero comm/fence terms;
    (3) re-fit with a held-out shape → a real R². Also still open: FD projection
    meter (rs_projection fires only on ES), and committing the context block +
    sweep script once you're happy.

- **2026-06-30 — PM-1 validated; PM-3 consumer verified; context block applied.**
  - **PM-1 confirmed working:** `p2.json` (Jun 23, h2o α, schema v1, 67 phases)
    proves the env-gated rank-0 emitter emits. Ran `refs/perf_model_fit.py` against
    it (first real breakdown): **`apply` = 89.7%** (`FunctionImpl::do_apply` — the
    phase to model), `sync` 7.3%, `rs_exchange_gamma` 0.3%, compress/reconstruct
    ~0 (folded inclusively into apply).
  - **Two gaps found.** (a) **`rs_projection` only fires on ES** — `rs::project`
    is called only from `es_solver.hpp`; the **FD (α/β) path projects inline**
    (`static.hpp:150`, `static.hpp:297/309`, `full.hpp:204/216/384`), so the
    `projection` phase is absent from every α/β profile. (b) **wall-est was wrong
    (1.96 vs 4.31 s)** because `context:null` hid the run's thread count — the
    offline fit fell back to a wrong cpu→wall divisor.
  - **Context block APPLIED (NOT yet built — core-lib, needs full `ninja`).**
    `WorldProfile::dump_json(world, path, context_json="null")`
    (`worldprofile.{h,cc}`); both v3 call sites (`main.cpp`,
    `test_calc_manager_run.cpp`) fill `{n_threads, k, thresh}` (the offline join
    still recovers n_occ/molecule via `--geometry`/`--metadata`). Default `"null"`
    preserves prior output exactly. Updated `refs/perf_model_fit.py` to **prefer
    `context.n_threads`** as the cpu→wall divisor (verified: clean fallback for the
    pre-context `p2.json`; context overrides a wrong `--threads`).
  - **NEXT (you run on alloc):** (1) `cm_rebuild` (full ninja for the worldprofile
    core-lib edit) + `MADQC_PROFILE_JSON=$PWD/p.r1.json cm_run h2o` ×2 → confirm
    `context` emits + structural counters (count/nmsg/nbyte) bit-identical across
    runs. (2) **FD projection meter** — `PROFILE_BLOCK(rs_projection)` at the FD
    inline-Q sites above so `projection` is populated for α/β (`cm_equiv` after —
    PROFILE_BLOCK is zero-numerics). (3) **The (k, n_occ) sweep** → `perf_model_fit.py
    --manifest` for the actual cost-model fit (≥3 runs).

- **2026-06-19 — thread bootstrapped.** Read the cross-thread board + brief +
  contracts + runtime/perf guide (`docs/parallel_runtime_guide/`, companion
  `parallel_runtime_and_performance_models.md` §9 model / §12 measurement plan).
  **Key finding:** the "per-phase timers/counters in the core" the brief asks for
  *already exist* — core MADNESS `WorldProfile` (`src/madness/world/worldprofile.{h,cc}`)
  is compile-gated (`ENABLE_WORLD_PROFILE`, default OFF; macros → no code when off
  = the zero-effect-when-off contract), already wraps apply/compress/reconstruct/
  project/truncate, and already captures CPU time + msg/byte counts parallel-reduced
  per call-site. **Decision (user-confirmed): extend `WorldProfile`** rather than
  build a v3-local meter. Gaps = this thread's deliverables: a machine-readable
  JSON emitter + pinned schema (PM-1), a canonical phase taxonomy + 2 named
  response-level meter blocks for exchange/projection (PM-2), and the cost-model fit
  (PM-3). Wrote design doc 29; pinned draft profile schema in `operator_contracts.md`.
  **PM-1 APPLIED (not yet built/committed):** `WorldProfile::dump_json(world,path)`
  added (`worldprofile.{h,cc}`; self-contained binary-tree gather mirroring `print`,
  hand-rolled JSON, no nlohmann dep in `world/`). Env-gated collective call
  (`MADQC_PROFILE_JSON`) wired into BOTH the app (`main.cpp`) and
  `tests/test_calc_manager_run.cpp` — the latter is what `cm_run` actually launches
  (cm_run runs `test_calc_manager_run`, NOT the app binary). Harness: added a
  `perf-model)` case to `cm_use` (`es_bench/cm.sh`) so the branch auto-builds with
  `-DENABLE_WORLD_PROFILE=ON` (mirrors the viz-branch VTK precedent). Schema v1
  pinned in `operator_contracts.md`.
  **NEXT (you run on alloc):** `cm_use perf-model; cm_rebuild` (force reconfigure to
  pick up the flag — turns on all PROFILE_ macros → full lib rebuild) then
  `MADQC_PROFILE_JSON=$PWD/p.r1.json cm_run h2o` ×2; confirm JSON emits + structural
  counters (count/nmsg/nbyte) bit-identical across runs (cpu within noise). Then
  PM-2 (phase taxonomy + `rs_exchange_gamma`/`rs_projection` blocks) and PM-3
  (`refs/perf_model_fit.py`).
