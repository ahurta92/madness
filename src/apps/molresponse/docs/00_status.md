
# molresponse_v3 — Active Workstream Status

**Purpose:** the single page an agent reads first so it does *not* re-derive
context. It says what is in flight, which files are hot (= conflict risk), how to
test each thing, and the standing contracts. Keep it short and current — update
the affected row in the same change that lands the work.

Last updated: 2026-06-19 (branch `madqc-refactor`) — fixed the v3 madqc adapter
ground-archive path resolution (`fs::proximate(work_dir, outdir)`, mirrors v2);
prior R3 entry below reconciled (climb fix is committed, not "pending"). Earlier:
2026-06-08 (`molresponse-feature-next`) — WS2 stabilized
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
~1149). **Real fix (COMMITTED + merged into `madqc-refactor`):** (1) madqc adapter
builds via `GroundState::from_archive` exactly like v2; (2) `from_archive` resets
`from_memory_=false` so `prepare()` reloads pristine MOs from disk on each climb
(restores pre-R3a behavior); (3) reverted the in-memory pristine-snapshot band-aid.
One mechanism (archive load) fixes climb + restart segfault. **Validation (SLURM job
`r4_resume_20260611`) surfaced a SECOND, separate adapter bug — now also fixed
(`madqc-refactor`):** the multi-node FRESH climb died with `could not find file:
resp/task_0/moldft/resp.restartdata` even though the archive existed. Root cause:
the v3 madqc adapter resolved the ground archive from `scf_calc->work_dir` **raw**,
but `work_dir` is stored relative to the top calc dir while `ResponseApplication::run`
has already chdir'd (ScopedCWD) into the response `outdir` → `ParallelInputArchive`'s
`access()` resolved the relative path against the wrong cwd. Fix mirrors v2's
`make_ground_context` (and CC2/TDHF/OEP): `fs::proximate(scf_calc->work_dir, outdir)`
(`madqc_adapter.hpp`). **Climb fix itself was never at fault** (archive existed,
prefix correct). **Root-cause follow-up FIXED:** the relative-vs-absolute work_dir
inconsistency (`SCFApplication` `Applications.hpp` ~158 stored a *relative* `pm.dir()`
while the Nemo path ~894 stored `current_path()` absolute) is resolved at the chokepoint
— `PathManager` now stores `std::filesystem::absolute(base)` (`PathManager.hpp`), so
`pm.dir()`/`work_dir` are absolute everywhere and a raw `work_dir` use is no longer
cwd-dependent (the original bug class is gone, not worked around). The adapter's
`proximate` is kept (now abs-vs-abs, still correct, and robust if PathManager reverts).
NB on coordinate systems: `proximate(p, base)` is only correct when `p` and `base`
share a system — making work_dir absolute alone (without absolute `outdir`) would have
broken the adapter fix; normalizing at PathManager makes both absolute together.
**Re-validation PASS 2026-06-19** (FRESH h2o climb via `cm_mq h2o` `PROTOCOL=1e-4,1e-6`,
k6→k8: archive opened, no "could not find file"; final `1e-06_k8` α_zz(static)=8.5328,
α_zz(ω=0.04)=8.5700; matches the R3a single-rung smoke 8.5346 to ~2e-3). On a 1TB
cache-mode node the HBM `numactl --preferred-many=8-15` had to be dropped (`cm.sh`
now auto-detects flat vs cache mode); run used `LAUNCHER="mpirun -n 2 --map-by=socket"`.
Caveat from Sweep 1: c2h4/c6h6/naphthalene
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

**Exchange branch (`exchange`, 2026-07-01):** FD exchange **tensor layer** (doc 28)
Inc-1 + Inc-2 **LANDED + VALIDATED** (`ad60c2673`). The `--fd-tensor` gate assembles
θ from shared Tx/g₀ convolution tensors (g₀ cached per protocol on
`ResponseGroundState::g0_alpha`) instead of the per-op `compute_V0x/compute_gamma`;
gate 0 (per-op reference) is untouched. Converged-α A/B vs the reference = rel 5.3e-6
(h2o, Static+Full, climb→1e-6/k8; `es_bench/fd_tensor_compare.py`). `exchange_ctx`
shape pinned in `operator_contracts.md` ("Tensor-layer exchange"). **Inc-3a/3b/3d
LANDED** (`ab0eb2ce1`/`d7f07c29d`/`9ef1004c0`): 3a fuses Tx+Ty into one Poisson wave;
3b tiles `build_pair_tensor` over the φ-row index (`--fd-tensor-tile=N`,
`ConvergencePolicy.exchange_tile`, default 0=off) to bound the n² peak; 3d adds coarse
`rs_ext_*` `PROFILE_BLOCK` meters (no-op unless `ENABLE_WORLD_PROFILE`) into
perf-model's WorldProfile. All bit-identical (`test_exchange_ctx`: Full 1.46e-6/1.61e-6;
tile=2 vs 0 → 6e-15; meters zero-effect). **Inc-3c slice 1 LANDED + VALIDATED**
(`0872818d8`): ES bundle γ via the cached g₀ tensor — `tda_batch::compute_gamma_flat`
(one Coulomb wave for M densities + per-root `contract_col(x_s,g0)`; TDA's only
exchange term is φ-only ⇒ per-iter exchange convolutions M·n²→0 amortized), gated
`--es-tensor` on `ESSolver::set_gamma_tensor` (SEPARATE from the bitwise `--es-batch`
gate; A/B-to-floor). Alloc A/B (h2o, 3 TDA roots, climb): converged roots max
|Δω|=4.0e-7 (tol 2e-5) and gate-1 wall 24% faster (1836s vs 2428s). Harness:
`es_bench/verify_es_tensor.sh`. **3c production wiring LANDED + VALIDATED** (`3410a45ec`):
`--es-tensor` → `ExecutorSettings.es_gamma_tensor` → `solve_es_tda_closed_shell`;
full-pipeline A/B on the alloc (cm_es h2o 3, climb, resonant-gradient) — both gates
ES 3/3 + derived FD 3/3 PASSED, persisted ω agree to max |Δω|=1.39e-08 across all 6
(pkey,root) pairs. NEXT = Full-ES/VBC reuse of the tensor layer (cross-thread w/
parallel-runtime); the M=1 FD path has no state-batch win. Cross-thread status on
the release board.

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

### WS-PR — parallel runtime / state-parallel  *(this branch: `parallel-runtime`)*
The runtime-architecture thread (cross-thread board: `madness_studies/RELEASE_STATUS.md`).
- **Fork RESOLVED — doc 31** (`docs/31_fd_multiworld_decomposition_decision.md`).
  `strong_scale_es` showed G=1 strong scaling walls (18% γ-eff at 32 ranks — *expected*;
  respect the ~5-occ/node floor). Parallelism comes from the **state** axis, not spatial.
  Decision: **two-axis decomposition (states × space), one solve per World, FD first**
  (FD states are coupling-free → pure fan-out/gather), first cut = one FD state per
  node-aligned subworld, simplest-first / measurement-driven. Sub-node packing (small)
  and multi-node-per-state (large) deferred; ES state-parallel reuses the World pool +
  keystone allreduce (doc 24 §2) later.
- **Proven primitives reused:** `make_node_aligned_subworld` (node_subworlds.hpp, Inc1),
  GS ship-in under `set_default_pmap(subworld)` (S1, doc 23), keystone A/S allreduce.
- **Seams:** `calc_executor.hpp` `run()` (:903; the single-group comment at :859-867 names
  this as the 15c STATE_PARALLEL design), `ExecutorContext::world` (:131), `solve_fd_protocol`
  (:254). Knob `--fd-subworlds=G` (0 = single-World reference).
- **F1 PASS** (`test_fd_subworld_fanout`): subworld fan-out α bit-identical to single-World
  (6.98e-12). Crash was a teardown bug (GS destructing after `finalize()`) — fixed by
  scoping World-bound objects in a block before `finalize()`.
- **F2a+F2b VALIDATED**: gated `run()` + `--fd-subworlds`. A/B on a real 3-node partition →
  α match 6.37e-12 vs G=0; per-group metadata shards merged by rank 0 (`merge_state_shards`).
  `fd_subworlds==0` byte-identical; `G≤1` short-circuits; ES stays single-World; cm_unit green.
- **F2d/F2f/F2g DONE + BENCHMARKED** (uncommitted): **F2d** tagged-stream logging (empty-default
  tag ⇒ G=0 byte-identical; suppresses redundant subworld GS/PROTOCOL banners). **F2f** per-node
  granularity — `--fd-subworlds=P` = subworlds PER NODE (`make_subworld_pool`: node-split →
  within-node contiguous split; NUMA via launch). **F2g** VBC fanned (fan subset = FD/NuclearFD/
  VBC; `merge_state_shards` unions vbc_states). **2026-06-30 bench (2 nodes×8, h2o β SHG, k6→k8):**
  ref=1646s · 2 subw=809s (2.04×) · **4 subw=644s (2.56×, BEST)** · 16 subw=710s (regress);
  **β bit-exact ≤1.1e-10 all 27 comps incl VBC.** Per-state k8 wall: 16rk=63s (over-decomposed),
  8rk=47s (sweet spot), 1rk=126s (serial). Two-axis S×R_state confirmed; sweet spot ~4 for h2o.
- **NEXT:** (a) **scheduler work-exposure** — SHG's 6 indep FD (3 axes × {ω,2ω}) serialized into
  3+3 waves → starves procs; batch all independent states per wave (Raman's Nocc×3 exposes more).
  (b) **F2e auto-selector** — cap P ≤ items/wave AND keep ranks/state ≥ spatial floor (~8). (c) size-
  aware partition (round-robin gave 10/8/7/5 at G=4). (d) **F2c** pre-flight mem abort (high P
  replicates φ per subworld → OOM risk on c6h6). Bigger systems = next test axis (c6h6/naphthalene).

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

## io-hdf5 log
*(thread `io-hdf5`, off trunk — IO survey + HDF5 prototype. Append, don't rewrite.)*

### 2026-06-19 — P0 audit + plan complete
- **Read-only IO audit done** + written plan → `docs/30_io_hdf5_survey_and_plan.md`.
- **Legacy stack mapped:** archive core (`world/archive.h`,
  `binary_fstream_archive.h`, `parallel_archive.h`), function persistence
  (`mra.h:2901/2907` free `save/load` → `ParallelOutputArchive` `nio=1`;
  `funcimpl.h`/`worlddc.h:2069` rank-0-gather single file), and the v3 call sites
  (`response_state.hpp`, `fd_/es_/vbc_save_load.hpp`; metadata = JSON, stays JSON).
- **Existing HDF5 prototype** (`src/examples/writecoeffs/`) is **not actually
  wired** — values-based text/JSON only; `writecoeff_hdf5.cc` includes the *text*
  `FunctionIO.h`; h5cpp demos aren't in the build list. `h5cpp`/`HighFive` not
  installed on the cluster.
- **Build reality:** no `find_package(HDF5)` in MADNESS CMake. GCC13/Milan
  **serial** HDF5 at `/gpfs/software/hdf5/gcc13/milan/1.14.3` matches the v3
  toolchain — enough for a rank-0-gather (`nio=1`-mirroring) prototype; parallel
  HDF5 deferred (would need a GCC13 parallel build).
- **`upstream/pr_writecoeff` checked** (user ask): it's the origin of the
  writecoeffs examples (already in trunk), and adds **no HDF5 build wiring** —
  furthest-along path is the validated values+JSON round-trip (`writecoeff3.cc` /
  `FunctionIO2.h`). Noted a NDIM==3 coords bug (`FunctionIO2.h:242-249`).
- **Decisions (doc 30 §6):** representation = **two layers** (A raw-coeffs for
  restart, B values+coords for plotting/MRChem, separate writer; build A first);
  concurrency = **rank-0 gather single file** (serial HDF5). Access layer =
  **native HDF5 C API** recommended, *pending user confirm*.
- **Access layer = native HDF5 C API (user-confirmed).**
- **P1 build wiring APPLIED** (isolated to molresponse_v3, zero core/root-CMake
  churn): `option(MADNESS_ENABLE_HDF5 OFF)` + gated `find_package(HDF5 COMPONENTS C)`
  + `test_hdf5_smoke` (pure HDF5 C-API round-trip, no MADNESS fn code) in
  `molresponse_v3/CMakeLists.txt` + `tests/test_hdf5_smoke.cpp`. No-op when OFF.
- **Test infra wired into cm.sh** (canonical `madness_studies/es_bench/cm.sh`):
  `cm_use io-hdf5` now injects `-DMADNESS_ENABLE_HDF5=ON -DHDF5_ROOT=…1.14.3`
  (override `HDF5_ROOT`); `test_hdf5_smoke` added to `cm_build`'s default targets
  (auto-skipped off-branch). Worktree verify script:
  `es_bench/verify_hdf5.sh` (self-sbatch → cm_use + cm_reconfigure + run →
  PASS/FAIL verdict), mirroring `verify_fd_tensor.sh` / `verify_es_build_subworld.sh`.
- **cm.sh bug fixed (affected ALL non-VTK verify scripts):** `cm_use` ended on a
  bare `[[ …VTK… ]] && echo`, returning 1 on non-VTK branches → `cm_use||exit`
  false-FATAL'd as "no worktree" (job 2031937 died in 3s despite resolving the
  worktree + flags correctly). Added explicit `return 0`. Confirmed `cm_use
  io-hdf5` → rc=0.
- **P1 COMPLETE — VERDICT PASS** (job 2032058, xeonmax node xm028, 2026-06-19):
  `find_package(HDF5)` resolved, C API compiled + linked + ran →
  `test_hdf5_smoke: PASS (HDF5 1.14.3, round-tripped 6 doubles)`. **Cross-arch
  confirmed:** the *milan*-built gcc13 HDF5 1.14.3 links AND runs on xeonmax — no
  separate xeonmax HDF5 build needed. Re-run anytime: `bash es_bench/verify_hdf5.sh`.
- **P2 Layer A WRITTEN** (first increment; NP=1, double; awaiting alloc build):
  `solvers/function_hdf5_io.hpp` (`save/load_function_hdf5` — structured `.h5`:
  `/meta` attrs + `/keys` int64 + chunked `/coeffs` double; stores raw coeff
  tensors + per-node has_children/has_coeff, mirroring `ar & coeffs`) +
  `tests/test_function_hdf5.cpp` (round-trip + size/time vs `madness::save/load`,
  prints VERDICT). Wired: CMakeLists HDF5 block + `cm_build` targets +
  `es_bench/verify_hdf5_function.sh`.
- **P2 Layer A NUMERICS PASS** (xeonmax, NP=1, 3D Gaussian k=8/1e-6): round-trip
  **bit-exact** — `(f-f_h5).norm2()=0.0`, `(f_h5-f_legacy).norm2()=0.0`. Size:
  HDF5 1,156,152 B vs legacy 1,167,512 B (**0.990×**). I/O: write 0.0129s/0.0063s,
  read 0.0048s/0.0022s (HDF5 ~2× slower at this size — chunked-per-leaf overhead,
  not a bottleneck). Two **test-harness** bugs fixed (not HDF5): (1) free
  `madness::load(f,name)` builds its archive from `f.world()` → aborts on a fresh
  (null-impl) Function; use `ParallelInputArchive(world,...) & f` instead (v3
  loaders already do). (2) `main()`-local Functions outlived `finalize()` →
  RecursiveMutex EINVAL abort at exit; wrap work in a scope (helper fn) so they
  destruct first.
- **P2 Layer A GREEN — sweep PASS** (xeonmax NP=1, off-center Gaussian, box
  [-50,50]³, k/thresh = 6/1e-4, 8/1e-6, 10/1e-8): all **bit-exact** (max_err 0.0),
  norm2=1.152477 steady (projection OK in big box; mid-process k-switch clean).
  **Efficiency:** HDF5 size 0.976→0.993× legacy (≈parity, slightly smaller); HDF5
  I/O ~2× slower — a *constant factor* ⇒ per-call overhead from **chunk-per-leaf**
  (`chunk=[1,k³]`), not data volume. (complex<T> dropped per user — efficiency
  focus.)
- **Contiguous `/coeffs` LANDED — speed parity** (re-ran sweep): write gap 2.3×→
  ~1.2×, read at k=8 now 0.92× (HDF5 *beats* archive), files shrank further (size
  0.949→0.988×). HDF5 now ≈parity on speed + consistently smaller; residual gap is
  small/noise (0.01 s ops). Bit-exact preserved.
- **Archive-backend path WRITTEN** (`save/load_function_archive_hdf5` in
  `function_hdf5_io.hpp`; awaiting build): the "vector-archive HDF5" — serialize
  via the optimized `ParallelOutputArchive<VectorOutputArchive>` gather (worlddc.h
  2067, thread-parallel + Gatherv), persist the byte vector to one HDF5 dataset,
  reverse via `ParallelInputArchive<VectorInputArchive>` (2328). **No custom
  archive class, no core edits** (Vector archives hold a *pointer* to the buffer,
  so the parallel wrapper's copy shares it; 2067/2328 are a symmetric pair). Key
  wins: **universal** (any type/tree-state, complex for free) + **multi-rank at
  nio=1** (unlike the NP=1-only structured path) + near drop-in for
  `madness::save/load`. Opaque ⇒ restart-only (structured stays for interop).
  `test_function_hdf5` rewritten as a 3-way A/B (legacy / structured / archive).
- **3-way A/B PASS** (xeonmax NP=1, box 100, k=6/8/10; all bit-exact 0.0):
  **archive reads FASTEST at every k** (~15-40% < legacy: e.g. k10 read 0.0046 vs
  legacy 0.0054 vs structured 0.0074) — one bulk `H5Dread` + in-RAM deserialize
  beats `BinaryFstream`'s interleaved read-during-deserialize. Archive size≈legacy
  (marginally smaller, one file vs header+.00000), archive write≈legacy. Structured
  = smallest files but slowest reads (multi-object open) ⇒ confirms it's the
  interop path, not restart. **Archive-backend is the restart winner.**
- **gzip (B) DONE — `deflate_level` opt-in, default OFF** (`save_function_archive_hdf5`,
  chunked + `H5Pset_deflate`; reads transparent; bit-exact). **Level-6 measured
  (box 100):** size 0.774×(k6)/0.832×(k8)/0.872×(k10) — win SHRINKS with k (dense
  high-entropy coeffs); CPU cost EXPLODES — write 2.9×/7.3×/**11.2×**, read
  3.3×/5.1×/4.7×. **Verdict:** bad trade for the hot restart path (memory, not
  disk, is the bottleneck here) ⇒ **restart = uncompressed archive**; gzip is an
  opt-in for disk-constrained/archival/transfer (legacy can't compress at all).
  **Level-1 measured:** ~identical size to L6 (0.779/0.836/0.875×) but still ~3-12×
  write cost ⇒ if compressing, use **level 1** (same size as 6, cheaper); gzip stays
  default-off. gzip question CLOSED.
- **NET P2 efficiency conclusion:** archive-backend (uncompressed) = the restart
  format — universal, multi-rank(nio=1), ≈legacy size/write, faster reads.
  Structured = interop only (Layer B). gzip = opt-in size knob.
- **MULTI-RANK NP=4 PASS** (xeonmax, box 100, k=6/8/10; all paths bit-exact 0.0):
  archive + archive+gz round-trip exact at 4 ranks (2067 gather → 1 file → 2328
  distribute); structured auto-skips (`skip(NP>1)`). Archive write/read ≈ legacy
  (read edges it at k10: 0.0289 vs 0.0326). **Archive file size byte-identical at
  NP=1 and NP=4** (legacy grew +57 B from per-client framing) ⇒ rank-count-stable,
  ideal for restart on a different node count (nio=1 ⇒ reader nio forced to 1).
- **P2 Layer A COMPLETE — COMMITTED `c8b8f994c`** (6 files, opt-in, no core edits).
  Archive-backend = validated restart format (universal, bit-exact, multi-rank,
  ≈legacy speed/faster reads, rank-stable). Structured = interop (Layer B). gzip =
  opt-in (level 1) size knob, default off.
- **PARALLEL disk I/O is serial (nio=1, both formats) — BUT archive gather scales:**
  no parallel disk I/O (rank-0 funnel; true MPI-IO is future, gated on a GCC13
  parallel HDF5 build — removes rank-0 gather + its memory buffer, the real
  large-system OOM lever). HOWEVER on real data at NP=4 the archive **write is 3×
  faster than legacy** (h2o 5 MOs: 0.419 vs 1.303 s) because archive uses one bulk
  `MPI_Gatherv` (2067) vs legacy's per-datum MPI streaming (2268) — gap grows with
  data/ranks. Reads ~tie (slight legacy edge from the distribute). Caveat:
  checkpoint I/O is a small fraction of wall time (solve dominates).
- **(1) REAL h2o MO round-trip PASS** (`--archive=<moldft restart>` mode in
  `test_function_hdf5`, +`GroundState.cpp`): 5 MOs, k=8, L=200, **bit-exact NP=1 &
  NP=4** (err 0.0). Archive total 28,698,479 B **byte-identical NP1≡NP4** (legacy
  28,765,466, ~0.2% larger + rank-varying); NP=4 write 3× faster (above).
- **(2) DONE — opt-in HDF5 restart wired into `response_state.hpp` + VALIDATED PASS**
  (job 2047076, xeonmax NP=1, 2026-06-29): `ResponseStateX<ClosedShell>::save/load`
  now writes/auto-detects `<file>.h5` when env `MADRESPONSE_IO_HDF5` is set (default
  = legacy archive, unchanged). New generic primitive `save/load_parallel_archive_hdf5`
  (callback over `ParallelOutputArchive<VectorOutputArchive>`; the single-`Function`
  `save/load_function_archive_hdf5` now delegate to it — no behavior change). New
  `tests/test_state_archive_hdf5.cpp`: 3-orbital state round-trip both paths
  **bit-exact (max_err 0.0)**, `.h5` auto-detected. CMake makes `molresponse_v3` +
  `test_calc_manager_run` HDF5-aware so the opt-in path is reachable from the real
  solver / `cm_run` / `cm_es`. **Closed-shell restart surface COMPLETE:**
  `ResponseStateXY<ClosedShell>` (dynamic-α / Full-ES, x_alpha+y_alpha) wired with the
  same opt-in pattern + validated PASS (job 2047082, X and XY both bit-exact 0.0,
  2026-06-29). Remaining: OpenShell X/XY (out of closed-shell scope) + cross-rank,
  Layer B interop, Parallel-HDF5/MPI-IO.
- **(3) END-TO-END VALIDATED in a live FD solve — PASS** (job 2049580, xeonmax NP=1,
  h2o static α_zz, 2026-06-30): the opt-in restart is reachable through the real
  solver because every FD/ES/VBC site already calls `ResponseStateX/XY::save/load`
  (`fd_save_load:115/239/328`, `es_save_load:107/368`, `vbc_save_load:57/105`), and
  the restart DECISION is metadata-driven (`response_metadata.json` stores the archive
  **basename**, format-agnostic) while the binary load auto-detects `.h5`. Proven with
  `MADRESPONSE_IO_HDF5=1`: (a) checkpoints write as `.h5` only (0 `.00000`); (b) a
  cross-protocol resume **loads** the `.h5` 1e-4 state to climb to 1e-6 (`action=
  restart`); (c) an idempotent rerun **skips** by recognizing the `.h5` via metadata
  (`nothing left to schedule`); (d) **α_zz = 8.532639 bit-identical to the legacy
  path**; (e) zero crashes. No code change needed — the wiring was already correct.
  `cm_run`/`cm_idem`/`cm_resume` drive `test_calc_manager_run` (HDF5-aware + in cm_build
  targets); just export `MADRESPONSE_IO_HDF5=1`.
- **⚠ GOTCHA — runtime HDF5 ABI conflict (cost ~5 jobs to find; cm.sh now pins it):**
  `load_<arch>.sh` auto-loads module `hdf5/parallel/intel2024.0/1.14.3`, which puts an
  **Intel parallel** `libhdf5.so.310` on `LD_LIBRARY_PATH`. We BUILD against the
  **gcc13 serial** HDF5 (`HDF5_ROOT`), but our binaries use **DT_RUNPATH** (searched
  *after* `LD_LIBRARY_PATH`), so the Intel `.so` won at load time → header≠runtime ABI
  → corrupt HDF5 ID tables → **infinite recursion in `H5E__push_stack`/`H5I_inc_ref`
  → SIGSEGV** on the first HDF5 call (looked like a write bug; wasn't). Diagnosed via
  `gdb` break on `H5E__push_stack` (backtrace named the wrong `.so` path). **Fix:**
  `cm_use io-hdf5` now prepends `${HDF5_ROOT}/lib` to `LD_LIBRARY_PATH` so the
  build-matched HDF5 loads first (confirmed: `ldd` resolves gcc13, both tests PASS).
  Defensive: `function_hdf5_io.hpp` now disables HDF5's (recursion-prone) auto-error
  printer + checks every H5 return code (`MADNESS_CHECK_THROW`) → a real future error
  aborts cleanly instead of stack-overflowing. **Permanent-fix follow-up:** link the
  HDF5 targets with `-Wl,--disable-new-dtags` (RUNPATH→RPATH, wins over LD_LIBRARY_PATH)
  so the binary self-protects even outside cm.sh.
- **Coordinate w/ `feat/amr-export`:** same system libhdf5 + bohr/cell + chunk
  conventions; pin in `operator_contracts.md` once P1/P2 land. **Directly relevant to
  the "one shared HDF5 stack" contract:** the cluster has ≥2 `libhdf5.so.310` (gcc13
  serial + intel parallel) with the same SONAME — viz must build+run against the SAME
  one we pin here, or hit the identical ABI crash.

## perf-model log
*(measurement-arm thread; design anchor: `docs/29_perf_model_design.md`. Append
newest-first; the status above is inherited from trunk — do not rewrite it.)*

- **2026-07-01 (pm) — slot-fixed sweep (job 2051561): FD meter ✓, rank axis = sync catastrophe.**
  7/7 profiles (rank axis runs now). **FD projection meter validated:** the
  `projection` phase appears in **all** FD breakdowns (was ES-only); small as
  expected (≤0.3%). **Rank axis result (important):** splitting a *single small*
  FD state (h2o, n_occ=5, k6) across ranks is catastrophic — NP=1 wall 60s →
  **NP=2 ~1900s / NP=4 ~1730s (~30×), 77–80% in fences** (millions of msgs, GBs).
  This **quantitatively confirms parallel-runtime doc-31**: parallelism must come
  from the *state* axis, not spatial single-state. So the comm/sync cost model
  belongs on the **state axis (subworlds)**, and the NP=1 compute model must NOT
  be pooled with these sync-bound points — doing so blew the fit up (LOO 1057%,
  negative fence coef). The good result stands: **NP=1 `apply·k³` compute model,
  LOO 13%.** ⚠ `total_wall_s` looks **over-accounted at NP>1** (the two multi-rank
  walls sum > whole-job time — likely spin-wait, no ENABLE_NEVER_SPIN); intrinsic
  sync%/msg counts are trustworthy, absolute wall needs a re-check.
  - **NEXT:** (1) split the fit — compute model on P==1 rows only; report multi-rank
    as a separate parallel-regime diagnostic (don't blend). (2) verify the NP>1
    `total_wall_s` accounting (spin vs real). (3) the meaningful multi-rank scaling
    is **state-parallel** (parallel-runtime `--fd-subworlds`), not single-state
    spatial — calibrate the comm/sync model there.

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
