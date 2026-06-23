<!-- THREAD BANNER (io-hdf5) — drop on merge to trunk -->
> **THREAD: io-hdf5** — branch `io-hdf5`, off `molresponse-feature-next`. IO survey + HDF5 prototype.
> Mandate, test recipe, do-not-touch, contracts: see the **io-hdf5** brief in `madness_studies/RELEASE_STATUS.md`.
> Start with a read-only IO audit + written plan before any core IO edits. Coordinate the HDF5 layer with `feat/amr-export` (VTKHDF) — one shared stack.
> Log progress under a `## io-hdf5 log` section in this file; the status below is inherited from trunk — append, don't rewrite.

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
- **PARALLEL I/O is NOT faster (and that's expected):** at NP=4 archive read got
  3-6× SLOWER than NP=1 (0.0052→0.0289 at k10) — both archive AND legacy use
  **nio=1**, so all disk I/O funnels through rank 0; more ranks just add gather
  (write) / scatter (read) MPI overhead, the disk transfer stays single-stream.
  HDF5 is not the bottleneck — the rank-0 funnel is, shared by both formats. TRUE
  parallel speedup needs **Parallel HDF5 / MPI-IO** (all ranks write 1 file via
  `H5Pset_fapl_mpio`, no rank-0 gather; also removes the rank-0 memory buffer — the
  relevant lever for large-system OOM). **Gated** on a GCC13 parallel HDF5 build
  (we have serial; cluster `-parallel` builds are Intel-ABI) + evidence the gather
  hurts at target sizes. Caveat: checkpoint I/O is a tiny fraction of wall time
  (solve dominates) — so this is a memory/scaling item, not a speed one.
- **NEXT (1→2):** (1) real h2o MO round-trip (true orbital tree) — extend
  `test_function_hdf5` with optional `--archive=<moldft restart>` mode reusing the
  FD-skeleton's molecule/header/`from_archive` scaffolding (h2o fixture exists at
  `…/v3_fd_skeleton/h2o/mad.restartdata`); then (2) wire the archive path into a
  real v3 restart site (opt-in). Later: cross-rank check (write NP1/load NP4),
  Layer B interop, Parallel-HDF5/MPI-IO (gated above).
- **Coordinate w/ `feat/amr-export`:** same system libhdf5 + bohr/cell + chunk
  conventions; pin in `operator_contracts.md` once P1/P2 land.
