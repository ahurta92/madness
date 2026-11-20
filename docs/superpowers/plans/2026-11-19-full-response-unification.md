# Full response unification implementation plan

**Execution status — November 20, 2026:** tasks 1–7 and the task 8 integration,
build, functional-test and review work are complete. Ten representative
scientific cases pass the retained numerical checks. An explicit virtual-AO
water 2PA variant also passes. The default cold solid-harmonic water case
converges to higher roots and fails its retained excitation-energy reference;
full milestone acceptance remains open for that case. No reference values or
guess defaults were changed. See the [test index](../../design/madqc-test-index.md#response-verification--november-20-2026).

> **For agentic workers:** Use superpowers:executing-plans to implement inline, task by task, with one independent final review. The execution status above and linked evidence distinguish completed work from outstanding acceptance gaps.

**Goal:** Bring every currently supported molresponse task and property through the agreed convergence, metrics, native checkpoint, restart, results, printing and timing conventions.

**Architecture:** Extend the adopted driven-state path and reuse the existing response DAG and numerical kernels. Applications retain functions, comparison history and numerical preparation; shared controllers consume measurements. Iterative and direct tasks share publication, source qualification and reporting, with direct tasks reporting completion rather than iterative metrics.

**Tech Stack:** Existing C++ templates, MADNESS World/native archives, JSON records, CMake/Ninja, CTest and Python workflow regressions.

**Spec:** [Response milestone](../../design/madqc-response-milestone.md), [design draft](../../../design_draft.md), [restart clarification](../../design/madqc-response-after-scf-refinement.md), [legacy convergence audit](../../design/madqc-legacy-convergence.md), [test index](../../design/madqc-test-index.md).

## What one pass means and what is realistic

One sustained implementation run can contain all eight tasks below, with focused interface checks at each boundary, followed by metric checks, representative convergence/value validation and one final review. Start from the existing implementation and the unfinished restart regression work; do not repeat completed SCF/alpha integration.

The realistic core target is tasks 1–6: same-model restart, common task/timing records, driven variants, coupled excited states, VBC/beta and the existing Raman/2PA workflows. Confidence is highest in restart/driven/reporting reuse; excited-state comparison history and derived-frequency replacement need the most careful integration around existing numerical work. Subworld ownership, imported seeds and inspection (task 7) need separate proof and remain required for the full milestone. Task 8 establishes functional coverage and representative convergence/value evidence across the integrated families. Full reference recalculation and detailed performance calibration belong to the later release-validation stage.

Aim to finish all integration tasks in the pass. If a supported feature is not connected or fails its functional/convergence checks, report it as unfinished. Do not equate this integration milestone with complete release validation or silently omit unsupported results.

The recorded costs of 15 response qctest decks in `src/examples/qc/CMakeLists.txt` total **3.90–4.09 hours in serial**, mostly at seven workers, with one ten-worker case. This is a scheduling estimate from historical comments, not a new benchmark. That full batch and detailed performance calibration are later validation work, not a mandatory before/after sweep for this first milestone. Representative checks, builds and interface tests still take time; there is no justified hour estimate for the coding yet.

## Integration boundary and validation order

This is component integration, not a solver redesign or performance-tuning project.
Preserve the existing equations, operator applications, subspace solve/rotation,
normalization, KAIN/plain update, restriction formulas and supported optimized
execution paths. Expose their existing work to the shared interfaces. Reuse
computed densities, kinetic quantities and raw BSH proposals wherever their
sampling conventions match; compute a missing required measurement once. Avoid
duplicate operator work, unnecessary function copies and new synchronization.
The goal is to preserve solver cost and convergence behavior at equivalent
settings, with minimal bookkeeping overhead.

The agreed measurement/stopping conventions still apply. A different stopping
request can change the number of updates without changing the numerical solver.
Record update/evaluation counts for analysis, but do not tune the solver to match
legacy iteration counts or treat a count difference alone as a regression. If satisfying the interface appears to require a numerical algorithm or
execution-path change beyond a demonstrated bug fix, document it separately and
resolve the scope before including it; do not hide it in integration work.

Tasks 1–7 use build and focused wiring checks: identities, source qualification,
archive round trips, decision dispatch, cached-work reuse, record structure,
collective errors and timing serialization. Small controlled calculation fixtures
are appropriate when needed to exercise those interfaces. Numerical metric checks
and representative converged-value checks run last, in task 8, after the response
paths are connected. Reasonable convergence to sensible values is the milestone
criterion, not an exact historical iteration count. Use existing numerical
oracles where useful and investigate clear incorrect values, nonconvergence or
false completion. Keep references unchanged; recomputing references and detailed
performance/calibration work come later. The full release matrix remains tracked
separately and is not a prerequisite for this integration milestone.

## Global constraints

- Reuse the current worktree/build and preserve unrelated edits. No commits, pushes or numerical-reference regeneration are implied by this plan.
- Keep one `response` convergence parameter set and `response.maxiter`; no excited-state aliases, named stages, new stall policies or automatic restart heuristics. Preserve supported CLI capabilities.
- Intermediate protocols use the approved legacy threshold-dependent floors. Final requested criteria and representation remain mandatory; maxiter gives a fresh update allowance on restart.
- Report maximum per-function BSH residual, model-defined density change, energy/excitation-energy change where applicable, and the positive response kinetic sum. Missing comparison history is unavailable, never zero.
- Retain SCF's existing orbital preparation. Do not add general occupied-orbital alignment or old-SCF retention for ordinary refinement.
- Native archives only for shared restart. Keep the existing safe publication/retention machinery during adoption; do not combine this work with a generation/attempt schema redesign or file-layout overhaul.
- Preserve existing solver formulas, operation ordering, supported tensor/batched/streamed algorithms and explicit imports; component integration must not silently retune or replace them. Do not add missing quadratic XC kernels, open-shell excited-state theory or open-shell property contractions.
- Optimizer redesign/restart and new geometry-workflow behavior remain deferred. Nuclear perturbation states and existing Raman analysis inside molresponse are included.
- Use shared steady-clock wall seconds, explicit participating-driver scope and current-invocation reporting. Keep legacy timers unchanged; add no synchronization solely for timing.
- Use SafeMPI/World operations, preserve native Function conversion rules and avoid shallow-copy mutations of comparison/checkpoint functions.

## Inventory: existing support versus shared adoption

Source code takes precedence over status badges. `FEATURE_STATUS.md` currently understates open-shell driven solving and overstates the scope of the new shared path; correct these descriptions during adoption.

| Path | Existing reusable implementation | Shared work remaining / evidence |
| --- | --- | --- |
| Closed-shell static/dynamic alpha | `FDSolver::evaluate/advance`, `iterate_shared`, per-state ladder/checkpoint session | Finish SCF-refinement restart; broaden HF/LDA/hybrid and molecular references; make shared path the default after migration. |
| Nuclear and open-shell driven states | `solve_fd_protocol<Type,Shell>`, X/XY alpha/beta storage, nuclear perturbation builders | Remove dipole/closed-shell adapter restrictions only as each variant passes shared measurements, persistence and timing tests. Open-shell properties remain explicitly unavailable. |
| TDA/full-RPA roots and transition properties | ESSolver subspace rotation, stable root labels, native roots/omega archives, kinetic sampling, analysis | Separate evaluation from KAIN in ordinary and streamed/batched variants; validate root correspondence; shared bundle decisions/checkpoints/results. |
| VBC and beta | Existing DAG, P/Q or equivalent VBC source, contraction and native storage | Exact qualified B/C inputs, direct completion/source record, selective rebuild, grouped tensors and timing. |
| 2PA | Root-derived frequency expansion, driven legs, source/contraction routines and regrouping tests | Same shared lifecycle for roots/legs/sources, refined-frequency continuation, exact property dependencies and partial results. |
| Raman | Nuclear states, mixed VBC, derivative tensors, existing Hessian/mode/spectrum analysis | Shared inputs/completion, partial/full derivative handling, grouped results and direct-step timings. No new optimizer work. |
| State-parallel execution | Node subworld scheduler and metadata shards | Bind native operations to the participating World, reconcile scalar records collectively, preserve task provenance and explicit timing scope. |
| DALTON seeds / native inspection | Import transformations, seed tools, tree/response inspection tools | Import as guesses followed by fresh checks; inspect shared records/native payloads without solving. Shared HDF5 remains deferred. |

## Interfaces and file boundaries

Reuse `IterationMetrics`, `MetricHelpers`, `ConvergenceSettings`, `LadderPlanner`, `CheckpointStore::publish`, `iterative_manifest` and `direct_manifest`; do not introduce another controller or scheduler.

- Generalize `solvers/shared_fd_checkpoint.hpp`, `shared_fd_session.hpp` and `shared_fd_input.hpp` into `shared_response_checkpoint.hpp`, `shared_response_session.hpp` and `shared_response_input.hpp`. Update callers and remove the superseded wrappers at cutover. FD-only property formulas stay in FD/property code.
- `response_shared::Checkpoint` keeps the existing `plan`, `activate`, `controls`, `expected`, `info`, `has_saved`, `reusable`, `complete`, `manifest` and native `load` responsibilities. Its iterative publication accepts `IterationMetrics`, `Decision`, `Progress` and an application archive-writer callback returning `checkpoint::Payload`; it does not inspect a `responses.size()==1` state shape.
- An iterative task definition supplies its expected subject IDs, model metric support/definitions, scientific settings, layout and input references. One FD state has one subject; one ES bundle has all required roots. Adapters own archive serialization and any conversion.
- `response_shared::DirectCheckpoint` uses the existing direct manifest/store primitives, a task definition, exact inputs and a writer callback. It has no maxiter, scalar convergence gate or invented residual.
- `response_shared::Session` owns the response request and task bindings, qualified-source access, grouped property publication and one response timing tree. The calculation manager owns dependency ordering and invalidation.
- Add `src/madness/chem/TimingRecords.h` for the approved small recursive timing node and collector. Expose a scoped collector with `begin(name,context)`, `end()` and `record()`; accept an injectable monotonic clock in tests. Keep one implementation of elapsed-clock arithmetic, reusing/moving `ElapsedTimer` as appropriate without changing old MADNESS timer APIs.
- Keep `solvers::iterate_shared` and application-owned `reset_shared`, `evaluate(state,info) -> IterationMetrics`, `advance(state)`. ES evaluation includes its normal preparation/root rotation and cached raw BSH work; advance alone performs the next update.

## Review focus

1. A changed model must not become a silent ordinary-restart guess; numerical refinement must not be rejected merely because a source identifier changed (tasks 1–3).
2. Root permutations/sign changes and genuine mixing must not fabricate zero changes or attach metrics to the wrong saved root (task 4).
3. Partial/capped states must never feed apparently completed VBC/properties; changing A alone must not unnecessarily rebuild VBC (tasks 5–6).
4. Subworld publications and error paths must be collective on the correct World; one failed rank must not strand the others (tasks 2 and 7).
5. Timing must describe current work, without double-counted bundle/root work, replayed history or misleading sums of concurrent tasks (tasks 2–8).

## Task 1: finish the ordinary restart rules

**Files:** `SCFScientificSettings.h/.cc`, `SCFCheckpoint.cc`, `shared_fd_checkpoint.hpp`, existing SCF tests, `test_shared_response_refinement.py` and madqc test registration.
**Interfaces:** Produce the model-only settings comparison consumed by both SCF and response; retain full settings for evidence qualification.

- [x] Before rebuilding, preserve the current executable with its checksum, build configuration, input/reference hashes and thread settings for later validation. Keep artifacts outside the source tree. This is preservation only; no baseline benchmark campaign is required before integration.
- [x] Reuse the unfinished [restart plan](2026-11-19-simple-restart-rules.md) and existing test edits. Its public geometry-rejection and changed-SCF-source tests already failed against the old behavior; record those logs, not a fresh claim of completion. Verify the rest of the new test's assertions after reaching them.
- [x] Implement ordinary model rejection and same-model numerical continuation. Preserve qualifying evidence for request-only changes, and per-state saved protocol for a refined SCF reference. Validate before replacing registered response inputs.
- [x] Exercise geometry/cell/XC/localization mismatches, smoothing/operator refinement, missing history, zero-update checks, repeated reuse and two-rank failures. Check saved payloads/provenance survive rejected requests.
- [x] Run existing settings/checkpoint/public-SCF tests and serial/two-rank response refinement tests. Expected: agreed behavior with unchanged SCF orbital preparation. Record any existing geometry-workflow interaction separately; no optimizer API is introduced here.

## Task 2: generalize task publication and implement shared timing records

**Files:** the three shared-response headers above; `IterationOutput.h`, new `TimingRecords.h`, chem CMake; response `iterate.hpp`, `calc_executor.hpp`, `response_workflow.hpp`; `test_shared_iteration.cpp`, new `test_timing_records.cc` and `test_shared_response_records.cpp`.
**Interfaces:** Produce common iterative/direct binding and session publication; preserve the existing FD evaluation API. Task 3–7 consume these interfaces.

- [x] Write tests showing publication supports a multi-subject bundle and a direct task; invalid/missing sources fail before publication and output-write failures propagate collectively. Watch the new behavior fail before changing the binding.
- [x] Generalize the existing FD implementation with archive callbacks and per-model subject/metric definitions. Preserve its current single-state behavior and all existing tests.
- [x] Add deterministic fake-clock tests for nested totals, missing/asynchronous splits, current-run reuse, and execution-scope overrides. Implement the design's recursive timing node and atomic `timings.json` writer at the existing calculation root.
- [x] Wire task totals and prepare/evaluate/decide/advance/checkpoint/property intervals. Measure shared ES bundle work once; record separately measured root-local work only where it exists. Use combined intervals when completion boundaries do not permit an honest split.
- [x] Feed console and JSON from the same records, preserving useful legacy diagnostics. Test a write failure on two ranks. Expected: no extra numerical evaluations, fences or archive writes solely for timing.
- [x] Run timing/shared-record/iteration wiring checks and the small public restart/ladder fixtures needed to verify integration. Reserve numerical metric and representative converged-value validation for task 8.

## Task 3: complete driven-state variants

**Files:** `fd_solver.hpp`, `fd_problem.hpp`, shared-response headers, `calc_executor.hpp`, `madqc_adapter.hpp`, model metric helpers as needed; new `test_shared_response_variants.py` and CMake entries.
**Interfaces:** All currently supported FD/nuclear FD tasks supply one-state measurements and native writer/reader callbacks to task 2.

- [x] Add failing public cases for nuclear perturbations, a one-electron open-shell solve and the existing multi-electron open-shell cases. Assert all stored spin/X/Y components participate, fixed frequency makes energy change inapplicable, and property limitations stay explicit.
- [x] Extend the existing evaluate/advance path to those variants using kernel-defined density and the agreed positive kinetic sum. Preserve plain/KAIN, exchange-tensor and step-restriction choices. Remove each adapter restriction only when its case passes.
- [x] Verify multiple frequencies/axes, partial state reuse, both changed-k and fixed-k transitions, interrupted intermediate work and new states. Continue existing per-state planning rather than imposing a global restart protocol.
- [x] Prepare the HF helium, LDA helium, hybrid helium and curated DALTON-reference water LDA alpha mappings for task 8. Do not infer supported meta-GGA response from SCF support; inventory the actual response kernels and queue their available numerical checks for final validation.

## Task 4: coupled excited-state adoption

**Files:** `es_solver.hpp`, `es_root_identity.hpp`, `es_kinetic_evidence.hpp`, `es_save_load.hpp`, `es_analysis.hpp`, `calc_executor.hpp`, shared-response definitions; new `test_shared_es_iteration.cpp`, `test_es_comparison.cpp`, `test_shared_excited_response.py`.
**Interfaces:** ESSolver provides `reset_shared`, `evaluate` and `advance`; one checkpoint contains the full roots/omega/identity bundle and one StateMetrics per required root.

- [x] Add focused interface checks for current-prepared-root record ownership, cached BSH work reaching advance, and no update after a stop decision. Cover TDA/full RPA and supported streamed/batched dispatch. Queue numerical residual reduction, density, omega and kinetic checks for task 8.
- [x] Expose the evaluation/update boundary around the existing preparation/subspace/rotation/BSH/KAIN sequence, preserving its numerical operations and reusing the raw proposal. Keep model equations and numerical kernels unchanged. No individual-root freezing or warmup is introduced in the shared path.
- [x] Keep a deep previous prepared sample per root. Reuse existing ordering and stable identities; validate the correspondence using normalized positive overlaps, then apply the same correspondence to density, omega and kinetic changes. Sign-align density comparisons. Ambiguity/replacement/rank deficiency makes changes unavailable and starts a new history.
- [x] Prepare pure permutation/sign, ambiguous 45-degree mixing, clearly matched and weak-overlap fixtures for final metric validation. Any internal comparison cutoffs govern measurement eligibility only; they must not rotate or alter the numerical solver. Record their implementation rationale and validate the boundary cases and H2 behavior in task 8. Do not inherit DALTON warning cutoffs silently.
- [x] Connect the coupled controller: all required roots and bundle checks must pass. Restart uses the whole bundle, fresh maxiter/KAIN/history, and saved protocol. Added roots retain usable functions as guesses but require fresh bundle checks.
- [x] Check bundle round trips, tightened-request dispatch, added-root membership, protocol transitions and partial-root record handling. Publish excitation energies and supported transition properties through grouped results with their qualified sources; verify numerical saved-state metrics in task 8.
- [x] Queue `response_h2_es_tda`, `response_h2_es_rpa`, `response_h2_es_rpa_nokain` and numerical matching/kinetic/native checks for task 8; use focused bundle/interface MPI checks during integration.

## Task 5: direct VBC and beta

**Files:** `vbc_save_load.hpp`, `calc_executor.hpp`, `calc_manager.hpp`, shared-response checkpoint/session, `response_property_metadata.hpp`; new `test_shared_response_dependencies.cpp`, `test_shared_quadratic_response.py`.
**Interfaces:** DirectCheckpoint consumes exact qualified B/C and reference inputs; the property publisher consumes exact A/B/C/VBC sources.

- [x] Add dependency tests for missing/incomplete inputs, tighter requests satisfied by saved evidence, A-only changes, B/C changes and damaged selected payloads. Expected behavior: wait/unavailable, reuse, or selective rebuild as specified; no fallback to a guess for property assembly.
- [x] Replace legacy VBC `converged=true` persistence with shared direct completion and exact input records, using the existing native serializer and source kernels. VBC runs after required inputs qualify at the final requested protocol; intermediate iterative advancement does not wait for VBC.
- [x] Keep readiness/invalidation in the existing manager. A-only changes reuse VBC; B/C/reference changes rebuild affected VBC and properties. Preserve earlier property source records.
- [x] Publish beta as grouped components/frequency combinations with units, availability and state metrics once. Capture direct construction, archive and contraction timings.
- [x] Run focused direct-task/restart checks on one and two ranks and preserve quadratic-DFT refusal/open-shell unavailability. Queue `response_lih_beta` and VBC source/spec/PQ-term numerical equivalence tests for task 8; the source and contraction equations stay unchanged.

## Task 6: 2PA and Raman through the same interfaces

**Files:** `ResponsePropertyPlanner.hpp`, `calc_manager.hpp`, `calc_executor.hpp`, `raman_analysis.hpp`, shared-response session; `test_tpa_regroup_identities.cpp`, `test_raman_spectrum.cpp`, new `test_shared_derived_response.py`.
**Interfaces:** Reuse task 3 driven, task 4 coupled-root, task 5 direct/property interfaces; retain symbolic derived requests and current task IDs.

- [x] Add restart tests for a refined excitation energy changing its derived FD frequency: reuse compatible saved response functions/protocol, rebuild the equation and perform fresh checks. Record replacement relations; old completed metrics/properties do not certify the changed frequency. Preserve unrelated independent states.
- [x] Cover both changed task IDs and a physical frequency change inside the existing rounded-ID bin without weakening ordinary independent-state identity checks. Retire obsolete tasks only after their live uses end, through existing retention primitives.
- [x] Route 2PA sources and contractions through qualified inputs, common completion, grouped tensors/observables and shared timings. Retain the existing preliminary scientific status.
- [x] Route nuclear Raman states, mixed VBC, derivative tensors and existing vibrational analysis through the same records. Test partial coordinate selections separately from full Hessian/mode/spectrum completion. Do not require optimizer work or add an outer geometry-restart system.
- [x] Run small dependency/restart cases during integration. Queue `response_lih_raman_full`, `response_h2o_raman_tpa` and `response_h2o_tpa` for task 8; the known lack of an independent full-Raman reference remains recorded.

## Task 7: supported parallel execution, seeds and inspection

**Files:** `node_subworlds.hpp`, `calc_executor.hpp`, `madqc_adapter.hpp`, `dalton_import.hpp`, `tools/dump_mra_trees.cpp`, `tools/seed_from_dalton.cpp`, `tools/tpa_from_dalton.cpp`; existing subworld/import tests and new `test_native_inspection.py`.
**Interfaces:** Same task definitions/native payloads, with participating World ownership and explicit execution scope. Full-run publication remains owned by the workflow/session.

- [x] Extend `test_fd_subworld_fanout` to shared publication/restart and inject an unavailable input/write failure. Keep a World-bound session out of another World; exchange serializable records at existing synchronization boundaries. Never dispatch universe collectives from a subset of ranks.
- [x] Merge task records/timing scopes without collisions or duplicate publication. Preserve per-state identity/progress and safe file publication. Queue numerical comparison with single-world references for task 8.
- [x] Preserve DALTON imported functions and explicit transformations as guesses. They receive fresh shared measurements at the requested protocol, no warmup/extra implicit stages. Check seed loading and shared dispatch; queue `response_h2o_tpa_seeded` for task 8 when its existing fixture/dependencies are available.
- [x] Migrate native inspection to shared task/source records. Inspect an FD state, ES root and VBC without numerical evaluation or checkpoint mutation; missing reclaimed payload reports unavailable. Keep optional HDF5 tools/builds separate from shared restart.
- [x] Expected: supported subworld/native seed/inspection capabilities remain available. Missing dependencies or failed cases remain explicit milestone gaps, not silently retired features.

## Task 8: default cutover, metric/convergence checks and review

**Files:** shared-response input adapter, `ResponseParameters.hpp`, `madqc_adapter.hpp`, `response_workflow.hpp`, response/madqc CMake, relevant result readers, `docs/design/madqc-test-index.md`, response guides/README/FEATURE_STATUS.
**Interfaces:** Ordinary `madqc --wf=response` selects the shared lifecycle for every migrated supported request; explicit metric keys are no longer an opt-in to a separate execution model.

- [x] Add default-input public tests proving all supported task families use shared records and decisions without an explicit convergence selector. Remove stale named-driver hooks and separate legacy stopping/publication branches only after all callers/test claims have migrated.
- [x] Curate each feature test by stable scientific intent, setup, physical quantity/units, source oracle, numerical tolerance and cost guard. Migrate format readers without replacing numerical references. Keep historical iteration/cost guards visible for later calibration; a difference in those counts alone does not block this milestone. Extend the existing test index rather than creating a competing inventory.
- [x] Verify grouped alpha/beta/2PA and excited/Raman outputs, numeric values plus availability maps, metrics reported once, exact sources, partial results, read-only inspection and consistent human output.
- [x] Verify `timings.json` for fresh/continued/reused/failed/direct/subworld paths. Terminal checks have no fictitious update; bundle work is not copied onto every root; concurrent totals are not added as elapsed time. Reconcile console and saved values within formatting precision.
- [x] Validate numerical metric definitions and sampling across all connected paths: maximum per-function BSH residual, density/omega/kinetic changes, relative scales, root correspondence and the actual functions loaded from archives. Exercise the numerical fixtures queued in tasks 3–7; test unavailable history and inapplicable quantities explicitly.
- [ ] Select representative converged calculations from the queued alpha, TDA/RPA, beta, Raman and 2PA cases, including seeded/parallel paths where supported. Confirm reasonable convergence, sensible values and consistent saved metrics; use existing numerical oracles and tolerances where available. Keep the remaining model/functional/reference matrix for later validation.
  Ten cases and the explicit virtual-AO water variant pass. The default cold water case fails the root-energy reference despite satisfying its selected convergence metrics; its assertion and default input remain intact. No representative run is still pending.
- [x] Build changed targets and run affected serial/MPI unit/functional tests plus the representative calculations. Use serial CTest with explicit MAD/BLAS thread settings. Do not regenerate reference data, change reference values or silently remove failing numerical assertions.
- [x] Save measured timings and update/evaluation counts for future study. Investigate obvious duplicate work or broken execution, but do not optimize toward historical iteration counts or require a full performance benchmark campaign for milestone acceptance.
- [x] Run clang-tidy on changed translation units. One independent final review covers this plan and the scoped increment; fix important findings with regressions and rerun affected tests. No repeated review loop is needed by default.
- [x] Record which response paths are integrated, functional checks passed, representative values/convergence, metric consistency and remaining gaps. Keep later release-matrix/reference/calibration work and existing scientific limitations explicitly separate.
- [ ] Accept the full first milestone once the outstanding default cold-water root-selection case is resolved or its scope is explicitly settled. Passing seeded and virtual-AO cases do not qualify the different default guess.

## Later validation stage

After this integration milestone, expand the full scientific/model matrix,
calibrate numerical defaults and performance where needed, and recompute the
appropriate MADNESS reference runs as a deliberate validation task. Preserve
independent external oracles, provenance and justified tolerances; new MADNESS
runs do not replace DALTON values. No reference recomputation is part of this
first milestone. Differences in iteration count are diagnostic data, not an
integration target.

## Execution checkpoints and boundaries

Task 1 -> task 2 -> task 3 -> task 4 -> task 5 -> task 6 -> task 7 -> task 8.
Keep each integration boundary reviewable with focused interface checks, then check metrics and representative convergence/values; proceed without repeated permission requests for the agreed work. If a physical ambiguity needs scientific judgment, complete independent work and ask a focused question. Do not lower standards to fit a session. At a deliberate interruption, state exactly which tasks are implemented and which are only planned.

The author's subsequent request to work on this plan authorized implementation.
The current execution status above supersedes the earlier planning-only pause.
Implementation and verification evidence are recorded in
`.superpowers/sdd/2026-11-19-full-response-unification/` in the active worktree.
