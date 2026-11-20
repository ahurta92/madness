# Drop 6 reuse audit: shared chemistry task lifecycle

- **Date:** 2026-11-10
- **Design:** [design_draft.md](design_draft.md)
- **Scope review:** accepted by the author on 2026-11-10; record details and remaining design choices are still under review.
- **Baseline:** `madqc-drop-6`, commit `9179f5d32076a2ba555906dd146dc330ffc040d7`
- **Current checkout:** `76255552e16e542f6b5ae2d0d7e2b6b034c7ed48`
- **Method:** source inspection of the pinned baseline, comparison with the current checkout, and inspection of existing test assertions. No builds or numerical tests were run. This is an implementation reuse audit, not an independent design review or release certification.

## Finding

Drop 6 supplies much of the infrastructure needed for the proposed functional
design: application callbacks, typed response states, iteration/protocol drivers,
metric operations, restart planning, archive readers/writers, and checked
checkpoint publication. The work is principally adapting those pieces to one
contract and connecting the retained workflow path to them.

The complete behavior is not already implemented. Current-state evaluation,
consistent stopping metrics, Nemo integration, excited-root correspondence,
direct-task dependency evidence, and the agreed output/timing records still
need work. The source inspection does not support an estimate of the percentage
complete or the implementation effort.

The current checkout adds three fixed-intermediate-work commits after drop 6.
Those changes and their tests must be adapted to the newly approved intermediate
heuristic policy. They are not evidence of behavior already present in drop 6.
For example, `test_iteration_protocol.cpp` exists in this checkout but not at
the pinned baseline. Preserve unrelated fixes and useful assertions when changing
that policy; this audit does not recommend resetting the checkout.

## Smallest shared boundary

Applications retain their numerical state. Operations can act on that state
through typed functions or callbacks; there is no requirement for a universal
orbital container, a new task class hierarchy, or a second scheduler.

| Boundary | Shared responsibility | Application-supplied operation |
| --- | --- | --- |
| Configure and identify | Resolve criteria, intermediate/final limits, protocol position and task identity; validate prerequisite evidence | Declare model, supported metrics, numerical settings and required inputs |
| Restore | Decide reuse, resume, converted guess or fresh work under the agreed rules | Read the existing payload format, check scientific compatibility, and prepare/project functions and operators |
| Iterate | Apply the common measurement, stopping, counting and recording conventions | Prepare the current state, evaluate its raw proposal and model quantities, then advance only if requested |
| Compute directly | Require eligible inputs and record completion | Build VBC or compute the requested property; no iterative stopping test |
| Publish | Validate and publish matching payload/evidence, retain required generations, and write common records | Serialize the complete numerical state and supply physical property values |

SCF and Nemo can keep their loops. Response can keep its shared drivers.
Workflows retain cross-application dependencies; molresponse retains its DAG.
The common evaluator handles reductions, availability and comparison rules;
model callbacks supply physical fields and quantities. Direct tasks use
completion and exact input identities rather than fabricated convergence metrics.

## Shared pieces: reuse, adapt, add, retire

Here, **reuse** means retain an existing mechanism, **adapt** means change its
contract or integration, **add** identifies missing behavior, and **retire**
means remove a superseded path only after migration and collaborator approval.
A reusable helper does not certify all of its callers.

| Area | Classification and source evidence | Remaining work |
| --- | --- | --- |
| Functional application binding | **Adapt.** `AppSpec` in [RunNamedStages.hpp](src/apps/madqc/RunNamedStages.hpp) already supplies settings, planning, observation, archive membership, publication and property callbacks. [Drivers.hpp](src/madness/chem/Drivers.hpp) already has consume/run/results/publish handoff. | Extract useful operations into the retained application lifecycle. Remove `Stage`/named-input dependencies; do not reproduce the whole registry in a new abstraction. |
| Workflow reference handoff | **Adapt.** [Applications.hpp](src/madness/chem/Applications.hpp) carries SCF and Nemo references, archive paths and model checks. Its old result-file reader can catch corruption and choose `Redo`. | Add the shared identity/completion checks and strict selected-checkpoint recovery. A successful result-file parse or a live reference pointer is insufficient evidence of final convergence. |
| Storage publication | **Reuse mechanism; adapt manifest and callers.** [StageCheckpoint.h](src/madness/chem/StageCheckpoint.h), `publish` and `read_published`, validate shards and collectively publish/select one generation. | Preserve checked close, integrity validation and replacement order. Adapt records and internal calculation IDs; supply complete payload membership and same-state scientific evidence. Storage validation itself makes no convergence claim. |
| Cadence and retention | **Reuse policy helpers; adapt ownership.** [CheckpointRetention.h](src/madness/chem/CheckpointRetention.h) has endpoint/optional interval selection and latest-plus-referenced generation retention. | Apply them to each independent state, coupled bundle or direct checkpoint. Derive the referenced set correctly and avoid retaining every numerical protocol merely to preserve history. |
| Restart planning | **Adapt.** [LadderPlanner.h](src/madness/chem/LadderPlanner.h) already separates fresh, seed, load and resume and records earlier iterations without subtracting them from the current cap. | Use the one selected-criteria contract, explicit stage completion and current reference/operator identities. Connect supported conversion to strict restart and the [per-task entry rule](design_draft.md#preserve-each-tasks-protocol-progress), subsequently approved on 2026-11-11: preserve matching progress and start new/unmatched guesses at the first requested protocol. |
| Criteria and limits | **Adapt.** [ConvergenceSettings.h](src/madness/chem/ConvergenceSettings.h) shares settings, outcome reporting and `rung_exit`; [ExplicitCriteria.h](src/madness/chem/ExplicitCriteria.h) validates finite scalar evidence and units. | Remove the engine-versus-requested decision split. Resolve selected defaults and common absolute/relative tests. Apply documented heuristic adjustments only to intermediate absolute tolerances. Defer stall exits. |
| Metric operations | **Reuse and extend.** [FunctionResidual.h](src/madness/chem/FunctionResidual.h) computes raw real-3D function norms; [FunctionMetricSummary.h](src/madness/chem/FunctionMetricSummary.h) provides maximum/trace reductions and availability; [KineticEnergy.h](src/madness/chem/KineticEnergy.h) provides gradient contractions. | Build the thin typed evaluation record around these operations. Add required model weighting, scales and non-real/non-3D support during the corresponding engine adoption. Do not duplicate the existing arithmetic. |
| Density differences | **Adapt.** [DensityChanges.h](src/madness/chem/DensityChanges.h) preserves inputs but currently returns charge/spin field-difference norms. | Support the agreed application-supplied fields and reductions, including SCF's maximum spin-channel change. The current helper's existence does not establish the target density convention. |
| Properties and timings | **Adapt existing writers; add common records.** [Results.h](src/madness/chem/Results.h), [ResponseProperties.hpp](src/apps/madqc/ResponseProperties.hpp), and response `StateMetrics` already collect useful values and elapsed times. | Implement grouped requested properties and compact state metrics once per calculation in `results.json` (author revision on 2026-11-11), with separate full runtime/analysis records and hierarchical `timings.json`. Existing total or state/protocol times do not supply iteration-step timings; missing timing must become unavailable rather than zero. |
| Superseded paths | **Retire after migration and approval.** Named parser/runner, duplicate policy decisions, active root-freezing/warmup/stall paths, and later forced-coarse-work behavior conflict with selected requirements. | Transfer numerical operations and useful tests first. Preserve the existing workflow CLI and algorithms. No code is retired by this audit. |

## Representative tasks

### SCF and Nemo

SCF already has a checkpoint observer, physical energy/kinetic calculations,
raw per-function residual support and archive serialization. Its
[solve/update boundary](src/madness/chem/SCF.cc) remains the main gap:
`solve` tests the stored BSH residual before `update_subspace` computes the next
one, and an observer can receive the previous residual alongside updated
orbitals. Separate equation evaluation from KAIN/update and publish the evaluated
state. Retain existing rotations, restrictions and normalization.

Nemo already has a workflow application and uses SCF orbital save/load machinery;
its integration should not start with a new workflow wrapper or archive format.
However, [Nemo::value and the solve loop](src/madness/chem/nemo.cc) reject some
explicit convergence controls. The loop computes an aggregate residual divided
by the square root of the function count, applies the solver update, then checks
density and convergence. This differs from the agreed maximum-per-function
residual and common sampling point. Adapt this loop and its observer/evidence
binding, and obtain scientific review of regularized versus physical density
and kinetic quantities.

**Acceptance evidence needed:** the residual, density, energy and kinetic
records describe the saved functions; stopping takes no extra KAIN update;
Nemo fresh/resumed/varying-k cases use the same policy with reviewed physical
fields. Existing SCF and Nemo QC decks remain the scientific baseline.

### Driven response

[FDSolver](src/apps/molresponse/solvers/fd_solver.hpp) already exposes a raw
fixed-point operation and `evaluate_saved_state`, which assesses output without
taking another solver step. [iterate_protocol](src/apps/molresponse/solvers/iterate_protocol.hpp)
already accepts preparation and endpoint callbacks. State save/load and
per-frequency planning already exist.

Adapt these operations into evaluate/decide/advance. The current `step` measures
input density/residual, applies KAIN, and measures output kinetics. Saved-state
reevaluation is conditional on the requested-criteria path; it is not yet the
uniform iteration boundary. Reuse the equation work when continuing, preserve
actual spin/X/Y blocks, and use the common max-per-function reduction and limits.
Use existing serialization with the shared state identity and publication rules.

**Acceptance evidence needed:** one unfinished frequency resumes while other
valid states remain reusable; each early exit and endpoint has current-state
metrics; amplitude scaling and final criteria use the agreed conventions.

### Excited-state bundles

[ESSolver](src/apps/molresponse/solvers/es_solver.hpp) already owns the coupled
bundle, root rotations, BSH/KAIN work and both retained execution variants.
[es_root_identity.hpp](src/apps/molresponse/solvers/es_root_identity.hpp) carries
labels through discrete permutations and save/load. These labels do not establish
physical correspondence through continuous root mixing.

[es_kinetic_evidence.hpp](src/apps/molresponse/solvers/es_kinetic_evidence.hpp)
explicitly marks kinetic changes unavailable because correspondence is
unvalidated. Add the agreed matching/phase/comparison checks and use them for
density, energy and kinetic changes. Align sampling across rotation and update
in both solver variants. The baseline's engine check uses density and omega
changes; it does not enforce the agreed selected BSH criterion. The executor
also enables root locking and has separate warmup paths. Adapt those choices to
the all-roots-active, no-warmup, common-criteria design.

**Acceptance evidence needed:** mixing/reordering cannot manufacture a small
change; the full bundle meets its criteria at the same evaluated state; extra
roots reuse usable saved roots and receive independent guesses without restoring
KAIN. Existing root IDs and bundle serializers are retained.

### VBC and direct properties

The [calculation manager](src/apps/molresponse/calc/calc_manager.hpp) already
owns prerequisites and scheduling; the [executor](src/apps/molresponse/calc/calc_executor.hpp)
builds VBC from loaded B/C states. [vbc_save_load.hpp](src/apps/molresponse/solvers/vbc_save_load.hpp)
already writes and loads VBC archives. `ResponseCheckpoint::members` in
[ResponseCheckpoint.hpp](src/apps/madqc/ResponseCheckpoint.hpp) includes those
archives in shared generation publication. There is no need for another VBC
storage subsystem.

Adapt readiness and saved records. Existing VBC metadata calls successful
construction `converged` and records zero iterations. Use direct completion
with exact B/C/reference identities and operator settings. Current prerequisite
checks can use a protocol-local convergence flag or a planner's `met` result;
the new boundary must require final requested convergence. Final properties
also require eligible A and VBC. Keep dependency knowledge in molresponse and
share the evidence/publication interface with SCF and response states.

**Acceptance evidence needed:** unchanged eligible VBC is reused; changed B, C,
reference or relevant operator settings invalidate it; a coarse or cap-exhausted
input cannot authorize a final property. Direct work has build timings and no
fabricated solver convergence criteria.

## Implementation plan changes justified by the audit

1. Retain the existing drivers, DAG, numerical serializers and storage publication
   mechanism. Extract the needed lifecycle operations from the named runner.
2. Make the new evaluation/timing record a small shared interface around existing
   metric operations. Keep state and temporary BSH work in application types.
3. Treat SCF/Nemo and response integration as concrete phase/metric changes, not
   merely adding a callback to otherwise complete solvers.
4. Reuse Nemo's existing workflow binding and VBC's existing archives. Add only
   their missing shared evidence, compatibility and publication integration.
5. Reuse current test assertions that still express the approved contract;
   replace assertions for superseded behavior. Change QC input/runner files only
   where their options or invocation require it, and update output assertions
   where the new records require it. Preserve numerical reference data.
6. Keep release-wide adoption for the approved engine scope. This representative
   audit establishes reusable boundaries, not adoption evidence for the other
   engines. On 2026-11-11 the author scoped the initial release to E01–E12;
   E13–E15 retain separate builds and their adoption is deferred.

## Validation evidence and remaining decisions

The baseline contains checkpoint corruption/missing-shard and failed-publication
assertions in [test_stage_checkpoint.cc](src/madness/chem/test_stage_checkpoint.cc),
plus [retention](src/madness/chem/test_checkpoint_retention.cc),
[ladder planning](src/madness/chem/test_plan_ladder.cc),
[response scheduling](src/apps/molresponse/tests/test_calc_manager.cpp),
metric and named-path integration tests. They are migration sources. Some
assertions intentionally describe old decision modes, so existing coverage must
be read and adapted rather than counted as proof of the new behavior.

No numerical defaults, root-match thresholds or scientific reference values
were validated in this audit. Remaining author decisions are already listed in
the draft: exact output record keys/units, option mappings that cannot preserve
intent, and unresolved model quantities for scientific-owner review. The
conversion-entry and research-program scope choices, open at the time of the
audit, were settled on 2026-11-11 in the draft. E13–E15 are deferred from the
initial release with their separate builds preserved. These later decisions do
not change the baseline source findings.

Source links point to the current checkout. Baseline claims were checked against
`git show 9179f5d32076a2ba555906dd146dc330ffc040d7:<path>`; later fixed-stage
amendments were inspected separately. The private audit record retains the
baseline source hashes and document validation results.
