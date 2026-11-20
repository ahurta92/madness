# Design review walkthrough

**Target:** [design_draft.md](design_draft.md)\
**Date:** 2026-11-09\
**Scope:** First read-through by the drafting assistant. The independent EGM
review stages have not started; the review-mode choice is pending. These notes
check how the written rules apply to example calculations. They do not verify
implementation or change the design.

## Current author-approved direction — 2026-11-12

The latest workflow-only revision below supersedes historical decisions about
named application stages, their executor and CLI targeting. The current target
keeps workflow selectors, existing input groups and command-line options;
workflow code manages cross-application dependencies. Shared convergence,
metrics, checkpointing and molresponse's internal DAG remain required. Older
entries are review history, not additional current requirements.

The approved intermediate-tolerance revision below also supersedes earlier
fixed-work decisions: use legacy threshold-dependent tolerances at intermediate
protocols, with `maxiter` as a cap. Shared metric definitions, selected criteria
and final convergence requirements remain in force. Stall-based stopping is
deferred.

Application stages are workflow-created calculations with requested properties
and resolved settings. Their internal tasks remain application-owned. The
approved record distinction on 2026-11-11 supersedes the earlier proposal for
one undifferentiated task record; it does not restore named-stage input.

Applications own numerical state, comparison history and equation/KAIN work.
They use shared numerical helpers and pass measurements to the shared controller,
which handles tolerance decisions and recording. The 2026-11-11 ownership
revision below supersedes wording that put state/history in a common evaluator.
Direct submission and callbacks use one measurement and decision implementation.

The evaluation return structure is also approved: one prepared sample and a
collection of subject measurements, including comparison identities, values or
availability reasons, relative scales and supporting scalar evidence. The
controller reply is also approved. Checkpoint metrics contain measured evidence;
current requests generate new checks. A compatible final checkpoint that already
satisfies tighter criteria is reused without another numerical evaluation or
archive write. This supersedes the earlier blanket fresh-check rule for changed
stopping tolerances. Exact persisted metric encodings remain proposed.

Protocol entry is now settled: preserve each task's compatible saved position;
new independent states and supported removed-protocol guesses enter the first
requested protocol. Response waves include only tasks still needing that level.
Excited roots retain one bundle entry, including the existing expansion rule.
This supersedes the earlier open choice about conversion entry.

The 2026-11-17 2PA refinement exception preserves saved functions and a requested
saved protocol for the same tracked root and perturbation at an updated photon
frequency. It requires fresh convergence evidence, empty KAIN and a fresh
allowance; the old equation's completion claims are not inherited. Genuinely
new/unmatched states and removed protocols keep their existing entry rules.

The initial release scope is now E01–E12. Periodic DFT `testperiodicdft`,
`tdse` and `tdse4` (E13–E15) remain standalone programs with separate builds;
their shared-layer migration is deferred. This supersedes historical wording
that left their applicability open or required adoption by every built engine.

Results now include compact measured state metrics once alongside grouped
properties in each application calculation. Alpha and beta retain one tensor
record per frequency/process context; states shared by properties appear once.
This supersedes the historical properties-only rule. Full histories, tolerance
assessments, dependency graphs and restart authority stay in supporting records.

Explicit convergence names are now required in the adopted shared interface.
`dconv`/`econv` are rejected with migration guidance, including mixed old/new
inputs and adopted CLI forms. Workflow selectors, group overrides and supported
key precedence remain. This resolves the tolerance-alias choice; distinct
inner-solver controls and other legacy mappings still need scoped review.

Read-only is now reserved for inspecting results and plotting saved functions,
including unfinished calculations, without new convergence tests or solver
iterations. Solver restarts use reuse-or-iterate. This supersedes the proposed
load-only solver refusal when a compatible checkpoint misses current tolerances.
"Requested final protocol" means final truncation threshold and polynomial order.

The response budget input is settled: require response `maxiter` and reject
explicit `excited.maxiter` with migration guidance. One configured cap applies
independently to each driven state or coupled excited bundle at each executed
protocol. This supersedes the proposed alias and the old excited-only override.

The `save` and `checkpoint.keep_all` inputs are removed entirely, including
values that previously matched the defaults. They receive only the normal
fatal unknown-input error, without migration messages or deprecated-key handling.
Endpoint saves, limited retention and optional checkpoint intervals remain.

## Start with the calculation behavior

The main unresolved behavior is what happens when a run fails or its requested
states change. The ordinary iteration and restart cases are already stated
clearly enough to walk through:

| Example | What the draft requires |
| --- | --- |
| Coarse stage allows up to 10 updates; all required intermediate targets pass after 2 | Finish this stage after 2 and still execute the requested final stage |
| A checkpoint records 1 iteration of an interrupted coarse stage | Resume that saved stage with a new allowance of 10; the old iteration stays in the history |
| A coarse stage reaches its update cap with targets unmet | Advance to the next protocol; record the cap exit and unmet targets without claiming convergence |
| Intermediate metrics plateau while targets remain unmet and updates remain | Continue; stall-based stopping is deferred |
| Final stage reaches its limit with criteria unmet | Report unconverged and retain usable work for another run |
| Final criteria pass on the last allowed iteration | Report converged |
| An unchanged final result is converged and its measurements remain valid | Load it without new solver iterations |
| One frequency state completed the current stage; another did not | Keep the completed work and resume the unfinished state before advancing |
| One root in a coupled excited-state bundle passes early | Keep updating it with the other roots |
| The user requests additional excited roots | Reuse usable roots, add independent guesses and iterate the enlarged bundle at the saved stage |
| Restart is disabled | Initialize fresh and continue writing checkpoints |

These examples follow the [iteration rules](design_draft.md#2-iteration-and-convergence)
and [restart table](design_draft.md#4-restart). The interrupted-stage example
assumes a checkpoint was actually written: unsaved iterations cannot be restored.

## Clarifications worth making

### 1. Distinguish a usable unconverged checkpoint from a diverged state

**Relevant sections:** [Iteration](design_draft.md#2-iteration-and-convergence),
[Restart](design_draft.md#4-restart), and
[Checkpoints](design_draft.md#5-checkpoints-and-retained-history).

The draft says failures stop a task, unconverged final work resumes, and controlled
stops write checkpoints. It does not explicitly say how those rules interact
when the latest numerical state has diverged.

**Example:** a valid checkpoint exists at iteration 5, but iteration 6 produces
an unusable state. Marking iteration 6 current and deleting iteration 5 would
violate the promise to retain the latest usable checkpoint. Treating every
unconverged state as resumable could also restart from the failed state.

**Clarification to consider:** retain the prior usable checkpoint on numerical
failure, record the failure separately, and limit the ordinary resume rule to
valid unconverged states. State whether a later run may resume the earlier
checkpoint or must use fresh initialization. This is a missing recovery rule,
not a reason to add a new convergence heuristic.

### 2. State where a newly added frequency starts

**Relevant section:** [Restart](design_draft.md#4-restart).

The draft specifies the saved stage for an enlarged excited-state bundle. For a
new independent frequency, it specifies a compatible guess but not its starting
protocol.

**Example:** existing frequencies have reached the final protocol and the user
adds another frequency. That state could start at the first requested protocol,
or use a compatible final-protocol guess and begin there. These choices change
cost and which coarse stages are executed.

**Clarification to consider:** add the starting-protocol rule for new independent
states, including what happens when no compatible guess exists. Existing valid
states should retain their completed work under either choice.

### 3. Define which work must finish before response advances — clarified

**Relevant sections:** [Iteration](design_draft.md#2-iteration-and-convergence)
and [State dependencies](design_draft.md#3-ground-state-response-states-and-properties).

**Author clarification — 2026-11-09:** VBC and property evaluation do not
iterate. The existing DAG builds VBC after B and C converge, then evaluates
⟨A;B,C⟩ after A converges and VBC is complete. They do not hold up advancement
of iterative states. This is now explicit in the draft. The original question
below is retained as review history.

For independent frequency states, the rule is clear: finish their configured
coarse work before advancing together. The draft also preserves dependencies
for quadratic response, Raman and two-photon calculations, without specifying
how those dependent equations participate in that stage boundary.

**Example:** a primary frequency state completes its coarse iteration count
without meeting its final criteria. A secondary equation needs that state as a
source. Does it run at this coarse stage, wait for a later stage, or block stage
advancement? A universal requirement for final convergence before any dependent
work could prevent completion of the coarse stage.

**Clarification to consider:** distinguish the ground-state prerequisite from
requirements for intermediate response sources, and name which solver tasks
participate in each shared stage boundary. The full scientific dependency graph
need not be repeated in the design; its rule for advancing must be unambiguous.

## Driver decision and remaining choices

**Author decision — 2026-11-09:** reuse the existing shared named-stage and
response iteration drivers. SCF keeps its numerical loop for the initial
implementation while adopting the common convergence, metrics and checkpoint
interfaces. Moving SCF into a common iteration driver is deferred; see
[the selected approach](design_draft.md#7-reuse-the-existing-drivers).

**Integration order — 2026-11-09:** molresponse and Nemo first, followed by MP2,
CC2 and the remaining chemistry engines. Nemo needs its own adapter for its
existing solve loop; retain the SCF reference support needed by molresponse.
See [integration order](design_draft.md#8-integration-order).

**Metric standard — 2026-11-09:** the derivative-sensitive quantity is kinetic
energy from orbital derivatives. Every iterative solver reports BSH residuals
and applicable energy/frequency, density and kinetic-energy changes. Applications
supply their functions and model-specific quantities to a shared evaluator that
enforces agreed norms, combination rules and comparison timing. Every convergence
decision must describe one evaluated state that matches the checkpoint.

**Measurement point — 2026-11-09:** evaluate the prepared current state before
the next update. Reuse its raw BSH proposal/residual for KAIN when continuing;
save that evaluated state when converged. SCF needs a split inside
`update_subspace`; driven response needs a pre-KAIN decision in the shared
driver. Excited states also need density, kinetic and frequency samples aligned
with the rotated bundle. The final produced state must be evaluated without an
extra update, with any terminal equation work recorded separately.

**SCF split approved — 2026-11-09:** the previous residual does not establish
convergence of the newer orbitals being saved. Separate evaluation and advancement
so the convergence decision and checkpoint describe the same state. Use common
prepare, evaluate, decide and advance phases across function-based solvers;
applications retain their equations and numerical operations. Evaluation must not
advance KAIN history. SCF can retain its own loop while adopting these phases.

**Residual reduction — 2026-11-09:** use the maximum per-function L2 norm of
`x - g(x)` over each state's actual blocks. Retain component and per-root/state
records; KAIN receives residual functions. Relative scaling and the treatment of
regularized functions still need review, and changing reduction requires a
review of tolerances. The loop carries an application-defined state and invokes
operations on it. The existing templated `Solver::State` pattern is a concrete
implementation approach; final C++ signatures remain open.

**Lifecycle and derived checkpoints approved — 2026-11-09:** refine the existing
outer lifecycle for settings, restart, protocol preparation, execution and
publication. Iterative tasks use the shared phases; VBC and properties run direct
operations with the same checkpoint interface. VBC records completion and its
B/C, reference and settings dependencies. MP2 and CC2 can adopt the phases
gradually while retaining coupled algorithms; PNO matrix residuals need their
own definition. The 3D residual helper needs extension for 6D functions. These
are implementation requirements, not claims of completed code changes.

**Density reporting approved — 2026-11-09:** keep one primary density-change metric per
state/root, with each model supplying its density construction operation. The
shared evaluator compares fields and applies a declared reduction; it need not
require separate charge/spin metrics. Explain the model quantity, factors, units
and comparison to the user. The interface must allow multiple fields to retain
SCF's spin-density comparison. Exact adapter definitions, SCF's target ambiguity
and excited-root comparison eligibility remain work to settle.

**SCF kinetic energy approved — 2026-11-09:** use the physical kinetic-energy
term from the total energy, including occupations and both spin contributions.
Report its value and absolute change between comparable prepared states. Migrate
the existing unweighted explicit kinetic-trace criterion and identify the changed
metric meaning. The response decision is recorded below.

**Response kinetic sum approved — 2026-11-09:** retain the positive sum of
per-function kinetic contributions over each state's declared stored blocks,
with unit weights and no extra occupation/spin factors or amplitude normalization.
Driven amplitudes keep their equation's scale; excited roots retain the solver's
TDA or signed full-response normalization. Evaluate before the next update and
report the value and absolute change between comparable states. Document the
static/full storage difference and amplitude-dependent units. Root correspondence
is still required for excited-state changes; relative tolerances remain undecided.

**Root comparisons approved — 2026-11-09:** retain the existing root identities,
ordering, sign fixes and near-degenerate alignment. Validate correspondence
between consecutive prepared roots using normalized positive overlaps and use
one accepted match for density, frequency and kinetic changes. Permutations carry
histories; density needs sign alignment. Ambiguous matches reset the comparison
sequence and leave affected changes unavailable. Required unavailable changes
prevent final convergence within the unchanged iteration allowance. Thresholds
and ambiguity margins need validation; subspace-based stopping is deferred.

**Criteria selection and initial defaults approved — 2026-11-09:** users choose
required metrics and tolerances for each named application stage. Without an
explicit selection, ground states require BSH residual, density and energy
changes; driven response requires BSH residual and density change; excited states
require BSH residual, density and frequency changes for every required root.
Kinetic change is reported and optional for stopping. These criteria determine
final convergence. The original fixed-work rule for intermediate MRA protocols
is superseded by the intermediate-tolerance revision below. Numerical tolerance
defaults remain to be specified; the syntax decision is recorded below.

**Shared input syntax approved — 2026-11-09:** use `convergence.criteria` and
`convergence.<metric>` tolerance fields inside each named application stage.
Tolerance overrides do not change the active set. All driven states and excited
roots within a molresponse stage share one parameter set; no separate driven or
excited groups. Older keyword compatibility remains a migration decision. The
energy naming and applicability decision below resolves mixed-stage behavior.

**Energy-change naming and applicability approved — 2026-11-09:** use one
`energy_change` metric: total energy for the ground state and excitation energy
for excited roots. Within molresponse it applies only to excited states; driven
states report it as not applicable. One shared response default list includes
BSH residual, density and energy changes, with each state requiring the applicable
subset. No separate frequency-change key or driven/excited parameter groups.
Every iterative state must retain at least one applicable selected criterion;
missing applicable measurements remain unavailable and cannot pass.

**Optional relative tolerances approved — 2026-11-09:** keep
`convergence.<metric>` as the absolute tolerance and add
`convergence.<metric>.relative`, defaulting to zero. Combine them as an absolute
limit plus the relative tolerance times the declared scale. Report the raw
measurement, both tolerances, scale and effective limit. Molresponse shares the
same settings across driven states and excited roots. Missing scales cannot
silently disable a requested relative test. The scale definitions are recorded
below; numerical absolute-tolerance defaults still need validation.

**Relative scales approved — 2026-11-09:** BSH uses the largest current prepared
function norm. Density uses the larger current/previous density size, applying
the model's field reduction to each sample. Energy and kinetic change use the
larger magnitude of their two comparable quantities. Compute scales per state
or root, with the agreed model definitions and root correspondence; never pool
states or rescale their functions. Use no division or artificial minimum scale.
The BSH rule does not impose equal fractional accuracy on each small component.
Total SCF energy scaling is optional and defaults to zero, since large total
energies can allow loose convergence of small energy differences. These limits
measure iteration convergence, not guaranteed final property accuracy.

**Absolute defaults deferred to implementation — 2026-11-09:** select and
validate absolute-tolerance defaults on representative SCF and molresponse
calculations using the agreed definitions, then Nemo once its model definitions
are settled. Compare candidate defaults with tighter calculations and record
the evidence, results and iteration costs. The implementation plan includes this
validation task; choosing numerical values does not block writing the plan.
The default criteria and zero relative tolerances remain selected.

**Checkpoint frequency approved — 2026-11-09:** default to saves at protocol
stage endpoints and controlled stops with a valid evaluated state. Archive
creation can be costly in memory and computation time. Saving every iteration
or every N iterations is optional; endpoint saves remain required. Independent
states save separately, coupled excited roots save as a bundle, and direct
tasks save on completion. Keep the latest usable checkpoint and continue
recording iteration metrics and timings independently. An unplanned interruption
can require repeating work since the last saved state, including an entire
stage under the default policy.

**Warmup disabled; separate future allowance — 2026-11-09:** disable separate
warmup iterations in the simplified implementation for fresh and restarted
solves. Preserve supported starting guesses and required preparation.
`maxiter` controls main solver iterations per protocol stage and run. Any future
warmup or iterative guess refinement must use its own allowance, counts and
timings, without consuming `maxiter`. Names such as `warmup.maxiter` or
`guess.maxiter` are illustrative future options, not selected input syntax.
Record warmup as disabled in new runs and preserve earlier warmup history
separately. The existing TDA warmup calls must be bypassed during implementation.

**Transition without backward compatibility — 2026-11-09:** preserve current
input and checkpoint paths until each application adopts the shared layer.
The new contract publishes numerical archives with metric records describing
those saved functions; direct tasks retain completion and dependency records.
Maintaining old input, checkpoint or output formats is out of scope. Existing
serialization primitives can be reused, but legacy aliases, archive importers
and schema-conversion support are not required. Keep users' old files intact
and report unsupported formats clearly. This replaces the compatibility proposal
in the inspection discussion below.

**All engines required for release — 2026-11-09:** the formal chemistry release
requires every chemistry engine to adopt the shared I/O, metrics, applicable
iteration/convergence rules and restart/checkpoint contract. Implement
molresponse and Nemo first, with required SCF support, then MP2, CC2 and the
remaining engines. Existing paths remain available during development; partial
adoption does not satisfy the release requirement. Preserve application-specific
algorithms and coupled loops. The implementation plan must enumerate all engines
and their required scientific and shared-layer tests, tracking each to completion.

**Next target:** freeze conventions between SCF and molresponse first, reviewing
unresolved scientific definitions with Florian and listing missing evaluator,
adapter or loop work. Apply those conventions to Nemo next. This specifies the
order of the convention review; the integration order above remains unchanged.
See [the functional interface and current gaps](design_draft.md#7-reuse-the-existing-drivers).

The remaining metric conventions and concrete formats remain open. Release
scope is settled; the complete engine inventory and coverage matrix belong in
the implementation plan.
Absolute defaults are assigned to implementation validation as described above.
These are planned work, not accidental omissions from a finished implementation
plan.

## Input and checkpoint compatibility discussion

**Inspection — 2026-11-09; author scope clarification incorporated.** These
findings describe this checkout, not a promise to retain its format support in
the new interface. No numerical calculations or tests were run.

| Area | What the code supports now |
| --- | --- |
| Older workflow inputs | `madqc --wf=...` uses application parameter groups and command-line overrides. This path includes Nemo and broader science than named stages. It writes `calc_info.json`; named stages write `results.json`. See [dispatch](src/apps/madqc/madqc.cpp). |
| Named inputs | Named `molecule`, `moldft` and `molresponse` blocks, an explicit parameter allowlist, and `--continue` to reopen an output root. Named moldft currently allows HF energy/dipole; named response requires a closed-shell reference. There is no named Nemo adapter yet. The proposed `convergence.criteria` and per-metric absolute/relative syntax are not implemented. See [parser](src/apps/madqc/NamedStages.hpp) and [execution](src/apps/madqc/RunNamedStages.hpp). |
| SCF and Nemo archives | Both use `restartdata`; the shared reader accepts versions 4 and 5, and current writes use version 5. Nemo marks regularized functions and its nuclear correlation factor; known SCF/Nemo representation mismatches are rejected. Version 4 lacks that representation field. SCF can adjust thresholds and project loaded functions to a different polynomial order, invalidating affected convergence evidence. Nemo uses the restartdata-only planner, not the SCF AO/NWChem fallback readers. See [header](src/madness/chem/Restart.h), [SCF loader](src/madness/chem/SCF.cc), and [Nemo restart](src/madness/chem/nemo.cc). |
| Shared named checkpoints | Version-1 manifests inventory archive shards, sizes, checksums, settings and evidence. App adapters check scientific consistency. State-record readers upgrade known v1 records to `madqc.state.v2` in memory without rewriting the old generation; unknown schemas are rejected. See [storage](src/madness/chem/StageCheckpoint.h) and [record upgrades](src/apps/madqc/NamedStages.hpp). |
| Response payloads | Driven states, excited bundles and VBC have dedicated readers with metadata; response supports native and optional HDF5 I/O. Named published checkpoints currently validate native payloads even when HDF5 mirrors are used. Existing readers tolerate some older missing fields, so read success alone is not evidence of convergence under the proposed rules. See [driven reader](src/apps/molresponse/solvers/fd_save_load.hpp), [bundle reader](src/apps/molresponse/solvers/es_save_load.hpp), [VBC reader](src/apps/molresponse/solvers/vbc_save_load.hpp), and [named adapter](src/apps/madqc/ResponseCheckpoint.hpp). |

The author selected a transition to the shared contract without maintaining
backward compatibility. Preserve existing paths until their applications adopt
that contract. After adoption, the supported checkpoint contains the numerical
archives and matching metric records, published and validated together. Existing
serialization code may be reused; old formats do not require continuing support.
Do not add legacy keyword mappings, automatic archive imports or schema converters
as compatibility deliverables. Keep existing files intact and explain unsupported
formats. The proposed preserve/convert/import policy is superseded by this scope.

The inspection still identifies implementation work. SCF's old `dconv` paths use
different density comparisons, and its old `accuracy.kinetic_trace_change` is an
unweighted trace rather than the physical kinetic-energy change selected here.
Implement the agreed definitions in the new interface; do not silently give old
measurements a new meaning. See [explicit targets](src/madness/chem/ExplicitCriteria.h).

The new shared metric records must retain missing measurements as unavailable;
the current driven loader's zero fallback for missing historical BSH diagnostics
must not become a passing criterion. VBC must write the B/C checkpoint identities
and completion evidence required by the new contract. Supporting imports of old
files missing those fields is not part of this implementation.

## A tradeoff to recognize, rather than a defect

Changing a tolerance alone currently requires fresh work even if saved
measurements appear to satisfy the new value. This is consistent with the
chosen conservative restart rule. It can cost extra iterations, but avoids
making saved-value reuse an additional policy in this first design. Changing
that behavior would be a new decision, not an editorial correction.

## Implementation drafting — 2026-11-09

At the author's request, the [implementation plan](design_draft.md#detailed-implementation-plan)
now starts with a source-backed inventory of all built chemistry solvers and
separately identifies direct tasks, workflow drivers and dormant sources. The
plan records the two distinct MP2 implementations, methods missing QC coverage,
and the applicability question for the periodic and TDSE research programs.

Eleven implementation steps list shared interfaces, existing/new files, engine
adoption and required tests. The approved molresponse/Nemo-first order and
all-engine release requirement remain in force. Existing numerical algorithms
and scientific coverage are retained; no backwards-compatibility work is added.

**Awaiting review:** the engine inventory and implementation sequence. Result
layout, recovery/removed-stage behavior, research-program release applicability
and unresolved model definitions are explicit decision gates. Numerical defaults
and root thresholds remain implementation validation work. The plan is not yet
ready for execution without those decisions. This update does not constitute a
new Goldfish review or numerical validation; only documentation was changed.

## Shared execution revision — 2026-11-10

The author requested revisions: named stages exist so users can connect
applications and dependencies in input, rather than select only predefined
workflows. If both input forms remain, they must share the technical execution
implementation. Matching results or common policy helpers alone are insufficient.
At this point preserving workflow shortcuts was optional; the later CLI
preservation revision below makes their continued support a requirement.

Source inspection found separate outer execution paths: `madqc.cpp` dispatches
to `run_if_named` or builds/runs `qcapp::Workflow`; `RunNamedStages.hpp` and
`Applications.hpp` contain different restart and publication paths despite
reusing numerical engines. The [revised design](design_draft.md#1-shared-responsibilities)
and [implementation step 5](design_draft.md#5-prepare-workflow-composition-and-the-shared-calculation-lifecycle)
now require one stage model, validator, application registry and executor.
Optional shortcuts construct requests for that executor. The plan requires
retiring duplicate execution code and testing the shared call path, including
restart through either input form when stage identities and settings agree.

This clarification preserves SCF's separate numerical loop, the existing
response iteration drivers and molresponse's internal dependency manager.
**Author approval — 2026-11-10:** the shared-execution revision is approved.
The subsequent decision below settles shortcut retirement. The engine inventory
and broader implementation plan remain under review; this approval does not
settle their other open decisions. No production
code was changed or numerical tests run.

## Workflow retirement decision — 2026-11-10

**Superseded for selector removal by the CLI preservation revision below.** The
collaborator-approval condition remains in force for obsolete executor cleanup.

**Author-approved target:** move to named-stage input and retire workflow
shortcuts after migration, with the amendment that the other project
collaborators must explicitly approve removal. Both conditions are required.
Passing tests and approval of this design are not substitutes for that later
approval. Record the approving collaborators and reviewed scope before removal;
no such approval has been claimed or requested from them in this conversation.

During the transition, migrated workflow input constructs stage requests for
the same executor as named input. The plan separates migration completion from
retirement completion, uses the existing entry point for temporary request
expansion, and removes the proposed permanent `WorkflowInput.hpp` deliverable.
Maintain named-stage examples for common calculations. This decision changes
neither the scientific engine inventory nor the agreed convergence/restart
conventions. No production code was changed and no workflow was retired.

**Next decision:** the common result-file layout. Other implementation choices
and full-plan review remain open.

## Requested results and supporting records — 2026-11-10

The author selected `results.json` for requested properties for each named
stage. Analysis-only data and runtime metadata belong separately. The same
scientific property type should use a common format across applications.

Section 6 now describes those roles and the common property-format requirement.
The implementation plan separates the current mixed output in
`RunNamedStages.hpp`, moves runtime reads away from property files and removes
embedded response solver evidence from requested-property records. Values keep
the units, scientific identity and validity needed for interpretation; detailed
convergence and restart evidence stays with its supporting records. A requested
quantity belongs in results even when an analysis routine computes it.
Checkpoint archives still carry matching measurements; this output separation
does not make analysis history the authority for restart.

The timing decision is recorded below. Exact property keys, numeric encodings,
units and physical supporting-file layout remain to be specified. No exact
schema or full-plan approval is inferred from this discussion. Only
documentation was changed; no numerical tests were run.

## Timing detail — 2026-11-10

The author requested a simple, understandable report showing the steps of each
named stage and the steps within each state's iteration. Section 6 now specifies
run/attempt → named stage → protocol → state or coupled bundle → iteration →
step, with the shared prepare/evaluate/decide/advance phases and measurable
application substeps. Terminal evaluation and direct VBC/property work retain
their actual meaning rather than acquiring artificial solver updates.

Independent-state overlap and shared excited-root work remain explicit: parent
elapsed times are measured, not obtained by adding nested or parallel timings;
per-root costs are only reported when directly measured. The code inspection
found optional excited-state phase timing that adds fences. Default timing must
use existing completion boundaries or report a combined interval; synchronized
fine-grained profiling remains explicitly diagnostic.

**Author-approved — 2026-11-10:** use a separate `timings.json` for the complete
breakdown, with readable timing tables generated from the same records. Keep
`results.json` focused on requested properties; do not copy timing summaries
into it. Exact record keys and directory placement remain to be specified.
The plan lists shared collector/record hooks and checks using a controlled
clock; no production code, timing file or benchmark was created.

**Next design item:** common property and timing schemas and the supporting
file layout. Other unresolved decisions and full-plan review remain open.

## Property-format discussion — 2026-11-10

At the author's request, section 6 now contains illustrative total-energy and
polarizability records for discussion. The examples retain the selected rule
that the same scientific property has one format across applications, while
keeping total, correlation and excitation energies distinct. They omit runtime
solver evidence and show only property data plus essential interpretation.

**Author-approved — 2026-11-10:** named components that align directly with
requested properties. Section 6 now retains the named polarizability example
and removes the indexed-array alternative. Requests and results share component
labels, with property, frequency and root context as applicable. Unrequested
components are absent; unavailable requested components retain explicit status
and a reason. A common reader can assemble an array when needed.

**Author-approved — 2026-11-10:** keep successful components as named numbers
and record unavailable requested components in a separate availability map.
The draft adds a partial-result example. The shared writer checks that each
requested name appears exactly once across the two maps, with no unrequested
names. The implementation plan includes these checks in the property fixtures.

**Author-approved — 2026-11-10:** complex values use two-number arrays in
`[real, imaginary]` order. Real quantities remain numbers; quantities declared
complex retain both entries even when the imaginary part is zero. Section 6
records the convention, and the implementation plan includes value-shape checks.

The frequency/root grouping decision is recorded below. Exact keys, units and
the overall stage container remain open. Example numbers are not numerical
reference data. No production serializer
or numerical output was changed; no broader implementation-plan approval is
inferred.

## Frequency and root layout discussion — 2026-11-10

The current named output groups results by property. Polarizability rows carry
a frequency and request index; excited-state records carry stable root IDs.
Source inspection also found that the two-photon formatter associates a row
with a bundle root by position. The shared output migration must preserve the
corresponding stable identity. IDs remain local to their source bundle, not
universal labels across calculations.

Further inspection clarified that the existing output already combines outer
property grouping with related excited-state quantities grouped per root. The
run-level results file indexes separate stage files. The earlier alternative
with a standalone excitation-energy list was a proposal, not the current shape.

**Author delegated the choice — 2026-11-10:** retain outer grouping by requested
property, frequency records for driven properties, and related requested
quantities together under each excited root. This follows the requests and
preserves useful existing structure with less migration work. Shared formats
still apply to every quantity and root collection across applications. Related
collections, including two-photon properties, use the same bundle-local root
IDs. Readers can join them; no second authoritative layout is stored.

Section 6 replaces the unselected alternatives with that decision and its
tradeoff: cross-collection views require a reader join. Frequency identity,
missing requests, archive-key checks and supported response processes retain
their stated rules. Exact keys, units, the overall stage container and timing
schema remain open. No production code was changed or numerical tests run, and
this delegated choice does not approve the complete implementation plan.

## Command-line preservation revision — 2026-11-10

The author requested preservation of existing command-line options and approved
keeping `--wf` and its aliases long-term as shorthand for shared stage requests.
This supersedes the earlier selector-retirement target. Named-stage input and
CLI commands use one validator, application registry and execution lifecycle;
only the obsolete separate workflow execution branches are retired, after
migration and the other collaborators' explicit approval.

Code inspection confirmed `--wf`/`--workflow`/`--w`, default SCF, geometry and
parameter-group options, optimization and inspection/run controls. Typed text
and JSON parameter readers apply command-line overrides after file settings.
The current named-stage runner instead rejects scientific CLI options. Section
1 and implementation step 5 now require preserving these capabilities and
precedence, removing that blanket restriction, and testing the shared call path.
Exact syntax for targeting an ambiguous named stage remains open. Inventory
parameter mappings before adoption and review any that cannot preserve intent
under the new scientific conventions.

The overview, scope, migration checks, file inventory and release evidence now
agree on retained CLI support. General compatibility for old input-file,
checkpoint and output schemas remains out of scope; CLI preservation is an
explicit requirement. This revision does not authorize removing internal code
before collaborator approval, contacting collaborators, or starting production
implementation. Only the design and review notes changed.

## Named-stage CLI targeting proposal — 2026-11-10

At the author's request for the simplest form, the recommendation is one
`--stage=<name>` selector with the existing parameter-group overrides. It
selects the owning application stage for settings, not which dependency graph
to execute. Unknown targets, unsupported options and repeated selectors fail;
CLI overrides without a selector require one unambiguous target. Multiple
stage edits remain expressible in the input file. Existing workflow commands
retain their parameter-group targeting.

The current tokenizer can read the spelling, and `value_raw` preserves stage
name case. Its map overwrites repeated keys, so duplicate-selector detection
must precede that loss of information. Section 1 records these implementation
constraints and examples. The spelling and omission rule remain proposed,
not approved or implemented. No source code or numerical tests were changed.

## Drop named stages; retain workflows and shared contracts — 2026-11-10

**Author-approved:** remove named application stages entirely from the target
design, keep existing CLI options and parameter groups, and let workflow code
manage application composition and dependencies. The proposed `--stage` syntax
is withdrawn. This also supersedes the plan to route workflow commands through
a retained named-stage model or runner. There is no input-defined application
graph or hidden named-stage adapter in the revised target.

The implementation builds on workflow builders, Drivers/Application interfaces
and their existing StepContext handoff. Shared lifecycle operations absorb
useful work from the named path; the common eligibility checks require the
right converged reference before releasing dependent calculations. Molresponse
continues to own its internal DAG. This workflow-only revision left numerical
policy unchanged; its fixed coarse-work assumption is superseded by the later
intermediate-tolerance revision below. Metric conventions and checkpoint
guarantees remain in force.

Results and timing records are grouped by application calculation, with
workflow-assigned internal identities instead of user-supplied stage names.
Requested-property grouping, per-root collections, named tensor components,
availability maps, complex pairs and separate timing/analysis/runtime records
remain as agreed. The exact outer schema and directory layout remain open.

The plan now migrates scientific checks and useful validation from named-stage
files into the retained workflow/shared-contract paths before retiring the
parser/runner and grammar-only tests. Collaborator approval still applies to
cleanup; no contact or production removal is claimed. The full engine/method
inventory and all-engine release requirement remain. Only the design and review
notes were changed; no production code, builds or numerical tests were run.

## Run-directory layout and checkpoint recovery — 2026-11-10

**Author-approved layout:** one run directory with root-level `results.json`
and `timings.json`, plus `analysis/` and `checkpoints/`. Per-calculation records
use workflow-assigned IDs. Runtime/dependency metadata is stored under
`checkpoints/` alongside saved-state records; only published numerical states
are restartable. Exact supporting record filenames and schemas remain open.
The layout does not change endpoint-first saving or latest-usable retention.

The shared checkpoint reader selects the published generation only, with
manifest and file validation and collective error reporting. Existing test
source covers corruption/missing shards and failed replacement writes. No new
tests were run. A failed replacement that preserves the previous publication
is distinct from corruption of the selected published checkpoint.

**Author-approved — 2026-11-10:** use strict restart and defer explicit
alternate-checkpoint selection. An expected selected checkpoint that is missing,
corrupt or in an unsupported format stops the affected calculation and its
dependents, with a clear error and existing files preserved. There is no silent
fresh start or older-generation fallback. The user can explicitly start fresh
or restore their own backup; normal validation still applies.

A new state without an expected checkpoint can initialize normally, and valid
unconverged checkpoints follow the existing resume rules. Failed replacement
writes retain the prior usable publication. The approved conversion rule below
handles a valid checkpoint whose saved protocol is absent from the revised
sequence; its entry-protocol rule remains open.

## Supported checkpoint conversion as a guess — 2026-11-10

**Author-approved:** automatically reuse a valid compatible checkpoint as a
guess when the engine supports conversion to the revised protocol. Load saved
representation defaults before projection, rebuild operators, start KAIN empty
and establish fresh metric comparisons and convergence checks. Preserve old
history without treating it as measurements of the converted state. Unsupported
or failed conversion stops with the checkpoint intact and an actionable reason.
This does not relax strict restart or add alternate-generation selection.

The remaining detail is which requested protocol receives that guess when the
saved representation is absent. Section 4 presents starting at the first
requested protocol versus entering later under an explicit compatibility rule.
Starting at the first is recommended for simplicity, but not yet approved.
Existing matching-stage and appended-stage rules remain unchanged. Only
documentation changed; no conversion or numerical tests were run.

## Intermediate tolerances and deferred stall detection — 2026-11-10

**Author-approved:** use only the legacy threshold-dependent tolerance heuristic
at intermediate MRA stages. Keep the shared metric definitions, selected
applicable criteria, measurement point and final convergence requirements.
This supersedes fixed intermediate iteration counts. `maxiter` is an update cap
for each executed stage and current invocation; historical updates are still
not subtracted from it.

An intermediate stage ends when all its required metrics pass their resolved
stage tolerances or its cap is reached. Record these as different stop reasons;
cap exhaustion does not claim convergence. Independent response states finish
that work before advancing together, and excited roots remain active within
their coupled bundle. The requested final threshold and polynomial order remain
mandatory, with no heuristic relaxation of the final requested tolerances.

The existing [response policy](src/apps/molresponse/solvers/convergence_policy.hpp)
and [SCF limit helper](src/madness/chem/SCF.h) supply the starting heuristic
formulas, including targets derived from `max(threshold, dconv)` and explicit
multipliers. Those formulas must be mapped to the agreed metric units and
reductions rather than copying each engine's old stopping-metric selection.
The shared resolver reports requested and intermediate absolute tolerances,
heuristic inputs and effective limits. Relative tolerances and scales retain
the agreed shared conventions. Metrics without an applicable legacy adjustment
keep their requested tolerances; record and validate these cases during adoption.

The group favors studying stall detection, but it will not stop or advance
solves in the first implementation. Metrics and timings provide the study data.
Warmup remains disabled, KAIN restarts empty, and the evaluation/update split,
checkpoint guarantees and strict recovery policy remain unchanged. The design,
implementation checklist and review examples now reflect this decision. No
production code or numerical tests were changed or run.

## Drop 6 implementation reuse audit — 2026-11-10

The author requested a reuse audit before finalizing the implementation plan.
[The audit](drop6_reuse_audit.md) examines pinned drop 6
`9179f5d32076a2ba555906dd146dc330ffc040d7`, with current checkout
`76255552e16e542f6b5ae2d0d7e2b6b034c7ed48` treated separately. Existing source
and test assertions were inspected; no builds or numerical tests were run.

The reusable foundation includes functional application bindings, response
drivers, metric operations, restart planning, checked checkpoint publication,
retention and numerical serializers. Missing behavior is concentrated in
consistent current-state evaluation, one criteria policy, Nemo integration,
excited-root comparisons, direct dependency evidence and common output/timing
records. VBC already participates in shared response checkpoint publication;
Nemo already has a workflow binding and orbital save/load machinery.

The plan now calls out those adaptations instead of implying replacement
infrastructure. It keeps state/equation work application-owned, limits the new
shared records to coordination and evidence, and preserves unaffected QC decks.
Current fixed-work amendments must be adapted to the approved heuristic policy,
with unrelated fixes and useful assertions retained. No production code was
changed or retired. Other engines still require individual adoption evidence;
all-engine release scope, open decisions and independent-review status remain
unchanged.

## Implementation scope approved; common record fields proposed — 2026-11-10

The author approved the revised scope following the drop 6 reuse audit. Proceed
with consolidation and adaptation of existing drivers, metric helpers,
serializers, restart planning and checkpoint publication. Keep the all-engine
release requirement and the scientific validation work identified by the audit.
This is scope approval, not independent review or implementation completion.

The next review item is the proposed common task record in section 6: identity,
settings, dependencies, outcome and checkpoint information shared by all tasks;
iteration/evaluation evidence only for iterative tasks; optional application
diagnostics. Task completion is distinct from finishing an intermediate
protocol. Direct work uses completion and input evidence without a fabricated
convergence flag. Property and timing data stay in their selected separate files.
The record grouping and field names remain proposed pending author review.
Other open choices have not been selected.

## Application stages and internal task records — approved 2026-11-11

The author clarified and approved three levels: workflow-managed application
stages, application-managed internal tasks, and MRA protocol stages within an
iterative solve. For the response workflow, moldft supplies the ground-state
reference to molresponse. Molresponse retains its internal linear/excited-state,
VBC and property dependencies, including any supported quadratic-response solve
that consumes VBC. No new method or workflow graph is introduced.

The application-stage record describes its requested properties, resolved
settings, upstream reference and overall completion. Requested values remain
in `results.json`, grouped by calculation. Internal task records hold numerical
state/evaluation evidence, metrics and checkpoints, with source identities
supplied by the application. Shared helpers handle metrics, validity/publication
conventions and I/O without owning the internal graph. A simple SCF stage can
use one orbital-solve task without adopting response's structure.

Finishing an intermediate protocol, completing an internal task, and completing
all requested stage properties remain distinct. Direct tasks use completion
and source evidence without manufactured iterations or residuals. Exact field
names, metric subfields, reason strings and supporting filenames remain open;
the previous common-task field list was a proposal and is superseded by this
approved separation. Other scientific and release decisions are unchanged.

## Minimal stage and task field proposal — 2026-11-11

At the author's request, section 6 now proposes a concrete small field set.
Shared context is `version`, `calculation_id` and `attempt_id`. Application-stage
records add application, resolved settings, requested properties, upstream
inputs and outcome. Internal-task records add task ID, kind, application-defined
definition, source inputs and outcome, with iteration evidence only for iterative
work. The existing checkpoint manifest continues to own generation and archive
inventory. Properties, timings and long analysis histories stay separate.

An exact saved settings context and source checkpoint references are required
for reuse; a mutable latest-stage input cannot reinterpret an older task.
Starting-guess provenance does not create a continuing equation dependency.
This is a field proposal, not approval of keys, reason vocabulary, nested metric
schemas or new filenames. The approved stage/task separation and scientific
requirements remain unchanged.

## Minimal record fields accepted — 2026-11-11

The author accepted the minimal stage/task field set and asked what still needs
consideration. The field names and their division are now recorded as approved.
Keep the existing checkpoint manifest for generations and archives; there is no
additional top-level field proposed here.

Clarify the interpretation during implementation: stable IDs do not prove
scientific compatibility; saved settings remain immutable, with changes treated
according to the existing restart table; and an incomplete application stage
retains its valid completed task checkpoints and available properties. Root
correspondence and source-generation checks remain required. The remaining
record decisions are nested metric/property/timing shapes, a consistent reason
vocabulary, units and supporting filenames. No new restart policy or independent
review approval is implied.

## Shared metric record proposed — 2026-11-11

Continue from the approved minimal stage/task fields with a metric subrecord.
The proposal separates raw measurement status/value from the optional selected
stopping test, with exact saved definition and resolved-criterion references.
Current sample context is inherited; change metrics identify a previous
comparable sample, and coupled bundles identify each root. Required supporting
quantities remain available without duplicating settings or numerical archives.

The test reports scale, effective limit and its result. Missing comparison data,
inapplicable quantities, invalid computations and values above tolerance remain
distinct. Absolute-only tests do not need a scale; a missing requested relative
scale cannot silently fall back to an absolute test. Existing intermediate/final
policy, root correspondence, restart and publication rules are unchanged.
Field names and grouping are proposed for author review. Exact unit identifiers
and the complete reason vocabulary remain open. No implementation or independent
review is claimed.

## Application state and metric ownership approved — 2026-11-11

The author approved applications owning current state, comparison history, raw
proposal/residual and numerical workspace. Applications call shared metric
helpers for agreed numerical conventions, then pass scalar measurements and
comparison evidence to the shared controller. The controller handles validity
checks, tolerance decisions and recording; it does not receive solver functions.

Existing application loops can submit measurements directly; shared drivers can
invoke callbacks that capture application-owned data. Both use the same typed
measurement interface and decision code. Cached quantities may be reused when
their definition and sampling agree. Scientific review and adapter tests still
establish physical meanings and root correspondence.

The design overview, metric descriptions, driver phases, implementation step 3,
validation requirements and functional walkthrough now follow this ownership.
Numerical metric definitions, stopping/restart rules, input knobs and checkpoint
cadence remain unchanged. The proposed persisted metric schema is not approved
by this decision. Next specify the minimal application-to-controller measurement
handoff and distinguish it from metadata filled in by the recorder. Exact C++
signatures and the other existing open choices remain under review.

## Measurement return walkthrough proposed — 2026-11-11

At the author's request, the functional walkthrough now shows a proposed
`app.evaluate()` result for one SCF subject and an excited-state root inside a
coupled bundle. Both carry the current sample and per-subject measurements,
with comparison identities for changes, values/availability and relative scales
when needed. Definitions resolve through the saved task context. Shared helper
results retain required scalar evidence; numerical state and history stay with
the application. The controller supplies tolerance tests and stopping decisions.

Examples cover unavailable first comparisons and ambiguous root correspondence;
a valid BSH residual does not require a previous iterate. All required roots
remain part of the bundle evaluation. Numbers and type notation are illustrative,
not scientific reference data or implemented interfaces. The return structure
and exact field names remain proposed for review. The main design's approved
ownership, numerical conventions and restart policies have not changed.

## Evaluation return approved; controller reply proposed — 2026-11-11

The author accepted the SCF/excited-root return structure and requested its
inclusion in the design. Section 7 now specifies one evaluated sample and a
collection of subjects, with subject identity, valid change comparisons and
metric values/availability, relative scales and supporting scalar evidence.
Definitions bind through the saved task context. The implementation step and
functional walkthrough now identify this structure as approved. Applications
retain functions, comparison data and numerical work; no convergence flags or
tolerances are added to their measurement return.

The next proposed interface is the controller's reply: sample, action, reason,
checkpoint request and per-subject checks. Actions distinguish an update,
intermediate protocol completion, final task completion and an incomplete stop.
The proposal applies existing policy, including intermediate cap advancement,
final cap nonconvergence, bundle-wide decisions and safe checkpoint publication.
Exact reply fields/action spellings remain proposed; numerical policies and the
existing approval status of persisted metric fields are unchanged.

## Saved measurements determine reuse — approved 2026-11-11

The author confirmed that a tighter request already satisfied by saved metrics
must reuse the checkpoint without another numerical evaluation. Changed criteria
or stopping tolerances are now distinct from changed scientific/operator settings.
The shared controller reassesses valid saved values and available scales/evidence
with the current request. Missing required evidence needs fresh measurements;
scientifically invalid evidence and corrupt checkpoints follow their existing
compatibility/strict-restart rules. A coarser representation cannot substitute
for the requested final threshold and polynomial order.

The author also accepted the controller's four-action reply in the preceding
exchange. Eligible saved final reuse returns complete with no checkpoint write;
keep original sample/attempt and generation identities. Current checks and
outcomes describe the current request. Old outcome labels do not control reuse.
Protocol progress, original settings and source metadata remain necessary for
interpretation and planning. Unchanged eligible sources preserve direct-product
reuse; updated sources invalidate dependent evidence as before.

Checkpoint metric objects now contain measured evidence only. The earlier
proposed test subrecord is removed from them; limits and pass/fail checks belong
in current run reports/analysis. Retain diagnostic measurements and available
relative scales, so later selections can use them. The design, walkthrough and
implementation/validation plan reflect this change. Exact metric encoding,
reason vocabulary and other open choices remain under review. No code was
implemented and no numerical tests were run for this documentation change.

## Per-task protocol entry approved — 2026-11-11

After examining mixed-progress response restarts, the author selected independent
entry planning. Preserve matching saved threshold/order positions: reuse a final
state when measured evidence qualifies, resume final unmet work there, resume
an interrupted intermediate stage there, and advance eligible completed work.
New independent states start at the first requested protocol. An otherwise
compatible checkpoint whose representation is absent from the sequence enters
there as a guess only when conversion is supported. A later seed transfers no
protocol progress. Other states are not reset.

Response retains its waves and internal dependency graph. Only tasks still
requiring a protocol participate in its barrier; later states sit out earlier
work and need no retained coarse archives. Keep explicit per-state plans during
migration and replace threshold-only/derived-k entry matching with the agreed
threshold/order identity. Excited roots remain one bundle; expansion keeps the
existing saved-protocol rule when matching, with the conversion fallback applying
to the whole bundle when necessary. Strict restart, scientific compatibility,
current-request metric checks and the mandatory final representation remain.

The main design, functional walkthrough and implementation/validation plan now
reflect this decision. The conversion-entry item is removed from open choices.
Remaining schema, scientific-owner, parameter-mapping and release-scope decisions
are unchanged. Existing source tests informed the discussion; no tests were run
and no production code was changed in this documentation update.

## Initial release scope approved — 2026-11-11

The author deferred E13–E15 (`testperiodicdft`, `tdse`, `tdse4`) from the initial
madqc release and required preservation of their separate builds. The design
retains all fifteen inventory entries, requires shared-layer adoption and
scientific validation for E01–E12, and replaces conditional research integration
with standalone-build preservation. New research adapters, TDSE-specific 4D
metric/archive work and research QC cases are deferred. Existing programs,
inputs and available coverage remain.

This resolves the release-boundary choice without changing convergence,
restart, protocol-entry or metric conventions. Remaining work concerns concrete
record schemas, existing option mappings and scientific model definitions.
This edit does not establish implementation, build success or independent review.

## Grouped tensors and shared result metrics approved — 2026-11-11

The author approved the alpha/beta layout with one compact state-metric
collection alongside the grouped tensors. A full request has nine alpha and
twenty-seven beta components; partial requests and availability use the same
format. Each actual contributing state has one measurement row, with sufficient
frequency/operator/root context to distinguish it. Shared ground measurements
appear in the ground-state calculation; shared response states are not repeated
under every property or component.

The result measurements are a snapshot of the states used to compute the
properties. Full definitions/comparison/scales/source evidence remain available
in retained records, while compact definitions, units and representation make
the public measurements readable. Updating a checkpoint does not silently change
previously published evidence. Direct VBC work has no invented residuals;
checkpoint-based restart and current-request tolerance decisions are unchanged.

The design, output examples and implementation/validation requirements now
reflect this approval. Exact encoding details and scientific unit/default
validation remain open. No numerical implementation, build or independent
review was performed by this document update.

## Explicit convergence names approved — 2026-11-11

The author selected explicit `convergence.criteria`, `convergence.<metric>`
and optional `.relative` inputs without `dconv`/`econv` aliases. The naming rule
applies to file input, parameter-group overrides and adopted standalone entry
points. Reject explicitly supplied old keys with a clear migration explanation;
do not silently ignore them, convert their numbers, or allow a mixed old/new
request to bypass the error. Use shared reviewed defaults without a fallback to
old coupled convergence controls.

Workflow selectors and CLI group forms remain supported. Engine-specific inner
controls retain their own reviewed scopes; the decision is not a blanket
source-code rename. The design, option table and implementation/validation plan
now record the rule. Numerical defaults remain subject to validation, and no
production code or builds changed. The next suggested decision is load-only
restart when saved evidence does not satisfy the current request.

## Read-only post-processing and solver restart approved — 2026-11-12

The author separated reading/plotting from solving. Read-only loads the data
needed for inspection or plotting, preserves its original scientific context
and completion status, and permits unfinished saved work. It checks readability
and operation requirements, not satisfaction of new solver tolerances. It may
write plot/export files but does not iterate, advance protocols or overwrite
solver checkpoints and convergence claims.

Solver restarts reuse compatible saved work if its measurements and requested
final protocol satisfy the request; otherwise they resume solving through the
existing progress/evaluation rules. Strict failures for expected corrupt,
missing or unsupported checkpoints remain. The earlier proposed load-only
solver exception is superseded. No new CLI spelling is selected; legacy entry
points must implement this split without silently invoking a solve from inspection.

The design, input mapping and functional walkthrough now express the split and
use "requested final protocol" for final threshold and k. The next suggested
input decision concerns the old ES-only `excited.maxiter` override. No production
implementation, numerical tests or independent review were performed here.

## Single response maxiter input approved — 2026-11-12

The author required response `maxiter` only and a clear migration error for
`excited.maxiter`. This is a rejection, not an alias or an excited-only override;
mixed old/new requests fail even with equal values. The standalone response
`--maxiter` entry supplies the same resolved setting. Remove the old override
from the response adapter and prevent any separate excited-state default from
replacing the shared value.

Each independent driven state and each coupled excited bundle retains its own
allowance per executed protocol/invocation. This does not establish a pooled
budget across states, change numerical defaults, or merge geometry/inner-solver
budgets. Warmup and terminal-evaluation rules remain unchanged. The design,
input table, walkthrough and planned validation now cover the decision. No
production implementation, numerical runs or independent review occurred.
The next suggested input choice concerns conflicting legacy checkpoint overrides.

## Checkpoint inputs removed with generic errors — 2026-11-12

The author chose to remove `save` and `checkpoint.keep_all` completely from the
new public input interface. Any supplied value is unknown input; no special
migration message, alias, deprecated registration, or warning-and-ignore path is
retained. This supersedes the proposal to handle only the conflicting values
with tailored migration errors. The decision is scoped to these checkpoint keys.

Mandatory endpoint/valid controlled-stop saves, latest-usable retention with
necessary temporary/referenced generations, and `checkpoint.every_iterations`
remain unchanged. The design, option table, walkthrough and planned validation
now record the decision. No production code or numerical tests changed.

## Consolidated contract for review — 2026-11-12

The author requested one consolidated revision and a small SCF/response review.
This revision selects mechanical record and routing details for that review; it
is not a new independent review or approval of unresolved scientific definitions.

- Main design section 1 fixes auto/none/iterate/continue/inspection routing, scoped
  response restart adoption, conflict checks and removal of deferred-policy knobs.
  Existing approved migration errors remain distinct from generic removed keys.
- Section 6 replaces repeated output proposals with concrete version-2 records,
  names, statuses, reasons, unit conventions, source associations and one timing
  node shape. Existing version-1 manifests are not silently reinterpreted.
- The [complete example](output_record_examples.md) follows the maintained LiH
  request through SCF, a shared static z state, direct VBC and alpha/beta output,
  then tests reuse, tighter requests, interruption and inspection. Numbers are
  illustrative. The grouped-tensor and functional examples now use this contract.
- Scientific definitions/normalization, scoped inner controls, property-specific
  scientific context and numerical validation remain explicit adapter prerequisites.
  E01–E12 adoption, E13–E15 standalone preservation and collaborator cleanup
  approval are unchanged.

The updated contract awaits the author's review as a whole, followed by fresh
independent design review. No production code, build or numerical test was run.

## Record layout accepted — 2026-11-12

The author accepted the example's execution/restart behavior, then the record
fields, with the qualification that repetition should stay limited. Shared
context is stored once at the applicable enclosing scope; states are reported
once across properties and detailed source associations remain separate. Each
checkpoint still retains the settings and definitions needed to validate it.
This accepts the layout, not unresolved model definitions, numerical defaults
or implementation. A fresh independent design review remains outstanding.

## Stage 6 comprehension clarifications — 2026-11-12

The author requested revisions and recommendations after the independent reader
reported C-01/C-02 (Major) and C-03 (Nit). The draft now makes the existing
contract explicit: bundle checks establish validity, while selected per-root
metrics alone supply convergence tests; unconverged observations remain in
analysis/human reports, while requested properties lacking eligible inputs are
unavailable without numeric values. Valid independent results remain available.
Stale record-approval wording is corrected. No new metric, threshold, result
status, dependency exception or production code was introduced. These revisions
are recorded for review; they do not certify scientific definitions or implementation.

## Stage 7 round 1 revisions applied — 2026-11-12

The author approved all six proposed edits and another critique round. The draft
now states the reused checkpoint trust/concurrent-writer limits, groups the
validation cases without dropping any acceptance text, consolidates protocol-entry
rules, removes repeated approval narration and effort labels, and drops two
hypothetical input spellings. The detailed implementation plan and collaborator
cleanup gate are unchanged. Round 1 remains not clean under the review-loop rule: its
Major findings were raised before these fixes. A new independent round follows.

## Stage 7 round 2 critique — 2026-11-12

Three fresh reviewers read the revised feature specification (lines 1–1896).
Design and structural reviewers returned no findings; the style reviewer raised
one Nit: remove two repeated durability caveats while retaining the canonical
storage limitation. The proposed deletions are prepared in the annotated report
and await the author’s decision. The draft and all linked context files remained
unchanged throughout review. This is the first consecutive clean round under the
EGM rule; Stage 7 remains open and another round needs explicit confirmation.
The detailed implementation plan was outside this round’s scope; the full final
review remains incomplete.

## Stage 7 round 2 cleanup applied — 2026-11-12

The author approved the remaining cleanup and continuation. Two repeated
sentences about MPI-process-loss and power-failure recovery were removed from
strict restart and output publication. The storage section retains the single
explicit limitation and its implementation link. This changes no restart,
publication, or scientific requirement. The detailed implementation plan is
unchanged. A third fresh critique round follows; round 2 remains the first clean
round.

## Stage 7 round 3 critique — 2026-11-12

Design and structural reviewers returned no findings. The style reviewer raised
two Nits: identify the baseline for the word “unchanged,” and consolidate two
old-format compatibility statements. Both passages predate the approved cleanup;
neither proposal reverses an earlier edit. Exact replacements are prepared but
await the author’s decision.

Rounds 2 and 3 are consecutive clean rounds under EGM's rule, so the critique has
converged without a fourth round. The target and all 72 linked context files
remained unchanged during review. Stage 8 awaits the agreed stage check-in;
Stage 9 and the full final review remain incomplete. The detailed implementation
plan was outside this round's scope.

## Stage 7 closed; Stage 8 authorized — 2026-11-12

The author approved both final wording fixes: direct submission and callbacks
now explicitly follow the criteria and policies defined above, and the archive
reuse paragraph states the old-format compatibility boundary once. The separate
collaborator-approval requirement is retained. The detailed implementation plan
is unchanged. Stage 7 is complete after two consecutive clean rounds, with all
proposed fixes applied and no open findings. Stage 8 now reviews the full draft
against the initial zoom-out baseline.

## Stage 8 intermediate zoom-out review — 2026-11-12

The fresh reviewer read the entire revised draft and compared it with the
initial review snapshot. No Blocker or Major findings: the additions clarify
existing boundaries without introducing another scheduler, recovery system or
convergence framework. Code searches support the principal reuse choices.

One Nit remains: the implementation plan says to reuse the existing wall-clock
source, but POSIX wall_time uses gettimeofday while the contract requires
monotonic elapsed timing. A one-line change to std::chrono::steady_clock in the
new collector is prepared, preserving the timing boundaries and no-extra-fence
rule. It awaits the author’s decision. The Stage 8 audit is complete; Stage 9
and full final review remain pending. No production code or tests were changed.

## Stage 8 closed; Stage 9 authorized — 2026-11-12

The author approved steady_clock inside the shared collector, with existing
MADNESS timers unchanged. Both the timing contract and implementation step now
state that boundary. No timer implementation was modified. Stage 8 has no open
findings. Stage 9 starts a fresh implementation-readiness loop on the detailed
plan with implementation, codebase and style reviewers; earlier feature-review
clean rounds do not count toward this separate loop.

## Stage 9 round 1 implementation-readiness review — 2026-11-12

All three fresh reviewers completed. The result is one Blocker, five Major
findings and two Nits. Required scientific contracts remain unresolved for
parts of E01–E12; correlated representation changes and restart calculation
matching also need explicit decisions. The implementation sequence has a
cutover-order problem, and shared test-helper/fixture migration plus old-oracle
comparison need specified transformations. The two Nits concern repeated
status wording and the term “controller tests.”

The primary agent verified the reported source constraints and prepared plan
repair proposals without applying them. The scientific Blocker is not treated
as resolved by adding a prerequisite label. No scientific definitions were
invented, release engines deferred, owner messages sent or source/tests changed.
Stage 9 remains open with zero clean rounds; the full review is incomplete.

## Stage 9 resolution — phase ordering approved

The author approved readiness-1-4 and asked to discuss the shared checkpoint-test
helpers next. The plan now prepares the lifecycle before engine activation,
requires each engine's input/readers/fixtures/tests to migrate with that engine,
and keeps global retirement after all callers and coverage have migrated and
collaborator approval is recorded. The detailed test-file inventory remains in
step 11, explicitly used during earlier engine adoption. E13–E15 preservation
checks run when affected shared changes occur. Scientific prerequisites and
release scope are unchanged. The repair awaits fresh review; other findings
remain open, and no checkpoint-test helper changes have been applied.

## Stage 9 resolution — checkpoint-test helpers approved

The author approved codebase-1-1 and asked to discuss numerical-reference
comparison next. The plan now inventories shared_records.py, its focused tests,
fixture generation/catalogue/decks, old-format rejection fixtures and all 24
existing helper consumers. It keeps numerical archives outside Git and preserves
source/input/file integrity checks and existing temporary-directory cleanup.
The reader exposes properties, immutable saved measurements and current-request
assessments separately. Engine activation includes the corresponding migrations.
No test helper, fixture, numerical reference or solver code was changed; this
updates the implementation plan. The repair awaits fresh review.

## Stage 9 resolution — durable test curation requested

The author requested a curation step and an index explaining why each test
exists and what it checks, independent of output format. Step 1 now plans
`docs/design/madqc-test-index.md`, linked from the QC and madqc READMEs, with
stable IDs, purpose, assertions, evidence/provenance, acceptance rationale and
current executable links. Engine adoption and final cleanup must account for
every migrated assertion. The index distinguishes scientific meaning from
file/key mappings and records changed expectations or intentional retirement.
This adds the curation requirement to the plan; it does not populate the index,
change tests/references or settle the pending reference-mapping proposal.
Stage 9 findings and fresh-review status are unchanged.

## Stage 9 resolution — numerical reference mapping approved

The author approved codebase-1-2 and asked to walk through the next finding.
The plan now specifies independent actual/reference filenames and field paths,
with defaults preserving unmigrated tests. The harness and comparison helper
change together with focused tests. Existing values, acceptance semantics,
near-zero safeguards and external-reference protection are retained; schema
migration cannot silently replace a reference. The curated index links every
old assertion to its destination, and current completion/refusal and lifecycle
checks use the new runtime records. The helium example retains its energy
tolerance and does not remove model or eigenvalue checks. No source, tests or
numerical reference files changed. This repair awaits fresh review. Five other
findings remain unaddressed; readiness-1-1, the scientific-definition Blocker,
is the next discussion. No new review round has been authorized.

## Stage 9 discussion — Nemo metric source review

The author asked to separate current Nemo definitions from questions for
Florian. The source review is recorded in
[the metric evidence document](docs/design/madqc-metric-validation.md). Nemo
already defines regularized orbitals, RMS equation residuals, a unit-weight
one-spin iteration density, total-energy change and a corrected physical
kinetic diagnostic (`energies[1]`). The source also shows different occupation
handling in the property density and a mixed pre-/post-update metric boundary.
The remaining scientific questions concern residual representation, intended
spin/occupation scope and confirmation of the physical kinetic estimator.
The shared maximum reduction and prepared-state/checkpoint boundary remain
approved. No scientific definition was approved by this audit; no code, tests,
references or design contract changed. No owner message or new review round
was sent or started. Readiness finding readiness-1-1 remains open.

## Stage 9 discussion — Nemo definitions recorded as proposals

The author selected maximum regularized residual, physical per-spin density
and the existing physical kinetic quantity as the proposed Nemo definitions,
pending Florian's confirmation. The design and metric evidence now record
those choices and keep occupation scope, derivative/smoothing conventions and
units explicit confirmation points. Total-energy change and the prepared-state
evaluation boundary remain unchanged. No scientific approval, numerical
validation or implementation is claimed. The broader readiness-1-1 Blocker
remains open, including other engines' definitions. The next discussion is
readiness-1-2: correlated representation transitions.

## Review scope update — CC2 and 6D engines deferred

The author asked not to consider CC2 or any 6D codes at the moment. Their
scientific-definition and representation-transition work is deferred from
the current review discussion. Readiness-1-2 is deferred, not resolved; its
original severity and evidence are retained. The current pass moves to
readiness-1-3, restart calculation matching, using SCF/Nemo and molresponse.
This temporary discussion boundary does not change builds, delete planned
work, or decide their eventual release inclusion. No full-scope readiness is
claimed. The design's release inventory remains unchanged.

## Stage 9 discussion — ground-state to response matching

The author requested the simplest checkpoint association for the ordinary
ground-state → response workflow. The proposed rule adds an internal workflow
role to existing run-inventory entries: ground_state → saved calculation ID and
response → saved calculation ID. Each application then resolves its own task
checkpoints. Identity remains separate from scientific compatibility, selected
criteria and exact dependency generations. The proposal covers tighter requests,
added frequencies and extending SCF to response. It is not yet approved; first
publication and repeated-calculation identity still need explicit rules.
[Worked proposal](/tmp/madqc-final-independent-review-yapnzvuv/stage9-round1/calculation-matching/proposal.md).
No draft contract, source, tests or numerical references changed.

## Stage 9 resolution — direct calculation locations approved, 2026-11-13

The author approved simplifying the ordinary ground-state → response case to
workflow-assigned names used directly as calculation IDs and checkpoint
locations: ground_state and response. The design and linked output/tensor
examples now use those names. No separate role-to-ID lookup or input naming
is introduced. Validation, exact reference-generation dependencies, per-task
response checkpoints and strict missing-checkpoint behavior remain unchanged.
Repeated-calculation occurrence names and first-publication ordering are still
open, so readiness-1-3 is only partially addressed. No implementation or tests
changed; CC2/6D work remains deferred from the current review.

## Stage 9 resolution — first-publication rule approved, 2026-11-13

The author approved null-before-first-publication in the existing task-level
latest pointer. The design now initializes that pointer before inventory
registration, permits numerical writes only after registration, and publishes
by atomically replacing the pointer after complete validation. Failed first
saves retain null; failed replacements retain the prior selection. Missing
registered pointers and missing selected payloads fail strictly. Direct tasks
without numerical payloads keep completion records. The implementation plan
includes collective selection/registration, test-reader handling and failure-
point checks. No source or tests changed. Repeated-calculation identity remains
open; CC2/6D work remains deferred from the current review.

## Stage 9 resolution — optimization composition, 2026-11-13

The author approved workflow assembly of a consistent energy and nuclear gradient
from SCF and response, while each application owns its computation, metrics and
checkpoints. Section 3 now defines total excited-state energy, explicit gradient
contribution semantics, common evaluated context and dependency evidence. The
plan prepares this boundary with deterministic providers and checks total-energy
assembly, gradient double counting and invalid dependencies. PySCF source and
optimization-guide references are included. This does not implement excited-state
forces or settle state following, gradient formulation or repeated-calculation
occurrence identities. Readiness-1-3 remains open; no new review round was run.

## Stage 9 resolution — saved optimization evaluations, 2026-11-13

The author approved persistent geometry-evaluation numbers shared by each
SCF/response pair, recorded before electronic work and restored on restart.
The design separates those numbers from optimizer steps and execution attempts,
keeps exact ground-state dependency checks and preserves ordinary singleton
names. It retains required optimizer/pending/reference archives, publishes a
replacement optimizer checkpoint before retiring obsolete work, and records
retirement before reclaiming unreferenced archives. The plan includes workflow
and interrupted-cleanup checks. Missing active checkpoints still fail strictly.
Optimizer payload/save boundaries and numerical-displacement occurrence naming
remain open. No production code, numerical references or review round changed.

## Stage 9 resolution — optimizer snapshot and save points, 2026-11-13

The author approved the minimum optimizer fields and publication boundaries.
The design now saves context, evaluation identities/coordinates, move count,
current/previous measurements and displacement, approximate Hessian, pending
line-search data and stopped outcome. The next operation remains either
evaluate_geometry or evaluate_trial, including across a controlled stop.
Snapshots publish before current/trial evaluations and after step selection;
phase, Hessian, history, coordinates and counters advance together. The plan
includes deterministic interruption checks for duplicate-update prevention.
Changed-settings history/reset, outer maxiter budgets and numerical-displacement
identities remain open. No production code or numerical references changed,
and no independent review round was run.

## Stage 9 resolution — optimizer history after settings changes, 2026-11-13

The author approved preserving valid optimizer history and resetting it when
changed settings invalidate its energy or gradient evidence. The design now
resumes at the saved current geometry, reuses eligible electronic work and
publishes evaluate_geometry with the initialized Hessian and cleared comparison/
trial history when required. Geometry/evaluation counters are preserved; discarded
trials are not promoted. Tighter requests already satisfied by valid evidence keep
history, while old/new precision samples do not create false geometry-change
measurements. The plan includes preservation/reset and publication checks.
Outer maxiter budgets and numerical-displacement identities remain open. No
production code or numerical references changed; no new review round was run.

## Stage 9 resolution — fresh optimization maxiter allowance, 2026-11-13

The author approved a fresh allowance of committed geometry steps on every
optimization restart. The design counts a step at publication of the selected
next geometry; historical steps remain separate. A resumed trial consumes one
new step when committed, while assessment of an already published geometry does
not count it again. Initial/trial/final evaluations and history resets do not
consume geometry steps; inner protocol changes do not renew the allowance.
The final geometry is assessed without another move, with incomplete electronic
work kept distinct from convergence. The implementation plan includes these
budget and interruption checks. Numerical-displacement occurrence identities
remain open. No production code, numerical references or review round changed.

## Stage 9 resolution — finite-difference restart, 2026-11-13

The author approved saved parent displacement plans and per-point evaluation IDs.
Changed geometries receive new IDs; unchanged geometries keep IDs and reassess
accuracy/dependency evidence. Eligible completed points can be reused individually,
but derivative assembly waits for the complete required set and records exact
sources. Plan/source changes invalidate the current derived result as appropriate;
publication precedes retirement and displacement work does not consume optimizer
moves. The plan includes deterministic reuse/invalidation/assembly checks, without
claiming numerical-gradient implementation. Readiness-1-3 is now recorded as an
approved fix applied, pending independent re-audit. Scientific definitions remain
open, and CC2/6D transition work remains deferred. No code, numerical references
or independent review round changed.

## Stage 9 preparation — Nemo confirmation checklist, 2026-11-13

Prepared docs/design/nemo-metric-confirmation.md with three source-backed items:
regularized residual representation, per-spin density/occupation weighting, and
the existing physical kinetic quantity with derivative/smoothing conventions.
Each separates current code from the proposed definition needing confirmation;
units and supported cases are explicit. The metric evidence links the checklist.
Already-agreed sampling/history work is identified separately. No message was
sent, scientific approval recorded, numerical validation run or design changed.
Readiness-1-1 remains open pending confirmation.

## Stage 9 resolution — wording cleanup, 2026-11-13

Applied both pending wording repairs: removed repeated scope/approval narration
from the reuse section and named runtime criterion assessments as `Decision.checks`.
Corrected `excited.maxiter` to the agreed generic unknown-input rejection in the
contract and planned checks; `dconv`/`econv` guidance is unchanged. Scientific
requirements, release scope and restart behavior are unchanged. Both Nit fixes
are applied pending independent review; scientific definitions remain open and
CC2/6D transition work remains deferred. No new review round or numerical tests ran.

Prepared [the acceptance walkthrough](docs/design/madqc-acceptance-checklist.md)
as an index of existing behavior, with stable case IDs and no passing-test claims.
The next discussion starts with SCF inputs, evaluation and restart.

## Acceptance checklist reviewed — 2026-11-13

The author reviewed and accepted coverage for all 16 cases in
[the acceptance checklist](docs/design/madqc-acceptance-checklist.md): SCF,
response, checkpoint recovery/storage, output/timings, and optimization/finite
differences. This records agreement on expected behavior, not passing tests.
Scientific confirmations remain pending; CC2/6D work remains deferred from the
current review without changing release scope. The six applied review repairs
still await independent re-review. No code or numerical references changed,
no tests ran and no further independent review round was started.

## Final consistency pass — 2026-11-13

Compared the examples, input mapping, record fields and implementation plan
against all 16 reviewed acceptance cases. Applied the corrections in
[the consistency report](docs/design/madqc-consistency-pass.md): obsolete-input
rejection, saved-evidence reuse, comparison reason names, atomic publication,
VBC timing placement, response timing/dependency order, root-comparison validity,
and finite-difference retention coverage. Linked acceptance IDs into planned test
curation. Scientific definitions, the reviewed checklist, numerical references
and release scope remain unchanged. Embedded JSON and local document links pass
the document checks; this is not a new independent review or numerical validation.
The two scientific/deferred findings remain open and the six previously applied
review repairs still await re-review. No code changed or owner message was sent.

## Independent Stage 9 review — round 2, 2026-11-13

All three fresh-context reviewers completed against design SHA-256
`226b64655bf0e110c20ebbe54ba6c9231bc9f4a04d048a994768a24edbd6eaf2`.
After merging one optimizer overlap, the round has eight Major findings and
one Nit. See [the readiness report](implementation_readiness_report.md) and
[annotated resolutions](goldfish_review_report.md). New work concerns optimizer
precision and caches, raw/collective CLI validation, existing SCF virtual and
Nemo CPHF paths, HDF5 publication, OEP coverage wording and prerequisite timing.
The prior scientific-definition Blocker and deferred CC2/6D finding remain open.
No design edits, builds or numerical tests were performed in this round.
The round is not clean; another review round requires author direction.

## Scope revision after review — 2026-11-13

The author deferred optimizer replacement until after workflow unification.
The [revised design](design_draft.md) preserves the existing optimizer and moves
new snapshots, Hessian/line-search recovery, geometry-step restart rules and
cross-application optimizer composition to
[follow-on requirements](docs/design/madqc-optimizer-follow-on.md). Acceptance
D01/D02 are deferred; standalone D03 planning remains in scope. No code changed.

- `readiness-2-3` / raw `codebase-2-4`: replacement stopping/precision policy is
  deferred, not fixed. The retained driver's `dconv`-derived defaults still need
  an agreed compatibility mapping before electronic input migration. The legacy
  forced-success rule remains documented; it does not establish electronic
  convergence or satisfy the replacement contract.
- `codebase-2-3`: the current plan now explicitly requires bypassing/migrating
  SCF/Nemo geometry-cache shortcuts and tests at the engine boundary. This is a
  documentation revision awaiting independent review, not an implemented fix.
- All other findings and scientific prerequisites retain their status. The
  review outcome and original severities/counts above describe the frozen round-2
  target. It has not become a clean round; no new review was run.


## Approved optimizer defaults — 2026-11-13

The author approved independent `optimization` defaults `gtol=1e-4` and
`gradient_precision=1e-5`, with explicit file/CLI overrides. The design and input
mapping now remove the retained driver's `dconv` dependency without coupling
these defaults to the new electronic criteria. Other optimizer defaults,
precision clamps, Hessian-update and stopping behavior remain unchanged.
The plan names the parameter/driver changes and checks for defaults, overrides,
changed electronic settings and the retained clamp through existing SCF/Nemo
optimization coverage.

The defaults/mapping portion of `readiness-2-3` (including raw `codebase-2-4`)
is revised in the plan and awaits independent review and implementation.
Replacement stopping-policy work remains deferred; no code or numerical
references changed. The next finding is `codebase-2-5`, the publication contract
for native and optional HDF5 archives. No backend decision is applied yet.


## Approved native checkpoint scope — 2026-11-13

The author selected native MADNESS archives as the sole shared restart format
and deferred HDF5 integration beyond parameter and I/O unification. The design,
input mapping, adapter/test inventory and release acceptance now agree. HDF5
remains important follow-on usability and analysis work. Existing HDF5 utilities,
build support and standalone tests are preserved separately; this revision does
not require an auxiliary HDF5 mirror, export path or backend-aware shared reader.
Adopted execution paths reject HDF5 selection explicitly, without silent fallback.

`codebase-2-5` is revised in the plan and awaits independent review and
implementation. Existing native immutability/reuse assertions from the HDF5
integration test move into native checkpoint coverage; HDF5-specific integration
assertions are recorded as deferred in test curation. Test retirement still
requires the agreed migration review and collaborator approval. No code or
numerical references changed, and no new independent review ran.

The next finding is `codebase-2-1`: SCF's frozen-occupied virtual solve and
virtual-count step-down preparation need explicit adoption coverage.


## SCF virtual paths and PR #798 — 2026-11-13

The author approved including frozen-occupied virtual solves and extra-virtual
preparation in the shared layer. [PR #798](https://github.com/m-a-d-n-e-s-s/madness/pull/798)
was inspected through its metadata, final SCF/restart patches and discussion;
its merged commit is already an ancestor of this checkout. The design now names
the task/reference/membership and preparation contracts, scoped `nv_its` budget,
metric applicability, padding/protocol preservation and canonical saved-state
boundary. The plan reuses the existing numerical corrections and KAIN helpers.
It does not reapply the PR or change KAIN algorithms.

The input mapping, E01 inventory, acceptance cases A02–A06 and
[metric evidence](docs/design/madqc-metric-validation.md) record the additions.
Validation distinguishes bound-virtual success from the documented unbound
plateau, covers both spin counts and interrupted preparation, and preserves
explicit occupations and supported seeding. PR-reported tests were not rerun.
`codebase-2-1` is revised in the plan, awaiting independent review and
implementation; scientific defaults/definitions remain subject to their normal
validation. No code, numerical reference or PR was modified.

The next finding is `codebase-2-2`, Nemo's existing Hessian/CPHF path.


## Nemo CPHF adoption — 2026-11-13

The author approved independent atom/axis CPHF tasks inside Nemo, using the
shared evaluation, iteration allowance, native checkpoints and reuse-or-iterate
assessment. The revised design and E02 inventory include exact ground-reference
dependencies, raw/corrected representation identity, sequential scheduling,
no separate LDA warmup, and eligible-response gates for direct Hessian/vibrational
properties. Existing CPHF restart inputs normalize to the shared policy.
Input mapping, planned tests and acceptance coverage have been updated.

[Florian's checklist](docs/design/nemo-metric-confirmation.md) now distinguishes
existing code from response residual/representation, density/kinetic and
Hessian-row acceptance questions, plus numerical validation of the warmup change.
No owner message was sent. `codebase-2-2` is revised in the plan and awaits
independent review/implementation; its scientific questions remain under the
existing definition prerequisite. No code or numerical references changed.

The next finding is `readiness-2-1`: repeated CLI groups can discard a retired
input before validation sees it. Collective CLI error propagation follows as
`readiness-2-2`.


## Validate every CLI occurrence — 2026-11-13

The author approved validating original command-line occurrences before overrides,
while preserving current precedence. The design and input mapping require an
ordered user-argument record alongside the existing final maps, shared owning
syntax/key/type/value validation without live-setting mutation, and normal
resolution only after those checks. Cross-setting checks apply to the resolved
request. Repeated groups retain last-payload replacement; they are not merged.

The plan includes `commandlineparser.h`, `QCCalculationParametersBase.h/.cc`,
existing reader tests and adopted workflow/standalone entry points. Acceptance
A01 now covers hidden retired/unknown/malformed inputs, valid repeated values,
disjoint payloads, case/provenance and failures before computation/publication.
`readiness-2-1` is revised in the plan and awaits independent review/implementation.
No code or numerical references changed.

The next finding is `readiness-2-2`: distribute CLI-validation failures to all MPI
ranks before any rank leaves the common input path. Listing the shared reader's
implementation file above does not itself resolve that error-propagation design.

## Collective CLI failures — 2026-11-16

The author approved reusing the file reader's collective error handling for both
original-occurrence validation and final CLI reading. The design requires all
ranks to receive the error category and message before any rank exits that input
path, with successful settings accepted only after validation. Existing input
precedence and diagnostics remain unchanged.

The plan covers `QCCalculationParametersBase.cc`, reader tests and their MPI
registration, plus the adopting application input tests. A bounded two-rank test
must reject retired and malformed inputs with the expected diagnostic and a
nonzero exit before numerical work or checkpoint publication. A timeout or
unrelated launcher failure is not a passing negative test. Serial/non-MPI
coverage remains required.

Acceptance A01 records these checks. `readiness-2-2` is revised in the plan and
awaits independent review/implementation; the MPI tests have not been run.
The next finding is `style-2-1`: missing OEP migration evidence must remain a
release gap, rather than permission to declare an existing model unsupported.

## OEP model migration evidence — 2026-11-16

The author approved migration evidence for every supported OEP model. E08 now
requires OAEP/OCEP/DCEP/MRKS coverage, preserving the existing OAEP case and
adding the planned missing cases. The model inventory, adapter plan and
acceptance checklist require numerical and applicable shared-layer evidence,
including fresh, resumed and refined calculations.

Pre-existing limitations are recorded separately with affected conditions and
supporting evidence. Missing migration tests remain release gaps; they cannot
be resolved by declaring an existing model unsupported. A supported-scope change
requires an explicit design revision. No numerical tests or scientific
validations were run. `style-2-1` is revised and awaits independent review.

The next finding is `style-2-2`: the prerequisite-table header says all evidence
is needed before implementation, while numerical calibration belongs during
implementation and must be reviewed before accepting the affected adapter.

## Prerequisite timing — 2026-11-16

The author approved distinguishing evidence needed before implementation from
acceptance evidence gathered during implementation. The table now uses
“Required evidence and acceptance point.” Scientific metric definitions and
scoped contracts must be settled before their dependent implementation;
absolute defaults, intermediate tolerance adjustments and root-match thresholds
are measured during implementation, with evidence recorded and reviewed before
accepting the affected adapter. Final document review applies to the latest draft.

`style-2-2` is revised and awaits independent review. All nine round-2 findings
now have approved revisions recorded; this is not a clean-review or scientific
approval claim. No numerical validation or new independent review was run.

## What remains before finalization

- Independently re-review the revised draft. All nine round-2 findings have
  approved revisions recorded, but those revisions have not been verified by a
  new review round.
- Resolve the model-definition prerequisite (`readiness-1-1`) before the affected
  adapter work. Florian's Nemo ground-state/CPHF confirmation and the other
  required model-specific units, normalization and metric contracts remain open.
- Keep CC2/6D definitions and representation transitions explicitly deferred
  from this discussion (`readiness-1-2` and the related part of `readiness-1-1`).
  Their full-release obligations are unchanged; deferral is not resolution.

Numerical calibration is implementation acceptance work. It does not need to
be completed to finish the design text, but its evidence must be reviewed before
accepting the affected adapter. Replacement optimizer work and HDF5 integration
remain deferred. The draft is not yet certified as implementation-ready across
the full engine inventory.

## Independent implementation-readiness review — round 3, 2026-11-16

Reviewed unchanged draft SHA256 `0a779c874cc3c17bfee42a86513d596156cfb632eb9c8be4f0d59fd4a2c48f3d` (3123 lines; detailed plan starts at 2318).
All three fresh reviewers completed: implementation 0 Blocker/2 Major/0 Nit,
codebase 0/3/0, clarity 0/0/0. Nine prior revised findings were not repeated at
their old scopes. Five new integration gaps remain: response seed-rung input,
active regression migrations, native visualization readers, retained optimizer
electronic identities, and standalone finite-difference plan storage/selection.

Shared policy, measurements, ordinary native task checkpoints and deterministic
workflow preparation can start while preserving existing callers. Affected engine
activation and geometry integrations wait for their specific decisions/migrations.
Scientific confirmations and deferred CC2/6D obligations remain open. Numerical
calibration remains implementation acceptance work. No design/source edits,
builds or numerical tests were performed by this review; no new round is implied.
See [readiness report](implementation_readiness_report.md) and
[annotated proposals](goldfish_review_report.md).

## Seed entry through the requested protocol — approved 2026-11-16

The author chose a simpler input contract: the seed supplies initial functions,
and `protocol` alone specifies the requested stages. Keep DALTON/basis-set
projection and refinement. Recommend a single final protocol for a well-converged
basis-set seed; users can explicitly request a ladder when they want one.
No automatic fine-entry rule or replacement seed-stage key is introduced.
Imported functions receive fresh MRA convergence checks, while eligible MADNESS
checkpoints keep their saved protocol progress.

The design, input mapping, acceptance checks and response migration plan now
remove `seed.start_rung`, the test-driver `--seed-start-rung` form and hidden
ladder collapse. The plan names the shared helper, adapter, registrations,
protocol-key tests and DALTON guide. `codebase-3-1` is revised in the plan and
awaits implementation and independent verification. Source code and numerical
references are unchanged. The other four round-3 findings remain open.

## Alpha curation and full response coverage — approved 2026-11-16

The author selected helium HF plus one DALTON-reference water case for the first
alpha set. The new test index curates `response_he_alpha` and
`response_h2o_lda_alpha`, with stable IDs, exact numerical oracles, reference
provenance, all 25 existing checks and their migration dispositions. Water keeps
its slower validation tier. No source, input deck, check or reference value was
changed, and no numerical calculation was run.

The author also required full molresponse functionality, numerical accuracy and
performance coverage: every supported property; one- and multiple-electron
systems; closed/open shells; static/dynamic response; exchange/correlation
variants and terms; functional components and full workflows. The design,
acceptance checklist and test index now record this requirement. The initial
alpha pair is only a starting set. Existing limitations remain separately tested;
missing supported coverage remains an acceptance gap.

The plan now explicitly lists the three active scripted consumers and two decks
from `codebase-3-2`, with owning-phase migration and preservation of their
numerical references and all assertions. That finding is revised in the plan,
awaiting implementation and independent verification. The remaining three
round-3 proposals and carried scientific/deferred prerequisites remain open.


## Native visualization migration — approved 2026-11-16

The design now explicitly migrates `VizManifest.hpp`, its madqc call site and
the native-discovery/loading portions of `dump_mra_trees.cpp` to shared records
and the shared native loading interface. Calculation/task identity, geometry,
representation and exact source generations come from that evidence, and the
visualization manifest points to the matching new results snapshot.

INSPECT-NATIVE records the planned unfinished-sample inspection check in C04:
saved functions remain discoverable with their original context and status,
without solver work, new convergence assessment or checkpoint publication.
Requested plot files are allowed. Prepare shared fixtures in step 4 and require
the relevant checks before SCF/response output activation in steps 6–7. Preserve
existing tool capabilities and limitations; shared HDF5 and new export features
remain deferred.

`codebase-3-3` is revised in the design, awaiting implementation and independent
verification. The two geometry-record proposals and carried scientific/deferred
prerequisites remain open. No source, numerical reference or test was changed;
no numerical run or independent review was performed for this revision.


## Retained optimizer electronic identities — approved 2026-11-16

The design assigns invocation-scoped internal IDs to distinct geometry/model
evaluations and reuses non-retired IDs for exact repeated calls. Saved evidence
must satisfy the current request; otherwise shared evaluation/iteration rules
apply. Separate electronic attempts preserve immutable settings and published
source maps. Old orbitals may seed a different geometry without proving its
convergence.

The run inventory selects one electronic calculation/attempt before work. Strict
selected-task recovery distinguishes null, published and missing/corrupt data;
no older-evaluation fallback or automatic optimizer geometry restoration is
introduced. Selected/actively needed and referenced archives survive retirement.
OPT-ELECTRONIC records deterministic and actual-binding tests, including geometry
changes/repeats, missing gradients, immutable sources and publication/interruption.
Optimizer history, Hessian/line-search recovery and D01/D02 remain deferred.

`readiness-3-1` is revised in the design, awaiting implementation and independent
verification. The standalone finite-difference plan-storage finding and carried
scientific/deferred prerequisites remain open. No source or numerical reference
was changed; no numerical tests or independent review ran for this revision.


## Finite-difference plan storage — approved 2026-11-16

The design selects one authoritative versioned file,
`checkpoints/finite_difference_plan.json`, with checked same-filesystem atomic
replacement. It holds the current revision, allocation counter, base geometry,
formula/displacements/request and point/source assignments. The run inventory
records only that this fixed file is required. Child checkpoints remain the
source of numerical progress and convergence evidence.

The contract specifies first registration, failed/stale updates, interruption
before/after replacement and child publication before source binding. Restart
and retention read the same plan; no alternate file or directory reconstruction
is permitted. FD-PLAN/D03 maps these cases to the shared record/retention tests
and deterministic workflow providers. New numerical-gradient science and
optimizer restart remain deferred.

`readiness-3-2` is revised in the design. All five round-3 proposals now have
approved revisions; they await independent verification and implementation.
The carried scientific-definition and deferred CC2/6D prerequisites remain open.
The additive shared measurement/controller/checkpoint work can start; this is
not full-engine acceptance or a clean independent review. No source, numerical
reference or test implementation was changed by this documentation revision.


## Independent verification — round 4, 2026-11-16

Three fresh reviewers read the complete unchanged 3,323-line draft (implementation
plan 2487–3323). Four revised areas passed: seed/protocol entry, active regression
migration and alpha curation, native inspection migration, and finite-difference
plan storage. The optimizer's electronic identity/reuse rules are specified, but
all reviewers identified the same missing result-to-electronic-attempt binding.
Three raw Major observations consolidate to one finding, readiness-4-1. One
editorial Nit about repeated finite-difference prose is recorded without fix.

A focused proposed edit defines a calculation-level source-map attempt override
with root inheritance and corresponding old-result/newer-attempt reader tests.
It has not been applied. The round is not clean; the new issue at an applied
revision triggers the review checkpoint. No additional review is authorized.
The additive measurement/controller foundation is startable, while this binding
is required before affected output readers/optimizer integration. Scientific
approvals and numerical acceptance evidence remain at their existing gates.

See [the readiness report](implementation_readiness_report.md) and
[annotated finding](goldfish_review_report.md). No source/design changes, builds
or numerical tests were made during verification.

## Result source-map attempt binding — approved 2026-11-16

The author approved readiness-4-1's focused clarification. The root results
attempt identifies the run writing the snapshot; a calculation-level attempt
identifies its immutable source-map owner, inheriting the root only when omitted
and identical. Numerical producing attempts remain in the map's state entries.
Missing/invalid/mismatched maps fail without a current-selector fallback.

The design, record examples and RESULT-SOURCES/C03/C04 cases now specify older
results alongside newer unfinished and concluded electronic attempts, mixed
owners, reused older samples and exact retention. The serializer/shared/native
reader plans reference these cases. Optimizer restart remains deferred. No
source, numerical reference or test implementation was changed.

The clarification is applied and awaits the separately authorized fresh
verification. This note does not claim that verification has already passed.


## Result-attempt clarification independently verified — round 5, 2026-11-16

The author approved applying the calculation-level attempt reference and then
verifying it independently. The design, examples and RESULT-SOURCES/C03/C04
planned tests were updated. Three fresh reviewers read the unchanged revised
draft (3348 lines; plan 2512–3348) and all verified the binding. readiness-4-1
is closed at the design level; implementation is still pending.

Two technical reviewers identified one separate Major, readiness-5-1: the design
must define which old results protect numerical archives and how that protection
ends. Clarity reported no findings. The two raw Major observations count once.
The lookup/serializer tests and shared core can start; retention-dependent cleanup
and old-snapshot native loading remain gated. No retention policy was applied.

The fifth-round cap is reached, with no automatic further review authorized.
Scientific approvals, numerical acceptance and deferred work retain their scopes.
No code, numerical reference, executable test, build or numerical run changed.
See [verification](implementation_readiness_report.md) and
[the remaining proposal](goldfish_review_report.md).


## 2026-11-16 — retention settled; shared-core implementation authorized

Applied readiness-5-1's bounded live-reference policy to design section 5,
source-map semantics, step-4 implementation work and RESULT-SOURCES/C02 tests.
Current results, latest non-retired tasks and active dependencies retain their
recursive sources. Historical metadata retains provenance without permanently
pinning payloads. Publication failures or incomplete reference checks prohibit
cleanup. No additional review round or optimizer restart work is introduced.
The user authorized starting the shared core; scientific definitions and engine
migration/numerical acceptance gates remain in place.

### First core increment

Started the scalar contract portions of implementation steps 2–3: typed
application-owned sample/measurement records; shared metric assessment; explicit
criteria resolution, intermediate floors and final limits; live decisions and
per-task restart planning. New APIs live in `madness::iteration` alongside legacy
callers until engine migrations. CORE-METRICS and CORE-RESTART curate these tests.
Record writers, native publication/retention integration, timing collection,
numerical helper adaptation and engine activation remain subsequent work. No
scientific definitions or absolute defaults have been chosen by this increment.

First-increment verification (2026-11-16): Ninja Release build in
`build/simple-policy`, GNU 13.3/C++20 on aarch64, MPI enabled and GENTENSOR off.
`madqc` and all five focused policy targets built; their five CTests passed.
Both new tests also compiled and ran without MPI/library linkage. Clang-tidy
reported no diagnostics in changed code; independent code review found no
Critical/Important issues or Minor requiring a change.

The unchanged `scf_he_hf` and `response_he_alpha` CTests passed in 8.87 s and
17.63 s with seven MADNESS workers and sequential BLAS/OpenMP settings. These
check existing application compatibility, not adoption of the new controller or
a performance comparison. All eleven protected alpha input/reference files
remain unchanged. No engine was switched to the new API. The additive interface
will coexist with the legacy path until each owning application migration.


### Native publication and retention increment — 2026-11-16

Implemented additive version-2 task evidence and a run-level checkpoint store.
Applications supply native archive callbacks and retain ownership of numerical
state and scientific eligibility. The store registers tasks before writes,
binds manifests to evaluated samples, validates native archives collectively,
and coordinates current-result source maps, dependency closure and retirement.
It preserves small manifests/publication receipts while reclaiming unreferenced
payloads. Missing scientific output disables cleanup without blocking numerical
restart. Required workflow plans supply validated live roots through callbacks.
CORE-CHECKPOINTS curates the storage tests and their limits.

Independent code review found two storage safety gaps, both fixed with regression
evidence: full checksums now prevent corrupt replacements from releasing valid
sources, and historical publication receipts prevent failed-save manifests from
authorizing reuse. Seven local assertions failed before the fixes; the repaired
suite and an explicit corrupt-B/valid-A preservation case pass. No Minor findings
or additional independent review round. The recorded scope exclusions preserve
the later application/workflow/scientific acceptance gates.

Verification: all eight requested build targets, including `madqc`, built in
`build/simple-policy`; 9/9 focused CTests passed, including native serial and
MPI2 cases. Clang-tidy completed with no reported warnings in changed code.
The unchanged SCF helium and response helium alpha smoke cases passed in 8.90 s
and 17.83 s with seven MADNESS threads and sequential BLAS/OpenMP settings.
All eleven protected alpha input/reference files remain unchanged. These runs
check compatibility, not production checkpoint performance or engine adoption.

This completes the publication/retention core increment, not all of step 4.
Stage and payload-free completion writers, full property formatting/readers,
workflow plan integration, timing, engine adapters and scientific validation
remain. Existing engine publication paths stay active until their owning
migrations. No optimizer restart redesign, HDF5 migration or numerical-reference
changes were introduced. Detailed evidence is in
`.superpowers/sdd/design_draft/checkpoint-verification.json` and adjacent logs.


### Record names and iteration identity — 2026-11-17

Renamed the shared records to `IterationInfo`, `MetricValue`, `StateMetrics`
and `IterationMetrics`. An iterative record uses task, producing attempt, actual
protocol and completed-update number; iteration zero is the prepared initial or
restarted state. Checkpoints store the number once as `updates_this_attempt`;
comparisons and result source maps carry it explicitly with the protocol.
There is no independent `sample_id`. Rechecking saved metrics preserves their
original identity. An attempt distinguishes solve sessions whose counters reset;
its ordinal has no numerical meaning, and detailed histories stay separate.

The VBC example uses the existing direct-task contract: exact B/C generations,
construction settings, producing attempt and native payload, without an iteration
block or artificial residual. Molresponse still owns dependency eligibility.
This updates the core and its examples; it does not activate new engine adapters.

The scoped review found a comparison-validation gap for unavailable metrics,
fixed with failing-then-passing regressions. Additional cases reject future
predecessors in the same session/protocol, stale live metrics and mismatched
source identities. Source-map tests now distinguish identity errors from
immutable-map rejection, and obsolete example identifiers were removed.

Verification: madqc and all focused targets built; 9/9 focused CTests passed,
including native serial/MPI2 checks. Clang-tidy reported no diagnostics in changed
code. Existing helium SCF and response alpha smoke cases passed in 8.89 s and
18.00 s. All eleven protected alpha input/reference hashes remain unchanged;
seven JSON examples parse. Detailed evidence is in
`.superpowers/sdd/design_draft/record-names-verification.json` and adjacent logs.
Broader completion writers, engine migration and workflow work remain deferred
while the author reviews the core.


### Refined 2PA frequency continuation — 2026-11-17

The author approved retaining saved response functions and their requested saved
protocol when refinement updates the photon frequency for the same tracked root
and perturbation. Other scientific settings, including the ground reference,
must match. Molresponse establishes the association; the shared planner receives
a qualified continuation source. The new equation requires fresh measurements
and comparisons, empty KAIN and a fresh allowance. Old completion claims do not
advance the updated task, and no frequency-distance input is introduced.

New/unmatched states and removed protocols retain the original fallback rules.
The new frequency and producing identity are recorded honestly; old checkpoints
and result sources are not relabeled. Affected 2PA sources/tensors are recomputed
and publication/retention protects current results and active continuation.

Updated the contract, examples, implementation-plan rows and acceptance/test
index (RESPONSE-2PA-CONTINUATION). This is a documentation change; no solver
code, numerical reference or engine-adoption status changed. The next code
walkthrough is checkpoint publication: matching functions and metrics, collective
validation, selecting a completed generation and interruption/replacement cases.


### Response replacement and task retirement — 2026-11-17

The author approved retaining current response task IDs, recording which saved
response each refined response replaces, and retiring the obsolete task once
nothing still needs it. The design uses optional `replaces` references with the
existing calculation/task/generation shape. Temporary continuation inputs keep
needed old functions alive; replacement history alone is not a numerical
dependency. Current results, other active consumers and requests must release
the old task before persisted retirement and archive cleanup. The containing
response calculation remains active.

Source review found that `fd_node_id` uses a five-decimal frequency key. The
contract therefore checks the actual frequency/context and handles same-ID
refinements with new generations, never by retiring the continuing task. No
task-ID format, user control or scientific default changes.

The publication walkthrough reran the existing serial and MPI2 checkpoint-store
tests successfully. Those tests cover existing publication behavior; the current
store still has only calculation-level retirement. Replacement metadata and
task-level retirement remain planned work, now explicit in the implementation
plan and acceptance/test index. No solver code changed in this design update.
Next review: scientific compatibility before reuse, including actual frequency,
ground-reference identity, protocol and operator settings.


## Simple restart checks — 2026-11-18

The author approved keeping restart to three checks: validate the selected
checkpoint, compare saved/current calculation settings, then assess valid saved
metrics against the current request. With unchanged settings, qualifying final
metrics permit reuse; unmet or missing evidence continues unfinished work while
preserving protocol progress. Changed settings require fresh measurements, using
saved functions where the application supports continuation. Stopping tolerances
and iteration limits are separate from the settings defining the calculation.

The contract, functional walkthrough, A04 acceptance coverage and CORE-RESTART
test plan now state this explicitly. No new policy, compatibility input, inferred
equivalence or frequency search is introduced. Existing protocol, matched 2PA,
strict-failure and retention rules remain in force.

The current planner still accepts `SavedEvidenceUse` and compares the supplied
context ID; application adapters must derive that outcome from settings rather
than assign it freely. This update changes documentation only. Next code review:
the first evaluation after restart, especially rebuilding BSH work and handling
density, energy and kinetic changes when previous comparison functions are not
saved. Valid saved scalar comparisons remain usable for unchanged calculations;
loading functions must never manufacture a zero change.


## SCF final orbital preparation — 2026-11-18

The author approved keeping final rotations, truncation and normalization, then
evaluating the resulting orbitals before the definitive convergence decision and
checkpoint publication. Changes retain the preceding comparable iteration as
the baseline. A post-rotation failure continues with the fresh residual when
updates remain, resetting KAIN if the orbital basis changed; at the final cap it
saves a valid unconverged endpoint. Intermediate caps keep their existing rule.
No further rotation occurs after the authoritative evaluation. The extra
evaluation contributes time but does not consume a main update.

Source review found a restricted-HF saved-state evaluator already used by the
legacy checkpoint adapter. Its measurements do not drive the SCF loop, and its
closed-shell all-electron HF restriction prevents treating it as general SCF
adoption. The implementation plan now requires generalizing that evaluation and
using the shared controller before publication. A02/A05 and CORE-METRICS record
post-rotation pass/failure, comparison-history and archive agreement coverage.

This turn updates documentation only. The previous restart walkthrough ran the
two focused metric/restart tests successfully; those scalar tests do not prove
the new SCF boundary. Next code review: protocol transitions, including function
conversion, operator rebuilding, cached BSH work and comparison/KAIN history.


## Protocol transition history — 2026-11-18

The author approved preserving the existing SCF transition behavior: carry
orbitals forward through projection/preparation, rebuild operators, discard old
raw BSH work and reset KAIN and comparison history at every protocol change.
This applies to both changed-k and fixed-k threshold transitions. Start a fresh
per-protocol update allowance; preparation and the first evaluation consume no
main update. Old metrics remain historical evidence, and the first new-protocol
change measurements remain unavailable until a valid new comparison exists.

The design records the automatic 1e-4/k=6 to 1e-6/k=8 example, explicit-k
precedence and independent smoothing settings. The functional walkthrough,
implementation plan and A03/A05/A06 test coverage agree. SCF already projects
orbitals, recreates operators and allocates history inside each solve; adoption
still must align the first current-state evaluation and common records.
No source or numerical-reference changes were made.

Next code review: the SCF-to-response prerequisite check. Trace how final
measurements and exact source generations qualify a ground-state reference,
while intermediate completion and an unconverged final checkpoint remain
resumable without releasing dependent response work. Then follow the same
source checks into molresponse's direct VBC tasks.


## Exact qualified inputs for VBC and properties — 2026-11-18

The author approved using exact published checkpoints qualified against current
final-protocol requirements for VBC and property computations. Molresponse owns
readiness, invalidation and selective rebuilding; the shared controller assesses
metrics and the store validates source bindings/publication and retains inputs.
VBC records completion plus B/C/ground-reference sources. Property records bind
all source generations their formulas use, preserving earlier result identities.

Source review found that load_fd_as_xy uses try_load_fd_state, which can select
coarser or unfinished guesses. Direct consumers must load their qualified exact
sources instead. The existing DAG and VBC invalidation remain useful; cap
acceptance or an old converged flag cannot satisfy the new requirements. Only-A
refinement can preserve VBC; B/C changes rebuild affected work; unchanged inputs
that already satisfy tighter limits permit reuse.

Updated the contract, functional walkthrough, response migration plan and B02
coverage; RESPONSE-DIRECT-INPUTS preserves the test purpose independently of
format changes. No source or numerical-reference changes were made. Next code
review: partial property publication, availability and exact state metrics in
results.json, keeping unconverged observations in analysis.

## Workflow completion gate — implementation, 2026-11-18

After reviewing partial results, tighter-request reuse and controlled stops, the
author approved implementing the workflow completion boundary. An application
can return a small complete/reason outcome; SinglePointDriver records it and
blocks context publication when incomplete. Workflow preserves that task's
summary before reporting a controlled stop. Numerical and publication exceptions
retain their failure path. SCF/Nemo use a transitional adapter to existing
engine stopping/restart rules; other engines adopt the interface in their phases.

The independent review resolved cached/reprojected references, stale Nemo
completion, same-instance geometry caching, changed-Hamiltonian read-only reuse,
invalid precision stamps, capped frozen-virtual continuation and report write
failures. These changes add no input controls, native archive fields or new
numerical stopping criteria. WORKFLOW-COMPLETION records the behavioral checks.
Final build/tests are recorded separately in the implementation ledger. The new
shared metric conventions, result schema, exact checkpoint handoff and molresponse
property-completion adapter remain migration work.

## SCF evaluation/update split — implementation, 2026-11-18

The SCF application owns the orbitals, comparison densities and KAIN history.
`evaluate_orbitals` computes densities, potentials, Fock matrices, energy and
kinetic quantities; `evaluate_residual` computes the raw BSH residual once for
the prepared state. `update_subspace` consumes that residual only when a main
update is needed, retaining existing KAIN, restriction and orthonormalization.

At a proposed stop, SCF performs final orbital preparation and evaluates the
result again against the preceding iteration. Failed checks resume with the new
residual when updates remain; exhaustion saves an unconverged evaluated endpoint.
No subsequent orbital rotation precedes native saving. Final checkpoint
publication copies the deciding evidence and its engine verdict.

This increment preserves existing convergence criteria, inputs and native archive
layout. The old named adapter remains until migration; its histories retain
original labels on read. New shared input/output adoption and Nemo scientific
adapters are still separate work. See SCF-EVALUATION in the test index and the
increment's verification record for executed coverage and limitations.

## SCF shared metric adapter — implementation, 2026-11-18

Ordinary SCF now has a programmatic shared-controller solve entry point. The
application keeps its density history and cached BSH work, uses shared scalar
helpers, and supplies one `ground` record containing residual, density, energy
and physical kinetic changes. Missing comparisons remain unavailable. The
collective observer receives the definitive `IterationMetrics`, `Decision` and
update counts for native checkpoint publication; it cannot coexist with the
legacy checkpoint observer in the same solve.

SCF supplies the existing threshold factors for intermediate residual and
per-spin density tolerances. Energy and kinetic changes have no invented floor.
All absolute defaults remain unset for this validation increment. Final requests
are not relaxed. Native headers retain conservative legacy precision stamps;
the typed metrics are the authority on the shared path, so an old reader cannot
mistake a residual-only completion for its own density-based convergence.

The CLI and its legacy workflow restart path remain unchanged. Their coordinated
input, saved-evidence and result migration follows this adapter validation;
frozen-only SCF and Nemo still need their adapters. SCF-SHARED-METRICS in the test
index records the purpose and limits of the new coverage.


### SCF default and input prerequisites — 2026-11-18

The increment validated the ordinary SCF adapter's initial absolute defaults
(residual/density `1e-4`, energy/optional kinetic `1e-5`) against curated decks,
including tighter comparisons and an intentionally incomplete coarse Li case.
The metric validation document and test index record the scientific scope.
This supersedes the earlier adapter entry's deliberately unset defaults.

The shared parser retains every human CLI occurrence and registered readers
validate discarded payloads without applying them. File/CLI error handling is
collective; last-whole-group precedence and path case are preserved. Standalone
owners and global retired-key cleanup still require migration. Optimizer gradient
defaults are now independent of dconv; explicit settings and MolOpt's precision
clamp remain. The public SCF input and native workflow-restart switch stays
coordinated with both SCF paths, avoiding a second stopping-policy option.

Final verification: 25 targeted tests pass, all 55 curated deck/reference files
remain unchanged, and the four Important code-review findings have failing and
passing regression evidence. Registered reads now treat file plus CLI input as
one collective transaction. Libxc-dependent test cases are conditional, and an
invalid-state protocol stop cannot be hidden by later success. A permissive
reader may still print an unknown-key warning twice; that minor diagnostic
cleanup remains for owner migration. Public input/restart activation is next.

### Frozen-occupied SCF shared adapter — 2026-11-18

The numerical adapter now passes one `virtuals` subject to the shared controller,
with the maximum raw residual and component evidence for every active alpha/beta
virtual function. It canonicalizes and reevaluates before recording metrics;
the same BSH work feeds a subsequent KAIN update. Each spin keeps its own KAIN
history. Fixed occupied density, energy and kinetic changes are inapplicable.
The occupied energy remains a quantity of the exact source reference.

The typed entry validates its protocol, subject, applicable criteria and budget.
Its observer uses the existing native checkpoint interface, including source
references. A caller must first qualify the occupied source; that workflow check
is not replaced by legacy archive stamps. The legacy entry still uses the same
numerical kernels with its existing policy until the coordinated public input
and restart cutover. SCF-FROZEN-METRICS records the new test's scientific limits;
virtual preparation metadata and broad spectrum/default validation remain open.

The independent review added three concrete fixes: detached occupied work and
preserved loaded prefixes for frozen orthogonalization; complete-array/nonempty
preflight before modifying prior outcome; and release of completed legacy spin
scratch before the next spin. Focused regression evidence is retained with
this increment. A direct comparison of every native round-trip
component and loaded virtual quantity before reevaluation remains a minor test
extension for the checkpoint workflow work.
