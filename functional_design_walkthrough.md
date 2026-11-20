# Functional design and user controls

Ownership approved; illustrative interfaces — 2026-11-11. Read alongside the
[design draft](design_draft.md#7-reuse-the-existing-drivers).
This illustrates the agreed responsibilities, not implemented C++ signatures
or a new solver framework. Persisted records follow the accepted design contract.

The application owns state and comparison history and calls shared numerical
helpers to compute measurements. The shared controller receives measurements,
decides and records. Users choose calculation settings; they do not configure callbacks,
record identities or serialization fields.

## What an application supplies

A callback is an application function passed to the common code. It can use the
application's existing state and capture its reference, operators and an
explicit `World&`. It does not require all applications to use one state type.

| Operation | Example for SCF | Example for molresponse |
| --- | --- | --- |
| Prepare the current state | Required orbital preparation and rotations | Prepare X/Y blocks or the coupled root bundle, including root alignment |
| Build the raw fixed-point proposal | Apply the SCF equation to the current orbitals | Apply the response equation to the current state or bundle |
| Supply density fields | Construct physical spin densities | Construct the induced or transition density for the model |
| Supply energy | Total energy of the current prepared state | Excitation energy for a root; explicitly inapplicable for a driven state |
| Supply kinetic inputs | Physical orbital blocks, occupations and derivative settings | Actual stored blocks, unit weights and derivative settings |
| Advance the state | Reuse proposal/residual for KAIN or a plain update, restriction and constraints | Reuse equation work in the existing response update, preserving its numerical ordering |
| Read/write numerical state | Existing orbital serialization | Existing state, bundle or VBC serialization |

Definitions accompany those operations once: which fields/blocks are supplied,
what they mean, and their weights, units and normalization. Applications call
shared numerical helpers for the agreed norms, reductions, differences and
scales. They pass scalar measurements and comparison evidence to the controller,
which applies tolerance and stopping rules without receiving solver functions.

A measurement callback can capture the application's state and history. An
application with its own loop can submit the same measurements directly; both
use the same controller. Existing computed quantities can be reused when their
definition and sampling match. Adapter tests and scientific review establish
that the physical inputs and comparisons follow the declared conventions.

## What happens during an iteration

At a protocol change, the application carries its functions forward, converts
them to the requested threshold/k and rebuilds the required operators. It clears
cached raw proposals/residuals, KAIN and convergence-comparison history. The new
protocol starts with a fresh update allowance and no measured changes until a
new baseline exists; old measurements remain historical records. This also
applies when k stays fixed. See the [protocol-change rule](design_draft.md#protocol-changes).

Conceptual pseudocode for a single protocol; names show responsibilities, not
new public API signatures. The outer lifecycle has already resolved settings,
validated inputs, planned restart, prepared the protocol and reset KAIN when
required. Each application retains its own state type and numerical workspace.

```text
repeat:
    app.prepare()
    measurements = app.evaluate()          # calls shared metric helpers
    decision = shared.decide(measurements, criteria, protocol, update_count)

    if app.needs_endpoint_preparation(decision):  # SCF final rotations
        app.prepare_endpoint()
        measurements = app.evaluate()      # keep the preceding iteration's baseline
        decision = shared.decide(measurements, criteria, protocol, update_count)

    shared.record(measurements, decision)

    if shared.checkpoint_due(decision, update_count):
        shared.publish(measurements, decision, app.write_evaluated_state)

    if decision ends this protocol or task:
        return decision

    app.retain_comparison_and_advance()
    update_count += 1
```

`app.evaluate` builds the raw proposal, forms or consumes the raw residual and
uses shared helpers to measure the prepared state. It retains numerical work
for the update; the controller receives scalar evidence and identities. A
successful convergence check needs no extra KAIN update. At the update cap,
the loop evaluates the produced state and exits without another update.

**SCF endpoint clarification — 2026-11-18.** The optional preparation branch
runs once at a prospective endpoint, before the definitive outcome is recorded.
Keep the existing final rotation/truncation/normalization, then evaluate those
orbitals and use the new decision. A failed check with allowance remaining
continues using the new residual; reset KAIN if final preparation changed the
orbital basis. An exhausted final allowance saves an unconverged endpoint.
No second rotation follows the definitive check, and evaluation does not consume
an update. This is conceptual application logic, not a new required callback
for every engine. See the [complete rule](design_draft.md#scf-final-orbital-preparation).

The application retains only the comparison data it needs: previous density
fields, energy/kinetic scalars, and root functions needed for correspondence.
The BSH residual needs no previous iterate. Missing valid comparison data gives
an unavailable change and a new baseline under the existing restart rules.

Evaluation preserves prepared functions and accelerator history. Retaining
comparison data requires deliberate ownership: copying a MADNESS function
handle is not a deep copy. Checkpoint callbacks must save the evaluated state
before it changes. State, history, raw work and KAIN remain in application scope.

SCF can keep its existing loop and call these shared operations at these points.
Response can adapt its existing shared driver. This example does not require
moving all engines into one loop immediately. The recorder creates the detailed
metric/checkpoint metadata; applications do not hand-build its JSON.

## What app.evaluate returns: SCF and an excited root

**Author-approved return structure — 2026-11-11.** Return an evaluation containing
one prepared sample identity and a collection of subject measurements. SCF has
one subject; a coupled excited-state solve returns all required roots together.
The following notation illustrates typed values, not JSON or final C++ syntax.
The numbers are invented examples, not tolerances or numerical reference results.

Use `measured(value, scale)` for a valid metric. Its scale is the agreed relative
scale, not an effective tolerance. It is needed for a selected relative test;
the examples show it for every metric to make the meaning visible. An
absolute-only test does not require a scale. Use `unavailable(reason)`,
`not_applicable(reason)` or `invalid(reason)` when there is no measured value.
A valid measurement may have a separately unavailable scale; preserve the
measurement, and prevent a requested relative test from passing. The exact
representation of those alternatives remains an interface detail.

Model definitions, units, normalization and reduction are declared in the saved
task context and bound to the evaluation. They need not be returned as repeated
strings on every call. Each measurement still resolves to its exact definition.
Iteration references carry task, producing attempt, actual protocol and the
number of completed main updates. Zero means the prepared initial/restarted
state. The `info` field also carries the application-validated scientific
context. There is no separate sample label. The common `comparison` below
applies to change metrics; BSH needs no preceding iteration. When metrics use
different comparable predecessors, preserve their individual references.

### SCF: one ground-state subject

SCF prepares its orbitals, builds the raw fixed-point proposal and retains the
residual functions for a possible update. It constructs the current physical
spin densities and obtains total and kinetic energies for the same prepared
orbitals. Shared helpers use SCF's retained previous fields/scalars to compute:

```text
IterationMetrics {
    info: {task_id: scf, attempt_id: a1, protocol: {threshold: 1e-6, k: 8}, iteration: 12, context_id: scf_context},
    states: [
        StateMetrics {
            id: ground,
            comparison: {task_id: scf, attempt_id: a1, protocol: {threshold: 1e-6, k: 8}, iteration: 11, subject: ground},
            metrics: {
                bsh_residual:   measured(8e-6, scale=1.0),
                density_change: measured(3e-6, scale=0.8),
                energy_change:  measured(5e-7, scale=75.1234565),
                kinetic_change: measured(2e-6, scale=74.500002)
            }
        }
    ]
}
```

Here the energy comparison is `-75.1234565` versus `-75.1234560` hartree;
the kinetic comparison is `74.500002` versus `74.500000` hartree. The helpers
retain these supporting scalars with the measurement evidence. Density change
is the maximum spin-channel L2 difference; its scale is the largest spin-density
norm over the two samples. The residual and its scale use maxima over the actual
orbital blocks. The controller does not need those functions.

### Excited states: one root inside a bundle

Molresponse prepares the coupled bundle and validates root correspondence.
For this TDA example, root 2 has a unique match in the previous bundle sample;
the transition density is phase-aligned before comparing fields. Energy means
excitation energy, and kinetic change uses the agreed positive response kinetic
sum. The displayed root is one entry in the bundle return:

```text
IterationMetrics {
    info: {task_id: excited_states, attempt_id: a1, protocol: {threshold: 1e-6, k: 8}, iteration: 12, context_id: response_context},
    states: [
        ...,
        StateMetrics {
            id: root_2,
            comparison: {task_id: excited_states, attempt_id: a1, protocol: {threshold: 1e-6, k: 8}, iteration: 11, subject: root_2},
            metrics: {
                bsh_residual:   measured(4e-6, scale=0.72),
                density_change: measured(2e-6, scale=0.12),
                energy_change:  measured(1e-6, scale=0.350001),
                kinetic_change: measured(3e-6, scale=1.400003)
            }
        },
        ...
    ]
}
```

The excitation energies are `0.350001` and `0.350000` hartree; the kinetic sums
are `1.400003` and `1.400000` in the declared model units. BSH still reports the
maximum raw per-function residual. Its scale is the largest current function
norm, not the excitation energy or the entire bundle norm. The density is the
model's transition density rather than an SCF ground-state density.

Molresponse returns measurements for all required roots at that bundle sample.
Application validity checks establish a usable bundle, and the shared layer
checks its expected subjects/context. Every root must pass its selected metrics;
there is no additional bundle convergence metric. A missing root cannot complete. A passing root does not authorize a separate
update, checkpoint or freeze; the bundle remains the task being iterated.

### Missing comparisons and supporting evidence

On the first SCF evaluation with no valid predecessor, BSH can be measured,
while density, energy and kinetic changes return
`unavailable(no_previous_sample)`. Current energies and kinetic quantities
can still be reported as supporting observations and retained for the next
comparison. There is no numerical zero standing in for an unavailable change.

If excited-root matching becomes ambiguous, retain a valid current BSH metric
but return `unavailable(root_correspondence_ambiguous)` for the affected changes.
Do not claim a previous comparison identity; begin a new baseline for that root.
An invalid current state can invalidate BSH too, and is a separate case from an
otherwise valid state with an ambiguous predecessor.

The compact displays omit supporting evidence already required by the draft:
per-function residual norms and the limiting component, underlying energy and
kinetic scalars, and valid comparison identities. Shared helper results can
carry this scalar evidence alongside the main values, so the recorder does not
need to access functions or rerun expensive computations. Detailed root-matching
diagnostics remain in analysis. The return contains no distributed functions,
raw proposal, KAIN subspace, checkpoint archive or previous numerical state.

The controller adds the selected criteria, requested/resolved tolerances,
effective limits and decisions from its settings and progress context. The
recorder adds the common identity/definition references to the stored record.
`app.evaluate()` does not return a `converged` flag or decide which metric is
required. Direct submission and callback invocation use the same handoff.
The evaluation/subject/measurement structure is approved. Concrete C++ types
and the representation of supporting evidence remain implementation details
within the accepted section 6 record contract. The controller reply below is
also approved.

## What the controller returns

**Author-approved reply — 2026-11-11.** Use a decision tied to the evaluated sample:
`info`, `action`, `reason`, `checkpoint` and per-subject `checks`.
The [main design](design_draft.md#controller-reply) defines the approved
fields and their relation to existing completion/checkpoint records.

| Situation | Action | Example reason | Checkpoint request for a valid sample |
| --- | --- | --- | --- |
| Criteria remain unmet and updates remain | `update` | `criteria_unmet` | Only if optional interval saving is due |
| Required changes lack a comparison and updates remain | `update` | `comparison_unavailable` | Only if optional interval saving is due |
| Intermediate targets pass | `finish_protocol` | `intermediate_criteria_met` | Yes |
| Intermediate cap reached with unmet targets | `finish_protocol` | `iteration_limit` | Yes |
| Newly evaluated final sample passes for all required subjects | `complete` | `final_criteria_met` | Yes |
| Valid saved final sample satisfies the current request | `complete` | `saved_metrics_satisfy_request` | No; retain the existing generation |
| Final cap reached with unmet targets | `stop` | `iteration_limit` | Yes |

At an intermediate endpoint, `finish_protocol` returns control to the existing
protocol manager; it does not bypass the barrier for other participating tasks.
For excited states, one decision covers the bundle. Root checks explain which
criteria are met or unavailable without separate root actions or freezing.

`checkpoint` requests publication before state mutation. It does not mean a
save succeeded; existing validation and publication rules still apply. Invalid
current numerical state cannot replace a usable checkpoint. Solver failures
stop affected work through the existing failure path and are distinct from
missing first comparisons or an iteration cap.

The application executes the action using its own state and numerical methods.
The controller's `checks` supply limits and reasons to current run reports and
analysis. Checkpoint metric records retain measured evidence without those tests.
The numerical examples above do not choose tolerances or declare convergence
until this policy evaluation takes place.

## Restart reassesses saved measurements

**Author-approved — 2026-11-11; simplified 2026-11-18.** Validate the selected
checkpoint, then let the application compare its saved and current calculation
settings. Unchanged settings allow the controller to check saved measurements
against the current request. If they pass at the requested final protocol
(threshold and k), reuse without calling `app.evaluate()`. Otherwise continue
the unfinished work while preserving its protocol progress.

Changed geometry, model/reference, actual response frequency or numerical
operator/smoothing settings require fresh measurements. Use the saved functions
where the application supports continuation; rebuild the affected work and
follow the existing protocol/conversion rules. Matching metric definitions are
also required for saved-evidence reuse. Stopping tolerances and iteration limits
are separate from these settings. These comparisons determine the action; no
discretionary compatibility policy, inferred equivalence or new input is needed.

A saved BSH residual of `3e-7` passes a new absolute tolerance of `1e-6` but fails
`1e-8`. Every other selected applicable criterion must also pass. A tighter
request does not change the original measured sample or its producing attempt.
The controller returns `complete` with `checkpoint: false` for eligible final
reuse, and records the new assessment without rewriting numerical archives.

Store diagnostic measurements and available scales with their scientific
context. A new relative criterion can use a saved scale or sufficient scalar
evidence to derive it. Missing required measurements/comparison evidence need
fresh numerical work; an unsupported requested metric is an input error.
An old pass/fail flag is never the authority for a new request.

Protocol completion and iteration counts remain progress metadata. Previous
comparison functions need not remain on disk merely to reassess valid saved
metric values. Unchanged source generations can keep their VBC/property records
eligible; changed source states trigger the existing dependency invalidation.

## Read-only inspection and plotting

**Author-approved — 2026-11-12.** Reading existing results or plotting saved
functions does not require satisfying a new convergence request. Check that the
data are readable and suitable for the requested operation. An unfinished
checkpoint is usable for plotting; display its saved protocol, metrics when
available and completion status. Keep the original scientific context.

Read-only can produce plot/export files, but does not call the solver controller
to satisfy new tolerances, iterate, advance protocols, rewrite solver checkpoints
or change convergence claims. A normal solver restart instead reuses qualifying
saved work or resumes the solve when the current request is unmet. "Requested
final protocol" means the final truncation threshold and polynomial order k.

## Each task keeps its protocol progress

**Author-approved — 2026-11-11.** Match saved and requested protocols by threshold
and polynomial order. Preserve each independent state's compatible saved
position, even when another state starts earlier. Tightening stopping tolerances
alone never sends a matching later state through coarse protocols again.

For P1 → P2 → P3, the restart plan can be:

| State | Remaining work |
| --- | --- |
| Saved at P3 and current criteria pass | Reuse without evaluation |
| Saved at P3 and current criteria fail | Resume P3 |
| Interrupted at P2 | Resume P2, then P3 |
| Eligible to advance after completed P2 work | P3 |
| New independent state | P1, then the remaining sequence |
| Matched 2PA response at P3, frequency refined | Keep its functions and P3; evaluate the updated equation with fresh convergence evidence |
| Compatible saved representation absent from the sequence | Supported conversion as a guess at P1, then the remaining sequence |

Only tasks needing a protocol participate in its scheduling wave/barrier. Later
states can remain on disk until needed; they are not projected down to join an
earlier wave. A seed from a later state supplies a guess, not progress or
convergence evidence for the receiving task. Unsupported conversion and invalid
expected checkpoints retain the strict restart behavior.

**2PA refinement exception — approved 2026-11-17.** Molresponse can match the
previous response for the same tracked excited root and perturbation when only
the derived photon frequency changes. Preserve its saved protocol if requested,
but discard old convergence/comparison claims, including completed-protocol
status. Rebuild equation work at the new frequency, start KAIN empty and use a
fresh allowance. The controller receives only measurements of the new equation.
Keep the original checkpoint identity as provenance; the new task records its
actual frequency. No nearest-frequency search or new input knob is introduced.
Keep existing frequency-based task IDs and record the predecessor in `replaces`.
That relationship preserves history; pending inputs protect the old functions
only while continuation needs them. Retire an obsolete old task after current
results and all other active work release it, without retiring the calculation.
If rounded frequency keys yield the same task ID, save a new generation and keep
the task active. Actual frequency/context still determine whether fresh checks
are required.
An unmatched root/frequency starts first, and an absent saved protocol follows
the existing supported-conversion rule. See the
[full continuation rule](design_draft.md#continue-response-functions-after-a-2pa-frequency-refinement).

An excited-state bundle has one entry plan. Adding roots follows the existing
expanded-bundle restart rule at a matching saved protocol; new roots do not
independently start at P1. If the saved bundle protocol is absent, conversion
and first-protocol entry apply to the bundle as a whole. This preserves the
all-roots-active rule and adds no new scheduling layer or input knob.

## A direct task uses the same outer lifecycle

Molresponse first checks B/C and ground-reference eligibility, using the shared
assessment of current requirements. It binds their exact published generations,
loads those functions and builds VBC. The shared lifecycle validates the source
records, records completion and publishes through the existing serializer.
VBC's `inputs` name B, C and the reference; it has no iteration loop or residual.

Property assembly similarly qualifies A and the other required inputs and binds
their exact sources. Direct consumers cannot fall back to a coarser or unfinished
starting guess. Molresponse owns readiness and invalidation: changing only A
preserves an otherwise valid VBC, while changing B or C requires rebuilding it.
Tighter requests that all unchanged inputs already satisfy need no rebuild.
Workflow code keeps the moldft-to-molresponse dependency; existing result-source
records preserve each published property's original inputs.

## The input knobs that matter

Keep settings in the existing application groups and preserve supported CLI
forms. The explicit convergence names below are approved target syntax; the
checkout does not yet implement them. `dconv`/`econv` aliases are not accepted
by the new interface: explicit old keys must produce a migration error rather
than an automatic numeric translation. Distinct inner-solver controls keep their
reviewed scopes. Numbers are examples, not defaults. This is a response-group
fragment, not a complete runnable calculation.

```text
response
    maxiter 30
    convergence.criteria [bsh_residual, density_change, energy_change]
    convergence.bsh_residual 1e-5
    convergence.density_change 1e-6
    convergence.energy_change 1e-6
end
```

| Control | What the user is choosing |
| --- | --- |
| Existing protocol and polynomial-order settings | Which representations to calculate, including the mandatory final threshold and k. Preserve reviewed group ownership; this example does not add a response-specific protocol key. |
| `maxiter` | Maximum main updates at each executed protocol in the current attempt. One response value applies to driven states and excited bundles, each with its own allowance. Explicit `excited.maxiter` produces the normal fatal unknown-input error. The solve can finish earlier; restart gets a fresh allowance. |
| `convergence.criteria` | Complete selected metric list. With no explicit list, use the agreed application defaults. |
| `convergence.<metric>` | Requested absolute tolerance. Changing it does not activate an unselected metric. |
| `convergence.<metric>.relative` | Optional relative term; default zero. The shared scale rule supplies its meaning. |
| `restart auto/none/iterate`, existing `--continue` | Scoped reuse/resume, explicit fresh start, or require-existing reuse/resume. `--continue` requires a recognized run and rejects fresh/seed/inspection conflicts. All solver paths use the same scientific checks. |
| `checkpoint.every_iterations N` | Optionally save inside a protocol every N updates; endpoints remain mandatory. Default is endpoint-only. `save` and `checkpoint.keep_all` are absent inputs and produce generic unknown-input errors for any value. |

**Response budget input approved — 2026-11-12:** set `maxiter` in the response
group, or use the standalone response `--maxiter`. The old `excited.maxiter`
key produces only the normal fatal unknown-input error, even if both keys have
equal values. There is no custom migration message. The same cap applies
independently to each driven state or coupled excited bundle; the states do not consume one shared total allowance.

The same response criteria and tolerances apply to driven and excited states.
Energy change is applicable only to excited roots. Each state must pass every
selected applicable criterion. Kinetic change is reported where supported and
becomes required only when selected. Required unsupported metrics fail input
validation; a metric without a first comparable sample is temporarily unavailable.

At intermediate protocols, the shared resolver applies the agreed legacy
threshold-based absolute-tolerance adjustment. It reports the resulting limits;
users do not need a separate tolerance table for every protocol. Final checks use
the requested tolerances without that adjustment. Representation settings and
iteration criteria have distinct meanings: satisfying the latter establishes
iteration convergence, not a general error bound on requested properties.

Existing KAIN, step restriction and model/operator controls still belong to the
engine. They affect how the solve proceeds. They do not create a second
convergence policy. Warmup remains disabled and stall-based advancement deferred. Their old policy
controls are removed, including disabled values, under the ordinary unknown-input
rule in [input routing](design_draft.md#input-routing-and-retired-controls).

## What the user sees

Before solving, report resolved settings and their file/default/CLI sources:
selected metrics, their model meanings and units, requested tolerances, protocol
threshold/k, effective intermediate limits, update allowance, restart action and
checkpoint frequency. Use the preserved parameter/restart inspection commands.
Definitions should name the physical quantity, for example maximum spin-density
L2 change, excitation-energy change, or the positive response kinetic sum.

During solving, use a compact table. These rows illustrate a driven state at its
final protocol; model-specific unit labels and surrounding state/frequency/k
context are omitted only from this excerpt. Real reports must include them.

| Metric | Value | Limit | Use | Result |
| --- | --- | --- | --- | --- |
| Maximum per-function BSH residual | 3e-6 | 1e-5 | Required | Met |
| Induced-density L2 change | 8e-6 | 1e-6 | Required | Unmet |
| Energy change | — | — | Inapplicable | Driving frequency is fixed |
| Response kinetic-sum change | 2e-5 | — | Diagnostic | Reported |

The stopping explanation is: “Continue: density change exceeds its limit.”
At a first comparison it is “Continue: no comparable previous density.” At an
intermediate cap it is “Advance: update limit reached; intermediate targets
remain unmet.” At an unconverged final cap it is “Stopped: update limit reached;
final criteria remain unmet.” Each message comes from the shared decision.

The console table, analysis output and checkpoint evidence are different views
of the same evaluation and decision. `results.json` keeps grouped requested
properties and compact measured state metrics once per application calculation.
Shared states are not repeated for every tensor component or property; see the
[alpha/beta example](tensor_result_examples.md). Full iteration histories and
restart evidence remain separate. Normal output explains what was measured and
why execution continues, advances or stops. The consolidated [record contract](design_draft.md#6-records-timings-and-output)
now supplies persisted encodings. The [complete example](output_record_examples.md)
connects inputs, decisions, results, checkpoints, timings and restart behavior.

The in-memory `IterationMetrics`/`Decision` examples above are API pseudocode. Persisted
metric scales use the explicit status/value object in design section 6; no new
state ownership or controller policy follows from that JSON encoding.
