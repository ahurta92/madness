# Geometry optimizer replacement — follow-on requirements

**Status:** Deferred until workflow unification is complete, by author direction
on 2026-11-13. This is retained design work, not an implementation-ready plan or
an acceptance gate for the current unification work.

The intended direction is to replace the existing optimizer after applications
use the shared convergence, metrics, restart and checkpoint layer. Keep the
existing optimizer and supported `--optimize` combinations during migration.
Retire it only after a validated replacement, migrated coverage and collaborator
approval. The [current design](../../design_draft.md) defines the electronic
compatibility work that remains in scope.

Current unification work does add internal IDs for the retained optimizer's
electronic evaluations and a selector for their saved SCF/Nemo work. These keep
geometry/model context and immutable request/source records correct. They do
not identify an accepted optimizer step or restore its history; the selected
electronic geometry may be a trial. The current design's electronic reuse and
strict checkpoint rules apply, while the optimizer still starts through its
existing entry path and supplied geometry.

The following notes preserve the agreed behavior for a replacement: consistent
energy/gradient composition, saved geometry evaluation identities, resumable
Hessian and line-search history, and a fresh geometry-step allowance on restart.
The snapshot fields and save points were worked out against the existing line
search. Revisit that mechanism for the chosen replacement algorithm while
preserving the requirement that restart does not apply a committed transition
twice. These notes do not require retrofitting the existing optimizer now.

## Decisions to revisit before implementation

- Choose the replacement algorithm and supported gradient providers; validate
  state selection across geometries before enabling excited-state optimization.
- Revisit geometry criteria and precision estimates for the replacement algorithm.
  The approved migration removes the retained driver's `dconv` dependency using
  independent defaults `gtol=1e-4` and `gradient_precision=1e-5`, with explicit
  optimization-group overrides. Its other defaults and precision clamp remain.
  This preserves the former default settings; it does not validate a measured
  gradient-error estimate or choose the replacement's stopping policy.
- Resolve precision-limited stopping. The existing `molopt.h` can report success
  when the gradient is below its assumed precision even if energy/displacement
  criteria have not passed. This remains a known legacy behavior, not evidence
  that the replacement satisfies the stricter contract below.
- Revalidate acceptance cases D01 and D02 with the chosen algorithm. Numerical
  finite-difference gradients and excited-state gradients still need their own
  scientific implementation and validation.

## Geometry optimization across applications

The workflow supplies the optimizer with a consistent total energy and nuclear
gradient for each geometry and selected electronic state. One evaluation may
combine several applications. Each application owns its numerical computation,
state, metrics and checkpoints through the shared interfaces; the workflow
owns their dependency order and assembles the optimizer's evaluation.

For a TDHF excited state I at nuclear coordinates R:

```text
E_I(R) = E_SCF(R) + omega_I(R)
grad_R E_I(R) = grad_R E_SCF(R) + grad_R omega_I(R)
```

SCF/Nemo supplies the ground-state energy and reference. Response supplies the
excitation energy and the excited-state gradient calculation. The gradient
interface must declare whether it supplies the full total-energy gradient or
only the excitation-energy contribution. In the latter case, the workflow also
obtains the ground-state gradient from its owning application and adds it once.
An energy-only line-search evaluation uses the same total-energy definition and
target-state selection as an energy-and-gradient evaluation. Pairing the SCF
energy alone with an excited-state gradient is invalid. Use explicit units and
the nuclear-gradient sign convention; forces are its negative.

At each geometry, finish the ground state at its required final protocol before
response work. Molresponse owns any additional iterative or direct tasks needed
for its gradient. Before releasing an energy-gradient pair to the optimizer,
the workflow checks the required convergence/completion evidence, common geometry
and compatible model. Response energy and gradient must describe the same selected
state and depend on the exact ground-state reference generation used for the
SCF energy. It keeps references to the contributing
application records rather than duplicating their numerical archives. Restart
reassesses those dependencies before reusing an assembled evaluation; a changed
geometry or reference invalidates the affected evidence under the main design’s restart rules (section 4).
Optimizer history and geometry stopping criteria remain workflow-owned and
separate from electronic solver history and criteria.

This specifies the composition contract. The current madqc entry point rejects
`--optimize --wf=response`; enabling it requires a validated excited-state
gradient provider. The gradient formulation and state-selection/following policy
remain scientific decisions. Existing root labels alone do not establish
correspondence between geometries.

PySCF provides concrete references: [TDBase.e_tot](https://pyscf.org/_modules/pyscf/tdscf/rhf.html#TDBase)
adds SCF energy and excitation energy; [TDSCF_GradScanner](https://pyscf.org/_modules/pyscf/grad/tdrhf.html#TDSCF_GradScanner)
returns the selected total energy and gradient; the [excited-state optimization guide](https://pyscf.org/user/geomopt.html#excited-states)
shows a TDHF gradient scanner passed to an optimizer.

## Optimization evaluation identities

For geometry optimization, the workflow assigns a saved evaluation number to
each geometry evaluation. It is separate from the optimizer step, solver
iteration, protocol and execution attempt: a line search can evaluate several
geometries during one optimizer step. Use that number directly in calculation
names, for example:

```text
checkpoints/opt_eval_0007_ground_state/
checkpoints/opt_eval_0007_response/
```

The application-owned tasks and generations beneath these locations keep the
shared layout. The optimizer checkpoint records the current and any pending
evaluation numbers, their exact coordinates with atom order and units, and the
next unused number. Publish a new pending evaluation and advanced counter before
starting its electronic calculations. Child task registration then follows
the main design’s checkpoint contract (section 5). Restart restores the saved numbers; it does not reconstruct them
from loop counters, directory order or rounded coordinates. Rejected trials
consume their numbers permanently.

An energy-only trial and a later gradient evaluation at that same selected
geometry can share an evaluation number. If the line search selects a different
geometry, allocate a new number before computing there. This is matching of
current and pending work, without a search through all historical geometries.
Protocol/tolerance changes preserve the saved evaluation number and trigger the
usual compatibility and convergence assessments. A new geometry gets a new
evaluation number; old functions may seed it, but old convergence evidence does
not certify it. Changing the target state still requires scientific validation.

For example, if evaluation 8 has a complete ground state and an unfinished
response bundle, restart opens opt_eval_0008_ground_state and
opt_eval_0008_response, reuses qualifying ground-state evidence and resumes the
bundle under the existing restart-mode rules. Response validates the exact
ground-state generation it used; a common evaluation number is insufficient.

## Optimizer checkpoint and save points

The workflow-owned optimizer saves a complete snapshot through the shared
checkpoint interface. Its next operation is either `evaluate_geometry` or
`evaluate_trial`. Publish that operation with the geometry, Hessian, history and
counters as one generation. Electronic application checkpoints remain separate.

| Part | Minimum saved contents |
| --- | --- |
| Context | Record version, optimizer algorithm/settings, model and target-selection context, atom ordering and coordinate/gradient units. Reuse the shared context and dependency fields. |
| Position | Next operation; current evaluation number and coordinates; pending trial number/coordinates when present; next unused evaluation number; number of committed geometry moves. |
| Measurements and history | Current energy/gradient when available, exact application source/sample references, previous comparison energy/gradient and the actual previous displacement. No previous sample is explicit at the initial geometry. |
| Approximate Hessian | Current value and update convention. This history is separate from the electronic KAIN subspace, which restarts empty. |
| Pending line search | For `evaluate_trial`, the base evaluation/geometry, search direction and signed trial scale; the base energy/gradient supplies the directional derivative. |
| Stopped outcome | Evaluated sample and available measurements, stop reason, and the next operation/history needed for compatible continuation. A stop reason is not convergence evidence. |

The current optimizer projects translations and rotations out of the gradient
before its convergence checks and Hessian updates. Preserve the previous gradient
in the convention the optimizer used, and identify the physical gradient from the
application separately or by its saved result reference. Rebuild projectors and
shifted Hessians from the snapshot. Retain current/previous samples and pending
trial data rather than every past optimizer state.

For an optimizer retaining this line-search structure, the proposed save points are:

1. **Before the current geometry evaluation:** publish `evaluate_geometry` with
   its incoming Hessian, previous comparison sample/displacement and evaluation
   number. Obtain the energy/gradient from eligible saved results or resume the
   child calculations, then assess geometry convergence.
2. **Before a trial evaluation:** if continuing, update or intentionally skip
   updating a separate candidate Hessian, compute the search direction and signed
   trial scale, and allocate the pending evaluation number. Publish `evaluate_trial`
   with that Hessian before starting the trial's electronic work. The base geometry
   and geometry-move count have not advanced. This phase never reapplies the base
   geometry's Hessian update.
3. **After selecting a geometry step:** calculate `R_next = R_base + alpha * direction`
   from the saved base geometry. Publish a replacement `evaluate_geometry` containing
   R_next, the updated Hessian, the old current sample as the previous sample,
   the actual displacement, and the incremented geometry-move count. Reuse the
   trial evaluation number if its exact geometry is selected; otherwise allocate
   a new number. Publish before evaluating R_next. Its coordinates do not establish
   a measured energy/gradient; a line-search energy prediction is not a result.

When stopping at an evaluated geometry, publish its measurements and reason
while retaining `evaluate_geometry` for a later compatible continuation. A
controlled stop with a pending trial retains `evaluate_trial` and its trial data.
Thus stopping does not erase whether the Hessian update has already happened.

Build replacement snapshots separately and use the shared complete-write,
validation and atomic-publication procedure. Publication failure stops advancement
and leaves the preceding snapshot selected. Before trial publication, restart may
recompute the Hessian update from the old saved Hessian; after publication, it
resumes the trial with the updated Hessian. Before publication of R_next, restart
may recalculate it from R_base; after publication, it uses R_next directly without
adding the old displacement again. Arithmetic may repeat, but a published state
transition is not applied twice. The existing storage failure guarantees apply.

These small optimizer saves do not request extra orbital archives. Retain child
generations required by the selected snapshot and pending work, and reclaim
obsolete evaluations only through the main design’s publication/retirement ordering (section 5).

## Optimization maxiter on restart

Keep the existing optimization-group `maxiter`. It caps additional committed
geometry steps in the current invocation. Each restart gets a fresh allowance;
previous steps stay in history and are not subtracted. Stop earlier when eligible
evidence satisfies the selected geometry criteria. The allowance spans all
geometry evaluations in that invocation and is not renewed by inner protocol
advancement, application changes or optimizer-history resets. SCF and response
retain their own separate per-protocol update allowances.

Count one step when publishing the selected next-geometry snapshot, together
with its history and geometry-move count. Count the committed optimizer step
regardless of displacement size. Trial evaluations, the initial energy/gradient
evaluation, checkpoint loading, history resets and final assessment do not consume
geometry steps. Keep historical totals and current-invocation usage separate,
with the attempt context identifying the new allowance.

A resumed `evaluate_trial` consumes one step when its selected geometry is
published. A resumed `evaluate_geometry` does not charge an already published
step again merely to assess that geometry. For example, after 12 committed steps,
a restart with `maxiter = 5` permits up to five more, for at most 17 total. A
pending trial's selected move is the first of those five.

After the last allowed step, obtain or reuse the final selected geometry's
energy and gradient and assess convergence without another geometry update.
Report convergence only if the requested criteria pass; otherwise report the
iteration limit and save the resumable state. An incomplete or failed electronic
evaluation retains that outcome and cannot certify geometry convergence. This
terminal assessment can require electronic work; time it separately from geometry
updates. The cap limits geometry steps, not energy/gradient calls or elapsed time.

## Optimizer history after changed settings

Preserve the optimizer Hessian, comparison history and pending trial when the
saved optimizer context and the energy/gradient evidence on which they depend
remain valid for the current request. Tightening requested tolerances alone does
not require a reset if the saved evidence already satisfies them. Reassess the
supporting application measurements and exact dependencies under the shared
restart rules; numerical agreement of two energy values is not a validity check.

If changed settings invalidate any energy or gradient used by that history,
continue from the saved current geometry with fresh optimizer history. This can
happen when a tighter protocol, operator threshold or smoothing setting requires
new electronic work for more accurate gradients. Reset the approximate Hessian
using its normal initialization, clear previous comparison samples and displacement
comparisons, and discard the pending search direction/trial. Publish a replacement
`evaluate_geometry` snapshot before new electronic work. Preserve the current
evaluation number, coordinates, historical move count and next unused evaluation
number. A reset is not a geometry move and does not recycle trial numbers.

The saved current geometry is the base geometry while a trial is pending, or
the selected next geometry once its step snapshot has been published. Do not
promote an unfinished or rejected trial to the restart geometry. Reuse eligible
SCF/response work at that geometry, or compatible functions as guesses, while
obtaining energy and gradient valid for the new request. Existing application
compatibility, exact-reference, conversion and strict-error rules still apply.
The optimizer then reassesses convergence and continues as required; it does
not return to the original input geometry solely because precision was tightened.

Retain the old samples as analysis history, but do not compare old-precision
and new-precision samples as consecutive geometry steps or invent zero changes.
Change measurements remain unavailable until comparable samples exist. Required
unavailable criteria cannot establish convergence. Publish the reset snapshot
before retiring now-unneeded trial archives through the existing retention rules.

## Finite-difference integration and retention

Reuse the [saved displaced-evaluation plan](../../design_draft.md#finite-difference-evaluations-and-restart)
when a validated numerical-gradient provider is available. Displacements inside
optimization share its saved evaluation sequence. A displacement is not a
committed optimizer move and does not consume or replenish its step allowance.
If changed point settings invalidate the optimizer gradient, reset history at
the saved current/base geometry under the rules above.

Protect generations required by the published optimizer snapshot, pending work,
retained results and any state-following reference. Publish a replacement snapshot
before recording retirement and deleting unreferenced archives. Do not retain
archives for every geometry solely because evaluation numbers exist. Reuse the
shared checkpoint publication and retirement operations.

Acceptance cases [D01 and D02](madqc-acceptance-checklist.md#4-optimization-follow-on-and-finite-differences)
remain reviewed behavioral requirements for this follow-on work. They are not
passing tests or current workflow-unification acceptance gates.
