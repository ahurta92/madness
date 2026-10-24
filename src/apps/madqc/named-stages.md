# Named moldft → molresponse calculations

Run the checked-in H2 example with a new output directory:

```
MAD_NUM_THREADS=2 madqc named_h2.in --prefix=H2-run
```

`named_h2.in` defines a molecule, `moldft gs`, and `molresponse optical`.
The ground stage requests HF energy and dipole; the response stage requests
`zz` polarizability at frequencies `0.0` and `0.02` atomic units. These are
example numerical settings, not new solver defaults.

Named stages run in file order. Each `moldft` names an earlier molecule, and
each `molresponse` names its `moldft` reference: an earlier stage of the deck,
or a stage outside the deck whose published generation already exists in the
output directory (run it there with `--continue`). Names are case-sensitive,
start with a letter, and otherwise contain letters, digits, `_`, or `-`.
All names must be unique. Each block ends with `end`; nested blocks are unsupported.
The complete deck is checked before any stage starts.

Every stage has its own typed parameters and directory. Two response stages
can reference the same ground stage while requesting different frequencies
and numerical settings. Geometry and the Hamiltonian come from that reference,
and only from its published generation: the orbital archive and the
`numerical.consumer_context` block of its manifest (the molecule as the SCF
oriented it, in bohr; charge, multiplicity, cell `L`, functional, spin
treatment, final `thresh` and `k`). Nothing is copied from the reference
stage's input, and no engine object passes between stages, so running the
reference alone and the response later gives the same results as one deck.
Response `protocol`, `dconv`, and `maxiter` belong to the response stage. Omitted
response protocol/maxiter use the existing parameter-group defaults; omitted
response dconv uses the existing `100 * finest protocol threshold` derivation.
No response setting silently tightens the upstream calculation.

The first implementation supports these controls:

| Block | Controls |
|---|---|
| `molecule NAME` | Ordinary atom coordinates, `units`, `eprec`, `no_orient` |
| `moldft NAME` | `molecule`, `xc hf`, `protocol`, `dconv`, `econv`, `maxiter` (one number or one per rung), `checkpoint.every_iterations`, `checkpoint.keep_all`, `convergence.use_default_criteria`, `accuracy.bsh_residual`, `accuracy.kinetic_trace_change`, `progress.*`, `L`, `k`, `charge`, `localize`, `print_level`, `maxsub`, `maxrotn`, `requested_properties [energy, dipole]` (either or both) |
| `molresponse NAME` | `reference`, `protocol`, `dconv`, `maxiter`, `kain`, `print_level`, `requested_properties [polarizability]`, `polarizability.frequencies`, `polarizability.frequency_units au`, `polarizability.components` |

Ground states must be closed shell. Response frequencies and components must
be explicit; components are Cartesian pairs such as `[zz]` or `[xx, xy]`.
The engine may calculate extra tensor elements needed to assemble the request;
`properties.polarizability.values` selects the requested components. The raw
engine evidence is retained alongside them. Response protocols finer than the
reference threshold/order are rejected (before any stage runs for a reference
in the deck, and against the published generation for every reference). A
generation published before the consumer context existed cannot be a reference;
run that stage again. A response `k` override is unsupported.
Frequencies and response protocol rungs that collide at the engine's existing
filename-key precision are rejected. Unknown controls and unsupported properties
fail rather than being ignored.
Command-line scientific overrides and the old `--wf` option do not apply to
named decks; put settings in their owning block. Existing unnamed decks retain
their workflow interface.

Outputs go to `<prefix>/<stage>/`, without changing the process working
directory. Each invocation writes engine files into a unique `attempt-...`
subdirectory; the stage-level `results.json` describes the latest attempt. `<prefix>/results.json` indexes the stage results; each stage also
has `results.json`. SCF's resolved input and engine files are under its stage,
and response metadata/state files are under the response stage. Results carry
method, units, effective reference cell, requested targets, and the existing
engine's convergence evidence. An unconverged stage is recorded and causes a
nonzero exit before descendants run.

Energy is in hartree, dipole in `e*a0`, and polarizability in `a0^3`; response
frequency is in atomic angular-frequency units. Property records identify their
state and, for response, the source ground state. A moldft state is identified
by its immutable checkpoint generation. Response states also identify immutable checkpoint generations and their exact
ground-state reference generation.

To continue an existing run, use the same output root explicitly:

```
MAD_NUM_THREADS=2 madqc named_h2.in --prefix=H2-run --continue
```

Without `--continue`, an existing root is refused. With it, each moldft stage
must have a published checkpoint. The adapter validates its complete orbital
inventory and checks the molecule, model, occupations implied by charge/spin,
localization and cell before the engine loads anything. Changed physical inputs
require a new output root; explicit seed conversion is not provided here.
Tighter protocol/density targets resume from saved orbitals. An unchanged
completed request loads without SCF iterations. Iteration limits and requested
energy/dipole properties can change. Results report `work` (`solved`, `loaded`,
`resumed`), `iterations_this_run`, and acceleration-history reset.

Where a moldft stage continues is decided from the saved generation's rung,
matched by value (`thresh`, `k`), not by position in the `protocol` list. A
finished rung is never redone: the stage continues at the next tighter
requested rung. An unfinished rung resumes where it stopped. A finished final
rung is loaded if it meets the current request and resumed otherwise.
`maxiter` (one number, or one per rung such as `maxiter [30, 10]`) is the
number of iterations a rung may take in this run; a `--continue` run starts
the count again, so a resumed rung gets `maxiter` more iterations. Each
generation records the rung (`position.rung`): its index, `thresh`, `k`, the
iterations in this run (`iterations_this_run`), the total across runs
(`iterations_at_rung`, recorded, never compared with `maxiter`), and how it
ended (`converged`, `stalled`, `iteration_limit`, or `in_progress`). The
stage's `results.json` `continuation` block says what happened and why, e.g.
"resumed at rung 2 (thresh 1e-06, k 8), iteration 3 of that rung".

The store is `<prefix>/<stage>/checkpoints/`. Each generation contains its own
`state.restartdata.*` physical shards and `manifest.json`; the `latest` pointer
selects exactly one generation. Occupations, orbital energies/order, molecular
metadata and functions are in the orbital archive. AO sidecars are unnecessary.
All writers use checked close; physical size/checksum and writer-count validation
precede publication and loading. A private attempt copy supplies the existing
SCF restart planner, keeping saved generations immutable. Never select archives
by modification time or combine shards from different generations.

A moldft stage publishes a generation at its decision points: every rung end
(the final rotation boundary) and every stop (converged, stalled, iteration
limit), plus the final properties. Saving the orbitals is not cheap for large
systems, so it does not publish after every update by default.
`checkpoint.every_iterations K` (K >= 1) also publishes every K-th iteration
at a rung, counted across `--continue` runs, so an interrupted run can resume
from inside the rung. After each publish, older mid-rung generations of the
stage are deleted (save and replace); `latest`, every rung end and stop, and
any generation another stage's manifest references are kept. A deletion first
renames the directory to `.deleting-<generation>`; one that fails leaves extra
files only and is reported in the stage's `checkpoint_warnings`.
`checkpoint.keep_all true` keeps every generation, for debugging. The settings
table shows both settings and who set them. molresponse rejects both keys: it
publishes at solve ends only, having no save boundary inside a solve yet.
A rung-end or stop generation carries the whole per-iteration history of this
run (`convergence_report.iteration_history`) and, after `--continue`, the
history of earlier runs (`earlier_runs_iteration_history`, a record only), so
the history never depends on which generations were kept.

An intermediate generation (with `checkpoint.every_iterations`) is resumable
but carries no convergence claim.
Its energy/density/residual measurements identify their input/previous-update
sampling phase; they are not presented as final properties of the saved output.
After the solver returns, energy and dipole are evaluated on the actual final
orbitals, then published together with that state. This adds a final evaluation,
not a new stopping test. Existing engine convergence evidence determines whether
a ground state is acceptable for response. An unfinished state remains saved
and causes a nonzero exit before any response stage runs.

The existing SCF gate uses `max(dconv, protocol threshold)`; both requested and
effective targets are recorded. `econv` is the target of the energy-change
metric, but only when the deck writes it; the engine's default `econv` is not a
request and is recorded as `set_by: engine_default_not_a_target`. A written
`econv` is judged (`relation: less_than_or_equal`) on the energy change
measured on the saved orbitals (`saved_state_validation.measurements.energy_change`),
never on the `econv` a generation recorded. Whether it is enforced depends on
`convergence.use_default_criteria` (below): with the engine check deciding, it
is reported only, and an unmet `econv` on `--continue` is reported "unmet, not
enforceable in engine mode" instead of making the stage resume.

On `--continue`, each response stage validates its own published generation.
A completed identical request reuses its saved properties with zero response
iterations (`work: loaded`). Saved functions without the requested completed
analysis are restored into a private attempt: the engine either analyzes those
prerequisites (`analyzed`), continues unfinished/tighter states (`resumed`), or
solves new frequencies/channels. Adding frequencies preserves compatible saved
states. `iterations_this_run` counts response iterations across channels and
protocol steps; `work_items` records the individual solves and restart seeds.
A previous invocation's exhausted iteration budget does not exhaust the new one.
Acceleration history is rebuilt on continuation.

Response generations contain the complete saved FD inventory and engine metadata.
Each static channel contains X; each dynamic channel contains X/Y together in one
native archive. Channel, model, shell, frequency, numerical protocol, per-state
convergence target and the exact ground generation govern reuse. Properties are
published only after analysis and identify the resulting response generation.
A result JSON file alone never establishes restart validity. Changing the ground
reference generation invalidates its response descendants; other named stages and
old generations remain intact. Supported coarser response work uses projected
working copies and preserves the fine source archive. Unsupported downstream
resolution is refused; upstream refinement must be requested explicitly.

Response publication occurs after each completed channel/protocol solve, including
an unfinished solve that reached its iteration limit. Interruption within a solve
loses work since that boundary. Every publication copies the current saved-state
inventory, trading disk space and I/O for a simple complete generation. Retention
and incremental storage are not implemented. Named response checkpoint payloads
use checked native archives even in HDF5-enabled builds; the optional ground-state
HDF5 mirror stays in the attempt directory. Standalone/legacy response backend
selection remains unchanged.

Beta, Raman, 2PA, excited states, Nemo, MP2, CC2, and optimization are not exposed
through this named frontend yet. The broader convergence contract is milestone 2d.

Checkpointing assumes a shared filesystem and supported archive reader counts
(at least the original writer count). An interrupted write leaves the prior
published generation selectable. This does not provide fault-tolerant MPI or
power-loss durability. Unpublished generations are ignored, and automatic
cleanup is deferred. A stale `.publish-lock` requires manual recovery after
confirming no publisher is active. Saving each iteration adds disk space and
I/O costs; retention/cadence controls are not included in this slice.

Tests: `test_named_stages` covers block parsing, ownership and errors in serial
and two ranks. `test_named_stages_e2e.py` compares fresh H2 energy/alpha with the
existing workflow, exercises two response instances, preserves an existing output
root on refusal, and checks the cell after automatic centering of translated H2.
Its numerical runs retain logs and JSON beneath a temporary test directory.
On OpenMPI, assign sufficient CPUs per rank (for example `--bind-to none` for
local testing); do not run threaded tests with one core assigned per rank.

`test_named_checkpoint_e2e.py` exercises complete reuse, unfinished continuation,
tighter targets, unpublished-generation isolation, publication failure, missing
or corrupt payloads, incomplete inventory and incompatible geometry in serial
and two MPI ranks. It compares final H2 energy and dipole and preserves the
original generation hashes across continuation.

## Convergence evidence (milestone 2d reporting slice)

Each stage's `results.json` now contains `convergence_report`, with
`gate_policy: existing_engine`. The numerical solvers keep their existing
stopping tests. An engine's completed verdict is distinct from a new,
all-criteria certification of the saved functions.

Ground stages report density change (maximum spin-density L2 change), raw BSH
residual (maximum orbital/spin norm), total-energy change, and the absolute
change of the occupation-weighted electronic kinetic energy summed over spins.
Energy and kinetic changes are diagnostic; no new stopping gate is enabled.
The energy record carries a user-written `econv` as its target and verdict,
explicitly marked `active_gate: false`; without one it has no target and no
verdict. `results.json` adds `convergence_report.energy_change_request`, the
current request judged on the saved state. A missing first comparison at each protocol rung
or restarted solve has `value: null`, never an invented zero.

Response reports contain separate channel/frequency/protocol records. Density
and BSH measurements include the actual absolute and relative engine gate
limits, the normalization scale and its sampling phase. The residual is from
the input to the last solver step; the response norm used by the existing
relative BSH gate is from the output after acceleration. These two phases are
reported separately. Units follow the native response amplitude and induced
density conventions, with perturbation/model/spin labels; no size scaling is
introduced. A response-local energy change is `not_applicable`; the exact
source ground generation's energy evidence is under `reference`. Response
kinetic values are reported per function, as described below; no response-wide
kinetic-change reduction is selected.

Every measured record identifies its sample and phase. The current engine
samples precede the archived update/final rotation, so
`valid_for_saved_state: false` prevents interpreting them as a residual
re-evaluation on the archived functions. Final energy/dipole/polarizability
properties retain their separate generation-bound evaluations. Reuse preserves
historical samples; old checkpoints without the new records report unavailable
evidence rather than acquiring new measurements from result files.

Existing response protocol progression remains: usable stalled or exhausted
nonfinal rungs can seed the next supplied rung; unmet final targets report an
unfinished stage and retain usable progress. No extra rung is invented. Named
adapters reject states flagged diverged or carrying invalid measurements before
load or publication, preserving the previous usable published generation.
This restriction does not alter the legacy response planner's recovery policy.

The additional `saved_state_validation` record reevaluates the archived output
without taking another solver step. Ground stages report density and raw BSH
residuals plus energy/kinetic changes. Changes compare against a genuine
iteration-input sample from the same rung; final diagonalization is not itself
treated as a convergence step. Response channels reevaluate raw BSH and induced
density after acceleration, using that saved state's norms in the existing
absolute-plus-relative formulas. The comparison is the last step's input.

These records carry `affects_stopping: false` and generation identity. A measured
record has `valid_for_saved_state: true`; this means the measurement applies to
those functions, while `satisfied` separately reports comparison to a limit.
A completed engine verdict can coexist with an unmet diagnostic comparison.
Identical reuse preserves the recorded evidence without new operator work.
Older checkpoints remain readable without inventing revalidation results.

Full convergence certification, new energy/kinetic gates, per-function kinetic
history policy, and new stall/default formulas remain later scientific work.


## Per-function response kinetic values

Named static/full response now records `kinetic_energy` for each stored
function: `K_i = 1/2 sum_d <D_d f_i | D_d f_i>`, the diagonal contribution used
by SCF before its separate occupation-weighted energy sum. Each X/Y function
retains its own scalar, indexed within its named block and the exact reference's
orbital ordering. No mean, RMS, maximum, spin factor, amplitude normalization,
or implicit Y duplication is applied.

`convergence_report.channels[].function_kinetic_history` holds the initial
iterate and every post-acceleration iterate within that solve/rung.
`measurements.kinetic_energy` describes the last step's input;
`saved_state_validation.measurements.kinetic_energy` independently reevaluates
the archived output. Values live under `blocks.x_alpha`, `blocks.y_alpha`, etc.,
as records containing `function_index`, `value`, and `availability`. The
containing record gives units, phase, iteration, threshold and k. The channel
and reference records identify perturbation, frequency and reference generation.

Only scalars persist between iterations. Reprojection, a new reference or a
restarted solve begins a new history segment. Exact reuse preserves existing
history; older checkpoints without it acquire no invented measurements.
Invalid individual scalars or history records prevent named publication/reuse.

The raw values remain unreduced. Each history record additionally contains
`change_from_previous`: matched per-function absolute changes and their maximum,
`max_i abs(T_i(n)-T_i(n-1))`, without occupation/spin weights. The record retains
component changes, the limiting function, and both sample identities. The initial
sample is unavailable, not zero. Only adjacent iterates of the same response
channel/reference/attempt/rung are compared; no SCF orbital matching is inferred.

`measurements.delta_kinetic_energy` describes the last step's input relative to
its genuine preceding iterate (unavailable for a one-step run). The corresponding
saved-state field compares the independently evaluated output with the last
input. It never differences two evaluations of the same state. These diagnostics
have `active_gate: false`, `target: null`, and `satisfied: null`: no kinetic
threshold or stopping criterion is selected. Derived evidence is validated
against its source records when publishing or loading a checkpoint.

The shared numerical definitions and current application-specific reductions are
recorded in [metric-definitions.md](../../madness/chem/metric-definitions.md).

## Per-function BSH residuals

Ground and response `measurements.bsh_residual.per_function` now records the
raw scalar `||f_i - F_i(f)||_L2` for each stored function, using the shared
`madness/chem/FunctionResidual.h` definition. Its `blocks` map contains arrays
of `{function_index, value, availability}`: `alpha` (and explicit `beta`) for
SCF, `x_alpha` and any stored Y/beta blocks for response. Values have no
occupation weights, normalization or per-function stopping gate. Units and
sample identity are inherited from the containing BSH measurement.

Ground components follow the same previous-update phase as its existing maximum
residual. Response components are measured before KAIN and retained in
`convergence_report.channels[].function_bsh_history`; there is one sample per
raw solver update, indexed by its input iteration. Independent saved-state
components live under `saved_state_validation.measurements.bsh_residual`.
Historical components describe the sampled working basis, not necessarily the
orbital ordering after a later rotation. Histories restart at each solve/rung.

The existing ground maximum and response quadrature aggregates and all stopping
rules are unchanged. Checkpoints predating these optional fields stay reusable
without fabricated component evidence. Present invalid components or response
history aggregates are rejected before reuse or publication.

## Diagnostic comparison of response BSH reductions

Every new `function_bsh_history` entry records that step's `absolute_limit`,
`relative_factor`, `normalization_scale`, `normalization_phase`, `effective_limit`
and existing scalar `satisfied` result. The raw residual is from the step input;
the existing relative normalization uses the output after KAIN. Limits are
captured during each step, never copied backward from the final iteration.

Historical and saved BSH measurements also contain `maximum_comparison`. It
reports the maximum raw component, its limiting block/index and whether it meets
**the same recorded channel limit**, using the existing inclusive comparison.
It changes only the residual reduction for this diagnostic; it retains the
existing normalization. Its `active_gate` is false. Passing this one comparison
is neither a complete convergence verdict nor a prediction of the trajectory
under a different stopping policy. Older records without it remain usable.

Response strict-convergence eligibility follows the existing solver rule: at
least two updates within the current solve/rung are needed to measure a density
change, and any configured minimum iteration count also applies. A stall or
divergence is a loop exit, not strict convergence. In particular, a one-step run
with residuals below deliberately loose targets remains unfinished and can be
continued. Restarted functions seed a new comparison history; loading them does
not supply a missing density sample. Named checkpoint validation rejects a
record claiming strict convergence after fewer than two updates or while
stalled. This corrects the executor's verdict; it does not change residual
formulas, thresholds, or enable the proposed SCF charge/spin or kinetic gates.

Ground-state measurement reports also include `charge_density_change` and
`spin_density_change`: L2 changes of rho_alpha+rho_beta and rho_alpha-rho_beta,
respectively, in a0^(-3/2). These use actual fields, not arithmetic on the
separate spin norms. Historical and saved-state samples remain distinct;
first differences without a comparable prior input are unavailable. In the
current restricted-HF frontend, measured spin change is zero and charge change
is twice the old single-spin norm. The original `density_change` and its gate
remain visible for comparison. The new fields have no target, effective limit
or pass/fail verdict. No new default thresholds or stopping rules are selected;
reference tests set their existing numerical criteria explicitly.

`diagnostic.iterations_per_rung N` is an explicit trajectory-study control on
moldft and molresponse stages. Positive N requests that many iterations at each
supplied rung even if the ordinary gates pass earlier; zero preserves the
ordinary loop. Invalid/divergent states still stop. It supplies no accuracy
threshold and does not turn an unmet final target into convergence. Reports keep
ordinary scientific verdicts and full-precision density/trace histories. This
experimental diagnostic is intended for fresh bounded protocol studies; it does
not establish a new stopping default.

### Who decides convergence: `convergence.use_default_criteria`

Every written convergence value of a named `moldft` stage is a **target** on
its own metric (DESIGN B8):

| Keyword | Metric | Quantity compared |
|---|---|---|
| `dconv` | `charge_density_change` | L2 change of the total charge density rho_alpha+rho_beta, `a0^(-3/2)` |
| `econv` | `energy_change` | absolute change of the total energy, hartree |
| `accuracy.bsh_residual` | `bsh_residual` | largest raw orbital BSH residual, dimensionless |
| `accuracy.kinetic_trace_change` | `delta_kinetic_energy` | absolute change of the unweighted kinetic trace (restricted alpha counted once), hartree |

Only a value the deck writes is a target; engine defaults are not (B3, D8).
Requested targets compare with `<=` and have no protocol floor or hidden
factor. `accuracy.density_change` and `accuracy.energy_change` are replaced by
`dconv` and `econv` and are rejected with that message.

One per-stage switch says which group stops the run:

- `convergence.use_default_criteria true` (the default): the **engine check**
  stops the run: both spin-density changes below `max(dconv, thresh) x max(5,
  natom)` and the largest BSH residual below `5 x max(dconv, thresh)`
  (`SCF::engine_check_limits`). The written targets are measured, recorded
  and reported, not enforced; each one gets a warning at stage start saying
  what the engine check does instead.
- `convergence.use_default_criteria false`: only the **written targets** stop
  the run; the engine check is still computed and recorded
  (`engine_check_passed`) but not applied. A request without a density or BSH
  target is accepted with a warning. A stage with no written target is
  rejected before solving.

There is no AND mode. Divergence and invalid states stop the run in both
modes. The engine records its own restart claim (`converged_for_dconv`) only
when its check passed.

At stage start the run prints the settings table (every metric: target, who
set it, enforced or reported only, relation, effective limit with its formula,
the per-rung `maxiter` and warnings); at stage end it prints each metric's
final value on the saved orbitals, the verdict and the limiting metric. Both
tables are rendered from the struct that `results.json` records
(`convergence_settings`, `convergence_report.convergence_outcome`), and the
manifest carries `targets.decided_by`, `convergence_report.decided_by` and the
settings.

`convergence_report.explicit_accuracy` records the last iteration's evaluation
of the written targets; `saved_state_validation.explicit_accuracy` evaluates the
archived state; `requested_accuracy` re-evaluates the saved measurements
against the current request, with `enforced` and, in engine mode, the outcome
"unmet, not enforceable in engine mode" for an unmet target.

On `--continue` (B6, B8): with the written targets deciding, a saved state is
loaded only when its saved measurements meet every written target, whichever
mode produced it; otherwise it resumes. With the engine check deciding, a
saved state that passes the engine check for the current request is loaded,
and written targets it doesn't meet are listed in
`continuation.unmet_not_enforceable` instead of forcing a resume. Historical
kinetic evidence with another quantity definition cannot certify a new
trace-change target, even when its units are also hartree.

Stop reasons (B5): `converged`, `stalled` (configured progress), and
`iteration_limit`; with the written targets deciding,
`unmet_requested_accuracy` means the final rung stopped on its targets but the
saved orbitals no longer meet them.

Nemo explicit accuracy policies are not supplied by this interface yet. A
forced-iteration diagnostic does not override invalid-state checks or certify
unmet explicit accuracy.

### Driven-response targets and who decides

A named `molresponse` FD stage follows the same rule as `moldft` (DESIGN B8):
written values are targets on their own metrics, and
`convergence.use_default_criteria` (default `true`) says who stops each solve.

| Keyword | Metric | Quantity compared |
|---|---|---|
| `dconv` | `charge_density_change` | L2 change of the kernel's induced charge density (`kernel_induced_density_L2_atomic_units`), absolute |
| `accuracy.bsh_residual` | `bsh_residual` | largest raw residual over all stored functions and X/Y blocks (`response_amplitude_L2_atomic_units`) |
| `accuracy.kinetic_trace_change` | `delta_kinetic_energy` | absolute change of the unweighted positive trace over those blocks (`hartree_times_response_amplitude_scale_squared`) |

These response quantities include the perturbation amplitude scale and must not
be read as a ground-state energy or density. A driven response has no local
energy, so there is no energy target; `accuracy.energy_change` is rejected
with that message, and `accuracy.density_change` is replaced by `dconv`.

- **Engine check decides** (default): each leg converges when `||dx|| <= 5
  max(dconv, thresh) (1 + ||x||)` and `||drho|| <= 5 max(dconv, thresh) (1 +
  ||rho||)` from the second iteration (`ConvergencePolicy::effective_for_thresh`),
  with the engine's plateau detector and best-effort acceptance as before. The
  written targets are evaluated on the solver's input and on the saved state
  and reported (`affects_stopping: false`), with a start-of-stage warning.
- **Written targets decide** (`false`): a solve converges when every written
  target is met on the last genuine input and on the independently reassessed
  saved output; the engine check is not applied, the plateau detector and
  best-effort acceptance are off, and only `progress.*` can declare a stall.
  A request without a target is rejected; one without a density or BSH
  target is warned about. An unmet target withholds the stage's properties.

Each saved entry records `checkpoint_decided_by`. On `--continue` a converged
entry is re-solved for an unmet written target only when the targets decide.
The settings table at stage start shows both groups (the engine's per-leg
limits as formulas, since they are relative); the outcome table at stage end
shows each written target's largest saved value over the channels of the
final rung. Forced diagnostics exhaust the specified iterations per supplied
rung without certifying unmet accuracy. Checkpoints preserve raw histories,
comparison identities and reconstructed decisions; missing or incompatible
historical evidence cannot certify a newly requested criterion.

### Experimental configured progress

Named SCF and driven-response FD stages may set `progress.window` (a positive iteration lookback) and
`progress.minimum_relative_improvement` (an explicitly supplied fraction strictly
between zero and one), together with at least one written target. For `moldft`
these stall rules act on the written targets, so they require
`convergence.use_default_criteria false`. Omitting
them preserves the existing iteration policy. There is no default scientific
stall threshold. A criterion still improving prevents a stall; unavailable or
invalid measurements cannot justify promotion. History resets at every rung.

A stalled nonfinal solve ends at its complete iteration boundary and the existing
protocol loop advances to the immediate next supplied rung. At the final rung,
the usable generation is retained with unmet accuracy. `configured_progress`
records the decision separately from engine convergence and explicit accuracy.
Forced iteration diagnostics record that decision but suppress its early exit;
`action_applied` distinguishes the recorded recommendation from execution.

For FD, the policy carries the entire supplied protocol ladder into each
single-channel solve. `configured_progress_history` records each input-side
decision separately from saved-state accuracy. Checkpoint readers reconstruct
the decisions from raw metric histories and the stored policy. A continued
solve resets its scheduler count and starts new histories; archived decisions
retain their original sample identity. No state marked invalid or divergent
can be promoted or replace the prior usable generation.
