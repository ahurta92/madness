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
each `molresponse` names an earlier `moldft` reference. Names are case-sensitive,
start with a letter, and otherwise contain letters, digits, `_`, or `-`.
All names must be unique. Each block ends with `end`; nested blocks are unsupported.
The complete deck is checked before any stage starts.

Every stage has its own typed parameters and directory. Two response stages
can reference the same ground stage while requesting different frequencies
and numerical settings. Geometry and the Hamiltonian come from that reference;
response `protocol`, `dconv`, and `maxiter` belong to the response stage. Omitted
response protocol/maxiter use the existing parameter-group defaults; omitted
response dconv uses the existing `100 * finest protocol threshold` derivation.
No response setting silently tightens the upstream calculation.

The first implementation supports these controls:

| Block | Controls |
|---|---|
| `molecule NAME` | Ordinary atom coordinates, `units`, `eprec`, `no_orient` |
| `moldft NAME` | `molecule`, `xc hf`, `protocol`, `dconv`, `econv`, `maxiter`, `L`, `k`, `charge`, `localize`, `print_level`, `maxsub`, `maxrotn`, `requested_properties [energy, dipole]` (either or both) |
| `molresponse NAME` | `reference`, `protocol`, `dconv`, `maxiter`, `kain`, `print_level`, `requested_properties [polarizability]`, `polarizability.frequencies`, `polarizability.frequency_units au`, `polarizability.components` |

Ground states must be closed shell. Response frequencies and components must
be explicit; components are Cartesian pairs such as `[zz]` or `[xx, xy]`.
The engine may calculate extra tensor elements needed to assemble the request;
`properties.polarizability.values` selects the requested components. The raw
engine evidence is retained alongside them. Response protocols finer than the
reference threshold/order are rejected. A response `k` override is unsupported.
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

The store is `<prefix>/<stage>/checkpoints/`. Each generation contains its own
`state.restartdata.*` physical shards and `manifest.json`; the `latest` pointer
selects exactly one generation. Occupations, orbital energies/order, molecular
metadata and functions are in the orbital archive. AO sidecars are unnecessary.
All writers use checked close; physical size/checksum and writer-count validation
precede publication and loading. A private attempt copy supplies the existing
SCF restart planner, keeping saved generations immutable. Never select archives
by modification time or combine shards from different generations.

An intermediate generation is published after each orbital update (and at the
final rotation boundary). It is resumable but carries no convergence claim.
Its energy/density/residual measurements identify their input/previous-update
sampling phase; they are not presented as final properties of the saved output.
After the solver returns, energy and dipole are evaluated on the actual final
orbitals, then published together with that state. This adds a final evaluation,
not a new stopping test. Existing engine convergence evidence determines whether
a ground state is acceptable for response. An unfinished state remains saved
and causes a nonzero exit before any response stage runs.

The existing SCF gate uses `max(dconv, protocol threshold)`; both requested and
effective targets are recorded. `econv` is recorded, and tightening it forces
iteration, but this slice does not introduce an energy-change stopping gate or
claim the full milestone-2d convergence contract.

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
The energy record includes the requested `econv` and its comparison, explicitly
marked `active_gate: false`. A missing first comparison at each protocol rung
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

These values are diagnostic (`active_gate: false`). The previous response-wide
`delta_kinetic_energy` placeholder becomes `not_applicable` when per-function
values are present: no channel-wide scalar reduction has been selected.
Per-function differences or trend-based stopping are separate future policy.

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
