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
| `molresponse NAME` | `reference`, `protocol`, `dconv`, `maxiter`, `kain`, `print_level` (moldft's scale: unwritten or 1 to 3, one row per solver iteration; 4 to 9, channel and save detail; 10, debug), `checkpoint.keep_all`, `seed.dalton_dir`, `seed.freq_tol` (see "Dalton seed"), `requested_properties` (any of `polarizability`, `hyperpolarizability`, `raman`, `excited_states`), and per optical property `<property>.frequencies`, `<property>.frequency_units au`, `<property>.components`; for Raman also `raman.atoms`, `raman.directions` (see "Hyperpolarizability and Raman"); for excited states `excited_states.roots`, `excited_states.method`, `excited_states.guess`, `excited_states.guess_basis`, `econv`, the excitation energy change target, and `stall.window`, `stall.ratio`, the bundle's plateau detector (see "Excited states") |

Ground states must be closed shell. Every requested property names its own
frequencies and `frequency_units au` (numbers are never defaulted); its
components are Cartesian tuples such as `[zz]` or `[xx, xy]` and mean every
tuple when omitted (9 pairs for the polarizability and Raman, 27 triples for
the hyperpolarizability). The engine takes axes and one frequency grid, not
components, so it may calculate extra tensor elements needed to assemble the
request; `properties.<property>.values` selects the requested ones. The raw
engine evidence is retained alongside them. The engine's own spellings
(`dipole.*`, `quadratic`, `beta.*`, `raman.nuc_*`, `nuclear.*`) are not stage
controls and are rejected. Response protocols finer than the
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
directory. Without `--prefix` the output root is `./<deck>/` in the current
directory, `<deck>` the input file's name without its extension
(`named_h2.in` writes `./named_h2/`). Each invocation writes engine files into a unique `attempt-...`
subdirectory; the stage-level `results.json` describes the latest attempt. `<prefix>/results.json` indexes the stage results; each stage also
has `results.json`. SCF's resolved input and engine files are under its stage,
and response metadata/state files are under the response stage. Results carry
method, units, effective reference cell, requested targets, and the existing
engine's convergence evidence. An unconverged stage is recorded and causes a
nonzero exit before descendants run. The root `results.json` names the binary
that wrote it in `build` (`commit`, `description` as `git describe --dirty`
gave it when the build was configured, and `built`); the startup banner prints
the same commit and description. Each generation's manifest keeps the
stage's scientific report in `metrics.convergence_report` (moldft's, and a
response stage's: its settings, decider, and for the analysed generation the
outcome rows and each needed state's judged targets); a moldft stage's
`results.json` gives it once, under `convergence_report`.

Four lines "Authorization required, but no authorization protocol specified"
at startup on a host that runs an X server come from the MPI runtime, not
from madqc: Open MPI's topology discovery in `MPI_Init` loads hwloc's GL
component, which probes the X display. They are harmless;
`export HWLOC_COMPONENTS=-gl` silences them.

Energy is in hartree, dipole in `e*a0`, polarizability in `a0^3`, the
hyperpolarizability in `e^3 a0^3 Eh^-2` and the polarizability gradient in
the engine's units, `a.u. (e^2 bohr^2 / Eh per bohr)` of displacement; response
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
localization and cell before the engine loads anything. A changed geometry or
other physical input is a different problem: the stage starts fresh in the same
stage directory (below); explicit seed conversion is not provided here.
Tighter protocol/density targets resume from saved orbitals. An unchanged
completed request loads without SCF iterations. Iteration limits and requested
energy/dipole properties can change. Results report `work` (`solved`, `loaded`,
`resumed`), `iterations_this_run`, and acceleration-history reset.
Before the work, a continued stage prints its plan, one line per state with
the reason it was planned that way (`stage <s> plan: <state> - <action> at
rung <i> (thresh, k): <why>`); after it, one line counts the states loaded,
resumed, seeded and started fresh, with the iterations this run. The same is
recorded under `continuation`. A response stage whose whole request is met by
its saved generation (request equal, stopped converged, properties present)
is loaded whole: its `continuation` is `{planner: request_complete,
generation, why, summary}`. A fresh run prints and records neither.

What a moldft stage does on `--continue` is decided by one planner
(`plan_ladder`, `src/madness/chem/LadderPlanner.h`) from the saved generation
and the current request. Saved values are re-judged against today's limits and
targets; a saved verdict is never trusted. Rungs are matched by value (`thresh`,
`k`), not by position in the `protocol` list. A finished rung is never redone:
the stage continues at the next tighter requested rung. An unfinished rung
resumes where it stopped. A finished final rung is loaded if it meets the
current request and resumed otherwise, except that a final rung that stalled
with the request unmet is loaded and reported unmet (resuming would stall
again; tighten the protocol or loosen a target). A changed problem (geometry,
model, charge, cell) starts fresh in the same stage directory, leaving the
old generations untouched until the new problem publishes: retention then
deletes them unless another stage's manifest references them or the deck
writes `checkpoint.keep_all true`, and `latest` names the new problem's. A generation with no recorded measurements (from
before the convergence report) cannot be re-judged and is resumed.
`maxiter` (one number, or one per rung such as `maxiter [30, 10]`) is the
number of iterations a rung may take in this run; a `--continue` run starts
the count again, so a resumed rung gets `maxiter` more iterations. Each
generation records the rung (`position.protocol_step`): its index, `thresh`, `k`, the
iterations in this run (`iterations_this_run`), the total across runs
(`iterations_at_protocol_step`, recorded, never compared with `maxiter`), and how it
ended (`converged`, `stalled`, `iteration_limit`, or `in_progress`). The
stage's `results.json` `continuation` block says what happened and why:
`action` (`loaded`, `resumed`, `fresh`), the start `protocol_step`,
`iterations_done_at_protocol_step`, `why`, `unmet` (e.g. a stalled final rung),
`unmet_not_enforceable` (engine mode) and a `summary`, e.g. "resumed at rung 2
(thresh 1e-06, k 8), iteration 3 of that rung".

A loaded state is one that meets the current request, judged on its saved
measurements. If its generation had stopped as something other than
`converged` (e.g. `iteration_limit` under tighter targets, now loosened), the
stage publishes a new generation recording the new verdict: the same orbital
shards, byte for byte, with `stop_reason: converged`, `position.reused_from`
naming the source generation, and `convergence_report.reuse_verdict` (the
stop reason before, who decided, why, and the settings it met). A later stage
then reads the same verdict whether it is in this deck or a later one. The
rung record keeps how the rung actually ended.

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
from inside the rung. molresponse publishes at solve ends only, having no
save boundary inside a solve yet, and rejects `checkpoint.every_iterations`.

Every stage keeps one copy of each saved function (D-070): after each publish,
every older generation of the stage is deleted except `latest` and any
generation another stage's manifest references. `checkpoint.keep_all true`
keeps every generation, for any application. A deletion first renames the
directory to `.deleting-<generation>`; one that fails leaves extra files only
and is reported in the stage's `checkpoint_warnings`. When another stage's
manifest exists but cannot be read (a half-written file, permissions), what
it references is unknown, so the stage deletes nothing that time and reports
the unread manifest there instead (D-087). When the stage ends, the
attempt directory's copies of its latest generation's archives (the engine's
working files, and the copies a `--continue` load makes) are deleted too
(D-072); the attempt keeps its text output and JSON records. The first
`--continue` on an output directory written before this rule deletes every
older generation of each stage that publishes, older-form names included,
unless the deck writes `checkpoint.keep_all true`: that is the rule applied
to an existing store, not data loss. The settings table shows both settings
and who set them.
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
iterations (`work: loaded`). Otherwise the moldft planner (`plan_ladder`)
decides each state the requested values are built from on its own (a dipole
or nuclear-displacement leg at one exact frequency: for `polarizability` the
dipole legs of each component's axes at each frequency; for
`hyperpolarizability` `beta_ABC(-2w; w, w)` the A leg at 2w and the B and C
legs at w; for `raman` the A and B legs at w and the nuclear leg at 0), from
that state's records in the generation (`numerical.states`), by the same
rules: saved values are re-judged against today's engine limits and written
targets; a finished rung is never redone; a state that meets the request is
loaded; a final rung that stalled with the request unmet is loaded and
reported unmet; a state with no record (a new frequency, or a changed exact
frequency) starts at the first rung from the same perturbation's saved state
at the nearest frequency, the lower one on a tie (`action: seeded`, with
`seed_from` and `seed_frequency`), and a new perturbation with nothing saved
is solved fresh. A state the request no longer names is left in the
generation, neither loaded nor deleted. The saved functions and the engine's
state entries are copied into a private attempt unchanged (only the derived
property rows and the Raman normal-mode analysis are left out), and the plans
are passed to the engine as input. Quadratic sources are not planned: the
engine builds each from its two legs at one rung, and drops and rebuilds one
whose leg is solved again at that rung.
From each state's start rung on, every rung is solved once in this run; after
that the engine's own scheduling continues. A loaded state that meets the
request counts as met in the analysis: its `row_accuracy` entry carries
`met: true` and `decided_by: "planner"`, and so does the tensor row built from
it (`decided_by: "planner"`; the same for hyperpolarizability and Raman rows),
while each leg's `converged` keeps the flag the
engine recorded, so the two verdicts stay distinguishable. Two kinds of saved
rung are not planned from, and the state is solved again: a rung whose engine
key equals a requested rung's but whose threshold differs (keys round it, e.g.
1e-6 and 1.01e-6; the next save redefines the key), and a rung tighter than
every requested rung (the analysis reads only states at or coarser than the
final rung). The stage's `results.json` `continuation` block has one record
per state: `state`, `perturbation`, `frequency`, `action` (`loaded`,
`resumed`, `seeded`, `fresh`), the start `protocol_step`, `iterations_done_at_protocol_step`,
`why`, `unmet`, `unmet_not_enforceable` and a `summary` (a seeded state also
`seed_from` and `seed_frequency`); the generation records the
same list as `numerical.plans`. The stage's `work` is `analyzed` when only the
analysis ran, else `resumed`. Adding frequencies preserves compatible saved
states. `iterations_this_run` counts solver iterations summed over every
response state (each driven leg and each excited-state bundle) and protocol
step, the excited-state warm-up steps excluded (a TDA stage's 12 is 7 for the
dipole leg plus 5 for the bundle); `iterations_this_run_counts` says so in
`results.json`, and moldft's says it counts SCF iterations over its protocol
steps. `work_items` records the individual solves and restart seeds.
A previous invocation's exhausted iteration budget does not exhaust the new one,
and forced diagnostic iterations (`diagnostic.iterations_per_rung`) on
`--continue` apply to the final rung only, since finished rungs are not redone.
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

## Hyperpolarizability and Raman

A `molresponse` stage may request the first hyperpolarizability and the
vibrational Raman polarizability gradient beside the polarizability:

```
molresponse quadratic
  reference gs
  protocol [1e-4]
  dconv 1e-4
  maxiter 60
  kain true
  requested_properties [polarizability, hyperpolarizability, raman]
  polarizability.frequencies [0.0]
  polarizability.frequency_units au
  hyperpolarizability.frequencies [0.0, 0.05]
  hyperpolarizability.frequency_units au
  hyperpolarizability.components [zzz, xzz]
  raman.frequencies [0.0]
  raman.frequency_units au
  raman.components [zz]
  raman.atoms [1]
  raman.directions [z]
end
```

**Hyperpolarizability.** A component `abc` is `beta_abc(-2w; w, w)` for each
driver frequency `w` in `hyperpolarizability.frequencies`: second-harmonic
generation, which at `w = 0` is the static hyperpolarizability. `a` is the
output dipole axis at `2w`, `b` and `c` the input axes at `w`. The engine
solves the dipole legs at `w` and `2w` and one quadratic source per `(b, c)`
pair; `properties.hyperpolarizability.values` holds one record per
frequency: `frequency`, `process` (`static` or `shg`), `output_frequency`,
`components`, `converged` (every leg of every requested element met its
target), and the engine rows as `engine_evidence`. An element the engine
wrote no row for (its quadratic source needs two converged legs) is published
with a null value under `components` and the engine's reason under
`components_unavailable` (from the run's `run_summary.dropped_work`, which
names the undelivered source); the record and the stage are unconverged, as
for the polarizability. Two rows for one element remain an error. Other
processes (optical rectification, the electro-optic Pockels effect) have no
quadratic source in the engine yet and no deck spelling here.

**Raman.** A component `ab` is `d(alpha_ab)/dQ` at each optical frequency in
`raman.frequencies`, for each nuclear coordinate `Q`. The coordinates are
`raman.atoms` (atoms as the reference's molecule numbers them, from 0;
omitted: every atom) times `raman.directions` (`x`, `y`, `z`; omitted: all
three). With both omitted the
engine computes the full `3N` set, builds the Hessian from the same nuclear
response states, and passes its normal-mode analysis through as
`properties.raman.vibrational_analysis` and `properties.raman.raman_spectra`
(`coordinates: all`, `spectrum_complete`); the stage is converged only when
every derivative row was. A written subset gives derivative rows only
(`coordinates: written`, `normal_mode_analysis: false` with the reason).
Each record in `properties.raman.values` names its `coordinate`
(`nuc_<atom>_<axis>`), `atom`, `axis`, `frequency`, `components`,
`converged` and the engine rows; a component without a row is published as
for the hyperpolarizability (null, `components_unavailable`, unconverged).

## Excited states

A `molresponse` stage may request the lowest excited states of its reference
(D-025), alone or beside the optical properties, whose dipole legs it shares:

```
molresponse excited
  reference gs
  protocol [1e-4]
  dconv 1e-4
  maxiter 20
  requested_properties [excited_states, polarizability]
  excited_states.roots 1
  excited_states.method tda
  polarizability.frequencies [0.0]
  polarizability.frequency_units au
  polarizability.components [zz]
end
```

`excited_states.roots` (how many roots, from the lowest) and
`excited_states.method` (`tda`, the Tamm-Dancoff approximation, or `rpa`,
the full TDHF response) are required: a root count is a number, never
defaulted, and the two methods are different calculations. Optional:
`excited_states.guess` (`solid_harmonics`, `virtual_ao` or `random`; omitted,
the engine's `solid_harmonics`) and, for `virtual_ao` only,
`excited_states.guess_basis` (omitted, the engine's `aug-cc-pvdz`). The
engine's `excited.*` spellings are not stage controls. The stage's `maxiter`
caps the bundle's rung as it caps the dipole legs' (one cap); the engine's
warmup length and subspace size are its defaults. Two-photon absorption is
not exposed.

The engine solves the roots together as one bundle per rung, each root
with its own verdict (a stalled root is never relabelled converged), and
the bundle is the unit the checkpoint saves and the planner judges on
`--continue` (on the worst root). `properties.excited_states` is one record:
`method`, `roots_requested`, `guess` and `guess_basis` (as the engine used
them, `guess_written` says whether the deck wrote them), `protocol_step_key` of
the final rung, `stop_reason`, `iterations_at_protocol_step`, `decided_by` (the
engine's check, or the planner for a met load), `converged` (the bundle
ended through `converged` and every root is converged and not stalled), and
`roots[]`: per root the engine's stable `root_id` (the `slot` and
`display_name` can change between rungs and are kept beside it), `omega`
(`omega_units au`), `converged`, `locked`, `stalled`, `excitation_energy_change`,
`density_residual`, `bsh_residual`, the targets, the `unavailable` reasons
the engine writes, and the transition analysis (`transition`: oscillator
strength, transition dipole and quadrupole) or `transition_unavailable` with
the reason (the engine writes the analysis of a converged bundle only). No
bundle at the final rung gives a record without roots and `unavailable`
with the reason; a bundle of another method or root count at that rung is
an error. The bundle entry travels as `engine_evidence`.

The named decks `named_h2_es_tda.in` and `named_h2_es_rpa.in` are the named
twins of the `response_h2_es_*` qctest decks.

The response checkpoint publishes the bundle at each rung as one state: the
per-root archives under `es__<key>/`, with the bundle's index (`roots.json`,
with `writer_nproc`) and, for a converged bundle, its transition analysis
kept in the manifest. On `--continue` the bundle is planned like a dipole
leg, judged on its worst root against the engine's own check (every root's
density change below `5 x max(dconv, thresh)` and its excitation energy
change below the same); the continuation record's entry for it has
`kind: excited_states`, `method` and `roots`. A loaded bundle that meets
the request is recorded `decided_by: planner`. A saved bundle of another
method or root count is not used: the bundle starts fresh, and the record's
`why` says so. A converged bundle loaded without its transition analysis (a restart
between the bundle's save and the analysis) gets the analysis computed
from the saved states on the load, with no solver iterations (D-089).

**The bundle's stall rule** (D-099, B3). A bundle stops as stalled when an
unconverged root has improved neither its density change nor its eigenvalue
change by `stall.ratio` over `stall.window` iterations while the other active
roots are converged or likewise flat; the check is the engine's and applies in
both modes (it is not a written target). A stalled bundle is done at its rung:
`--continue` loads it and does not iterate it further. The keys are the
engine's own: `stall.window <iterations>` and `stall.ratio <fraction>` on a
stage that requests `excited_states`; unwritten, the engine's defaults 6 and
0.1 apply, and the settings table shows them as engine rows ("bundle stall
window", "bundle stall ratio"). `stall.window 0` disables the detector, so the
bundle runs to `maxiter` and `--continue` resumes it. On a stage without
`excited_states` the keys do nothing and are refused: the dipole legs' stall
rules are `progress.window` and `progress.minimum_relative_improvement`,
which act on the legs only and never on a bundle. A root that stalls at the
edge of a cluster of states (the warmup lists them) is usually told apart by
requesting more roots.
A bundle written on another number of ranks is refused
("ES bundle in es__<key> was written with N process(es) but loaded with M").

`econv` on a stage that requests `excited_states` is the target of the
excitation energy change per root (DESIGN D8: a written value only; the
engine's default `econv` is not a request): the eigenvalue change of each
root across its last step, judged on the worst root with `<=`, in hartree.
It is never the loop's target: the bundle's solve stops on the engine
check in both modes (D-092), whose own gate on the same quantity is `5 x
max(dconv, thresh)`, strict, per root (the `dconv` engine row of the settings
table), so no bundle rung runs to its cap for a target the loop does not
watch. With `convergence.use_default_criteria true` (the default) the written
`econv` is reported; with `false` it is enforced as the stage verdict (an
unmet `econv` leaves the stage unconverged), by the shared judge, not by the
loop. The settings table shows it as the requested row "excitation energy
change per root" (keyword `econv`, metric `excitation_energy_change`), and a stage that
wrote it in engine mode starts with the warning "econv 1.0e-06 is reported,
not enforced: the engine check stops the bundle when each root's eigenvalue
change is below 5.0e-04 (5 x max(dconv, thresh)), strictly, in both modes; a
written econv is judged by the stage verdict". In criteria mode the
deck check requires a target for every loop that would otherwise have none:
a stage that requests only `excited_states` needs `econv`; one that requests
a driven-response property needs one of `dconv`, `accuracy.bsh_residual`,
`accuracy.kinetic_trace_change` for its legs' loop, and may add `econv`. The
bundle's state record carries the worst root's
`excitation_energy_change` (the engine row the shared verdict reads, with its units,
as the saved measurement it judges the target on) (one judge, B1: `requested_accuracy` under the stage's
convergence report holds the status and, in engine mode, the "unmet, not
enforceable in engine mode" outcome). The record, `properties.excited_states`,
describes the request beside the engine's gate: `requested_omega_target` with
the `target` (null when the deck wrote none), `set_by`, `units`, `relation`,
the `worst` root (`value` and `root_id`; null when a root has no eigenvalue
step yet) and `verdict_in`, the bundle's state key under
`convergence_report.requested_accuracy` where the status lives; with a
written target each root carries `excitation_energy_change_within_requested`. The
record holds no status of its own. On `--continue` a saved bundle that
passes the engine check but not `econv` is loaded and reported unmet (B8).
On a `molresponse` stage without `excited_states`, `econv` is rejected: a
driven response has no local energy.

The named deck `named_lih_beta.in` requests the static `beta_zzz` of LiH on
the benchmark geometry. The response checkpoint publishes the nuclear legs
and the quadratic sources with the dipole legs, one property record per
requested property, and `--continue` judges every leg the requested values
are built from (above).

2PA, Nemo, MP2, CC2, and optimization are not exposed
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

## Dalton seed

A `molresponse` stage can take its starting guesses from a Dalton
linear-response output directory (S5, DESIGN B3):

```
  seed.dalton_dir dalton/
  seed.freq_tol 0.01
```

- `seed.dalton_dir` names the directory, relative to where madqc runs. It must hold
  `molden.inp` and `RSPVEC` loose (or a unique `*.molden.inp`/`*.molden` and
  `*.RSPVEC`), and may hold the Dalton `.out`; a directory with only a `.tar.gz` is
  refused on the named path. A `.out` naming another method than HF is refused.
- The excited-state bundle starts from the `RSPVEC` eigenvectors (the seed guard's
  overlaps are recorded); a driven leg starts from the `RSPVEC` record within
  `seed.freq_tol` au of its frequency. `seed.freq_tol` is required when the stage
  requests a driven-response property and refused otherwise; there is no default.
- A seed is a starting guess only: a state with a saved state starts from it
  instead (`--continue`), and a fresh state's plan says so in its `why`. The
  geometry in `molden.inp` must match the reference molecule; the engine refuses
  the stage before the response solve otherwise.
- Every generation the stage publishes holds a copy of the seed files under
  `seed/`, listed in the manifest's `seed` block by size and SHA-256 with what the
  engine recorded seeding from. A generation whose copies are missing or altered is
  not published, and is refused when read.
- The engine spellings `dalton.dir`, `io.dalton.dir` and `seed.es_warmup` are refused
  in a named stage. The ground state's seed from `molden.inp` is not available on the
  named path yet.

## Record names

The records name what they hold (D-106). Generations are written in
`madqc.state.v2`; a generation written before, in `madqc.state.v1` or the
response's earlier `madqc.response.state.v1`, is still read, as the current
keys, and is never rewritten. The engine's own metadata
(`response_metadata.json`, the bundle and leg records the engine writes)
keeps its spellings: it is raw evidence, and the records above it name it.

| Old key (v1) | Current key | Where |
|---|---|---|
| `omega_residual` | `excitation_energy_change` | the excited-state bundle's measurement and root rows, `properties.excited_states.roots[]`, the requested-accuracy criterion |
| `omega_residual_within_requested` | `excitation_energy_change_within_requested` | each root under a written `econv` |
| `delta_kinetic_energy` | `absolute_kinetic_trace_change` | measurement and target name in both applications (the deck key `accuracy.kinetic_trace_change` is unchanged) |
| `kind: fd` | `kind: frequency_dependent_response` | a driven response state's record |
| `engine_converged` | `converged` | a response stage's evidence about its reference (the reference stage's verdict, not the engine check) |
| `rung`, `iterations_at_rung`, `iterations_done_at_rung` | `protocol_step`, `iterations_at_protocol_step`, `iterations_done_at_protocol_step` | state records, moldft's `position`, the continuation record |
| `at_rung`, `beyond_rung` | `at_protocol_step`, `beyond_protocol_step` | the verdict's per-state entries |
| `protocol_key` | `protocol_step_key` | `properties.excited_states` |
| `channels[].protocol` | `channels[].protocol_step` | the response report's leg entries (the stage's sequence stays `protocol`) |
| `stop` | `stop_reason` | state records and the verdict's entries; values unchanged |

Unchanged on purpose: `requested_omega_target` (it names the omega target;
only its contents rename), `stop_reason` of the manifest, and every deck
keyword. A driven state's `identity` prints its frequency as the shortest
decimal that reads back to the same number (`dipole_z|0.04`, not
`dipole_z|0.040000000000000001`); an older generation's identity is read in
this form, so its states still match today's request.

Where the report is: a stage's `results.json` keeps this run's judgement
(`convergence_report`: `published_in`, `decided_by`, `requested_accuracy`,
`unmeasured_targets`); the scientific report is the published generation's
`metrics.convergence_report`, which `published_in` names (the manifest path and
the JSON pointer). A stage loaded from a generation that holds no report
keeps the report it assembled. The root `results.json` lists each stage with
its `status`, `work` and `path`.

What some names mean:

- `protocol`: the deck's sequence of thresholds; a `protocol_step` is one
  threshold and its polynomial order `k`.
- `thresh`, `k`: the representation threshold and the polynomial order of a
  protocol step. `maxiter` is the iteration cap per protocol step, per run.
- `dconv`: on a moldft stage, the total charge-density change target; on a
  response stage, the induced-density change target. Both are written
  targets only when the deck writes them, and they are not the engine's
  density gates, which derive from them.
- `econv`: the energy-change target of a ground state; on an excited-state
  stage, the excitation energy change target. A driven response has no
  energy target.
- `bsh_residual`: the bound-state Helmholtz fixed-point residual; the
  written target is its maximum over functions, which differs from the
  engine's aggregate per leg.
- `energy_change`: the absolute change of the total energy; a response stage
  records it as `not_applicable`.
- `density_residual`: the excited-state engine's per-root density change
  check, distinct from `charge_density_change`.
- `omega`: a root's excitation energy; the optical properties'
  `frequencies` are driving angular frequencies. Both are in atomic units.
- `reference`: the ground-state stage a response stage reads, not a
  starting guess; `seed` is a starting guess.
- `key`: a state record's lookup key (a ground state's fixed label; a
  response state's channel, protocol step and frequency). `identity`: the
  problem a state solves (a ground state's physical description; a driven
  state's perturbation and frequency; a bundle's method and root count).
- `kind`: `ground_state`, `frequency_dependent_response`, or
  `excited_states` (a saved bundle of roots).
- `decided_by`: in a state record, `engine_check` or `planner` (loaded under
  a met request); in a report, which rule decided the stage (`engine_check`
  or `requested_criteria`).
- `engine_measurements`, `saved_measurements`: the engine's deciding sample
  and the values measured on the saved state; their phases can differ.
- `phase`: when in the solve a value was taken.
- `availability`: `measured`, `unavailable` (no sample), `invalid`, or
  `not_applicable` (no such quantity here); never a pass or fail.
- `stop_reason`: `converged`, `iteration_limit`, `stalled` or
  `in_progress` for a state; a stage's own `stop_reason` also has
  `unfinished` and `unmet_requested_accuracy`.
- `requested_accuracy`: the written targets judged, with statuses
  `not_requested`, `satisfied`, `unmet`, `unavailable`, `invalid`; not an
  uncertainty of a property.
- `position`: where in the run a generation was written (phase, protocol
  step), not molecular coordinates; `numerical`: settings and
  implementation records; `metrics`: measurements and the report;
  `properties`: the requested results with units.
- `state_id`, `checkpoint_generation`, `generation`: after a publish, the
  same checkpoint name; `reference_generation`: the reference's.
- `latest`: the stage's published generation, not necessarily a converged
  one. `g-<stage>-<run>-<n>` names are opaque.
- `consumer_context`: what a dependent stage needs to rebuild the
  reference's problem.
- `root_id`, `slot`: a root's label and its place in the bundle, not a
  physical state identity.

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
| `accuracy.kinetic_trace_change` | `absolute_kinetic_trace_change` | absolute change of the unweighted kinetic trace (restricted alpha counted once), hartree |

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
when its check passed. One exception to "the written targets stop the run":
an excited-state bundle's solve stops on the engine check in both modes
(D-092); its written target, `econv`, is enforced as the stage verdict (see
"Excited states").

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

One verdict for every stage (D-046, D-077): a stage is converged when every
state its request needs (moldft: the ground state; molresponse: each
perturbation and frequency the properties need, and the excited-state bundle)
has a record at the request's final rung (or finer, when a finished state
was loaded under R3, D-085) that is not diverged and either was
loaded under a met request (`decided_by: planner`) or stopped `converged`
and, with the written targets deciding, its saved measurements meet them. The
records are the generation's `numerical.states`. A moldft state loaded
after its final rung stalled (R2) is not published again: the stalled
generation stays latest, the stage is unconverged and the next stage does
not run (D-091). The written targets are
judged once per state; a molresponse stage records the result at
`convergence_report.requested_accuracy` too, one entry per state key, in the
same form as moldft's.

Each state is judged on the written targets its kind has rows for (D-082),
the same set the planner reads on `--continue`: a perturbation-and-frequency
state on `dconv`, `accuracy.bsh_residual` and `accuracy.kinetic_trace_change`;
the excited-state bundle on `econv` alone; the moldft ground state on every
moldft target. A row of its kind that a record lacks is `unavailable` and
fails as before. A state none of whose kind's rows was written is judged on
its stop alone, and its `requested_accuracy` says `not_requested`. A written
target that no state the stage needs has a row for (`dconv` on a stage that
requests only `excited_states`) judges nothing: it is listed once at
`convergence_report.unmeasured_targets` and in `checkpoint_warnings`, not as
a failed state.

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
| `accuracy.kinetic_trace_change` | `absolute_kinetic_trace_change` | absolute change of the unweighted positive trace over those blocks (`hartree_times_response_amplitude_scale_squared`) |

These response quantities include the perturbation amplitude scale and must not
be read as a ground-state energy or density. A driven response has no local
energy, so there is no energy target; `accuracy.energy_change` is rejected
with that message, and `accuracy.density_change` is replaced by `dconv`. The
one exception is `econv` on a stage that requests `excited_states`: the
excitation energy change target of its bundle, reported and never enforced
(see "Excited states"); on any other `molresponse` stage `econv` is rejected
with the same message.

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
