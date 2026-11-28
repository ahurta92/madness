# MADQC — Multiresolution Quantum Chemistry

`madqc` is the public executable for the supported chemistry engines. It reads
one input file plus command-line overrides, builds a workflow, and calls the
numerical applications. Dependencies between applications belong to workflow
code. Dependencies between response states and properties belong to molresponse.

There is no named-stage input language or separate `moldft`, `nemo`,
`molresponse`, `cis`, `oep` or `cc2` executable. Engines not represented by madqc
retain their separate builds; see the [application inventory](../README.md).

## Running a calculation

```sh
madqc --wf=scf h2o.in
madqc --wf=response response_he_alpha.in
madqc --help=response
madqc --print_parameters=dft
```

Prefer the tested decks in [src/examples/qc](../../examples/qc/README.md).
`--workflow`, `--wf` and `--w` select the same workflows, in that precedence
order. Supported workflows are `scf`, `nemo`, `response`, `cis`, `oep`, `mp2` and
`cc2`; the last two require a GENTENSOR build. `--optimize` applies to `scf` and
`nemo`. Optimizer restart/replacement is deferred.

Input supports plain parameter blocks and JSON. Command-line group overrides,
for example `--dft="maxiter=30;convergence.energy_change=1e-6"`, take precedence
over the file. `--geometry=he`, XYZ geometry files, `--prefix=RunName`, help,
parameter listing and inspection remain available. Invalid named blocks or
retired selectors fail before calculations create output.

## SCF convergence and restart

The SCF application owns its functions, history and KAIN update. It evaluates
the functions it will save, then passes measurements to the shared controller.
For example:

```text
dft
    xc hf
    protocol [1e-4, 1e-6]
    maxiter 30
    convergence.criteria [bsh_residual, density_change, energy_change]
    convergence.bsh_residual 5e-4
    convergence.density_change 1e-4
    convergence.energy_change 1e-6
end
```

Tolerances are explicit metric requests, not aliases for the retired SCF
`dconv`/`econv` inputs. Protocol changes carry orbitals forward, rebuild the
operators and reset KAIN and comparison history. Missing history is reported
as unavailable, never as a zero change. The final requested representation must
be reached even if an earlier protocol meets the requested tolerances.

`restart auto` reuses a usable checkpoint whose measured metrics satisfy the
request; otherwise it continues with fresh checks. Each restart receives a
fresh `maxiter` update allowance. `--continue` additionally requires an existing
shared SCF checkpoint. `restart none` starts fresh. `restart read_only` and
`--restart_info=<run-directory>` are for inspection and plotting.
Checkpoint saves default to endpoints; `checkpoint.every_iterations N` enables
additional saves. Selected-checkpoint damage is an error, not a silent fallback
to an analysis file. Shared restart archives use the native format.

The response workflow first obtains a qualified SCF checkpoint. The response
engine loads that exact source. DALTON ground-state seeding uses the same SCF
setup for `--wf=scf` and `--wf=response`; an existing shared restart wins over a
seed. Configure `io.dalton.dir` to use an import directory.

## Printed output and inputs

`--print_parameters=dft` and `--print_parameters=response` use aligned columns,
including long convergence names. Default and derived settings appear as
comments; explicitly supplied settings remain active input lines, so copying
the listing preserves which controls were requested. Long values remain on one
input line; their explanations continue on comment lines. Retired options are omitted from this
listing and rejected when supplied.

The response input no longer accepts named-stage accuracy/progress controls,
`diagnostic.iterations_per_rung`, `stall.*`, `seed.start_rung`,
`seed.es_warmup`, `excited.guess_max_iter`, or `excited.maxiter`. Choose the protocol explicitly and use
`response.maxiter` for the response iteration allowance. `dconv` is also retired;
all response families use the explicit shared convergence names. The
[remaining input inventory](../../../docs/design/madqc-input-cleanup-inventory.md)
records the active controls that still need a migration decision.

For the shared SCF and driven-response paths, `print_level 1` includes aligned
metric tables, current limits, decisions and iteration step times. Both applications
use the same printer, with model-specific density and kinetic labels. Measured
metrics that are not selected for stopping appear as `reported only`; missing
comparison history appears as `unavailable`, never a zero change. Timings use
driver-rank elapsed seconds at existing completion boundaries; they add no
synchronization or equation evaluations. SCF reports preparation/evaluation/
decision together, checkpoint/report, and advance. Response reports evaluate,
decide, checkpoint/report, and advance. An endpoint has no advance time; a
reused checkpoint has no new iteration history; its current native load/check
interval is reported separately. These broad intervals include
associated diagnostics and are not fine-grained CPU profiles.

The final summary shows actual SCF completion and representation instead of
legacy `dconv` stamps, current calculation wall time and execution size, plus
response planning/solve/property-assembly times. `print_level 0` suppresses
iteration lines but keeps this summary. Detailed SCF operation timers remain
available at `print_level 3`, and energy/orbital detail at `print_level 2`.
Restart output shows saved measurements against current tolerances when they are
comparable, and states whether the functions are reused or need fresh checks. A
reuse prints only this invocation's load/check time. At level 3, SCF also shows
final-orbital overlap, eigensolve and transformation wall intervals; these are
parts of the broad iteration interval, not additional times to sum into it.
The [stdout preservation audit](../../../docs/design/madqc-stdout-preservation.md)
records the retained information and gaps from the SCF evaluation split.
Both SCF and response publish their measured intervals in the shared
`timings.json`, including reuse and failed invocations: one entry per launch,
identified by its start time, holding each application's steps.

The full engine output goes to stdout (your terminal or captured job log).
`<prefix>.out` receives an appended final summary; it contains the full run only
if the launcher also captures stdout there. A separate stdout/stderr job log
avoids mixing that stream with the separately written summary.

## Current run directory

For `--prefix=he-alpha --wf=response`, single-point SCF and response share:

```text
he-alpha.out                         appended human-readable summary
he-alpha.calc_info.json              existing combined workflow report
he-alpha/
  results.json                      requested properties and state metrics
  summary.txt                       results.json formatted for people
  timings.json                      SCF and response interval records, one entry per launch
  analysis/
    ground_state/                   SCF diagnostics, plots and import scratch
    response/                       excited-state transition analysis (es_analysis__<protocol>.json)
  states/
    ground_state/                   the saved ground state: description.json, history.jsonl, archive/
    virtuals/                       frozen-occupied virtual orbitals, when requested
    response/<node>/                one folder per driven state, excited-state set or second-order source
```

The same ground-state location is used by `--wf=scf`, so adding response keeps
its existing SCF work. `results.json` reflects the calculations published so
far. Restart qualification reads each node's description.json; analysis files
are diagnostics.
Optimizer and engines awaiting adoption retain their existing workflow layout.

Previous development runs under `task_0/moldft` are rejected before writing.
A run folder written by an older version, one with a `checkpoints/` folder, is
also refused before anything is written: use a new output prefix. It is never read or
migrated.

## Output and adoption status

All supported response requests use the shared lifecycle: driven states,
coupled TDA/RPA roots, VBC and alpha/beta/Raman/2PA properties. Explicit metric
inputs customize its convergence request. For example, with the helium HF SCF input at
`protocol [1e-4]` and `k 6`:

```text
response
    dipole.directions z
    dipole.frequencies [0.05]
    kain true
    maxiter 30
    convergence.criteria [bsh_residual, density_change, kinetic_change]
    convergence.bsh_residual 1e-5
    convergence.density_change 1e-5
    convergence.kinetic_change 1e-5
end
```

These are test-case tolerances, not calibrated defaults. The default BSH, density and applicable energy absolute tolerances are
provisionally `5e-4`; selecting kinetic change requires an explicit tolerance.
Each metric also accepts `.relative` (default zero). Full calibration is later work.
The residual is the maximum over the stored X/Y functions before KAIN. Density
change uses the kernel's induced density. Kinetic change uses the positive,
unit-weighted sum over those functions. Energy change is inapplicable to a
driven state because its frequency is an input.

Restart is automatic: satisfying saved measurements are reused without another
evaluation; otherwise the saved functions resume with fresh KAIN, comparison
history and `maxiter`. Checkpoints default to endpoints, with
`checkpoint.every_iterations N` available. A capped state is saved, but its
alpha row is unavailable and the command fails; qualified independent rows remain
available. Each state keeps its own saved progress: new states start first,
unfinished states resume their saved protocol, and completed intermediate states
advance. Qualified final states are reused even when other states need earlier
protocols. Protocol changes convert functions and reset KAIN/comparison history.
Intermediate BSH and density limits use the legacy floor `5 * threshold`; final
limits are the requested tolerances. The scheduler uses automatic polynomial
orders. Same-model SCF refinement keeps saved response functions as guesses
with fresh checks; geometry, box, XC and localization mismatches reject. General
conversion between occupied orbital frames remains deferred. Subworlds, imported
DALTON guesses and supported open-shell HF states use the same native layer.
Independent frequencies must have distinct task identities. Refined 2PA frequencies
reuse functions only through their explicit root provenance, with fresh checks
and recorded replacement of the old task.

See the [runnable helium example](../../examples/qc/response_shared_ladder/README.md)
for frequency scans, full tensors, tighter requests and adding states on restart.
Polarizability `values` is an array grouped by frequency, with `components` and
`availability` maps; each state's metrics appear once alongside these tensors.

SCF and response share the run-root records shown above. Measurements describe
the saved native functions; each property records its exact qualified sources.
An unfinished prerequisite is recorded and blocks dependent computation with a
nonzero exit status. Nemo and remaining adapters retain their existing records
and controls until adoption; their release gates remain in the
[design](../../../design_draft.md).
Existing optional HDF5 analysis/import code remains; shared HDF5 restart is a
future feature.

## Implementation and tests

- [madqc.cpp](madqc.cpp): runtime, help, inspection and exit handling.
- [WorkflowInput.h](WorkflowInput.h): collective outer-input validation.
- [WorkflowDriver.hpp](WorkflowDriver.hpp): workflow assembly and common SCF seed setup.
- [WorkflowBuilders.hpp](../../madness/chem/WorkflowBuilders.hpp): engine wiring.
- [SCFRun.cc](../../madness/chem/SCFRun.cc): shared SCF iteration/restart session.
- [madqc_adapter.hpp](../molresponse/madqc_adapter.hpp): response workflow adapter.

Build `madqc`, then run its tests with `ctest -R madness/test/madqc` and the
scientific decks with `ctest -L qctest`. Follow the repository thread-budget
instructions when running MPI or concurrent tests. The
[test migration index](../../../docs/design/madqc-driver-test-migration.md)
records every retired named test, its scientific replacement and remaining gaps.
See the [workflow guide](WORKFLOWS.md) for adding an engine.
