# Shared driven-response restart example

Helium HF, alpha at frequencies 0 and 0.05 au, over protocols 1e-4, 1e-5,
1e-6 (automatic polynomial orders 6, 8, 8). SCF completes its full protocol
before response starts. The selected response metrics are BSH residual, density
change and kinetic change. These tolerances are example settings, not defaults.

From this directory, with `madqc` on PATH:

```sh
# Solve two independent z-directed response states.
./run.sh
# Repeating the same request reuses both checkpoints, without solver evaluation.
./run.sh
# Tighten the response tolerances; continue only states that fail the new request.
./run.sh '--response=convergence.bsh_residual=1e-6;convergence.density_change=1e-6;convergence.kinetic_change=1e-6'
# Add a frequency. Existing qualified states are reused; 0.1 starts at 1e-4.
./run.sh '--response=dipole.frequencies=[0,0.05,0.1]'
# Request every tensor component at all three frequencies: nine driven states.
./run.sh '--response=dipole.directions=xyz;dipole.frequencies=[0,0.05,0.1]'
```

CLI overrides apply to that invocation only. Repeat tighter tolerance overrides
in later commands if they should remain in force, or edit the response block.
`response.maxiter` provides a fresh update allowance per state and protocol on
each restart. An interrupted state resumes its saved protocol; a completed
intermediate protocol advances to the next. States saved later in the ladder
skip earlier protocols even when newly requested states start there.

The final representation is required regardless of early convergence. Intermediate
BSH and density tolerances have the legacy floor `5 * protocol threshold`;
final tolerances are exactly those requested. Each protocol change converts the
functions, rebuilds operators, and discards KAIN and comparison history. Missing
comparison measurements are unavailable, never reported as zero.

Read `he-alpha/results.json`: `response` contains polarizability
`values` grouped by frequency, with named components and a separate availability
map. Each state's metrics appear once in `state_metrics`. Shared native
checkpoints and exact source maps live beside this file; response iteration
history remains under `he-alpha/analysis/response`. Only qualified final
states contribute property values; an unmet request exits unsuccessfully after
saving its endpoint and availability information.

This deck illustrates independent driven states; the shared response layer also
covers coupled excited states and direct/property tasks. Compatible SCF
refinement can reuse response functions as guesses with fresh checks; model
changes are rejected. Frequencies must have distinct checkpoint IDs (currently
rounded to five decimal places). A different physical frequency that collides
with an independent state's ID is rejected. Refined 2PA frequencies use the
explicit replacement path instead. See the
[response milestone](../../../../docs/design/madqc-response-milestone.md) for
coverage and remaining scientific acceptance gaps.

For a previous development run still under `task_0/moldft`, stop its solver and
use the [offline relocation procedure](../../../apps/madqc/README.md#current-run-directory)
before restarting. The updated solver will not silently start fresh beside it.

`test_shared_response_ladder.py` runs this deck, including selective reuse,
interrupted publication, new frequencies, all tensor components and serial/MPI
checks. Its numerical references remain the existing helium alpha reference;
this directory does not introduce a regenerated scientific reference.

In the recorded helium run, alpha_zz was approximately 1.322329 bohr^3 at zero
frequency and 1.325806 bohr^3 at 0.05 au. Tightening all three response tolerances
to 1e-6 reused the static checkpoint with no evaluation; the dynamic state needed
one update. The full xyz tensors have equal diagonals within 3e-6 and negligible
off-diagonal components. These observations document the example, not new
reference data or recommended tolerances for other molecules.
