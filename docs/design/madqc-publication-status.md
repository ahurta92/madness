# madqc development snapshot — November 20, 2026

This branch is an implementation snapshot. Scientific release acceptance is
still open. The active scope is workflow unification: convergence measurements,
native checkpoints and restart, results, console output and timings.

## Implemented

- SCF evaluation is separate from KAIN advancement. Shared decisions use metrics
  of the prepared orbitals, including the final orbitals saved to checkpoints.
  Ordinary and frozen-occupied SCF paths use the shared lifecycle.
- Molresponse integrates independent driven states, coupled excited-state roots,
  VBC and supported quadratic/derived properties. Applications own their state
  and numerical work; the shared layer handles measurements, decisions and
  publication. Dependency logic stays in molresponse.
- Restart preserves independent-state protocol progress, reassesses saved
  measurements against current requests, and performs fresh checks when needed.
  Selected native archives are checked strictly; incomplete prerequisites block
  downstream calculations.
- SCF/response single-point runs publish results, timings and checkpoints at the
  run root. Application analysis lives in separate subdirectories. The explicit
  offline migration tool relocates older development layouts.
- The public driver is madqc. Named stages and duplicate executables for engines
  already routed through madqc are retired. Other engines retain separate builds.

See the [response milestone](madqc-response-milestone.md),
[test index](madqc-test-index.md) and
[run layout](../../src/apps/madqc/README.md#current-run-directory).

## Remaining work

1. **Investigate the tight-protocol helium SCF result below.** Do not count this
   as a successful response refinement or relax the criteria to hide the failure.
2. **Resolve the default cold-water excited-state acceptance case.** The current
   solid-harmonic guess converges to higher roots and fails the retained energy
   reference. DALTON-seeded water and an explicitly selected virtual-AO variant
   pass; neither establishes acceptance of the default guess.
3. **Adopt the remaining in-scope engines.** Nemo ground-state and CPHF work needs
   Florian's [metric confirmation](nemo-metric-confirmation.md), then shared
   evaluation, restart, properties and timing integration. Legacy TDHF/CIS and
   OEP also need their adapters and scientific coverage. Require evidence for
   each supported OEP model, not only an OAEP smoke case.
4. **Finish release validation.** Cover supported properties, electron counts,
   spin cases, functional/kernel variants, seeds and parallel modes. Validate
   default tolerances and compare costs at equivalent settings. Keep independent
   numerical references and documented scientific limitations explicit.
5. **Finish input/documentation cleanup.** Reconcile historical inventories with
   the implemented interfaces and remove remaining redundant public controls
   after checking their active callers. The run-root migration is complete;
   older notes describing task_0/moldft are historical.

Optimizer redesign and restart/history remain deferred until workflow
unification is finished. Shared HDF5 and new 6D/correlated-engine adoption are
outside this increment. Existing workflow availability does not certify an
engine's adoption of the shared contract.

## Helium refinement at protocol 1e-8

Manual continuation of the helium HF alpha example from its converged 1e-7
checkpoint, with one MPI rank, two MADNESS threads and sequential BLAS. Both
protocols selected k=10. Existing selected limits were tightened tenfold;
SCF kinetic change was additionally selected at 1e-7, below the preceding
measured change of 2.62e-7. Response kinetic change was already selected.

| SCF metric | Requested absolute limit | Final saved measurement | Outcome |
| --- | ---: | ---: | --- |
| Maximum BSH residual | 1e-7 | 3.417663921418733e-7 | Unmet |
| Maximum spin-density change | 1e-8 | 1.3450371398942312e-7 | Unmet |
| Energy change (Ha) | 1e-9 | 7.823846903676213e-9 | Unmet |
| Kinetic change (Ha) | 1e-7 | 1.965210438825693e-8 | Met |

The run exhausted 30 updates, saved its incomplete endpoint and exited with
status 1. Total reported wall time was 2959.9 seconds (49.3 minutes). Response
did not run, so there are no new 1e-8 polarizabilities. The current results mark
SCF energy unavailable with reason `criteria_unmet`. The residual fluctuated
above its target throughout; this test alone does not establish the cause.

Response was configured with absolute residual, density-change and kinetic-change
limits of 1e-7. These limits were not exercised because SCF did not qualify.
The earlier 1e-7 run completed both response frequencies, yielding alpha_zz
1.3223287206351033 at omega=0 and 1.3258057376316152 at omega=0.05 au (bohr^3).

Local evidence is under
`~/testing_madqc_design/he-alpha-layout.IZdqoY/comparison/protocol-1e-8/`:
`invocation.json` records the command, environment and exit status;
`before-1e-7/` and `after-1e-8/` preserve the full-precision result/timing records.
The sibling `tighter-1e-8.log` contains the full iteration history.
These are manual experiment records, not replacement regression references.

## Publication verification

The current madqc executable and selected native test targets built successfully
with Ninja. Thirteen focused CTest checks passed in 115.59 seconds, run serially
with `MAD_NUM_THREADS=2`, `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1` and
`OMPI_MCA_hwloc_base_binding_policy=none`:

- CLI occurrence validation and native checkpoint store, each serial and MPI2;
- iteration metrics, metric helpers, timing records and iteration restart;
- SCFRun and the shared response iteration driver;
- offline run-layout migration;
- public driver contract, serial and MPI2.

The staged diff passes `git diff --check`. Numerical reference files are
unchanged. Generated calculation output, Python bytecode, build directories and
local verification scratch files are excluded from the commit. Previous
integration evidence, including its numerical gaps, is indexed in the test
index. This snapshot does not claim a clean full release suite.
