# Phase 4 Validation - Iteration and Convergence Logic

Date:

- 2026-03-16

Related roadmap phase:

- Phase 4 — Iteration and Convergence Logic

Related notes:

- [`../phase_notes/phase_4_iteration_and_convergence_notes.md`](../phase_notes/phase_4_iteration_and_convergence_notes.md)

## Validation Scope

- Validate that the current Phase 4 implementation builds under the `xeonmax` environment.
- Attempt the focused excited-state scripted slice against the same build to determine whether
  runtime validation is possible inside the sandbox.
- Record the current support boundary if runtime validation is blocked by environment rather than
  solver code.

Not covered until a runtime-capable MPI environment is used:

- legacy per-iteration numerical parity
- convergence-decision parity against legacy
- subgroup/runtime compatibility after the Phase 4 changes
- unrestricted accelerated behavior

## Environment

- Repository root: `/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next`
- Build directory: `/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/build-xeonmax-sandbox-mpi`
- Execution mode: sandboxed `xeonmax` configure/build plus attempted scripted test slice
- Loader / environment script: `/gpfs/home/ahurtado/load_xeonmax.sh`
- Slurm job id if used: none
- Additional runtime baseline: `analysis/results/h2_excited/v4`

## Commands Run

```bash
source /gpfs/home/ahurtado/load_xeonmax.sh
cmake -S /gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next \
  -B /gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/build-xeonmax-sandbox-mpi \
  -G "Unix Makefiles" \
  -DCMAKE_BUILD_TYPE=Debug \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_TESTING=ON \
  -DENABLE_MPI=ON \
  -DENABLE_MKL=ON \
  -DENABLE_LIBXC=ON \
  -DENABLE_PCM=ON \
  -DENABLE_INTEGRATORXX=ON \
  -DCMAKE_C_COMPILER=/gpfs/software/openmpi/xeonmax/gcc13.2/4.1.6/bin/mpicc \
  -DCMAKE_CXX_COMPILER=/gpfs/software/openmpi/xeonmax/gcc13.2/4.1.6/bin/mpicxx \
  -DMPI_C_COMPILER=/gpfs/software/openmpi/xeonmax/gcc13.2/4.1.6/bin/mpicc \
  -DMPI_CXX_COMPILER=/gpfs/software/openmpi/xeonmax/gcc13.2/4.1.6/bin/mpicxx \
  -DMPIEXEC_EXECUTABLE=/gpfs/software/openmpi/xeonmax/gcc13.2/4.1.6/bin/mpiexec
```

```bash
source /gpfs/home/ahurtado/load_xeonmax.sh
cmake --build /gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/build-xeonmax-sandbox-mpi \
  --target madqc -j 96
```

```bash
source /gpfs/home/ahurtado/load_xeonmax.sh
export MAD_NUM_THREADS=4
ctest --test-dir /gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/build-xeonmax-sandbox-mpi \
  --output-on-failure -j 1 \
  -R "molresponse_excited_(metadata_smoke|restart_reuse|protocol_projection|guess_fallback|restart_mode_compat)"
```

## Tests / Checks

| Check | Result | Notes |
| --- | --- | --- |
| `xeonmax` configure in writable sandbox build tree | Passed | Required explicit MPI wrapper paths because `load_xeonmax.sh` did not expose `mpicc` / `mpicxx` cleanly to CMake discovery |
| `madqc` build under `xeonmax` | Passed | Build completed successfully after fixing the `IterationContract` declaration order in `ExcitedStateBundleSolver.cpp` |
| Focused excited-state scripted slice launch | Blocked by environment | OpenMPI singleton startup failed before `MPI_Init_thread` completed; solver runtime was not reached |
| Metadata smoke runtime validation | Blocked by environment | `test_molresponse_excited_metadata_smoke.py` aborted during OpenMPI startup |
| Same-protocol restart runtime validation | Blocked by environment | `test_molresponse_excited_restart_reuse.py` aborted during OpenMPI startup |
| Lower-protocol reprojection runtime validation | Blocked by environment | `test_molresponse_excited_protocol_projection.py` aborted during OpenMPI startup |
| Guess-fallback runtime validation | Blocked by environment | `test_molresponse_excited_guess_fallback.py` aborted during OpenMPI startup |
| Serial/subgroup artifact compatibility runtime validation | Blocked by environment | `test_molresponse_excited_restart_mode_compat.py` also lacked `mpirun` on `PATH` in the sandboxed loader environment |
| Real MPI dynamic restricted baseline (`analysis/results/h2_excited/v4`) | Failed to converge | Runtime reached Stage 2c, wrote restart/metadata artifacts, but both protocols exhausted the iteration budget |
| Workflow metadata parity in `v4` | Passed | `v4.calc_info.json` and `task_1/molresponse/response_metadata.json` carry the same final `excited_states` subtree |
| Bundle rotation stability in `v4` | Failed | `EXCITED_ROTATE_FAIL reason=sygvp_exception` appears repeatedly after the first rotation in both protocols |
| Lower-protocol snapshot reuse in `v4` | Not exercised | The stalled `1e-04` snapshot was rejected, so the `1e-06` protocol reused the guess archive after bundle reprojection instead |

## Variant Coverage

| Variant | Covered in this validation? | Notes |
| --- | --- | --- |
| `static_restricted` | Compile only | Runtime validation blocked before solver execution |
| `static_unrestricted` | No | Not targeted in this phase validation |
| `dynamic_restricted` | Yes, runtime baseline only | `v4` reaches the solver and persistence layers, but does not converge |
| `dynamic_unrestricted` | No | Not targeted in this phase validation |

## Metadata / Restart Checks

- `response_metadata.json` parity checked? `no`
- Workflow metadata parity checked? `no`
- Restart source verified? `no`
- Restart support mode verified? `no`

Reason:

- sandboxed OpenMPI failed during startup, so the scripted runs never reached Stage 2c output.

Additional runtime baseline:

- `response_metadata.json` parity checked? `yes`
- Workflow metadata parity checked? `yes`
- Restart source verified? `partial`
- Restart support mode verified? `partial`

Notes:

- In `v4`, both protocols record `restart_source=guess_archive`.
- The `1e-06` protocol logs `EXCITED_RESTART_SKIP ... reason=stalled_snapshot` for the
  `1e-04` protocol snapshot, then reprojects the guess bundle to the tighter protocol.
- The same runtime baseline shows that the current restart/metadata contract is functioning,
  but the numerical kernel is not converging.

## Results Summary

- Phase 4 compile validation passed under `load_xeonmax.sh` using a writable sandbox build tree.
- The current solver code is therefore compile-valid for the `xeonmax` configuration.
- The sandbox still does not allow the scripted MPI slice to execute, so job-backed validation is
  still required for closure.
- The `v4` dynamic restricted runtime baseline shows that the actual Phase 4 blocker is now
  numerical parity:
  - both protocols fail to converge
  - bundle rotation succeeds at iteration 1, then repeatedly fails with `sygvp_exception`
  - step restriction binds frequently and may be materially altering the update path
  - metadata persistence and restart sidecars are no longer the limiting issue

## Remaining Validation Gaps

- Run the focused excited-state scripted slice in an MPI-capable job environment, preferably under
  the previously used interactive workflow.
- Compare restricted-shell per-iteration energies, density changes, and relative residuals against
  the legacy solver on at least one small closed-shell case.
- Verify whether restoring the legacy overlap-subspace reduction ahead of `sygvp` eliminates the
  repeated post-iteration-1 bundle-rotation failures seen in `v4`.
- Reassess the explicit restricted-shell step-restriction policy against legacy behavior before
  treating remaining convergence differences as purely accelerator-related.
- Re-run metadata parity, restart reuse, lower-protocol reprojection, and restart-mode compatibility
  after that runtime validation environment is available.
