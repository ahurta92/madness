# Phase 3 Validation - Solver Initialization Path

Date:

- 2026-03-13

Related roadmap phase:

- Phase 3 — Solver Initialization Path

Related notes:

- [`../phase_notes/phase_3_solver_initialization_notes.md`](../phase_notes/phase_3_solver_initialization_notes.md)

## Validation Scope

- Validate fresh-start initialization behavior after the Phase 3 implementation lands.
- Validate consistency between fresh-start initialization and guess-archive reuse for the tested cases.

Not covered until implementation exists:

- full restricted-shell numerical parity
- unrestricted full-restart support
- owner-group execution semantics

## Environment

- Repository root: `/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next`
- Build directory: `/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/build-40core-tests`
- Execution mode: MPI scripted slice under `load_40core.sh`; serial smoke under `load_xeonmax.sh`
- Loader / environment script: `/gpfs/home/ahurtado/load_40core.sh`, `/gpfs/home/ahurtado/load_xeonmax.sh`
- Slurm job id if used: `1750034`

## Commands Run

```bash
srun --jobid=1750034 --overlap --nodes=1 --ntasks=1 bash -lc 'source /gpfs/home/ahurtado/load_40core.sh && gmake -C /gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/build-40core-tests CMAKE_COMMAND=/usr/bin/cmake madqc -j 4 && cd /gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/build-40core-tests && ctest --output-on-failure -R "molresponse_excited_(metadata_smoke|restart_reuse|protocol_projection|guess_fallback|restart_mode_compat)"'
```

```bash
srun --jobid=1750034 --overlap --nodes=1 --ntasks=1 bash -lc 'source /gpfs/home/ahurtado/load_xeonmax.sh && unset MADQC_LAUNCHER && cd /gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/build-40core-tests/src/apps/madqc_v2 && ./test_molresponse_excited_metadata_smoke.py'
```

## Tests / Checks

| Check | Result | Notes |
| --- | --- | --- |
| Fresh-start initialization behavior | Passed | `test_molresponse_excited_guess_fallback.py` verifies `EXCITED_FRESH_GUESS_SELECT`, `guess_bundle`, and stable ids/names after forced guess-archive reuse |
| Guess-archive consistency | Passed | `test_molresponse_excited_guess_fallback.py` compares guess-archive roots/names with first-run metadata and with fallback rerun results |
| Metadata parity | Passed | Covered by `test_molresponse_excited_metadata_smoke.py`, `test_molresponse_excited_restart_reuse.py`, `test_molresponse_excited_protocol_projection.py`, `test_molresponse_excited_guess_fallback.py`, and `test_molresponse_excited_restart_mode_compat.py` |
| Same-protocol restart | Passed | `test_molresponse_excited_restart_reuse.py` |
| Lower-protocol reprojection | Passed | `test_molresponse_excited_protocol_projection.py` |
| Serial/subgroup artifact compatibility | Passed | `test_molresponse_excited_restart_mode_compat.py` |

## Variant Coverage

| Variant | Covered in this validation? | Notes |
| --- | --- | --- |
| `static_restricted` | Yes | Covered by `test_molresponse_excited_restart_reuse.py`, `test_molresponse_excited_guess_fallback.py`, and `test_molresponse_excited_restart_mode_compat.py` |
| `static_unrestricted` | No | Not touched in Phase 3; still outside validated support |
| `dynamic_restricted` | Yes | Covered by `test_molresponse_excited_metadata_smoke.py` and `test_molresponse_excited_protocol_projection.py`; serial smoke under `load_xeonmax.sh` also used this path |
| `dynamic_unrestricted` | No | Not touched in Phase 3; still outside validated support |

## Metadata / Restart Checks

- `response_metadata.json` parity checked? `yes`
- Workflow metadata parity checked? `yes`
- Restart source verified? `yes`
  - `fresh_guess`
  - `current_protocol_snapshot`
  - `lower_protocol_snapshot`
  - `guess_archive`
- Restart support mode verified? `yes`
  - `full_bundle_resume` for the covered restricted variants

## Results Summary

- Phase 3 validation passed for the scoped initialization-path work.
- The fresh-start path now produces an explicit, restart-visible `guess_bundle` seed when the selected
  trial bundle can be promoted into the active typed response bundle.
- The forced guess-archive fallback path preserves root ids and state names for the covered restricted case.
- The `dynamic_restricted` and `static_restricted` paths remain restart-capable and metadata-consistent
  after the Phase 3 initialization changes.

## Remaining Validation Gaps

- No unrestricted variant validation was added in this phase.
- No legacy numerical parity benchmark was run; Phase 3 only validates initialization-path consistency
  and metadata/restart correctness.
- `load_xeonmax.sh` produced missing-module warnings in the interactive allocation, so only the serial smoke
  was validated there; the MPI scripted slice still relies on `load_40core.sh`.
