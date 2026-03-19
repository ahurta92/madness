# Phase 4 Notes - Iteration and Convergence Logic

Date:

- 2026-03-16

Related roadmap phase:

- Phase 4 — Iteration and Convergence Logic

Related contract:

- [`../contracts/phase_4_iteration_and_convergence_contract.md`](../contracts/phase_4_iteration_and_convergence_contract.md)

Detailed comparison note:

- [`phase_4_legacy_excitedresponse_comparison.md`](./phase_4_legacy_excitedresponse_comparison.md)

## What Was Implemented

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp :: ExcitedBundleProtocolInput`
  now carries `dconv`.
- `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage(...)`
  now wires `response_params.dconv()` into the excited-state protocol input.
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` now defines a solver-local
  `IterationContract` that derives:
  - `density_target`
  - `relative_target`
  - `max_rotation`
  - `convergence_mode`
  - `accelerator_mode`
  - `accelerator_subspace`
  from the active protocol input and `ResponseVector` variant.
- Restricted-shell iteration now:
  - computes and stores per-root density-change and relative-residual diagnostics
  - creates one local `response_solver` accelerator per root when `maxsub > 1`
  - keeps explicit step-restriction plumbing available through
    `ResponseSolverUtils::do_step_restriction(...)` but leaves it disabled by
    default to match the legacy excited-state main update path
  - records accelerator and update metadata during restricted-shell updates
- Restricted-shell convergence now uses a dual gate:
  `max_density_change <= density_target` and
  `max_relative_residual <= relative_target`.
- `src/apps/molresponse_v2/ResponseRecord.hpp` and
  `src/madness/chem/MolresponseLib.hpp` now persist the new convergence and accelerator fields
  through the authoritative Stage 2c metadata path.
- `src/apps/madqc_v2/test_molresponse_excited_metadata_smoke.py` and
  `src/apps/madqc_v2/test_molresponse_excited_restart_reuse.py`
  now assert the new convergence/accelerator metadata contract.

## Files Changed

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/apps/molresponse_v2/ResponseRecord.hpp`
- `src/madness/chem/MolresponseLib.hpp`
- `src/apps/madqc_v2/test_molresponse_excited_metadata_smoke.py`
- `src/apps/madqc_v2/test_molresponse_excited_restart_reuse.py`

## Key Symbols

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: IterationContract`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: make_iteration_contract`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: configure_iteration_contract`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_typed_bundle_legacy_sequence`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_state_from_potentials`
- `src/apps/molresponse_v2/ResponseRecord.hpp :: record_excited_protocol_result`
- `src/madness/chem/MolresponseLib.hpp :: protocol_result_to_json`
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate`
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::update_response`

## Planned Focus

- Translate the restricted-shell runtime kernel and convergence gates far enough to support
  parity debugging against the legacy solver.
- Make `maxsub` behavior explicit instead of leaving it as inert input.
- Publish machine-readable convergence diagnostics without changing the Stage 2c orchestration
  model or the Phase 1/2 contracts.

## Actual Design Decisions

- Kept the current typed-bundle architecture and layered the translated convergence contract
  into `ExcitedStateBundleSolver.cpp` rather than trying to reproduce legacy `X_space`
  structure directly.
- Used a solver-local `IterationContract` instead of scattering threshold/accelerator policy
  across the loop body so the active convergence semantics are explicit and serializable.
- Activated restricted-shell `maxsub` behavior through per-root local `response_solver`
  objects rather than a global bundle-level accelerator object.
- Kept unrestricted variants explicitly out of the accelerated path:
  they still advertise scaffold/fallback modes instead of silently claiming parity.

## Deviations From The Original Plan

- The phase now exposes explicit convergence and accelerator metadata, but it has not yet been
  reconciled against legacy iteration-by-iteration output in a runtime-capable MPI environment.
- Restricted `maxsub` behavior is live in the current code path, but accelerator history is still
  solver-local to one run and has not been introduced as restart-visible state.
- This pass completed compile validation under `load_xeonmax.sh` in the sandbox, but the
  contract-required runtime validation still needs a job environment where OpenMPI can start.
- Real MPI runtime evidence now exists in
  `analysis/results/h2_excited/v4`, and it shows that the main Phase 4 blocker is numerical,
  not orchestration or metadata. The dynamic restricted solver path runs, writes restart artifacts,
  and persists metadata correctly, but both protocols fail to converge.

## New Constraints Discovered

- The shared external 96-core build tree is not writable from the sandbox, so `xeonmax` compile
  validation needs a writable local build tree inside the worktree.
- In the sandbox, `load_xeonmax.sh` does not leave `mpicc` / `mpicxx` on `PATH` reliably enough
  for CMake discovery, so explicit MPI compiler wrapper paths were required for local configure.
- Sandboxed OpenMPI singleton startup fails before `MPI_Init_thread` completes because the runtime
  cannot create the local OOB daemon sockets it expects. As a result, the scripted excited-state
  runtime slice still needs an interactive MPI-capable job environment for Phase 4 validation.
- The current bundle rotation path is less robust than legacy. In the `v4` dynamic restricted run,
  `diagonalize_excited_bundle(...)` begins throwing `sygvp_exception` after the first rotation,
  so most later iterations skip bundle rotation entirely and fall back to the non-rotated branch.
  See
  `analysis/results/h2_excited/v4.log.0`,
  `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: diagonalize_excited_bundle`,
  and
  `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::reduce_subspace`.
- Early `v4` debugging showed that explicit step restriction was a concrete
  parity deviation from legacy and was binding frequently. The current
  restricted-shell path now leaves that restriction disabled by default, which
  matches the legacy `if (false)` branch in
  `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::update_response`.

## Remaining Follow-Up

- Run the focused excited-state scripted slice in an MPI-capable environment and compare the
  restricted-shell per-iteration diagnostics against legacy behavior.
- Restore the missing legacy overlap-subspace conditioning before `sygvp` so bundle rotation stays
  available after iteration 1. The current new-code diagonalizer floors only the overlap diagonal;
  the legacy path also had a `reduce_subspace(...)` SVD screen for near-singular overlap modes.
- Revisit the current restricted-shell step-restriction policy. `v4` suggests it is binding often
  enough to materially change the iteration path relative to legacy.
- Decide whether any restart-visible accelerator history is needed beyond the current bundle-state
  restart contract.
- Keep unrestricted variants outside the claimed accelerated/parity boundary until they are
  implemented and validated explicitly.
- Reconcile the master roadmap from `In Progress` to `Completed` only after runtime validation.

## Validation Summary

- `madqc` now builds successfully under `load_xeonmax.sh` in a writable sandbox build tree.
- The focused excited-state scripted slice could be launched from that build tree, but all runtime
  cases failed before solver execution because sandboxed OpenMPI could not complete singleton
  startup.
- A real MPI baseline now exists in `analysis/results/h2_excited/v4` for the dynamic restricted
  path. That run confirms:
  - the Stage 2c runtime path is active
  - workflow metadata and `response_metadata.json` stay in sync
  - guess/protocol restart artifacts are written
  - both protocols still fail on the iteration budget
  - repeated `EXCITED_ROTATE_FAIL reason=sygvp_exception` was the main observed
    blocker before the later overlap-conditioning and step-restriction fixes

## Reconciliation Checklist

- [x] roadmap updated
- [x] status dashboard updated
- [x] active contracts unchanged and still point to Phase 4
- [x] validation record written
