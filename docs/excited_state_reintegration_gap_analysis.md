# Excited-State Reintegration Gap Analysis

Line references below were verified against the checkout at `/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next` on 2026-03-13. Statements are facts unless explicitly labeled `Inference:` or `Uncertainty:`.

## 1 Scope and Inputs

This document compares:

- the legacy excited-state solver in `src/apps/molresponse`
- the current partially reintegrated excited-state implementation in `src/apps/molresponse_v2` and `src/madness/chem/MolresponseLib.hpp`

Primary synthesized inputs:

- `docs/legacy_excited_state_report.md`
- `docs/current_excited_state_implementation_report.md`

Repository files re-checked directly for this gap analysis:

- `src/apps/molresponse/molresponse.cc`
- `src/apps/molresponse/ResponseBase.cpp`
- `src/apps/molresponse/ExcitedResponse.cpp`
- `src/apps/molresponse/x_space.h`
- `src/apps/molresponse/response_functions.h`
- `src/apps/madqc_v2/madqc.cpp`
- `src/apps/molresponse_v2/molresponse2.cpp`
- `src/madness/chem/WorkflowBuilders.hpp`
- `src/madness/chem/Applications.hpp`
- `src/madness/chem/Drivers.hpp`
- `src/madness/chem/MolresponseLib.hpp`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/apps/molresponse_v2/ResponseRecord.hpp`
- `src/apps/molresponse_v2/ResponseVector.hpp`
- `src/apps/molresponse_v2/EXCITED_STATE_EXECUTION_CHECKLIST.md`
- `src/apps/molresponse_v2/EXCITED_STATE_REINTEGRATION_PLAN.md`

This report is intentionally implementation-oriented. It does not propose a general redesign of `molresponse_v2`, and it does not recommend reintroducing the legacy `X_space` container as a public architectural model.

## 2 Executive Summary

The legacy solver is a protocol-looped bundle solver. It constructs all requested roots together inside one `X_space`, refines an oversized initial trial subspace, then iterates a coupled excited-state update that repeatedly:

- computes response potentials,
- diagonalizes a generalized bundle eigenproblem,
- rotates the bundle,
- applies BSH updates,
- evaluates residuals and density changes,
- optionally applies KAIN,
- and stops on density plus relative residual convergence.

The new code already has the workflow plumbing required to host that behavior in `molresponse_v2`:

- an explicit Stage 2c excited-state stage in `molresponse_lib::solve_all_states`
- protocol planning and protocol metadata
- stable root descriptors and slot-permutation tracking
- typed restart snapshots and lower-protocol restart reuse
- root naming and root manifests
- a restart-aware scaffold solver with restricted-shell bundle diagonalization and per-iteration restart saves

Phases 1, 2, and 3 of the reintegration plan are now complete. The foundational and initialization-path gaps identified in the original plan are resolved in the current code:

- stable excited-root identity is now represented by `ExcitedRootDescriptor`
- Stage 2c metadata writes are unified through `ResponseRecord2::record_excited_protocol_result(...)`
- final excited metadata is written both to workflow results and to on-disk `response_metadata.json`
- restart snapshots now encode typed bundle state, restart support mode, and explicit full-snapshot versus guess-only semantics
- fresh-start Stage 2c initialization now explicitly refines the trial space, sorts candidate roots, selects the lowest requested roots, and seeds a typed `guess_bundle` when the selected bundle is restart-capable

The remaining work is now concentrated in the numerical and downstream integration layers. The highest-value remaining pieces are:

- legacy-equivalent restricted-shell update-kernel, convergence, and KAIN behavior
- unrestricted/full-response completion beyond the current explicit guess-only restart boundary
- owner-group execution semantics
- property-stage consumption of finalized excited roots
- regression coverage against legacy numerical behavior

The recommended implementation strategy is to finish reintegration in place inside `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` and Stage 2c orchestration, while translating legacy bundle algorithms into the existing `ResponseVector`-based model. The bundle math should remain local to the excited solver; the public state model should remain manifest-driven and restart-aware.

## 3 Legacy Solver Call Graph and Control Flow

### 3.1 Narrative Control Flow

1. Driver entry:
   `src/apps/molresponse/molresponse.cc :: main` reads input via `initialize_calc_params(...)`, checks `response_parameters.excited_state()`, constructs `ExcitedResponse`, then calls `calc.solve(world)` and `calc.output_json()` (lines 60-101).

2. Object construction:
   `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::ResponseBase` prepares shared ground-state context, copies orbitals and orbital energies, initializes the XC functional, and allocates the excited-state bundle `Chi` as `X_space(world, num_states, num_orbitals)` (`ResponseBase.cpp`, constructor near lines 12-29; `x_space.h` line 33).

3. Protocol loop:
   `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::solve` loops over `r_params.protocol()`. At the first protocol it branches:
   - restart path: `load(world, r_params.restart_file())`
   - fresh path: `ExcitedResponse::initialize(world)`
   Each protocol then calls `check_k(...)` and `ExcitedResponse::iterate(world)` (`ResponseBase.cpp` lines 1741-1780).

4. Fresh initialization:
   `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::initialize` creates an oversized trial bundle with `2 * num_states` rows, chooses one guess generator, projects against occupied space, orthonormalizes, refines guesses via `iterate_trial`, sorts by `omega`, selects the lowest roots into `Chi.x`, zero-initializes `Chi.y`, and writes `guess_restart` (`ExcitedResponse.cpp` lines 7-115).

5. Trial refinement:
   `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate_trial` repeatedly computes TDA transition densities, calls `ResponseBase::compute_response_potentials(..., "tda")`, diagonalizes the current trial bundle with `rotate_excited_space`, forms `theta_X`, builds BSH operators, applies them, reprojects, and re-orthonormalizes (`ExcitedResponse.cpp` lines 479-677; `ResponseBase.cpp` line 1116).

6. Main production iteration:
   `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate` runs the main iteration loop. After the first iteration it checks density and relative BSH residual convergence; otherwise it calls `update_response(...)`, updates densities and residuals, and either continues or saves/breaks (`ExcitedResponse.cpp` lines 2032-2373).

7. Numerical update kernel:
   `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::update_response` is the legacy numerical core. It:
   - computes `Lambda`, `V0`, and `gamma` with `ResponseBase::compute_response_potentials(...)`
   - rotates the bundle with `rotate_excited_space(...)`
   - forms `theta_X = rotated_v_x - rotated_E0X + rotated_gamma_x`
   - applies `bsh_update_excited(...)`
   - computes residuals with `ResponseBase::update_residual(...)`
   - optionally applies `ResponseBase::kain_x_space_update(...)`
   - normalizes and truncates
   (`ExcitedResponse.cpp` lines 2375-2477; `ResponseBase.cpp` lines 1116, 1440, 1570).

8. Bundle diagonalization:
   `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::rotate_excited_space` and `ExcitedResponse::excited_eig` build and solve the generalized bundle eigenproblem, rotate the state bundle, and sort the roots (`ExcitedResponse.cpp` lines 814-1084 and 964-1178).

9. Save/load and finalization:
   `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::save` and `load` serialize the full excited bundle (`omega`, `Chi.x`, optional `Chi.y`) using MADNESS archives (`ExcitedResponse.cpp` lines 2746-2825). Final JSON goes through `ResponseBase::output_json()` after `main` returns from `solve`.

### 3.2 Compact Call Graph

```text
src/apps/molresponse/molresponse.cc :: main
  -> src/apps/molresponse/global_functions.cc :: initialize_calc_params
  -> src/apps/molresponse/ExcitedResponse.hpp :: ExcitedResponse::ExcitedResponse
     -> src/apps/molresponse/ResponseBase.cpp :: ResponseBase::ResponseBase
  -> src/apps/molresponse/ResponseBase.cpp :: ResponseBase::solve
     -> src/apps/molresponse/ResponseBase.hpp :: ResponseBase::set_protocol
     -> [branch] src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::load
     -> [branch] src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::initialize
        -> make_random_trial / make_nwchem_trial / create_trial_functions2 / create_virtual_ao_guess
        -> src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate_trial
           -> src/apps/molresponse/ResponseBase.cpp :: ResponseBase::compute_response_potentials
           -> src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::rotate_excited_space
              -> src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::excited_eig
           -> create_shift
           -> create_bsh_operators
     -> src/apps/molresponse/ResponseBase.cpp :: check_k
     -> src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate
        -> [loop] src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::update_response
           -> src/apps/molresponse/ResponseBase.cpp :: ResponseBase::compute_response_potentials
           -> src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::rotate_excited_space
              -> src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::excited_eig
           -> src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::bsh_update_excited
           -> src/apps/molresponse/ResponseBase.cpp :: ResponseBase::update_residual
           -> [branch] src/apps/molresponse/ResponseBase.cpp :: ResponseBase::kain_x_space_update
        -> [branch] src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::save
  -> src/apps/molresponse/ResponseBase.cpp :: ResponseBase::output_json
```

### 3.3 Legacy Branch Points That Matter For Translation

- First protocol:
  `ResponseBase::solve` branches between `load(...)` and `initialize(...)`.
- Guess creation:
  `ExcitedResponse::initialize` branches among random, NWChem, Cartesian, and virtual-AO guesses.
- Protocol boundaries:
  `check_k(...)` reprojects stored functions when polynomial order changes.
- TDA versus full response:
  `ExcitedResponse::update_response` and `bsh_update_excited` branch on `r_params.tda()`.
- Convergence and save:
  `ExcitedResponse::iterate` branches on density and relative residual targets and saves on convergence or iteration exhaustion.
- KAIN:
  `update_response` only applies `kain_x_space_update` when `r_params.kain()` and `iter > 0`.

## 4 Current Solver Call Graph and Control Flow

### 4.1 Narrative Control Flow

1. Driver entry:
   `src/apps/madqc_v2/madqc.cpp :: main` and `src/apps/molresponse_v2/molresponse2.cpp :: main` both route response execution through `workflow_builders::add_response_workflow_drivers(...)`, which installs `ResponseApplication<molresponse_lib>` (`WorkflowBuilders.hpp` lines 75-84).

2. Response workflow entry:
   `src/madness/chem/Applications.hpp :: ResponseApplication::run` creates the task directory, enters `task_<n>/molresponse`, and calls `molresponse_lib::run_response(...)` (lines 252-265).

3. Planning:
   `src/madness/chem/MolresponseLib.hpp :: molresponse_lib::run_response` builds `GroundContext`, then calls `plan_required_states(...)`. `plan_required_states(...)` builds:
   - linear states via `StateGenerator`
   - derived-state requests via `DerivedStatePlanner`
   - linear scheduling via `StateParallelPlanner`
   - an `ExcitedStateBundlePlan` via `build_excited_state_bundle_plan(...)`
   (`MolresponseLib.hpp` lines 977-1049, 4302-4354).

4. Linear stage first:
   `solve_all_states(...)` executes the linear solver path first, either through `execute_subgroup_state_solve(...)` or `execute_serial_state_solve(...)`, then validates the final protocol state (`MolresponseLib.hpp` lines 3521-3568).

5. Excited Stage 2c:
   After linear completion, `solve_all_states(...)` calls `execute_excited_state_bundle_stage(...)` (`MolresponseLib.hpp` line 3552). This is the top-level excited-state entry point in the new architecture.

6. Stage 2c orchestration:
   `execute_excited_state_bundle_stage(...)`:
   - constructs `ExcitedBundleSolverConfig`
   - creates the adapter with `make_excited_state_bundle_solver_adapter(...)`
   - ensures `state_metadata_json["excited_states"]` and per-protocol placeholder nodes exist
   - skips protocols only when metadata indicates a restart-capable saved bundle (`saved && converged && bundle_state_present && restart_capable`)
   - builds `ExcitedBundleProtocolInput` and dispatches `solve_protocol(...)`
   - writes normalized energies, names, residuals, timings, root manifests, restart provenance, and slot permutations through `ResponseRecord2::record_excited_protocol_result(...)`
   (`MolresponseLib.hpp` lines 2784-3092).

7. Concrete excited solver:
   `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: RestartAwareExcitedScaffoldSolver` contains the current implementation. Fact: despite the retained scaffold-oriented class name, Phases 1-2 have already made it the authoritative implementation for stable root metadata and typed restart/archive handling.

8. Protocol execution:
   `ExcitedProtocolWorkflow::solve_protocol(...)`:
   - reads the current protocol restart snapshot
   - skips only if metadata already describes a restart-capable saved protocol bundle
   - calls `prepare_protocol(...)`
   - calls `initialize_protocol_guess(...)`
   - calls `iterate(...)`
   (`ExcitedStateBundleSolver.cpp` lines 357-399).

9. Guess initialization and restart choice:
   `initialize_protocol_guess(...)` chooses restart input in this order:
   - current protocol snapshot
   - nearest lower protocol snapshot
   - generic guess archive
   - carryover in-memory guess
   - fresh guess from `build_fresh_guess(...)`
   It restores root descriptors, names, energies, residual history, and any persisted active bundle state via `load_restart_seed(...)`. For full-bundle same-protocol and lower-protocol restart, it then reprojects the restored bundle with `align_active_bundle_protocol(...)` before synchronizing the trial space (`ExcitedStateBundleSolver.cpp` lines 2805-2998).

10. Fresh guess path:
   `build_fresh_guess(...)` now chooses a generator based on `tda` and `num_states`, runs `iterate_trial(...)`, explicitly sorts and truncates the trial space through `sort_trial_space(...)` and `select_lowest_trial_roots(...)`, retries with a larger random trial when the initial generator under-fills the request, then seeds `active_response_bundle_` through `seed_active_bundle_from_trial_space(...)`. The path emits `EXCITED_FRESH_GUESS_SELECT` so the fresh-start selection and `guess_bundle` / `guess_trial_seed` distinction are visible in logs (`src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_fresh_guess`, `sort_trial_space`, `select_lowest_trial_roots`, `seed_active_bundle_from_trial_space`).

11. Main excited iteration:
    `iterate(...)` now starts from `build_response_bundle_seed(...)`, which prefers a restored active bundle when it is variant- and protocol-compatible and otherwise rebuilds from `trial_space_.x_states`. It then iterates up to `input.maxiter`. For restricted-shell variants it executes a translated bundle sequence:
    - `compute_bundle_response_potentials(...)`
    - `rotate_excited_bundle_states(...)`
    - `iterate_state_from_potentials(...)`
    For unrestricted-shell variants it skips bundle rotation and iterates each state independently with `iterate_state(...)`.
    It records per-state residuals, applies a single max-residual convergence gate, stores the active bundle in `active_response_bundle_`, then writes both a protocol snapshot and a generic guess archive. The snapshot now records `response_variant`, `restart_support_mode`, `snapshot_kind`, `bundle_state_present`, and `restart_capable` (`ExcitedStateBundleSolver.cpp` lines 3033-3228).

12. Finalization:
    `solve_all_states(...)` continues with derived-state execution and property assembly. Final metadata is returned through `ResponseApplication::results` and written into `*.calc_info.json` by the workflow layer, and Stage 2c also writes the same excited subtree into on-disk `response_metadata.json` through `ResponseRecord2` (`Applications.hpp` lines 271-275; `Drivers.hpp :: qcapp::Workflow::run`; `ResponseRecord.hpp :: ResponseRecord2::record_excited_protocol_result`).

### 4.2 Compact Call Graph

```text
src/apps/madqc_v2/madqc.cpp :: main
  -> src/madness/chem/WorkflowBuilders.hpp :: add_response_workflow_drivers
  -> src/madness/chem/Applications.hpp :: ResponseApplication<molresponse_lib>::run
     -> src/madness/chem/MolresponseLib.hpp :: molresponse_lib::run_response
        -> src/madness/chem/MolresponseLib.hpp :: plan_required_states
           -> StateGenerator::generateStates
           -> DerivedStatePlanner::build_vbc_driven_quadratic_plan
           -> StateParallelPlanner::build
           -> build_excited_state_bundle_plan
        -> src/madness/chem/MolresponseLib.hpp :: solve_all_states
           -> [branch] execute_subgroup_state_solve
           -> [branch] execute_serial_state_solve
           -> prepare_and_validate_final_protocol_state
           -> execute_excited_state_bundle_stage
              -> make_excited_state_bundle_solver_adapter
              -> ExcitedStateBundleSolver::solve_protocol
                 -> read_restart_snapshot
                 -> [branch] restart_ready_skip (only for restart-capable saved bundle)
                 -> prepare_protocol
                 -> initialize_protocol_guess
                    -> select_restart_seed
                    -> load_restart_seed
                    -> [branch] align_active_bundle_protocol
                    -> ensure_trial_space_matches_guess
                    -> [branch] build_fresh_guess
                       -> iterate_trial
                 -> iterate
                    -> build_response_bundle_seed
                    -> [restricted] iterate_typed_bundle_legacy_sequence
                       -> compute_bundle_response_potentials
                       -> rotate_excited_bundle_states
                          -> build_rotation_matrices
                          -> diagonalize_excited_bundle
                       -> iterate_state_from_potentials
                    -> [unrestricted] iterate_state
                    -> write_restart_snapshot
                    -> write_guess_archive
              -> ResponseRecord2::record_excited_protocol_result
           -> execute_derived_state_requests
        -> compute_requested_properties_with_property_group
```

### 4.3 Current Branch Points That Matter For Reintegration

- Serial versus subgroup linear solve:
  `solve_all_states(...)` branches into `execute_subgroup_state_solve(...)` or `execute_serial_state_solve(...)` before Stage 2c.
- Stage 2c protocol skip:
  `execute_excited_state_bundle_stage(...)` skips protocols only when the metadata says `saved && converged && bundle_state_present && restart_capable`.
- Restart seed choice:
  `select_restart_seed(...)` chooses current protocol, lower protocol, or guess archive; stalled snapshots are explicitly rejected, and same/lower snapshots must satisfy full-restart-capable bundle checks.
- Fresh versus carried guess:
  `initialize_protocol_guess(...)` branches between full bundle restart reuse, guess-archive reuse, carryover, and fresh guess generation.
- Restricted versus unrestricted:
  `iterate(...)` branches on the active `ResponseVector` variant; unrestricted bundle rotation is explicitly skipped with `EXCITED_ROTATE_SKIP reason=unrestricted_wip`, and restart support is advertised as `guess_only_unrestricted_variant`.
- Bundle rotation versus energy-estimate fallback:
  If bundle rotation fails or is skipped, `iterate(...)` falls back to `estimate_state_energies(...)`.

### 4.4 Current Manifest and Restart Read/Write Points

- Metadata plan seeding:
  `JsonStateSolvePersistence::initialize_excited_bundle(...)` -> `ResponseRecord2::initialize_excited_bundle(...)` (`MolresponseLib.hpp` lines 209-212; `ResponseRecord.hpp` lines 84-108).
- Per-protocol placeholder creation and authoritative result writeback:
  `MolresponseLib.hpp :: ensure_excited_protocol_placeholder_node(...)` and `ResponseRecord2::record_excited_protocol_result(...)` from `execute_excited_state_bundle_stage(...)` (`MolresponseLib.hpp` lines 2716-2752 and 2784-3092; `ResponseRecord.hpp` lines 269-297).
- Restart snapshot I/O:
  `ExcitedStateBundleSolver.cpp :: read_restart_snapshot(...)` and `write_restart_snapshot(...)`.
- Trial-state archive I/O:
  `read_restart_trial_states(...)` and `write_restart_trial_states(...)`.
- Active bundle archive I/O:
  `read_restart_response_bundle(...)` and `write_restart_response_bundle(...)`.

### 4.5 Resolved Stage 2c Persistence Split

Fact:
Stage 2c now normalizes per-protocol excited results and writes them through `ResponseRecord2::record_excited_protocol_result(...)`, which updates both the in-memory metadata tree and the on-disk `response_metadata.json` used by restart/persistence (`ResponseRecord.hpp` lines 269-297; `MolresponseLib.hpp` line 3136).

Implication:
the detailed excited metadata written by Stage 2c is now present in both workflow result metadata / `*.calc_info.json` and in the on-disk `response_metadata.json`. The remaining persistence concern is schema evolution and variant coverage, not a split write path.

## 5 Legacy vs Current Comparison Matrix

| Component | Legacy Behavior | Current Implementation | Status | Evidence | Notes |
| --- | --- | --- | --- | --- | --- |
| Solver entry point | Direct excited-state driver in `molresponse.cc`; `main` constructs `ExcitedResponse` and calls `solve()` when `excited_state()` is enabled. | Excited work is Stage 2c inside `molresponse_lib::solve_all_states`, reached through `madqc --wf=response` or `molresponse2`. | Partially Implemented | Legacy: `src/apps/molresponse/molresponse.cc :: main`; Current: `src/madness/chem/MolresponseLib.hpp :: run_response`, `solve_all_states`, `execute_excited_state_bundle_stage` | Plumbing exists, but the runtime model is now a workflow stage rather than a standalone solver object. |
| Initialization logic | `ExcitedResponse::initialize` creates an oversized trial bundle, projects, orthonormalizes, refines, sorts, and selects lowest roots into `Chi`. | `build_fresh_guess` plus `initialize_protocol_guess` now create `ExcitedTrialSpace`, run `iterate_trial`, explicitly sort/truncate to the lowest requested roots, and seed a typed active bundle before writing the guess archive. | Partially Implemented | Legacy: `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::initialize`; Current: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_fresh_guess`, `sort_trial_space`, `select_lowest_trial_roots`, `seed_active_bundle_from_trial_space`, `initialize_protocol_guess` | Phase 3 made fresh-start ordering and guess-archive seeding explicit and reproducible for the covered cases. The exact oversized-trial heuristic still differs from legacy. |
| Trial-space refinement | `iterate_trial` is a TDA bundle refinement loop using `compute_response_potentials`, subspace diagonalization, BSH updates, and re-orthonormalization. | `iterate_trial` exists in the new solver, now leaves consistent trial-space counts, and its refined roots are explicitly sorted, truncated, and preserved as the seed for `guess_bundle` restart state. | Partially Implemented | Legacy: `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate_trial`; Current: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_trial`, `sort_trial_space`, `select_lowest_trial_roots`, `seed_active_bundle_from_trial_space` | Phase 3 improved the initialization contract, but the retained state is still thinner than the full legacy bundle and exact numerical parity remains open. |
| Main iteration/update kernel | `update_response` computes bundle potentials, rotates the bundle, forms `theta_X`, applies BSH, updates residuals, and optionally applies KAIN. | Restricted-shell path uses `iterate_typed_bundle_legacy_sequence`, adds per-root `response_solver` accelerators when `maxsub > 1`, records step norms / relative residuals, and applies explicit step restriction through `iterate_state_from_potentials`; unrestricted path still skips bundle rotation and iterates roots independently. | Partially Implemented | Legacy: `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::update_response`; Current: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_typed_bundle_legacy_sequence`, `iterate_state_from_potentials`, `iterate_state`, `make_iteration_contract` | Phase 4 translated more of the restricted runtime kernel, but runtime parity against legacy trajectories is still unverified and unrestricted variants remain scaffold-style. |
| Convergence checking | Density residual and relative BSH residual are both checked; save/break occurs on convergence or iteration budget. | Restricted variants now use a dual gate on `max_density_change <= density_target` and `max_relative_residual <= relative_target`; unrestricted variants keep a fallback max-residual threshold. Convergence diagnostics are persisted in protocol metadata. | Partially Implemented | Legacy: `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate`; Current: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate`, `make_iteration_contract`; `src/apps/molresponse_v2/ResponseRecord.hpp :: record_excited_protocol_result`; `src/madness/chem/MolresponseLib.hpp :: protocol_result_to_json` | Phase 4 replaced the old scalar restricted-shell gate with an explicit dual-gate contract, but runtime parity against legacy convergence decisions still needs MPI-capable validation. |
| Bundle diagonalization / excitation energies | `rotate_excited_space` and `excited_eig` solve the generalized subspace eigenproblem and reorder roots. | `build_rotation_matrices`, `diagonalize_excited_bundle`, and `rotate_excited_bundle_states` implement the restricted-shell bundle rotation and sorting sequence. | Partially Implemented | Legacy: `src/apps/molresponse/ExcitedResponse.cpp :: rotate_excited_space`, `excited_eig`; Current: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: diagonalize_excited_bundle`, `rotate_excited_bundle_states` | The new restricted path is clearly derived from legacy logic. |
| Response-vector representation | All roots live together inside `X_space`, with `response_space x` and `response_space y` storing state-major bundles. | Runtime solving uses `std::vector<ResponseVector>` plus `ExcitedTrialSpace::x_states`, explicit variant helpers, and persisted `active_response_bundle_` state for restart-capable variants. | Partially Implemented | Legacy: `src/apps/molresponse/x_space.h :: X_space`, `src/apps/molresponse/response_functions.h :: response_space`; Current: `src/apps/molresponse_v2/ResponseVector.hpp :: response_variant_name`, `make_response_vector_from_variant`; `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: active_response_bundle_` | The replacement model is now variant-aware and restartable for restricted variants, but unrestricted numerical parity is still incomplete. |
| Excited-state container / root identity model | Root identity is implicit: each root is a row in the shared bundle and can move when the bundle is rotated or sorted. | Root identity is represented explicitly by `ExcitedRootDescriptor` (`root_id`, `stable_index`, `slot_index`, `energy`, `display_name`), and protocol results persist `slot_permutation` separately from display names. | Implemented | Current: `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp :: ExcitedRootDescriptor`; `src/apps/molresponse_v2/ResponseRecord.hpp :: record_excited_protocol_result` | Phase 1 resolved the foundational identity problem. Remaining work is downstream adoption of the descriptor contract. |
| State ordering | Sorting is bundle-global and directly mutates the state bundle. | Ordering remains energy-driven inside the runtime bundle, but the mapping from current slot to stable root is persisted explicitly through `slot_permutation` and normalized root descriptors. | Implemented | Legacy: `src/apps/molresponse/ResponseBase.cpp :: sort`; Current: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: rotate_excited_bundle_states`; `src/apps/molresponse_v2/ResponseRecord.hpp :: normalize_excited_protocol_result` | Phase 1 translated ordering into a restart-safe manifest contract without exposing bundle slot as identity. |
| Naming scheme | Mostly positional; analysis labels and JSON are not a formal root identity system. | `assign_excited_state_names` still derives `es1`, `es2`, `es3a`, ... from energy order, but names are now persisted as `display_name` on stable root descriptors and are not treated as identity keys. | Implemented | Current: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: assign_excited_state_names`; `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp :: ExcitedRootDescriptor` | Names are now derived labels layered on top of stable root ids. Degeneracy policy still matters for future property/UI consumers. |
| Save/load archive structure | `save` and `load` archive the full bundle: archive name, TDA flag, dimensions, `omega`, all `Chi.x`, and optional `Chi.y`. | Restart snapshot JSON now stores root descriptors, restart support mode, snapshot kind, protocol metadata, and trial-state metadata; active bundle data are archived separately through `read_restart_response_bundle(...)` / `write_restart_response_bundle(...)`. | Partially Implemented | Legacy: `src/apps/molresponse/ExcitedResponse.cpp :: save`, `load`; Current: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: RestartSnapshot`, `read_restart_snapshot`, `write_restart_snapshot`, `read_restart_response_bundle`, `write_restart_response_bundle` | Phase 2 implemented a typed full-bundle archive for supported variants. Unrestricted variants still do not claim full resume support. |
| Restart behavior | First protocol can restart from a saved bundle archive; later protocols rely on in-memory continuation plus `check_k`. | Restart still chooses current-protocol snapshot, lower-protocol snapshot, generic guess archive, or carryover guess, but same/lower protocol reuse now requires a restart-capable full bundle snapshot. Lower-protocol restart reprojects the restored bundle with `align_active_bundle_protocol(...)`. | Partially Implemented | Legacy: `src/apps/molresponse/ResponseBase.cpp :: solve`, `check_k`; Current: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: select_restart_seed`, `initialize_protocol_guess`, `align_active_bundle_protocol`, `load_restart_seed` | Phase 2 made restart semantics reliable for restricted variants while explicitly leaving unrestricted restart in guess-only mode. |
| Manifest / metadata system | JSON is mainly solver diagnostics (`response_base.json`). No formal excited-root manifest exists. | `excited_states.plan`, protocol nodes, timings, names, energies, residuals, `roots`, `slot_permutation`, and restart provenance fields are present and normalized through `ResponseRecord2`. | Implemented | Current: `src/apps/molresponse_v2/ResponseRecord.hpp :: initialize_excited_bundle`, `record_excited_protocol_result`; `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage` | The Stage 2c manifest contract is now explicit and restart-visible. Remaining work is consumer adoption, not schema absence. |
| Metadata write path | `ResponseBase::output_json` writes final JSON after solve. | Stage 2c now writes authoritative per-protocol results through `ResponseRecord2::record_excited_protocol_result(...)`, and the same subtree is returned in workflow metadata. | Implemented | Current: `src/apps/molresponse_v2/ResponseRecord.hpp :: record_excited_protocol_result`; `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage` | Phase 1 eliminated the split ad hoc write path. |
| Protocol handling | `ResponseBase::solve` loops protocols and calls `check_k(...)` between thresholds. | `ExcitedStateBundlePlan.protocols` drives Stage 2c; `align_trial_space_protocol(...)` reprojects trial states, and `align_active_bundle_protocol(...)` now reprojects restored active bundles for lower-protocol restart. | Partially Implemented | Legacy: `src/apps/molresponse/ResponseBase.cpp :: solve`, `check_k`; Current: `src/madness/chem/MolresponseLib.hpp :: build_excited_state_bundle_plan`, `execute_excited_state_bundle_stage`; `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: align_trial_space_protocol`, `align_active_bundle_protocol` | The protocol ladder and bundle reprojection are in place, but owner-group execution and numerical parity are still unfinished. |
| Owner-group scheduling | Legacy solver assumes the full world owns the excited bundle solve. | `owner_group` is planned and recorded in metadata, but the excited solver still runs on the full `World` and does not enforce owner-group execution. | Missing | Current: `src/madness/chem/MolresponseLib.hpp :: build_excited_state_bundle_plan`, `execute_excited_state_bundle_stage`; `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp :: ExcitedBundleProtocolInput` | This is a direct architecture gap relative to the new scheduler model. |
| KAIN / `maxsub` | Bundle KAIN is available through `ResponseBase::kain_x_space_update`; legacy `iterate` sets `set_maxsub(10)`. | Restricted-shell iteration now creates one local `response_solver` accelerator per root when `maxsub > 1`, records accelerator mode/subspace in metadata, and applies the accelerator during per-root updates. Unrestricted variants remain explicitly non-accelerated. | Partially Implemented | Legacy: `src/apps/molresponse/ExcitedResponse.cpp :: iterate`, `update_response`; Current: `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp :: ExcitedBundleProtocolInput`; `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: make_iteration_contract`, `iterate`, `iterate_typed_bundle_legacy_sequence` | Phase 4 made restricted `maxsub` behavior live, but it still needs runtime parity validation and does not widen unrestricted support. |
| Unrestricted / full-response support | Legacy code has `Y`-channel update paths, though the old report documents internal inconsistencies around `omega()`-gated `Y` handling. | Response variants are now explicit (`static_restricted`, `static_unrestricted`, `dynamic_restricted`, `dynamic_unrestricted`). Restricted variants advertise `full_bundle_resume`; unrestricted variants are persisted with typed metadata but explicitly marked `guess_only_unrestricted_variant`. | Partially Implemented | Current: `src/apps/molresponse_v2/ResponseVector.hpp :: response_variant_name`; `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: restart_support_mode_for_variant`, `make_state_response_vector`, `iterate` | Phase 2 made the support boundary explicit instead of silently pretending all variants restart the same way. |
| Property integration | Legacy has optional `analysis(...)`, but it is commented out from the main iteration path. | No property stage consumes excited-state roots, names, or archives. Existing property assembly uses linear and derived states only. | Missing | Current: `src/madness/chem/MolresponseLib.hpp :: compute_polarizability`, `compute_hyperpolarizability`, `compute_raman`; no excited-root consumer found | Downstream consumers need a stable excited-root contract first. |
| Output / finalization | Save archive plus `response_base.json` are produced directly by the solver path. | Excited results are returned through workflow metadata and restart files, and the same final excited subtree is now written into `response_metadata.json`. | Implemented | Legacy: `src/apps/molresponse/ExcitedResponse.cpp :: save`, `ResponseBase::output_json`; Current: `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage`; `src/apps/molresponse_v2/ResponseRecord.hpp :: record_excited_protocol_result` | Phase 1 resolved metadata parity for Stage 2c finalization. |
| Validation coverage | No dedicated regression suite was identified in the legacy code. | Focused regression tests now cover metadata parity, same-protocol restart reuse, lower-protocol reprojection, guess-archive fallback, and serial/subgroup artifact compatibility. Phase 4 also added assertions for convergence metadata, but sandboxed `xeonmax` runtime validation is blocked by OpenMPI singleton startup limits. | Partially Implemented | Current: `src/apps/madqc_v2/test_molresponse_excited_metadata_smoke.py`, `test_molresponse_excited_restart_reuse.py`, `test_molresponse_excited_protocol_projection.py`, `test_molresponse_excited_guess_fallback.py`, `test_molresponse_excited_restart_mode_compat.py` | The validation surface now checks the new convergence/accelerator metadata contract. Full runtime validation still requires an MPI-capable job environment. |

## 6 Architectural Translation Plan

The goal is to translate legacy solver concepts into the existing `molresponse_v2` architecture, not to recreate the old class layout.

| Legacy Concept | Legacy Assumption | Translation Into `molresponse_v2` | Why This Fits The New Design |
| --- | --- | --- | --- |
| `X_space` bundle | All roots, all channels, and bundle ordering live in one mutable container. | Keep bundle-local math inside `ExcitedStateBundleSolver` as `std::vector<ResponseVector>` plus bundle helper utilities; expose only manifest/root descriptors externally. | Preserves the new per-state and metadata-driven architecture while still allowing coupled bundle operations internally. |
| `response_space` row = one root | Root identity is the row index. | Represent one root as one `ResponseVector`, but attach explicit root metadata (`root_index`, stable name/id, current energy, ordering metadata) alongside it. | Array position can still be used for the current bundle, but restart/property interfaces stop depending on it. |
| Legacy positional ordering after `sort` / `transform` | Root identity changes whenever the bundle is rotated or re-sorted. | Track permutations explicitly at each diagonalization/sort boundary and write them into the excited manifest. | `molresponse_v2` needs restart- and property-safe root identity across protocol transitions. |
| `ResponseBase::solve` protocol loop | Excited solving owns its own top-level protocol loop. | Keep Stage 2c as the protocol loop owner for excited bundles; do not fold excited roots into the linear point scheduler. | This matches the current workflow separation: linear states first, excited bundle second, derived states third. |
| `rotate_excited_space` / `excited_eig` | Generalized subspace diagonalization is a bundle operation over all roots together. | Continue using `rotate_excited_bundle_states(...)` and `diagonalize_excited_bundle(...)` as the bundle-rotation layer; extend them to cover all supported response variants. | The restricted-shell translation already lives here; this is the right place to finish the port. |
| `bsh_update_excited` | Bundle rotation produces root-specific `omega`; BSH then updates each root with those shifts. | Finish the translation inside `iterate_state_from_potentials(...)` and related typed helpers, keeping bundle rotation separate from per-root operator application. | This maps well onto the current typed `ResponseVector` helpers. |
| Legacy restart archive | Restart writes the full active bundle (`omega`, `X`, optional `Y`) in solver-owned archive format. | Define a new excited bundle restart contract for `molresponse_v2`: manifest metadata plus enough active response data to resume exactly or deterministically. | The new architecture already has per-protocol snapshots and generic guess archives; they need to be made complete, not replaced. |
| `check_k(...)` | Protocol changes require explicit reprojection of stored functions. | Keep `align_trial_space_protocol(...)` and extend the same idea to any persisted full response bundle state. | This preserves protocol-aware restart in the new threshold ladder model. |
| Legacy solver-owned JSON | The solver directly owns its diagnostic JSON file. | Use `ResponseRecord2` / Stage 2c metadata as the single source of truth, with final writeback into both workflow results and restart-visible metadata files. | Fits the workflow-based persistence model already used by linear and derived stages. |

### 6.1 Recommended Translation Principles

- Do not restore `X_space` as a public or persistent data model.
- Do keep coupled bundle math local to the excited solver implementation.
- Separate "root identity" from "bundle slot index".
- Persist the minimum state needed for exact restart of the chosen algorithm, not just enough state to rebuild a guess.
- Keep Stage 2c as one bundle solve per protocol, owned by one scheduler lane or subgroup.

## 7 Gap Inventory

Phases 1, 2, and 3 are complete. The inventory below tracks the remaining gaps plus the residual follow-up exposed by those completed phases.

| Gap Name | Description | Affected Files / Classes | Dependency Relationships | Priority | Suggested Implementation Strategy |
| --- | --- | --- | --- | --- | --- |
| Stable root descriptor adoption by downstream consumers | Phase 1 introduced `ExcitedRootDescriptor` plus `slot_permutation`, but downstream property code and any future state-selection logic still need to treat `root_id` / `stable_index` as the authoritative identity instead of silently falling back to slot order. | `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp`, `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`, `src/madness/chem/MolresponseLib.hpp` | Foundation for property consumers and any future excited-state UI/API surface | Core | Keep the descriptor local to Stage 2c internals, but require all downstream consumers to read `roots` manifests and `slot_permutation` instead of positional bundle order. |
| Metadata schema hardening | Phase 1 resolved the split write path, but the now-authoritative excited metadata schema needs to remain stable as more consumers are added. The risk is no longer missing persistence; it is schema drift between protocol JSON, workflow metadata, and restart readers. | `src/madness/chem/MolresponseLib.hpp`, `src/apps/molresponse_v2/ResponseRecord.hpp` | Should stay aligned with restart/property work | Downstream | Extend the excited subtree only through `ResponseRecord2::record_excited_protocol_result(...)`, and version any future schema changes deliberately. |
| Variant-complete restart coverage | Phase 2 now persists typed active bundle state and supports full-bundle resume for restricted variants, but unrestricted variants are still explicitly `guess_only_unrestricted_variant`. | `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`, `src/apps/molresponse_v2/ResponseVector.hpp` | Depends on unrestricted numerical completion | Core | Either extend the bundle replay/reprojection path to unrestricted variants or gate those variants out of full-restart workflows until validated. |
| Restricted-shell numerical parity | The current restricted-shell path now has a translated iteration contract, explicit step restriction, and per-root accelerator hook, but the full update trajectory still has not been validated against legacy output iteration-by-iteration. | `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`; legacy references: `src/apps/molresponse/ExcitedResponse.cpp`, `src/apps/molresponse/ResponseBase.cpp` | Depends on restart/data model and initialization seeding being stable enough to debug and compare | Core | Compare the translated restricted-shell iteration path against legacy per-iteration energies/residuals in an MPI-capable validation environment and close any remaining ordering/update mismatches. |
| Convergence and iteration diagnostics parity | Legacy solve used density plus relative residual gates and recorded per-iteration solver diagnostics. Phase 4 now publishes density/relative residual diagnostics and a dual-gate contract for restricted variants, but those decisions are not yet runtime-validated against legacy behavior. | `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`, `src/apps/molresponse_v2/ResponseRecord.hpp`, `src/madness/chem/MolresponseLib.hpp` | Depends on restricted-shell numerical parity | Core | Validate the new `density_relative_dual_gate` metadata and convergence outcomes against legacy runs, then tighten any remaining contract differences. |
| Excited KAIN / `maxsub` translation | Restricted-shell iteration now creates per-root `response_solver` accelerators and records the active accelerator contract, but the behavior is only compile-validated so far and unrestricted variants remain intentionally non-accelerated. | `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp`, `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` | Depends on the translated residual/update kernel and MPI-capable runtime validation | Core | Validate restricted `maxsub` behavior against legacy expectations, decide whether any restart-visible accelerator state is needed later, and keep unrestricted variants explicitly outside the accelerated support boundary until implemented. |
| Unrestricted and non-TDA completion | Unrestricted bundle rotation is currently skipped, and rebuilt dynamic/unrestricted states start with zero `y` / beta blocks. | `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`, `src/apps/molresponse_v2/ResponseVector.hpp` | Depends on full restart state and restricted-shell parity | Core | Finish typed bundle rotation/update for unrestricted variants and persist all active channels needed for restart. |
| Owner-group execution semantics | The plan carries `owner_group`, but Stage 2c still executes on the full world and does not bind the excited bundle to one group/lane. | `src/madness/chem/MolresponseLib.hpp`, `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp`, `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` | Depends on metadata and restart behavior because subgroup ownership changes output paths and merge behavior | Core | Execute one excited bundle task per protocol on the designated owner subgroup or serial lane, and record that ownership in metadata and restart artifacts. |
| Root naming and protocol-state matching stability | Current names are protocol-local energy-order names. Phase 3 now preserves them across no-op guess-archive reuse for the covered restricted case, but cross-protocol mapping and degeneracy policy are still not fully defined. | `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`, `src/madness/chem/MolresponseLib.hpp` | Depends on stable root descriptor/permutation tracking | Core | Make naming a derived view of stable root identity. Preserve names across restarts and protocol transitions unless an explicit re-mapping step says otherwise. |
| Property-stage consumption of excited roots | No property code currently reads excited roots, root manifests, or excited archives. | `src/madness/chem/MolresponseLib.hpp`, downstream property code, and any future property manager layer | Depends on stable root descriptors plus complete persistence | Downstream | Define a clear property input contract: where finalized excited roots live, how they are named, and how response data are loaded. Then add the first consumer. |
| Numerical validation and regression harness | Current tests cover metadata/restart/projection plumbing plus fresh-start/guess-archive consistency for restricted variants, but not solver parity. | `src/apps/madqc_v2/test_molresponse_excited_*.py`, future benchmark fixtures | Depends on all core numerical work | Downstream | Add small parity fixtures first, then benchmark fixtures that compare against legacy/reference outputs across protocol ladders and restart paths. |

## 8 Ordered Implementation Roadmap

### Phase 1 - Foundational Data Structures

Status: Completed

Purpose:
establish stable excited-root identity and a single metadata contract before changing the numerical kernel.

#### Implementation Summary

- `ExcitedRootDescriptor` now lives in `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp` and separates `root_id` / `stable_index` from `slot_index` and `display_name`.
- Stage 2c results persist `roots` plus `slot_permutation`, so bundle reordering no longer changes the public identity contract.
- `ResponseRecord2::record_excited_protocol_result(...)` is now the authoritative Stage 2c write path.
- The same final excited subtree is now written to both workflow metadata and `response_metadata.json`.

Actual design decisions / deviations from the original proposal:

- The descriptor stayed local to the excited-state solver header rather than moving into a broader shared state type.
- Human-readable names remain derived labels (`display_name`), not the identity key.

Follow-up work:

- Downstream consumers still need to adopt `root_id` / `stable_index` instead of relying on slot order.
- Degeneracy and presentation policy still need hardening for future property/UI use.

Implemented validation:

- no-op restart stability checks for root ids and names
- metadata parity checks between workflow metadata and `response_metadata.json`

Original planned tasks (completed in this phase):

- Introduce a lightweight excited-root descriptor for Stage 2c results and manifests.
- Record bundle-slot permutations when roots are rotated or re-sorted.
- Make state names a view over stable root identity rather than the identity itself.
- Decide where the descriptor lives:
  likely `ExcitedStateBundleSolver.hpp` result structures plus metadata helpers in `MolresponseLib.hpp`; optionally a shared type in `ResponseState.hpp` if broader reuse is needed.
- Unify Stage 2c metadata writes so there is one authoritative path for `excited_states`.

Likely affected files/classes:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/madness/chem/MolresponseLib.hpp`
- `src/apps/molresponse_v2/ResponseRecord.hpp`

Deliverables:

- stable root descriptor schema
- stable `roots` manifest contents
- deterministic root naming behavior across a no-op restart
- one metadata write path for excited protocols

Validation steps:

- add a test that runs one protocol, reloads metadata, and confirms root ids and names are unchanged
- confirm `response_metadata.json` and workflow metadata contain the same excited subtree

### Phase 2 - Archive and Restart System

Status: Completed

Purpose:
make Stage 2c restart semantics reliable enough for iterative numerical development.

#### Implementation Summary

- `RestartSnapshot` now records explicit restart contract fields, including `response_variant`, `restart_support_mode`, `snapshot_kind`, `bundle_state_present`, `restart_capable`, protocol metadata, root descriptors, and slot permutations.
- Active bundle state is round-tripped through typed sidecar archives via `read_restart_response_bundle(...)` / `write_restart_response_bundle(...)`.
- Restart precedence remains `current protocol -> lower protocol -> guess archive -> carryover -> fresh guess`, but same/lower protocol reuse now requires a full restart-capable bundle snapshot.
- Lower-protocol restart reprojection now acts on the restored active response bundle through `align_active_bundle_protocol(...)`, then synchronizes the trial space.
- Variant support is explicit: restricted variants support full-bundle resume; unrestricted variants are persisted with typed metadata but remain `guess_only_unrestricted_variant`.

Actual design decisions / deviations from the original proposal:

- The restart artifact is split between JSON snapshot metadata and a typed sidecar response-bundle archive rather than one monolithic archive file.
- Phase 2 stopped short of claiming unrestricted full restart support; the code records that boundary explicitly instead of silently reconstructing missing channels.

Follow-up work:

- Extend exact restart semantics to unrestricted variants if those modes are to remain exposed as full workflow features.
- Continue using the typed bundle snapshot as the basis for any future convergence/KAIN debugging so restart semantics stay exact for the supported paths.

Implemented validation:

- same-protocol restricted restart reuse
- lower-protocol restart with bundle reprojection
- guess-archive fallback separation
- serial/subgroup artifact compatibility

Original planned tasks (completed in this phase):

- Extend the excited restart snapshot so it captures the full state required to resume the chosen algorithm.
- Persist dynamic and unrestricted channels if the resumed algorithm depends on them.
- Keep current restart precedence:
  current protocol -> nearest lower protocol -> generic guess archive -> carryover -> fresh guess.
- Write final Stage 2c metadata back to restart-visible files, not only workflow metadata.
- Ensure subgroup and serial paths produce compatible excited restart artifacts.

Likely affected files/classes:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/apps/molresponse_v2/ResponseVector.hpp`
- `src/madness/chem/MolresponseLib.hpp`
- `src/apps/molresponse_v2/ResponseRecord.hpp`

Deliverables:

- revised restart schema
- restart round-trip support for the active excited-state bundle
- final writeback of excited metadata to on-disk metadata files

Validation steps:

- restart from a converged same-protocol snapshot
- restart from a lower-protocol snapshot after protocol reprojection
- restart in both serial and subgroup-enabled runs

### Phase 3 - Solver Initialization Path

Status: Completed

Purpose:
match the legacy initialization and guess-refinement behavior closely enough that the production solver starts from comparable subspaces.

#### Implementation Summary

- `iterate_trial(...)` now always writes back consistent `ExcitedTrialSpace` counts after refinement.
- `sort_trial_space(...)` and `select_lowest_trial_roots(...)` make fresh-start root ordering and truncation explicit before the solver commits to a protocol-local active bundle.
- `build_fresh_guess(...)` now retries with a larger random trial when needed, then seeds `active_response_bundle_` from the selected trial roots and emits `EXCITED_FRESH_GUESS_SELECT` to make the chosen fresh-start seed observable.
- `initialize_protocol_guess(...)` now reports fresh-start `snapshot_kind`, `bundle_state_present`, and `restart_capable` consistently with the seeded bundle state.
- `test_molresponse_excited_guess_fallback.py` now verifies the fresh-start marker, guess-archive snapshot contract, and stable ids/names across forced guess-archive reuse.

Actual design decisions / deviations from the original proposal:

- The implementation kept the existing guess generators and layered explicit refinement / sort / select behavior on top of them instead of enforcing a literal legacy `2 * num_states` bundle size in all cases.
- The fresh-start seed is now a real typed `ResponseVector` bundle for restart-capable variants, built from the selected trial-space roots through `make_state_response_vector(...)`.
- Phase 3 stopped at initialization-path behavior. It did not change the main numerical kernel, convergence logic, or `maxsub` / KAIN handling.

Follow-up work:

- Phase 4 still needs restricted-shell numerical-kernel parity, convergence semantics, and `maxsub` / KAIN behavior.
- Phase 5 still needs cross-protocol naming and degeneracy policy hardening.
- If later parity work shows the exact legacy oversized-trial heuristic matters numerically, revisit the generator sizing strategy then.

Implemented validation:

- interactive 40-core scripted slice covering metadata parity, same-protocol restart reuse, lower-protocol reprojection, guess-archive fallback, and serial/subgroup restart-artifact compatibility
- serial `load_xeonmax.sh` metadata smoke for the dynamic restricted path

Specific tasks:

- Port the legacy initialization sequence into the new model:
  oversized trial bundle, occupied-space projection, Gram-Schmidt/normalization, trial refinement, energy sort, lowest-root selection.
- Review the current guess selection heuristics and align them with the intended legacy cases.
- Reconcile fresh-start seeding with the now-fixed bundle snapshot contract (`guess_bundle` / `guess_trial_seed` during initialization and `protocol_bundle` after iteration).
- Keep `align_trial_space_protocol(...)` as the protocol-boundary projection mechanism and extend it to any additional persisted state.

Likely affected files/classes:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/apps/madqc_v2/test_molresponse_excited_guess_fallback.py`

Deliverables:

- documented fresh-start behavior
- protocol-consistent refined trial bundle
- reproducible initial root ordering

Validation steps:

- compare initial guessed root count/order with the legacy solver on a small closed-shell case
- confirm a fresh run and a restart from the generic guess archive seed the same root names/order

### Phase 4 - Iteration and Convergence Logic

Status: In Progress

Purpose:
restore the core legacy solver behavior inside the `ResponseVector`-based implementation.

Implementation summary:

- `ExcitedBundleProtocolInput` now carries `dconv`, and Stage 2c wires that value from
  `response_params.dconv()` into the excited-state solver path.
- `ExcitedStateBundleSolver.cpp` now defines a solver-local `IterationContract` that derives
  `density_convergence_target`, `relative_convergence_target`, `max_rotation`,
  `convergence_mode`, `accelerator_mode`, and `accelerator_subspace` from the active protocol
  threshold, `dconv`, `maxsub`, and `ResponseVector` variant.
- Restricted variants now publish a machine-readable dual-gate convergence contract
  (`density_relative_dual_gate`) and per-iteration diagnostics:
  `density_change_norms`, `relative_residual_norms`,
  `iteration_max_density_changes`, and `iteration_max_relative_residuals`.
- Restricted variants also create one local `response_solver` accelerator per root when
  `maxsub > 1`, making `accelerator_mode=kain_per_root` live in the current solver path.
  Unrestricted variants remain explicitly scaffold/fallback mode and do not claim equivalent
  accelerator support.
- `iterate_state_from_potentials(...)` now records relative residuals, step norms, optional
  accelerator application, and explicit step restriction before the updated state is committed.
- `ResponseRecord2` and the Stage 2c protocol-result JSON now persist the new convergence and
  accelerator fields through the same authoritative metadata path introduced earlier.
- Current validation state:
  the `madqc` target now builds successfully under `load_xeonmax.sh` in a writable sandbox
  build tree, and a real MPI baseline exists in `analysis/results/h2_excited/v4` for the
  dynamic restricted path. That run confirms the runtime solver path, restart sidecars, and
  metadata persistence are live, but both protocols fail on the iteration budget.
- The strongest current numerical blocker is bundle diagonalization stability. In `v4`,
  `diagonalize_excited_bundle(...)` succeeds at iteration 1 and then repeatedly throws
  `sygvp_exception`, so most later iterations bypass bundle rotation and fall back to the
  non-rotated update branch.
- A second concrete parity deviation is that the current restricted path applies explicit step
  restriction on nearly every large update, while the legacy excited-state update loop leaves
  `x_space_step_restriction(...)` disabled in the main path.

Specific tasks:

- Finish the restricted-shell translation of the legacy update kernel:
  `compute_response_potentials`, bundle diagonalization, `theta` assembly, BSH application, residual update, normalization, and any step restrictions.
- Port legacy convergence gates:
  density residual plus relative response residual.
- Translate KAIN or an equivalent accelerator and make `maxsub` live.
- Emit per-iteration machine-readable diagnostics needed for parity debugging.
- Treat unrestricted/non-TDA as a follow-on inside the same phase:
  do not claim the feature complete until typed bundle updates cover the supported variants.

Likely affected files/classes:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- possibly shared low-level helpers in `src/apps/molresponse_v2/ResponseSolver.hpp/.cpp`

Deliverables:

- restricted-shell translated iteration path with documented remaining parity gaps
- explicit convergence metadata
- live restricted-shell `maxsub` / accelerator behavior or an intentional removal of the input

Validation steps:

- compare per-iteration energies and residual trends against the legacy solver for at least one restricted TDA case
- confirm convergence decisions match legacy expectations at the same protocol thresholds
- rerun the focused excited-state scripted slice in an MPI-capable environment; sandboxed `xeonmax`
  currently only provides compile validation
- use `analysis/results/h2_excited/v4` as the current dynamic restricted baseline and close the
  remaining gap by addressing the repeated post-iteration-1 `sygvp_exception` failures and the
  current step-restriction behavior

### Phase 5 - Naming and Protocol-State Mapping

Status: Planned

Purpose:
make root ordering, naming, and restart continuity stable across protocol changes and repeated solves.

Specific tasks:

- Harden the already-implemented stable root descriptors for downstream use:
  confirm every protocol transition, skip path, and consumer reads `root_id` / `stable_index` instead of slot order.
- Define the long-term naming policy on top of stable ids, including degenerate and near-degenerate roots.
- Confirm that lower-protocol restart reuse and later property consumers preserve root identity instead of silently renumbering everything.
- Add explicit edge-case diagnostics when reordered roots would otherwise cause ambiguous presentation.

Likely affected files/classes:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/madness/chem/MolresponseLib.hpp`
- `src/apps/molresponse_v2/ResponseRecord.hpp`

Deliverables:

- stable naming/order policy
- protocol-aware root-mapping metadata

Validation steps:

- protocol ladder test with restart at each threshold
- degenerate or near-degenerate root naming test

### Phase 6 - Property Integration

Status: Planned

Purpose:
make finalized excited roots usable by downstream property code without exposing solver-internal bundle slots.

Specific tasks:

- Define the property-stage contract for excited roots:
  where root descriptors live, how energies are read, and how response data are loaded.
- Add the first concrete consumer of excited-state output in the property stage chosen by project priority.
- Ensure property consumers do not depend on internal bundle ordering.

Likely affected files/classes:

- `src/madness/chem/MolresponseLib.hpp`
- any downstream property manager or property helper used by the selected feature
- excited metadata / archive code as needed

Deliverables:

- one stable excited-root consumption API
- one downstream property path using finalized excited-state results

Validation steps:

- end-to-end workflow test proving that an excited-state run produces property-stage consumable data

### Phase 7 - Verification and Cleanup

Status: Planned

Purpose:
lock in behavior and remove scaffolding assumptions that are no longer correct.

Specific tasks:

- Preserve and extend the existing regression tests for:
  metadata persistence, same-protocol restart, lower-protocol restart, guess fallback, serial/subgroup artifact compatibility, subgroup owner-group execution, restricted-shell numerical parity, and any enabled unrestricted path.
- Compare against the legacy solver and/or trusted references on a small benchmark set.
- Clean up scaffold-only names and comments once numerical parity is achieved.
- Revisit compile-time options such as `MADNESS_ENABLE_LEGACY_EXCITED_BUNDLE_ADAPTER` and either implement or remove dead configuration paths.

Likely affected files/classes:

- `src/apps/madqc_v2/test_molresponse_excited_*.py`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/apps/molresponse_v2/CMakeLists.txt`

Deliverables:

- regression suite covering architecture and numerics
- cleaned-up excited-stage implementation

Validation steps:

- automated CI-style runs for the excited-state test matrix
- documented benchmark comparison results

## 9 Risks and Design Tensions

### 9.1 Legacy Bundle Math Versus New Public State Model

Problem:
legacy math assumes one mutable `X_space` bundle; `molresponse_v2` organizes runtime state around typed per-root vectors and metadata-driven workflow stages.

Why it conflicts:
directly restoring `X_space` would cut across the new planning, persistence, and property interfaces.

Mitigation:
keep bundle math internal to the excited solver as `std::vector<ResponseVector>` plus dense rotation helpers, but keep external interfaces descriptor- and manifest-based.

### 9.2 Positional Root Identity Versus Restart-Safe Root Tracking

Problem:
the core Phase 1 fix is in place, but new consumers can still regress to treating bundle slot order or `state_names` as identity.

Why it conflicts:
`molresponse_v2` now has a stable identity contract; violating it in downstream code would reintroduce the legacy positional-coupling problem by accident.

Mitigation:
keep `ExcitedRootDescriptor` and `slot_permutation` authoritative. Treat array order and `display_name` as transient presentation details only.

### 9.3 Partial Restart State Versus Exact Resumption

Problem:
Phase 2 resolved exact bundle restart for restricted variants, but unrestricted variants still do not advertise full resume support.

Why it conflicts:
the workflow now exposes explicit variant metadata, so pretending all four `ResponseVector` variants are equally restart-capable would be incorrect and would hide real numerical/state gaps.

Mitigation:
keep `restart_support_mode` and `restart_capable` explicit in metadata and dispatch logic. Extend full-bundle replay only when the unrestricted numerical path and required channels are validated.

### 9.4 Workflow Metadata Split Versus One Source Of Truth

Problem:
the previous split write path has been resolved, but Stage 2c now spans two persistence forms: normalized JSON metadata and typed sidecar bundle archives.

Why it conflicts:
schema drift between JSON metadata and bundle archive readers would make restart/debug behavior harder to reason about even if the write path remains unified.

Mitigation:
keep `ResponseRecord2` as the single JSON authority and evolve `RestartSnapshot` deliberately, with explicit versioning and compatibility checks.

### 9.5 Bundle Coupling Versus Owner-Group Scheduling

Problem:
the excited solver is a coupled bundle solve, but `molresponse_v2` scheduling is group- and lane-aware.

Why it conflicts:
splitting one root bundle across groups would break the bundle-coupled diagonalization/update steps.

Mitigation:
keep the excited bundle as one protocol task pinned to one owner group. Parallelize within that task only where mathematically safe.

### 9.6 Restricted Progress Versus Unrestricted Feature Surface

Problem:
the restricted-shell path is meaningfully ahead of the unrestricted path.

Why it conflicts:
claiming broad excited-state support before unrestricted completion will produce inconsistent behavior across spin cases.

Mitigation:
gate feature completeness by variant. Land restricted-shell parity first, then explicitly finish unrestricted/non-TDA support before enabling those modes in production workflows.

Phase 3 note:
fresh-start bundle seeding now promotes selected trial roots into typed `ResponseVector` state for all variants, but only restricted variants still advertise full restart-capable semantics. The unrestricted support boundary therefore remains explicit and should not be widened accidentally in later phases.

### 9.7 Restart-Capable Snapshot Semantics Versus Skip Logic

Problem:
Stage 2c now has multiple artifact strengths: full restart-capable protocol snapshots, guess-only archives, carryover state, and fresh guesses.

Why it conflicts:
if protocol-skip or restart-selection logic treats these artifacts as equivalent, the solver can silently skip required work or report misleading restart provenance.

Mitigation:
keep the selection rules explicit in metadata and code:
same-protocol skip requires `saved && converged && bundle_state_present && restart_capable`, while weaker artifacts remain seed-only paths.

### 9.8 Live Accelerator State Versus Restart-Visible State

Problem:
Phase 4 made restricted-shell `maxsub` behavior live through per-root `response_solver`
accelerators, but that accelerator history still lives only inside one solver run.

Why it conflicts:
if later phases start to assume restart should reproduce the same accelerator subspace history,
the current bundle-only snapshot contract will be insufficient even though the numerical state
itself is restart-capable.

Mitigation:
keep the current claim narrow:
restart snapshots are bundle-state restart-capable for the supported variants, not full
accelerator-history checkpoints. Only widen that contract if future parity work proves the
accelerator history must also be serialized.

## 10 Acceptance Criteria

### 10.1 Architecture

- The public `molresponse_v2` excited-state interface remains Stage-2c-based and manifest-driven.
- No public or persistent API depends on legacy `X_space`.
- Each excited root has a stable descriptor independent of its current bundle slot.

### 10.2 Solver Functionality

- Restricted-shell excited-state runs execute the full translated sequence:
  initialization, trial refinement, bundle rotation, per-root update, residual update, and convergence.
- If non-TDA or unrestricted modes remain exposed, their bundle update and restart paths are fully implemented.
- `maxsub` is either operational in the excited solver or removed from the public excited-state inputs.

### 10.3 Restart / Save-Load Behavior

- Same-protocol restart resumes from the saved excited-state bundle without silently degrading to a fresh guess.
- Lower-protocol restart plus reprojection works across a protocol ladder.
- Excited restart artifacts are sufficient for exact or documented deterministic resume of the supported algorithms.
- `response_metadata.json` and workflow result metadata contain the same final excited-state protocol data.

### 10.4 State Naming and Ordering

- Root identity survives bundle rotations and protocol transitions.
- Root names are stable across restart and protocol reuse.
- Degenerate or near-degenerate roots have deterministic naming behavior.

### 10.5 Protocol Compatibility

- Stage 2c runs once per protocol threshold in the excited plan.
- Cross-protocol reprojection is performed when required by `k` / threshold changes.
- `owner_group` is honored in subgroup-capable execution modes.

### 10.6 Property Integration

- Finalized excited roots are available to downstream property code through a documented API or manifest contract.
- At least one end-to-end property path consumes excited-state output without depending on internal bundle ordering.

### 10.7 Numerical Validation

- A restricted-shell benchmark set shows parity with the legacy solver or trusted references within documented tolerances.
- Regression tests cover metadata, restart, protocol reprojection, subgroup execution, and numerical convergence behavior.

## 11 Codex Implementation Brief

Objective:
complete the excited-state reintegration inside `molresponse_v2` by translating the legacy bundle algorithm into the current Stage 2c, `ResponseVector`-based architecture.

Ordered tasks:

1. Port the restricted-shell legacy numerical kernel, including convergence diagnostics and `maxsub` / KAIN behavior, on top of the now-explicit fresh-start seed.
2. Extend the typed bundle update/restart path to unrestricted variants, or explicitly gate those variants out of full-restart workflows.
3. Harden naming and protocol-state mapping policy on top of the already-implemented stable root descriptors.
4. Enforce `owner_group` execution semantics for the excited bundle stage.
5. Add downstream property consumption and numerical regression coverage.

Constraints:

- Do not reintroduce legacy `X_space` as the new public state model.
- Keep excited solving as a Stage 2c bundle task, not a linear state-point task.
- Preserve the existing `molresponse_v2` workflow layering and metadata conventions.
- Treat bundle-slot order as transient; stable identity must live in descriptors/manifests.

Expected deliverables:

- updated excited-state solver and metadata code
- restart-complete archive format for all supported variants
- stable root manifest and downstream consumer contract
- owner-group-aware Stage 2c execution
- regression tests for restart, protocol behavior, and numerical parity

Validation procedure:

- run the existing excited metadata/restart/projection tests
- add parity tests against the legacy solver on a small restricted-shell case
- validate same-protocol and lower-protocol restarts
- validate subgroup execution with `owner_group`
- validate at least one downstream property consumer

## 12 Code Reference Index

### Legacy Solver

- `src/apps/molresponse/molresponse.cc :: main`
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::ResponseBase`
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::solve`
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::compute_response_potentials`
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::update_residual`
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::kain_x_space_update`
- `src/apps/molresponse/ResponseBase.cpp :: sort`
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::initialize`
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate_trial`
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::rotate_excited_space`
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::excited_eig`
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate`
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::update_response`
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::bsh_update_excited`
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::save`
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::load`
- `src/apps/molresponse/x_space.h :: X_space`
- `src/apps/molresponse/response_functions.h :: response_space`

### Current Solver

- `src/madness/chem/MolresponseLib.hpp :: build_excited_state_bundle_plan`
- `src/madness/chem/MolresponseLib.hpp :: plan_required_states`
- `src/madness/chem/MolresponseLib.hpp :: solve_all_states`
- `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage`
- `src/madness/chem/MolresponseLib.hpp :: build_excited_root_manifest`
- `src/madness/chem/MolresponseLib.hpp :: JsonStateSolvePersistence::initialize_excited_bundle`
- `src/madness/chem/MolresponseLib.hpp :: ensure_excited_protocol_placeholder_node`
- `src/madness/chem/MolresponseLib.hpp :: compute_polarizability`
- `src/madness/chem/MolresponseLib.hpp :: compute_hyperpolarizability`
- `src/madness/chem/MolresponseLib.hpp :: compute_raman`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp :: ExcitedRootDescriptor`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp :: ExcitedBundleSolverConfig`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp :: ExcitedBundleProtocolInput`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp :: ExcitedBundleProtocolResult`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: RestartSnapshot`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: ExcitedTrialSpace`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: assign_excited_state_names`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: sort_trial_space`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: select_lowest_trial_roots`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: seed_active_bundle_from_trial_space`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_fresh_guess`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: RestartAwareExcitedScaffoldSolver`
- `src/apps/madqc_v2/test_molresponse_excited_guess_fallback.py`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: ExcitedProtocolWorkflow::solve_protocol`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: load_restart_seed`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: select_restart_seed`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: initialize_protocol_guess`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_fresh_guess`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_response_bundle_seed`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: align_active_bundle_protocol`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_trial`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: compute_bundle_response_potentials`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_rotation_matrices`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: diagonalize_excited_bundle`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: rotate_excited_bundle_states`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: make_state_response_vector`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_state_from_potentials`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_typed_bundle_legacy_sequence`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: read_restart_snapshot`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: write_restart_snapshot`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: read_restart_response_bundle`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: write_restart_response_bundle`
- `src/apps/molresponse_v2/ResponseVector.hpp :: ResponseVector`
- `src/apps/molresponse_v2/ResponseVector.hpp :: response_variant_name`
- `src/apps/molresponse_v2/ResponseVector.hpp :: make_response_vector_from_variant`
- `src/apps/molresponse_v2/ResponseSolverUtils.hpp :: align_response_bundle_protocol`
- `src/apps/molresponse_v2/ResponseRecord.hpp :: ResponseRecord2::initialize_excited_bundle`
- `src/apps/molresponse_v2/ResponseRecord.hpp :: ResponseRecord2::record_excited_protocol_result`

### Validation / Planning References

- `src/apps/madqc_v2/test_molresponse_excited_metadata_smoke.py`
- `src/apps/madqc_v2/test_molresponse_excited_restart_reuse.py`
- `src/apps/madqc_v2/test_molresponse_excited_protocol_projection.py`
- `src/apps/madqc_v2/test_molresponse_excited_guess_fallback.py`
- `src/apps/madqc_v2/test_molresponse_excited_restart_mode_compat.py`
- `src/apps/molresponse_v2/EXCITED_STATE_EXECUTION_CHECKLIST.md`
- `src/apps/molresponse_v2/EXCITED_STATE_REINTEGRATION_PLAN.md`
