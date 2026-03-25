# Current `molresponse` Inventory

## Scope

- Primary workflow entry: `src/madness/chem/MolresponseLib.hpp`
- Workflow dispatch also touches `src/madness/chem/WorkflowBuilders.hpp` and `src/madness/chem/ResponseParameters.hpp`.
- Production `molresponse_v2` library sources per `src/apps/molresponse_v2/CMakeLists.txt`:
  - `ResponseManager.cpp`, `GroundStateData.cpp`, `ResponseSolver.cpp`, `ResponseState.cpp`, `ResponseInitializer.cpp`, `FrequencyLoop.cpp`, `ExcitedStateBundleSolver.cpp`, `ExcitedResponse.cpp`
- Production `molresponse_v2` headers per `CMakeLists.txt`:
  - `GroundStateData.hpp`, `MolecularProperty.hpp`, `ResponseBundle.hpp`, `ResponseManager.hpp`, `ResponseState.hpp`, `ResponseIO.hpp`, `ResponseRecord.hpp`, `DerivedStatePlanner.hpp`, `StateParallelPlanner.hpp`, `FrequencyLoop.hpp`, `PropertyManager.hpp`, `ResponseSolver.hpp`, `ResponseSolverUtils.hpp`, `ResponseKernels.hpp`, `ResponseVectorKernels.hpp`, `excited_ops/*.hpp`, `ResponseInitializer.hpp`, `ResponseDebugLogger.hpp`, `ResponseDebugLoggerMacros.hpp`, `broadcast_json.hpp`
- Executable in scope: `src/apps/molresponse_v2/molresponse2.cpp`
- Also inventoried because they define current runtime behavior even though they live outside `src/apps/molresponse_v2/`:
  - `src/madness/chem/MolresponseLib.hpp`
  - `src/madness/chem/WorkflowBuilders.hpp`
  - `src/madness/chem/ResponseParameters.hpp`
- Explicitly excluded from the main inventory tables:
  - tests: `test_parameter_manager.cpp`, `test_tda_h2.cpp`, `test_response_vector_alias.cpp`, `test_preliminaries.cpp`
  - docs/notes/checklists/markdown files
  - validation helper script `validate_property_component_claims.py`
- This tree is the active response implementation. Unlike the legacy tree, solver orchestration, restart metadata, property staging, and subgroup scheduling are all first-class current features.

## Top-Level Algorithmic Flow

### Main workflow path

1. Workflow selection routes the response task into `molresponse_lib::run_response`.
   - Implemented by `add_response_workflow_drivers` in `src/madness/chem/WorkflowBuilders.hpp:75` and `molresponse_lib::run_response` in `src/madness/chem/MolresponseLib.hpp:4874`.
2. Serialize the resolved response input block, rebuild `ResponseParameters`, and construct the ground/runtime context from the SCF checkpoint.
   - Implemented by `run_response` in `MolresponseLib.hpp:4874`, `ResponseParameters::ResponseParameters(World&, parser)` in `ResponseParameters.hpp:8`, and `make_ground_context` in `MolresponseLib.hpp:1146`.
3. Generate all required linear states, derived-state requests, excited-state bundle plan, and the state-parallel ownership plan.
   - Implemented by `plan_required_states` in `MolresponseLib.hpp:1223`, `StateGenerator::generateStates` in `StateGenerator.hpp`, `DerivedStatePlanner::build_vbc_driven_quadratic_plan` in `DerivedStatePlanner.hpp`, and `StateParallelPlanner::build` in `StateParallelPlanner.hpp`.
4. Build a restart-aware runtime schedule that decides, per protocol threshold, whether to run channel-series ownership or independent channel-frequency point ownership.
   - Implemented by `build_state_solve_schedule_context` in `MolresponseLib.hpp:2080`, with helpers `compute_runtime_point_ownership_policy` at `:1448` and `build_protocol_execution_policy` at `:1341`.
5. For each active protocol threshold, set MADNESS numerical defaults, reprepare the ground orbitals/Fock context, and solve every pending linear response point, either serially or in subgroup subworlds.
   - Serial path: `execute_serial_state_solve` in `MolresponseLib.hpp:2269`
   - Subgroup path: `execute_subgroup_state_solve` in `MolresponseLib.hpp:2414`
   - Per-state/channel-series solve: `computeFrequencyLoop` in `FrequencyLoop.cpp:566`
   - Per-point solve: `computeFrequencyPoint` in `FrequencyLoop.cpp:620`
6. During each point solve, choose a restart seed, possibly promote static archives to dynamic layout, iterate the coupled response equations to convergence, and record metadata/debug information.
   - Implemented by `seed_response_guess_for_point` in `FrequencyLoop.cpp:187`, `promote_response_vector` in `FrequencyLoop.cpp:509`, `iterate` in `FrequencyLoop.hpp:57`, `solve_and_record_point` in `FrequencyLoop.cpp:373`, and `JsonStateSolvePersistence` in `MolresponseLib.hpp:230`.
7. Run the protocol-indexed excited-state bundle stage, reusing restart-ready bundle snapshots when possible and recording per-protocol excited metadata.
   - Implemented by `execute_excited_state_bundle_stage` in `MolresponseLib.hpp:3029`, `ExcitedResponse::solve_protocol` in `ExcitedResponse.cpp`, and `ResponseRecord2` excited-state helpers in `ResponseRecord.hpp`.
8. Evaluate the derived-state dependency gate at the final protocol and execute ready VBC requests, either in subgroup lanes or a deterministic serial fallback.
   - Implemented by `execute_derived_state_requests` in `MolresponseLib.hpp:3531`, `DerivedStatePlanner::evaluate_dependency_gate` in `DerivedStatePlanner.hpp`, and `run_derived_request` in `MolresponseLib.hpp:2948`.
9. Rebuild the final protocol context and compute the requested properties, optionally precomputing components in subgroups and assembling on one designated property subgroup.
   - Implemented by `prepare_and_validate_final_protocol_state` in `MolresponseLib.hpp:3854`, `compute_requested_properties` in `MolresponseLib.hpp:4475`, and `compute_requested_properties_with_property_group` in `MolresponseLib.hpp:4510`.
10. Return structured JSON fragments for metadata, properties, Raman outputs, and debug logs.
   - Implemented by `run_response` in `MolresponseLib.hpp:4874`.

### Linear-response solve flow

1. Build the perturbation descriptor and file key for one channel/frequency/protocol target.
   - `LinearResponseDescriptor` / `LinearResponsePoint` in `ResponseState.hpp` and `ResponseState.cpp`
2. Prepare the protocol context by setting wavelet order/threshold, reprojection rules, Coulomb operator, and ground-state preliminaries.
   - `ResponseManager::setProtocol` in `ResponseManager.cpp:14`, `GroundStateData::prepareOrbitals` in `GroundStateData.cpp:100`, `GroundStateData::computePreliminaries` in `GroundStateData.cpp:166`, and `prepare_protocol_context` in `MolresponseLib.hpp:2197`
3. Decide whether the point still needs work from persistence metadata.
   - `point_needs_solving` in `MolresponseLib.hpp:1692`, `point_needs_solving_from_metadata` in `MolresponseLib.hpp:1706`, and `should_solve_point` in `FrequencyLoop.cpp:33`
4. Seed the guess from the current protocol archive, a coarser protocol archive, previous-frequency carryover, previous-frequency archive, or a fresh initializer.
   - `seed_response_guess_for_point` in `FrequencyLoop.cpp:187`, `seed_retry_guess_from_previous_frequency_archive` in `FrequencyLoop.cpp:264`, and `initialize_guess_vector` in `ResponseInitializer.cpp:3`
5. Run the KAIN/fixed-point iteration for the chosen response-vector variant.
   - `solve_response_vector` in `FrequencyLoop.cpp:485`, `iterate` in `FrequencyLoop.hpp:57`, and the type-dispatched operator kernels in `ops/*.hpp` and `ResponseKernels.hpp`
6. Record convergence/failure diagnostics, timing, restart provenance, and saved/converged flags.
   - `solve_and_record_point` in `FrequencyLoop.cpp:373`, `JsonStateSolvePersistence` in `MolresponseLib.hpp:230`, and `ResponseRecord2` in `ResponseRecord.hpp`

### Property flow

1. Filter property work to frequencies/states that are actually converged at the final protocol.
   - `compute_polarizability` in `MolresponseLib.hpp:4167`
2. Load final response vectors from archive and contract them with perturbation vectors or derived VBC states.
   - `load_response_vector` in `ResponseIO.hpp`, `compute_alpha` / `compute_hyperpolarizability` / `compute_Raman_components` in `PropertyManager.hpp`
3. Save flat property rows to `properties.json`, then build Raman Hessian/post-processing tables if requested.
   - `PropertyManager::save` in `PropertyManager.hpp`, `compute_hessian` / `compute_Raman` in `PropertyManager.hpp`, and `print_raman_table` in `MolresponseLib.hpp:4271`

## Missing or Partial Features

- Saved-state naming convention:
  - Present and much more structured than legacy.
  - Linear-response keys are centralized in `LinearResponseDescriptor::make_key` (`ResponseState.cpp:51`) and `LinearResponsePoint::response_filename` (`ResponseState.cpp:96`).
  - Metadata keys are also centralized in `ResponseRecord2::freq_key`, `ResponseRecord2::protocol_key`, and excited-state naming helpers in `ResponseRecord.hpp`.
- Restart logic:
  - Present and layered.
  - Linear-response restart precedence is implemented in `seed_response_guess_for_point` (`FrequencyLoop.cpp:187`) and metadata gating in `point_needs_solving*` (`MolresponseLib.hpp:1692`, `:1706`).
  - Excited-state restart metadata and normalization live in `ResponseRecord.hpp`; protocol-stage execution uses them in `execute_excited_state_bundle_stage` (`MolresponseLib.hpp:3029`).
- Protocol handling / threshold ramping:
  - Present and explicit.
  - Protocol thresholds come from `CalculationParameters::protocol()` and are consumed in planning, stage-2 loops, response archive naming, and excited metadata.
  - Numerical setup is centralized in `ResponseManager::setProtocol` (`ResponseManager.cpp:14`) and `prepare_protocol_context` (`MolresponseLib.hpp:2197`).
- Parallelization logic:
  - Present and substantial.
  - `StateParallelPlanner` computes static ownership; `MolresponseLib.hpp` computes restart-aware runtime ownership and uses `MacroTaskQ::create_worlds(...)` for subgroup subworlds.
  - Property assembly also has dedicated property-group execution and component precompute.
- Partial / not-yet-complete areas:
  - Unrestricted response kernels are stubbed or placeholder in `ops/StaticUnrestrictedOps.hpp` and `ops/DynamicUnrestrictedOps.hpp`.
  - The legacy-bundle adapter pathway is optional at build time; current excited-state bundle execution is still under active reintegration and uses metadata scaffolding plus `ExcitedResponse`.
  - There is still some duplication of naming helpers between `ResponseState.cpp`, `ResponseRecord.hpp`, and `ExcitedStateBundleSolver.cpp` for excited-root/state-name formatting.

## Mixed-Concern Findings

### SOLVER functions flagged `MIXED`

- `solve_all_states` (`MolresponseLib.hpp:4063`): `SOLVER (MIXED: protocol handling, restart logic, state-parallel logic)`
- `execute_serial_state_solve` (`MolresponseLib.hpp:2269`): `SOLVER (MIXED: restart logic, protocol handling, state-parallel logic, IO via persistence)`
- `execute_subgroup_state_solve` (`MolresponseLib.hpp:2414`): `SOLVER (MIXED: restart logic, protocol handling, state-parallel logic, IO/log sharding)`
- `execute_excited_state_bundle_stage` (`MolresponseLib.hpp:3029`): `SOLVER (MIXED: restart logic, naming logic, protocol handling, model-specific branching)`
- `execute_derived_state_requests` (`MolresponseLib.hpp:3531`): `SOLVER (MIXED: restart logic, IO claim files, state-parallel logic)`
- `computeFrequencyLoop` (`FrequencyLoop.cpp:566`): `SOLVER (MIXED: restart logic, protocol progression, carryover logic)`
- `computeFrequencyPoint` (`FrequencyLoop.cpp:620`): `SOLVER (MIXED: restart logic, protocol-specific point solve)`
- `solve_and_record_point` (`FrequencyLoop.cpp:373`): `SOLVER (MIXED: restart logic, IO persistence, failure policy)`
- `seed_response_guess_for_point` (`FrequencyLoop.cpp:187`): `SOLVER (MIXED: restart logic, protocol handling, model-specific static/dynamic promotion)`
- `solve_response_vector` (`FrequencyLoop.cpp:485`): `SOLVER (MIXED: model-specific branching)`
- `iterate` (`FrequencyLoop.hpp:57`): `SOLVER (MIXED: model-specific type dispatch via templates, convergence policy)`
- `ExcitedResponse::solve_protocol` (`ExcitedResponse.cpp`): `SOLVER (MIXED: restart logic, protocol handling, model-specific TDA/full branching)`
- `RestartAwareExcitedScaffoldSolver::*` in `ExcitedStateBundleSolver.cpp`: `SOLVER (MIXED: naming logic, restart logic, protocol handling)`

### PROPERTY functions flagged `MIXED`

- `compute_polarizability` (`MolresponseLib.hpp:4167`): `PROPERTY (MIXED: IO archive loads, metadata/restart readiness filtering)`
- `compute_hyperpolarizability` (`MolresponseLib.hpp:4255`): `PROPERTY (MIXED: model-specific frequency-triplet policy)`
- `compute_raman` (`MolresponseLib.hpp:4346`): `PROPERTY (MIXED: IO/property saves, model-specific vibrational workflow)`
- `compute_requested_properties_with_property_group` (`MolresponseLib.hpp:4510`): `PROPERTY (MIXED: IO shard merge, state-parallel logic)`
- `PropertyManager::save` (`PropertyManager.hpp:282`): `PROPERTY (MIXED: IO)`
- `compute_alpha` (`PropertyManager.hpp:450`): `PROPERTY (MIXED: IO archive loads, naming/state lookup)`
- `compute_beta` overloads (`PropertyManager.hpp:789`, `1021`): `PROPERTY (MIXED: naming logic, IO archive loads, task-claim files)`
- `compute_Raman_components` / `compute_Raman` (`PropertyManager.hpp:1109`, `1147`): `PROPERTY (MIXED: naming logic, IO saves, model-specific Raman assembly)`

## Naming Convention

- Linear response states are named from the perturbation descriptor plus encoded frequency and protocol threshold.
- The canonical linear filename/key format is:
  - `<PerturbationDescription>_f<fixed-3-decimal-frequency>_p<scientific-threshold>`
  - Implemented by `LinearResponseDescriptor::make_key` in `src/apps/molresponse_v2/ResponseState.cpp:51`
- `LinearResponsePoint::response_filename` returns that key and is what `ResponseIO.hpp` uses for archives.
  - Implemented at `ResponseState.cpp:96`
- Frequency normalization is centralized by `canonicalize_response_frequency` in `ResponseState.hpp:23`, which zeroes tiny values and rounds to a stable 3-decimal representation before naming/planning.
- Second-order / derived response states use a distinct combined key:
  - `<prefix><pertB>_<pertC>_f<fB>_<fC>_p<thresh>`
  - Implemented by `SecondOrderResponseDescriptor::make_key` in `ResponseState.cpp:148`
- Metadata naming is separately centralized in `ResponseRecord2`:
  - `freq_key(double)` for frequency JSON keys
  - `protocol_key(double)` for protocol JSON keys
  - `make_excited_root_id(size_t)` for stable excited-root ids like `es_root_0000`
  - `assign_excited_state_names(...)` for display names like `es1`, `es2a`, `es2b`
- Excited-state naming is mostly centralized in `ResponseRecord.hpp`, but `ExcitedStateBundleSolver.cpp` also carries local helper duplicates (`protocol_key`, `alpha_suffix`, `assign_excited_state_names`, `make_excited_root_id`) for solver-side restart/scaffold work, so this part is not fully singular yet.

## Restart Logic

- Linear-response restart detection uses two layers:
  - archive existence/loading in `ResponseIO.hpp`
  - metadata status checks through `StateSolvePersistence` / `ResponseRecord2`
- The solver decides whether to reload vs recompute with:
  - `point_needs_solving` in `MolresponseLib.hpp:1692`
  - `point_needs_solving_from_metadata` in `MolresponseLib.hpp:1706`
  - `should_solve_point` in `FrequencyLoop.cpp:33`
- Linear-response seed precedence in `seed_response_guess_for_point` (`FrequencyLoop.cpp:187`) is:
  1. exact current protocol archive
  2. nearest lower protocol archive
  3. previous-frequency in-memory carryover
  4. previous-frequency archive fallback
  5. fresh initializer guess
- If the target is dynamic but only a static archive exists, `promote_response_vector` (`FrequencyLoop.cpp:509`) expands the static vector into the dynamic layout instead of forcing a fresh solve.
- Final-protocol convergence matters: a point can be saved but still require recomputation at the last threshold if it is not marked converged.
- Removed-frequency policy is explicitly tracked through solver diagnostics and can be overridden with `force_retry_removed_frequencies`.
- Excited-state restart logic uses protocol-indexed metadata plus restart-capability flags in `ExcitedProtocolResult`, normalized through `ResponseRecord2::normalize_excited_protocol_result(...)` and consumed by `execute_excited_state_bundle_stage`.

## Protocol Handling

- Protocol thresholds come from `CalculationParameters::protocol()` and drive stage planning, runtime execution, metadata indexing, and archive keys.
- Numerical protocol setup is centralized in `ResponseManager::setProtocol` (`ResponseManager.cpp:14`), which sets:
  - MADNESS threshold/order defaults
  - `vtol`
  - Coulomb operator
  - derivative operators
- Ground-state reprojection and protocol-consistent orbital preparation happen in `GroundStateData::prepareOrbitals` (`GroundStateData.cpp:100`).
- Ground-state preliminary operators/Fock data are rebuilt per protocol in `GroundStateData::computePreliminaries` (`GroundStateData.cpp:166`).
- `prepare_protocol_context` (`MolresponseLib.hpp:2197`) is the stage-2 orchestration wrapper that applies the current protocol consistently before any linear solve work.
- Protocol ramping is implemented as explicit loops over the ordered threshold schedule:
  - generic loop helper: `run_protocol_threshold_loop` (`MolresponseLib.hpp:1749`)
  - linear solve paths: `execute_serial_state_solve` / `execute_subgroup_state_solve`
  - excited-state bundle path: `execute_excited_state_bundle_stage`
- The runtime can skip fully complete protocol prefixes on restart by setting `restart_start_protocol_index` from existing metadata.

## State-Parallel Logic

- Static planning lives in `StateParallelPlanner.hpp`.
  - `StateParallelPlanner::build` decides whether state-parallel mode is off/auto/on, how many mapping groups to use, how many groups can participate in point mode, and which protocol index may start point-level fanout.
  - `PointOwnershipScheduler` provides deterministic `(channel, frequency) -> owner group` mapping.
- Restart-aware runtime scheduling lives in `MolresponseLib.hpp`.
  - `compute_runtime_point_ownership_policy` (`:1448`) reads existing metadata and can advance the effective point-parallel start protocol to the earliest pending protocol.
  - `build_protocol_execution_policy` (`:1341`) computes per-protocol active groups and whether each threshold should use channel-series or channel-point scheduling.
  - `build_pending_work_manifest` (`:1800`) turns those policies into actual per-lane worklists.
- Subworld creation/management is done with `MacroTaskQ::create_worlds(...)` in `execute_subgroup_state_solve` and `compute_requested_properties_with_property_group`.
- Metadata and debug logs are sharded per subgroup with helpers:
  - `group_shard_file`
  - `group_console_file`
  - `merge_state_metadata_json`
  - `merge_debug_log_json`
- Property work also has a state-parallel mode:
  - component precompute can happen across all execution groups
  - final property assembly runs on one chosen `state_parallel_property_group`
  - outputs are broadcast back to the whole world

## Structure and Code-Quality Notes

- The current code is much more staged than legacy: Stage 1 planning, Stage 2 solves, and Stage 3 properties are explicit in `MolresponseLib.hpp`.
- `ResponseRecord2`, `JsonStateSolvePersistence`, and `ResponseDebugLogger` are clean abstractions that keep most restart/status bookkeeping out of the low-level iteration kernel.
- `StateParallelPlanner` is notably readable: it separates static planning from restart-aware runtime policy, which makes subgroup behavior auditable.
- `ResponseState`, `Perturbation`, and `StateGenerator` form a good naming/descriptor layer; they make the top-level orchestration substantially easier to follow than the legacy filename-and-mode plumbing.
- The main structural weakness is that `MolresponseLib.hpp` is now the new god-file: it is much cleaner than legacy `TDDFT`, but it still mixes planning, scheduling, restart policy, JSON merge helpers, subgroup logging, excited-state staging, derived-state execution, and property orchestration in one very large header.
- A second weak spot is duplicated naming/restart helper logic for excited states between `ResponseRecord.hpp` and `ExcitedStateBundleSolver.cpp`.

## File-by-File Inventory

### `src/madness/chem/WorkflowBuilders.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `workflow_kind_from_name` | function | 31 | UTILITY | Maps user workflow text to the internal workflow enum. | `std::string_view user_workflow` | `WorkflowKind` | string comparisons | `add_workflow_drivers` |
| `runnable_workflow_list` | function | 49 | UTILITY | Returns the printable list of supported workflow names. | none | `const std::string&` | none | workflow error/reporting |
| `add_scf_workflow_drivers`, `add_nemo_workflow_drivers`, `add_cc2_workflow_drivers`, `add_cis_workflow_drivers`, `add_oep_workflow_drivers` | functions | 62, 68, 86, 109, 118 | UNCLEAR | Register non-response workflow handlers into the shared driver table. | `World&`, `Params&`, workflow container refs | `void` | workflow-app registration helpers | `add_workflow_drivers` |
| `add_response_workflow_drivers` | function | 75 | IO | Registers the response workflow so `madqc --wf=response` dispatches SCF first and then `molresponse_lib`. | `World&`, `Params&`, workflow container refs | `void` | `ResponseApplication<molresponse_lib>` | `add_workflow_drivers` |
| `add_workflow_drivers` | function | 136 | IO | Top-level workflow-driver registration entrypoint. | `World&`, `Params&`, workflow container refs | `void` | `workflow_kind_from_name`, all `add_*_workflow_drivers` helpers | MADNESS workflow bootstrap |

### `src/madness/chem/ResponseParameters.hpp`

This file is the current response input contract and mostly consists of field registration plus typed getters.

| Symbol(s) | Kind | Line(s) | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `ResponseParameters` | struct | 6 | IO | Owns all current response input knobs, including state-parallel, restart policy, excited-stage, and property options. | none | object | `initialize<T>` for many fields | current response code |
| constructors | methods | 8, 17 | IO | Read input/CLI data, derive default state-parallel groups, derive property list, and validate user-specified property requests. | `World&`, `commandlineparser` or none | object | `read_input_and_commandline_options`, `set_derived_properties`, `validate_user_specified_properties` | `molresponse_lib::run_response`, tests |
| typed getters (`prefix`, `fock_json_file`, `archive`, `print_level`, `state_parallel*`, `excited_*`, `dipole_*`, `nuclear_*`, etc.) | methods | 80-217 | IO | Expose typed access to parsed response settings. | none | values of configured types | `get<T>` | nearly every planning/solve/property function |
| `validate_user_specified_properties` | method | 220 | IO | Enforces that explicitly requested properties have the supporting dipole/nuclear knobs defined. | none | `void` | `requested_properties`, `is_user_defined` | constructor |
| `set_derived_properties` | method | 241 | IO | Backfills `requested_properties` from older `dipole`/`nuclear`/`quadratic` flags when the newer explicit property list was not used. | none | `void` | `set_derived_value`, property booleans | constructor |

### `src/madness/chem/MolresponseLib.hpp`

This is the main orchestration file. The inventory below groups related helpers where they form one cohesive stage.

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `molresponse_lib` | struct | 137 | UNCLEAR | Owns the end-to-end response workflow orchestration and result packaging. | none | type | stage helpers below | workflow driver |
| `Results`, `GroundContext`, `ExcitedStateBundlePlan`, `PlannedStates`, `SolvedStates`, `PropertyStageOutput`, `ProtocolExecutionPolicy`, `RuntimePointOwnershipPolicy`, `StateSolveScheduleContext`, `PendingPointWorkItem`, `PendingProtocolManifest`, `DerivedRequestTiming`, `DerivedExecutionResult`, `ExcitedExecutionResult`, `FinalProtocolState`, `PropertyContext` | structs | 147-223, 1260-1798, 2939-3027, 3846-4165 | UTILITY | Value objects that carry state between planning, solve, metadata, and property stages. | field data | objects/JSON | trivial methods like `to_json` | stage orchestration |
| `JsonStateSolvePersistence` | class | 230 | IO | Implements the solver persistence interface by writing status, timing, diagnostics, restart provenance, and debug logs into stable JSON files. | metadata/debug filenames, baseline metadata, flags | object | `ResponseRecord2`, `ResponseDebugLogger`, `point_ready_in_metadata` | serial/subgroup stage-2 solve paths |
| JSON/file/log helpers: `group_shard_file`, `group_console_file`, `group_derived_timing_file`, `sanitize_manifest_token`, `derived_request_manifest_scope`, `derived_request_done_file`, `derived_request_claim_file`, `derived_request_done_record_exists`, `write_derived_request_done_record`, `FilteredLineStreambuf`, `ScopedRankLogRedirect`, `write_json_file`, `read_json_file_or_object`, `broadcast_json_object`, `broadcast_json`, `with_subworld`, `merge_state_metadata_json`, `merge_debug_log_json`, `point_ready_in_metadata`, `point_marked_for_frequency_removal_in_metadata` | functions/classes | 362-1130 | IO | Provide shard naming, JSON broadcast/merge, manifest claim-file helpers, and subgroup console redirection used by the solver and property stages. | filenames, JSON payloads, subgroup counts, points | filenames, JSON, booleans, side effects | filesystem I/O, MADNESS broadcast/fence | many orchestration helpers |
| `make_ground_context` | function | 1146 | IO | Builds the reusable ground-state/runtime context from SCF checkpoint files and resolved archive paths. | `World&`, `CalculationParameters`, `shared_ptr<SCF>`, `outdir` | `GroundContext` | checkpoint JSON read, `GroundStateData`, `ResponseManager` | `run_response` |
| `build_excited_state_bundle_plan` | function | 1201 | SOLVER | Translates response knobs and subgroup limits into the protocol-indexed excited-bundle execution plan. | `CalculationParameters`, `ResponseParameters`, `StateParallelPlan` | `ExcitedStateBundlePlan` | response getters | `plan_required_states` |
| `plan_required_states` | function | 1223 | SOLVER | Performs stage-1 planning for generated linear states, derived requests, excited-bundle metadata, and static subgroup ownership. | `World&`, `CalculationParameters`, `GroundContext`, `ResponseParameters` | `PlannedStates` | `StateGenerator`, `DerivedStatePlanner`, `StateParallelPlanner`, `build_excited_state_bundle_plan` | `run_response` |
| Scheduling helpers: `build_owner_by_channel_index`, `build_protocol_execution_policy`, `compute_runtime_point_ownership_policy`, `build_local_channel_workset`, `use_channel_series_ownership_for_protocol_runtime`, `active_owner_groups_for_protocol_runtime`, `point_needs_solving`, `point_needs_solving_from_metadata`, `any_state_point_needs_solving`, `run_protocol_threshold_loop`, `run_frequency_loop_with_flush`, `build_pending_work_manifest`, `build_pending_manifest_from_metadata`, `print_state_solve_execution_mode`, `build_state_solve_schedule_context` | functions | 1328-2175 | SOLVER (MIXED: restart logic, protocol handling, state-parallel logic) | Convert the static state-parallel plan plus restart metadata into per-protocol runtime execution policy and concrete pending-work manifests. | schedule context, metadata, state lists, protocol index | policies, manifests, booleans, side effects | `read_json_file_or_object`, `point_ready_in_metadata`, `build_pending_work_manifest`, `compute_runtime_point_ownership_policy` | `solve_all_states`, subgroup progress polling |
| Linear-stage helpers: `prepare_protocol_context`, `log_pending_manifest`, `execute_manifest_work`, `cached_or_built_manifest`, `execute_serial_state_solve`, `execute_subgroup_state_solve` | functions | 2197-2922 | SOLVER (MIXED: restart logic, protocol handling, state-parallel logic, IO) | Apply the active protocol, dispatch stage-2 linear work either on the whole world or in subgroup shards, and emit merged metadata/debug output. | worlds, contexts, schedule info, metadata/debug refs | `void` / `bool` | `ResponseManager::setProtocol`, `GroundStateData`, `computeFrequencyLoop`, `computeFrequencyPoint`, JSON merge helpers, `MacroTaskQ::create_worlds` | `solve_all_states` |
| `run_derived_request` | function | 2948 | SOLVER (MIXED: IO) | Runs one VBC derived-state request and returns timing/success information. | `World&`, `GroundStateData`, `DerivedStateRequest`, `SimpleVBCComputer`, counters | `DerivedRequestTiming` | `DerivedStatePlanner::make_vbc_state`, `SimpleVBCComputer::compute_and_save` | serial and subgroup derived execution |
| `ensure_excited_protocol_placeholder_node`, `build_excited_root_manifest`, `execute_excited_state_bundle_stage` | functions | 3016, 3024, 3029 | SOLVER (MIXED: restart logic, naming logic, protocol handling, model-specific branching) | Manage per-protocol excited-state placeholder metadata, dispatch the excited solver when needed, normalize results, and record them into metadata. | `World&`, plans, contexts, metadata | `ExcitedExecutionResult` / JSON | `ResponseRecord2`, `ExcitedResponse::solve_protocol`, JSON broadcast/normalization helpers | `solve_all_states` |
| `execute_derived_state_requests` | function | 3531 | SOLVER (MIXED: restart logic, IO claim files, state-parallel logic) | Evaluates the dependency gate, skips precompleted derived requests, and executes ready requests in subgroup or serial mode. | `World&`, contexts, plans, final metadata info | `DerivedExecutionResult` | `DerivedStatePlanner::evaluate_dependency_gate`, `run_derived_request`, `try_claim_property_component_task`, `MacroTaskQ::create_worlds` | `solve_all_states` |
| `prepare_and_validate_final_protocol_state`, `build_state_stage_metadata`, `k_from_thresh`, `preflight_memory_estimate` | functions | 3854, 3892, 3950, 3967 | UTILITY | Rebuild final protocol state, package runtime metadata, map threshold to `k`, and estimate memory before subgroup allocation. | world/context/protocol data | `FinalProtocolState`, JSON, `int`, side effects | `ResponseManager::setProtocol`, `GroundStateData`, metadata helpers | `solve_all_states` |
| `solve_all_states` | function | 4063 | SOLVER (MIXED: protocol handling, restart logic, state-parallel logic) | Master stage-2 orchestrator for linear states, excited bundle, derived states, and combined metadata/debug output. | `World&`, `CalculationParameters`, `GroundContext`, `ResponseParameters`, `PlannedStates` | `SolvedStates` | `build_state_solve_schedule_context`, stage-2 serial/subgroup solves, excited and derived stages, `build_state_stage_metadata` | `run_response` |
| `parse_property_name`, `compute_polarizability`, `compute_hyperpolarizability`, `print_raman_table`, `compute_raman`, `compute_requested_properties`, `compute_requested_properties_with_property_group` | functions | 4144, 4167, 4255, 4271, 4346, 4475, 4510 | PROPERTY / PROPERTY (MIXED where noted above) | Stage-3 property dispatch, filtering, Raman post-processing, subgroup property execution, and property-output packaging. | property context, solved metadata, subgroup settings | property JSON and vibrational/Raman outputs | property helpers from `PropertyManager.hpp`, metadata readiness helpers, `MacroTaskQ::create_worlds` | `run_response` |
| `run_response` | function | 4874 | UNCLEAR | Public workflow entrypoint that drives stages 1 through 3 and returns structured results. | `World&`, `Params`, `shared_ptr<SCF>`, `outdir` | `Results` | `ResponseParameters`, `make_ground_context`, `plan_required_states`, `solve_all_states`, `compute_requested_properties_with_property_group` | `ResponseApplication<molresponse_lib>` |

### `src/apps/molresponse_v2/molresponse2.cpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `print_help` | function | 12 | IO | Prints CLI usage for the standalone `molresponse2` executable. | none | `void` | printing only | `main` |
| `main` | function | 29 | UNCLEAR | Standalone executable wrapper that parses CLI input and launches the current response codepath. | `int argc`, `char** argv` | `int` | parser/bootstrap helpers, current response library | process entry point |

### `src/apps/molresponse_v2/GroundStateData.hpp` / `GroundStateData.cpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `GroundStateData` | class | `GroundStateData.hpp` | IO | Owns ground-state orbitals, energies, occupations, molecule data, and protocol-dependent preliminaries used by response stages. | archive path or `SCF` handle | object | checkpoint load and compute helpers | stage 1/2/3 |
| constructors | methods | `GroundStateData.cpp:5`, `:16` | IO | Build ground-state data either from a restart archive/molecule or from an SCF object. | `World&`, archive/molecule or `shared_ptr<SCF>` | object | `load` and member initialization | `make_ground_context`, subgroup property/solve setup |
| `load` | method | 29 | IO | Reads the archived orbitals and basic ground-state quantities from disk. | `World&` | `void` | archive input helpers | constructor, restart setup |
| `print_info` | method | 87 | IO | Prints a summary of the loaded ground-state information. | none | `void` | printing only | diagnostics |
| `prepareOrbitals` | method | 100 | OPERATOR | Reprojects/truncates the occupied orbitals to the active protocol `k` and threshold before response work begins. | `World&`, `int current_k`, `double thresh` | `void` | `project`, `truncate`, protocol helpers | `prepare_protocol_context`, property subgroup setup |
| `tryLoadHamiltonianFromJson` | method | 147 | IO | Attempts to reuse a serialized Fock/Hamiltonian tensor from JSON for the active protocol. | `World&`, `json`, `double thresh`, `int k` | `Tensor<double>` | JSON parsing | `computePreliminaries` |
| `computePreliminaries` | method | 166 | OPERATOR (MIXED: IO, model-specific branching) | Builds the density, nuclear/Coulomb/XC/HF-exchange intermediates, optionally loading a saved Hamiltonian instead of recomputing it. | `World&`, `operatorT`, `double vtol`, `std::string fock_json_file` | `void` | `computeDensity`, `computeNuclearPotential`, `computeXCPotential`, `computeHFExchangeEnergy`, `tryLoadHamiltonianFromJson` | `prepare_protocol_context`, property setup |
| `computeDensity`, `computeNuclearPotential`, `computeXCPotential`, `computeHFExchangeEnergy`, `computeKineticEnergy` | methods | 220, 229, 241, 247, 260 | OPERATOR | Compute the ground-state density and one-electron/two-electron operator pieces reused by response kernels. | `World&` or scalar args | density/potential/tensor values | MADNESS density/operator helpers | `computePreliminaries`, solver kernels |
| trivial getters (`getOrbitals`, `getEps`, `getOcc`, `getNumOrbitals`, `getL`, `getK`, `getMolecule`, `isSpinRestricted`, etc.) | methods | `GroundStateData.hpp` inline | IO | Expose stored ground-state data to the solver/property code. | none | member refs/copies | none | throughout current code |

### `src/apps/molresponse_v2/ResponseManager.hpp` / `ResponseManager.cpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `ResponseManager` | class | `ResponseManager.hpp:17` | UTILITY | Holds protocol-dependent numerical defaults and operator handles shared by the response solver. | `World&`, `CalculationParameters` | object | member initialization | `make_ground_context`, subgroup setup |
| constructor | method | `ResponseManager.cpp:10` | UTILITY | Initializes the manager from the global calculation parameters. | `World&`, `CalculationParameters` | object | none | stage setup |
| `setProtocol` | method | 14 | OPERATOR | Centralizes protocol-dependent MADNESS defaults, operator creation, and `vtol` selection for one threshold. | `World&`, `double L`, `double thresh`, `int override_k=-1` | `void` | `FunctionDefaults`, Coulomb/gradient operator factories | `prepare_protocol_context`, final property rebuild |
| getters (`getCoulombOp`, `getGradOps`, `getVtol`, `params`) | methods | `ResponseManager.hpp` inline | UTILITY | Expose the active response numerical context to kernels and property rebuild paths. | none | refs/values | none | solver/property setup |

### `src/apps/molresponse_v2/Perturbation.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `DipolePerturbation`, `NuclearDisplacementPerturbation`, `MagneticPerturbation` | structs | 10, 15, 20 | PERTURBATION | Strongly typed perturbation descriptors for linear and derived response states. | field values | objects | none | `ResponseState`, `StateGenerator`, property code |
| `describe_perturbation` | function | 27 | PERTURBATION | Converts a perturbation descriptor into the canonical state-name fragment such as `Dipole_x` or `Nuclear_0_x`. | `const Perturbation&` | `std::string` | variant visitation | naming, state generation, properties |
| `perturbation_type_string` | function | 39 | PERTURBATION | Returns the perturbation family name used in diagnostics and branching. | `const Perturbation&` | `std::string` | variant visitation | planning and logging |

### `src/apps/molresponse_v2/ResponseState.hpp` / `ResponseState.cpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `canonicalize_response_frequency` | function | `ResponseState.hpp:23` | UTILITY | Normalizes frequencies to stable naming/lookup values by zeroing tiny values and rounding. | `double raw_frequency` | `double` | rounding helpers | `StateGenerator`, naming |
| `AbstractResponseDescriptor` | struct | `ResponseState.hpp:35` | UTILITY | Common descriptor interface for response states. | descriptor fields | object | virtual helpers | descriptor hierarchy |
| `LinearResponseDescriptor`, `LinearResponsePoint`, `SecondOrderResponseDescriptor`, `VBCResponseState`, `XBCResponseState` | structs | `ResponseState.hpp:44`, `:81`, `:104`, `:150`, `:155` | PERTURBATION | Represent linear points and second-order/derived state requests, including their perturbation, frequency, protocol, and spin restriction metadata. | perturbation/protocol/frequency fields | objects | inline accessors and helpers | planning, IO, properties, derived stage |
| `LinearResponseDescriptor` ctor and accessors (`num_thresholds`, `num_frequencies`, `threshold`, `frequency`, `is_static`, `is_spin_restricted`, `make_key`, `response_filename`, `perturbationDescription`) | methods | `ResponseState.cpp:8-72` | PERTURBATION | Own the canonical per-channel state description and archive key generation for linear response. | perturbation, threshold/frequency indices | values/strings | `describe_perturbation`, formatting helpers | planning, IO, restart, properties |
| `LinearResponsePoint` accessors (`threshold`, `frequency`, `is_static`, `is_spin_restricted`, `response_filename`, `perturbationDescription`) | methods | `ResponseState.cpp:80-100` | PERTURBATION | Provide per-point views and filenames on top of a channel descriptor. | none | values/strings | `LinearResponseDescriptor` | solver, restart, IO |
| `SecondOrderResponseDescriptor` methods (`B_state`, `C_state`, `get_states`, `current_threshold`, `current_frequency`, `perturbationDescription`, `is_spin_restricted`, `is_static`, `make_key`, `response_filename`, `make_vector`) | methods | `ResponseState.cpp:108-181` | PERTURBATION | Build the combined naming and state decomposition for VBC/XBC-like second-order states. | second-order descriptor fields | descriptors, strings, vectors | `LinearResponseDescriptor`, formatting helpers | derived state execution, property code |
| `dir_index` | function | 184 | UTILITY | Maps Cartesian direction characters to tensor indices. | `char c` | `int` | none | perturbation operator builders |
| `make_perturbation_operator` overloads, `raw_perturbation_operator`, `project_perturbation_onto_orbitals`, `perturbation_vector` overloads | functions | 197-273 | PERTURBATION | Build raw perturbation operators and project them into orbital-space RHS vectors for dipole, nuclear, magnetic, linear, and second-order states. | `World&`, `GroundStateData`, descriptor/perturbation values | functions or orbital-space vectors | ground-state data, `describe_perturbation`, MADNESS projection helpers | initializer, properties, RHS construction |

### `src/apps/molresponse_v2/StateGenerator.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `GeneratedStateData` | struct | 15 | PERTURBATION | Holds the generated linear-state list and the string-to-descriptor lookup map. | state collections | object | `print_generated_state_map` | planning |
| `GeneratedStateData::print_generated_state_map` | method | ~20 | IO | Prints the planned state map for inspection. | state map | `void` | printing only | `plan_required_states` |
| `StateGenerator` | class | 35 | PERTURBATION | Expands user property requests into the concrete linear response channels/frequencies required by stage 2. | molecule, protocols, spin restriction, response params | object | none | `plan_required_states` |
| `StateGenerator::generateStates` | method | ~45 | PERTURBATION | Generates dipole, nuclear, and quadratic-driven linear states and the state lookup map, canonicalizing frequencies as it does so. | none | `GeneratedStateData` | `canonicalize_response_frequency`, `describe_perturbation`, response getters | `plan_required_states` |
| `PropertyComponentPlan` | struct | 167 | PERTURBATION | Small helper value type describing planned property components. | fields | object | none | planning internals |

### `src/apps/molresponse_v2/StateParallelPlanner.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `PerturbationChannelAssignment`, `StateParallelPlan` | structs | 34, 49 | UTILITY | Hold static subgroup ownership assignments and planning diagnostics for state-parallel execution. | fields | objects | `to_json` | planning/runtime metadata |
| `PointOwnershipScheduler` | class | 110 | UTILITY | Provides deterministic `(channel, frequency)` owner-group mapping in point mode. | linear-state list, owner-group count | object | assignment math | runtime scheduling |
| `StateParallelPlanner` | class | 165 | UTILITY | Builds the static state-parallel plan from user knobs, world size, and generated states. | response params, world size, states | object | owner-group assignment helpers | `plan_required_states` |
| `StateParallelPlanner::build` | method | ~170 | UTILITY | Decides whether subgroup mode is active, how many groups to use, when point-mode can start, and which channels each group owns initially. | `ResponseParameters`, `size_t world_size`, `vector<LinearResponseDescriptor>` | `StateParallelPlan` | `PointOwnershipScheduler`-style logic, policy heuristics | `plan_required_states` |

### `src/apps/molresponse_v2/ResponseIO.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `save_response_vector` overloads | functions | 14, 121 | IO | Save a response vector archive under the canonical state filename. | `World&`, descriptor/point, response vector | `void` | archive output, `response_filename` | solver and derived/property code |
| `load_response_vector` overloads | functions | 40, 104 | IO | Load and protocol-align a response vector archive if it exists. | `World&`, orbital count, descriptor/point indices, output vector | `bool` success | archive input, `align_response_vector_protocol`, `response_filename` | restart seeding, property assembly |

### `src/apps/molresponse_v2/ResponseRecord.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `ExcitedRootDescriptor`, `ExcitedProtocolInput`, `ExcitedProtocolResult` | structs | 29, 56, 71 | IO | Carry excited-state naming, protocol-input, and protocol-result metadata in JSON-serializable form. | fields | objects | `to_json`/`from_json` | excited-stage metadata |
| `to_json` / `from_json` overloads for excited metadata structs | functions | 37, 45, 105, 140 | IO | Serialize/deserialize excited-state metadata records. | metadata structs / JSON | JSON / structs | JSON utilities | `ResponseRecord2`, excited stage |
| `ResponseRecord2` | class | 187 | IO | Central metadata authority for linear and excited response state status, timings, diagnostics, restart provenance, and naming helpers. | world, metadata filename | object | JSON file I/O | persistence layer and excited stage |
| key helpers: `freq_key`, `protocol_key`, `protocol_numeric`, `make_excited_root_id`, `excited_alpha_suffix`, `assign_excited_state_names` | methods | ~863, ~868, ~1068 | IO | Generate stable metadata keys and stable excited display/root names. | frequencies, thresholds, stable indices, energies | strings / vectors | formatting helpers | naming, restart, excited metadata |
| metadata builders: `initialize_states`, `initialize_excited_bundle`, `ensure_root`, `ensure_excited_root`, `ensure_excited_protocol`, `ensure_state`, `ensure_protocol` | methods | throughout class | IO | Create the canonical metadata skeleton before solver stages write status into it. | descriptors/plans | `void` / JSON refs | key helpers | `JsonStateSolvePersistence`, excited stage |
| state/diagnostic writers: `record_status`, `record_timing`, `record_solver_diagnostics`, `record_restart_provenance`, `record_excited_protocol_result`, `normalize_excited_protocol_result`, `is_saved`, `is_converged`, `is_marked_for_frequency_removal`, `to_json`, `print_summary` | methods | throughout class | IO | Read/write the canonical persisted status and normalize excited-stage results before storage. | points/protocols/diagnostics | booleans, JSON, side effects | JSON tree updates | `JsonStateSolvePersistence`, excited stage, restart gating |

### `src/apps/molresponse_v2/FrequencyLoop.hpp` / `FrequencyLoop.cpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `ResponseSolveDiagnostics` | struct | `FrequencyLoop.hpp:23` | SOLVER | Holds convergence/failure diagnostics for one point solve. | fields | object | none | persistence/logging |
| `response_point_stall_timeout_seconds` | function | `FrequencyLoop.hpp:36` | UTILITY | Reads the optional environment-based stall timeout for progress polling. | none | `double` | `getenv` | `iterate` |
| `iterate<R>` | function template | 57 | SOLVER (MIXED: model-specific branching) | Performs the actual KAIN/fixed-point iteration for one typed response-vector variant. | world, response manager, point, ground state, params, initial guess | `ResponseVector` plus diagnostics | response kernels, BSH helpers, debug logger | `solve_response_vector` |
| `solve_response_vector`, `promote_response_vector`, `computeFrequencyLoop`, `computeFrequencyPoint` declarations | functions | 299, 308, 358, 368 | SOLVER | Public frequency-solver entrypoints used by stage-2 orchestration. | see `.cpp` definitions | response or side effects | `iterate`, restart helpers | `MolresponseLib` |
| `StateSolvePersistence` | abstract class | 315 | IO | Abstracts saved/converged/restart-provenance/debug-log persistence away from the iteration kernel. | interface methods | polymorphic API | none | `JsonStateSolvePersistence`, `FrequencyLoop.cpp` |
| local structs `RestartGuessProvenance`, `GuessSeedResult`, `PointSolveResult`, `SolverAttemptResult` | structs | `FrequencyLoop.cpp:8-27` | SOLVER | Internal carriers for seed-selection and solve-attempt results. | fields | objects | none | point-solve helpers |
| `should_solve_point`, `apply_failure_policy`, `print_point_solve_report`, `print_solve_start` | functions | 33, 53, 104, 163 | SOLVER | Control whether a point needs work, apply failure/removal policy, and emit point-level diagnostics. | persistence, points, diagnostics | bool / `void` | persistence API | point solve drivers |
| `seed_response_guess_for_point`, `seed_retry_guess_from_previous_frequency_archive`, `should_retry_from_previous_frequency_seed` | functions | 187, 264, 299 | SOLVER (MIXED: restart logic, protocol handling, model-specific promotion) | Implement the layered restart/carryover policy for one linear-response point. | world, manager, point, ground state, persistence | seed structs / bool | `load_response_vector`, `initialize_guess_vector`, `promote_response_vector` | `solve_and_record_point`, `computeFrequencyLoop`, `computeFrequencyPoint` |
| `run_solver_attempt`, `solve_and_record_point` | functions | 312, 373 | SOLVER (MIXED: restart logic, IO persistence) | Run one solver attempt, optionally retry, and write final status/timing/provenance. | world, manager, point, ground state, persistence, seed | attempt/point result structs | `solve_response_vector`, persistence writers, `save_response_vector` | `computeFrequencyLoop`, `computeFrequencyPoint` |
| `solve_response_vector` | function | 485 | SOLVER (MIXED: model-specific branching) | Dispatches the iteration kernel to the correct response-vector variant. | world, manager, point, ground state, initial vector | `ResponseVector` | `iterate<StaticRestrictedResponse>`, `iterate<DynamicRestrictedResponse>` | `run_solver_attempt` |
| `promote_response_vector` | function | 509 | UTILITY | Converts a static restricted archive layout into a dynamic layout when restart reuse needs it. | world, input vector, output vector | `void` | flat-layout helpers | restart seeding |
| `computeFrequencyLoop` | function | 566 | SOLVER (MIXED: restart logic, protocol handling) | Solves one channel across all its frequencies at one protocol threshold, carrying guesses forward along the frequency series. | world, manager, state descriptor, threshold index, ground state, persistence, final-protocol flag | `void` | `solve_and_record_point` | serial/subgroup stage-2 solve |
| `computeFrequencyPoint` | function | 620 | SOLVER (MIXED: restart logic, protocol handling) | Solves one independent channel-frequency point at one protocol threshold. | same plus `freq_index` | `void` | `solve_and_record_point` | point-mode stage-2 solve |

### `src/apps/molresponse_v2/ResponseInitializer.hpp` / `ResponseInitializer.cpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `initialize_guess_vector` | function | `ResponseInitializer.cpp:3` | PERTURBATION | Builds the fresh initial response guess from the perturbation RHS and active response-vector type. | `World&`, `GroundStateData`, `LinearResponsePoint` | `ResponseVector` | `perturbation_vector`, `make_response_vector` helpers | `seed_response_guess_for_point` |

### `src/apps/molresponse_v2/ResponseVector.hpp`

This header is inline container and variant-dispatch logic for the current linear-response vector representation.

| Symbol(s) | Kind | Line(s) | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| response-vector variant structs (`StaticRestrictedResponse`, `DynamicRestrictedResponse`, `StaticUnrestrictedResponse`, `DynamicUnrestrictedResponse`) | structs | header top | UTILITY | Represent the typed channel layouts for current response solves. | orbital count / function vectors | objects | flatten/sync helpers | solver, IO, properties |
| `ResponseVector` | variant alias | header | UTILITY | Type-erased wrapper over the supported typed response-vector layouts. | typed response | variant | `std::variant` | throughout current code |
| creation/query helpers (`make_response_vector`, `response_num_orbitals`, `response_has_y_channel`, `response_is_unrestricted`, `response_alpha_factor`, `response_num_flat_slots`, `response_num_channels`) | inline functions | 543-660 area | UTILITY | Create/query the runtime response-vector layout and physical prefactors. | typed or variant response objects | vectors, sizes, booleans, doubles | variant visitation | initializer, solver, properties |
| access/shape helpers (`response_all`, `response_x`, `response_y`, `flatten_response`, `sync_response`, `get_flat`, `assign_all_and_sync`, `assign_flat_and_sync`) | inline functions | 577-698 area | UTILITY | Provide safe access and synchronization between typed channels and the flattened storage used by solvers and archives. | typed or variant response objects | refs / `void` | flat-channel bookkeeping | solver, IO, properties |

### `src/apps/molresponse_v2/ResponseBundle.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `ResponseBundle<R>` | class | 34 | UTILITY | Container for bundles of typed response states used by excited-state bundle code. | world/sizes/state vectors | object | flatten/sync/inner helpers | excited-state solver |
| constructors, iterators, `states`, `flatten_all`, `sync_all`, `to_flat`, `from_flat` | methods | 34-152 | UTILITY | Manage ownership and flatten/sync operations for a bundle of typed response states. | bundle/state data | objects/refs/`void` | response-vector helpers | excited-state solver |
| `inner` | function | 153 | PROPERTY | Computes bundle-wise inner products. | two bundles | `double` | state inner-product helpers | excited kernels |
| `ResponseBundleAllocator`, `make_response_bundle` | struct/function | 173, 208 | UTILITY | Allocate zero-initialized response bundles for excited-state iterations. | world, state/orbital counts | bundle | constructors | excited solver |

### `src/apps/molresponse_v2/ResponseSolver.hpp` / `ResponseSolver.cpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `response_vector_allocator` | struct | `ResponseSolver.hpp:26` | UTILITY | Allocator adapter for nonlinear/KAIN solvers working on `ResponseVector`. | world, orbital count, flags | allocator | `make_response_vector` | `iterate` |
| response-solver wrapper types/functions in `ResponseSolver.hpp` / `ResponseSolver.cpp` | types/functions | file scope | SOLVER | Provide the nonlinear solve wrapper used by `iterate` to apply KAIN acceleration over `ResponseVector`. | response-vector state, residuals | updated vectors | MADNESS nonlinear solver helpers | `iterate` |

### `src/apps/molresponse_v2/ResponseSolverUtils.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `print_iteration_table_border`, `print_iteration_line` | functions | 20, 26 | UTILITY | Print the standard per-iteration diagnostics table. | iteration diagnostics | `void` | printing only | solver diagnostics |
| `make_bsh_operators_response`, `compute_bsh_x_shift` | functions | 61, 84 | OPERATOR | Build BSH operators and stable energy shifts for linear-response updates. | world, manager, energies, frequency | operator list / `double` | operator factories | solver kernels |
| `inner`, `infer_function_k`, `infer_state_bundle_k`, `align_function_vector_protocol`, `align_state_bundle_protocol`, `align_response_vector_protocol`, `align_response_bundle_protocol`, `do_step_restriction` | functions | 90-169 | UTILITY | Provide protocol-alignment, inner-product, and step-limiting helpers shared by IO and solver code. | vectors/bundles/functions | values / bool / `void` | projection/truncate/math helpers | `ResponseIO`, excited restart, `iterate` |

### `src/apps/molresponse_v2/ResponseDebugLogger.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `ResponseDebugLogger` | class | 14 | IO | Collects per-state/per-iteration debug traces and can emit them as JSON or on-disk logs. | debug filename | object | internal JSON helpers | `JsonStateSolvePersistence`, solver |
| `ResponseDebugLogger` methods (`set_enabled`, `start_state`, `begin_iteration`, `log_value`, `end_iteration`, `finalize_state`, `write_to_disk`, `to_json`, timing/value table printers) | methods | 14-340 | IO | Record and serialize solver debug instrumentation for one run. | state ids, values, timing data | JSON / `void` | formatting helpers | persistence/solver |
| `TimedValueLogger` | class | 343 | UTILITY | Lightweight console-oriented value logger used by the solver. | static/iteration context | object | formatting helpers | solver/debug output |
| `TimedValueLogger` methods (`set_iteration_context`, `clear_iteration_context`, `set_console_enabled`, `log`, formatting helpers) | methods | 343-end | UTILITY | Emit compact timed value diagnostics without baking console policy into the solver. | labels, values, iteration metadata | `void` | printing helpers | solver |

### `src/apps/molresponse_v2/MolecularProperty.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `MolecularProperty` | struct | 11 | PROPERTY | Small typed descriptor for requested properties. | property type/value fields | object | none | response parameter/property plumbing |
| `property_type_to_string` | function | 25 | PROPERTY | Converts the property enum to its canonical string form. | `MolecularPropertyType` | `std::string` | switch | printing/logging |
| `print_requested_properties` | function | 39 | IO | Prints the requested-property list for diagnostics. | vector of properties | `void` | printing only | setup/debug |

### `src/apps/molresponse_v2/PropertyManager.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `components` stream overload | function | 36 | UTILITY | Formats a property-component string vector for printing. | `std::vector<std::string>` | stream | none | diagnostics |
| `PropRow`, `PropKey` | structs | 47, 55 | PROPERTY | Flat row/key representation used for persisted property tables. | fields | objects | comparison operators | `PropertyManager` |
| `to_json`, `from_json` for `PropRow` | functions | 89, 102 | IO | Serialize flat property rows to/from JSON. | row / JSON | JSON / row | JSON access | `PropertyManager` |
| `checked_size_to_long`, `checked_long_product`, `iso_timestamp` | functions | 116, 126, 138 | UTILITY | Safe tensor-extent conversion and stable timestamp formatting for property bookkeeping. | sizes or none | `long` / `std::string` | std helpers | property assembly |
| `try_claim_property_component_task` | function | 147 | IO | Claims a component task by atomically creating a lock file. | claim filename, optional payload | `bool` | `open(O_CREAT|O_EXCL)` | subgroup property/derived execution |
| `compute_response_inner_product_tensor` | function | 177 | PROPERTY | Contracts sets of response vectors into a dense tensor, optionally recording per-orbital contributions. | world, vector sets, save flags | `Tensor<double>` | `matrix_inner`, JSON contribution store | alpha/beta/Raman assembly |
| `PropertyManager` | class | 258 | PROPERTY | Owns the flat persisted property table and read/write/query/update operations on it. | world, filename | object | JSON file load | property stage |
| ctor, `to_json`, `save`, `has_alpha`, `has_beta`, `has_raman`, `get_alpha`, `get_beta`, `get_raman`, `set_alpha`, `set_beta`, `set_raman`, `print_table`, `freq_str` | methods | 260-446 | PROPERTY / PROPERTY (MIXED: IO for ctor/save) | Load, query, update, and print property rows. | frequencies, components, values | JSON / bool / optional values / `void` | row/key helpers | stage-3 property functions |
| `compute_alpha` | function | 450 | PROPERTY (MIXED: IO) | Loads final linear dipole response states and contracts them with dipole perturbations to produce the polarizability tensor. | world, state map, ground state, frequencies, directions, property manager | `void` | `perturbation_vector`, `load_response_vector`, `compute_response_inner_product_tensor`, `response_alpha_factor` | `molresponse_lib::compute_polarizability` |
| `compute_hessian` | function | ~560 | PROPERTY | Computes the nuclear Hessian and vibrational data from converged nuclear-response states and the SCF density. | world, state map, ground state, directions, `SCF` handle | `VibrationalResults` | `load_response_vector`, derivative functors, frequency analysis helpers | `molresponse_lib::compute_raman` |
| `compute_beta` template overloads | functions | 789, 1021 | PROPERTY (MIXED: naming logic, IO, task claims) | Compute beta or Raman components by building/loading the needed A/B/C response states and VBC derived state. | world, ground state, perturbation lists, frequency pairs, property manager, property type | `void` | `SimpleVBCComputer`, `load_response_vector`, `compute_response_inner_product_tensor`, property setters | hyperpolarizability/Raman assembly |
| `make_frequency_pairs_cartesian`, `select_beta_frequency_triplets` | functions | 1008, 1035 | UTILITY | Build the `(ωB, ωC)` grids used by beta assembly from the user-selected policy. | frequency vectors, policy flags | vector of pairs | simple loops | `compute_hyperpolarizability` |
| `compute_hyperpolarizability` | function | 1075 | PROPERTY (MIXED: model-specific policy) | Convenience wrapper that builds dipole A/BC lists and dispatches beta computation for the selected frequency-triplet policy. | world, ground state, frequencies, directions, property manager, flags | `void` | `select_beta_frequency_triplets`, `compute_beta` | `molresponse_lib::compute_hyperpolarizability` |
| `compute_Raman_components`, `compute_Raman` | functions | 1109, 1147 | PROPERTY (MIXED: naming logic, IO) | Build Raman component rows and then assemble them into Raman tensors over nuclear modes. | world, ground state, frequencies, dipole/nuclear directions, property manager | `void` / vector of tensors | `compute_beta`, property getters/setters | `molresponse_lib::compute_raman` |

### `src/apps/molresponse_v2/DerivedStatePlanner.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `DerivedStateDependency`, `DerivedStateRequest`, `DerivedStatePlan`, `DerivedStateGateEntry`, `DerivedStateGateReport` | structs | 38, 55, 93, 106, 125 | PERTURBATION | Represent derived-state requests, their linear-response dependencies, and dependency-gate results. | fields | objects | `to_json` methods | planning/runtime metadata |
| `DerivedStatePlanner` | class | 154 | PERTURBATION | Builds and evaluates derived-state plans, currently centered on VBC-driven quadratic response. | none | type | helper parsing functions | planning and stage 2d |
| helpers like `canonical_direction`, `parse_perturbation`, `canonical_property` | functions | file inline | UTILITY | Normalize property/perturbation tokens used by derived planning. | strings/chars | normalized values | parsing helpers | planner internals |
| `build_vbc_driven_quadratic_plan` | method | file inline | PERTURBATION | Expand quadratic requests into VBC-derived requests and their required linear dependencies. | `ResponseParameters`, `Molecule`, spin restriction, protocols | `DerivedStatePlan` | naming/descriptor helpers | `plan_required_states` |
| `evaluate_dependency_gate` | method | file inline | SOLVER | Determine which derived requests are ready from final linear-state convergence metadata. | plan, generated states, final protocol index, readiness callback | `DerivedStateGateReport` | dependency scans | `execute_derived_state_requests`, subgroup tail polling |
| `make_vbc_state` | method | file inline | PERTURBATION | Convert a derived request into the concrete `VBCResponseState` descriptor to execute. | `DerivedStateRequest`, spin restriction | `VBCResponseState` | `ResponseState` naming helpers | `run_derived_request` |

### `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp` / `ExcitedStateBundleSolver.cpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `ExcitedBundleSolverConfig` | struct | `ExcitedStateBundleSolver.hpp:15` | SOLVER | Config object for the excited-state bundle solver adapter. | fields | object | none | excited-stage setup |
| `ExcitedStateBundleSolver`, `ExcitedStateBundleNoopSolver` | classes | 22, 31 | SOLVER | Abstract interface and no-op implementation for bundle-solver backends. | config / virtual API | solver objects | none | excited-stage plumbing |
| `make_excited_state_bundle_solver_adapter` | function | `ExcitedStateBundleSolver.hpp:49`, `.cpp:3773` | SOLVER | Factory that chooses the active bundle-solver adapter implementation. | config | solver pointer | build flags / solver classes | excited-stage setup |
| naming/restart helpers: `protocol_key`, `alpha_suffix`, `assign_excited_state_names`, `make_excited_root_id`, `file_exists_world`, `RestartSnapshot`, `trial_space_archive_path`, `response_bundle_archive_path`, `archive_exists_world`, `read_restart_*`, `write_restart_*`, `read_restart_snapshot`, `write_restart_snapshot` | functions/structs | `.cpp:36-558` | IO | Handle excited-state protocol/root naming and on-disk restart snapshot read/write. | thresholds, paths, bundle/trial-space data | strings, structs, side effects | archive I/O | excited scaffold solver |
| `LocalizedGaussianGuess` | class | `.cpp:150` | PERTURBATION | Functor used to build localized excited-state trial guesses. | guess metadata | callable functor | `operator()` | excited solver |
| `RestartAwareExcitedScaffoldSolver` | class | `.cpp:615` | SOLVER (MIXED: naming, restart, protocol handling) | Current restart-aware excited-state scaffold/bundle adapter implementation. | solver config | solver object | restart helpers, bundle kernels | excited stage |
| key scaffold methods (`protocol_k`, `prepare_protocol`, `load_restart_seed`, `normalize_root_descriptors`, `sync_root_descriptors_for_protocol`, `ensure_ground_data`, `sort_trial_space`, `align_trial_space_protocol`, `ensure_state_names_for_protocol`) | methods/functions | `.cpp:980-1258+` | SOLVER (MIXED: naming, restart, protocol handling) | Prepare protocol context, align restart seeds, and keep stable root identity across protocol stages. | thresholds, world, snapshot data | vectors/metadata/`void` | `ResponseRecord` naming helpers, protocol alignment helpers | `ExcitedResponse` / bundle adapter |

### `src/apps/molresponse_v2/ExcitedResponse.hpp` / `ExcitedResponse.cpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `ExcitedSolverParams`, `ExcitedIterDiagnostics`, `ExcitedSolverConfig` | structs | `ExcitedResponse.hpp` | SOLVER | Hold configuration and diagnostics for the current excited-state protocol solver. | fields | objects | none | excited stage |
| `ExcitedResponse` | class | `ExcitedResponse.hpp` | SOLVER | Owns the current protocol-aware excited-state solve path used by `MolresponseLib`. | solver config | object | bundle solver helpers | excited stage |
| helper functions like `read_gs_energy_from_calc_info`, `ensure_ground_state` | functions | `.cpp` file scope | IO | Load the SCF energy/ground-state context needed by excited-state stages. | paths/world | energy / side effects | checkpoint JSON/archive helpers | `ExcitedResponse` |
| `ExcitedResponse::solve_protocol` and supporting methods | methods | `.cpp` main body | SOLVER (MIXED: restart logic, protocol handling, model-specific TDA/full branching) | Solve one excited-state protocol stage or skip it from restart-ready metadata. | `World&`, `ExcitedProtocolInput` | `ExcitedProtocolResult` | bundle solver adapter, restart helpers, protocol setup | `execute_excited_state_bundle_stage` |

### `src/apps/molresponse_v2/ResponseKernels.hpp` / `ResponseVectorKernels.hpp` / `src/apps/molresponse_v2/ops/*.hpp`

These files are inline kernel code for the current typed response-vector solver.

| Symbol(s) | Kind | Line(s) | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `ResponseStatePotentials` and helper kernels in `ResponseKernels.hpp` | inline functions/types | header | OPERATOR | Build the local, gamma, and lambda/source pieces for current linear-response equations. | world, response manager, ground state, response vector | operator/source bundles | ground-state/operator helpers | `iterate` |
| `K`, `compute_ground_exchange`, `compute_gamma_response` in `ResponseVectorKernels.hpp` | functions | 32, 57, 91 | OPERATOR | Low-level vector-kernel building blocks for exchange and response-density operator application. | world, ground data, vectors | functions/vectors | MADNESS operators | `ResponseKernels.hpp`, `iterate` |
| per-model kernels in `ops/StaticRestrictedOps.hpp`, `ops/DynamicRestrictedOps.hpp`, `ops/TDARestrictedOps.hpp`, `ops/StaticUnrestrictedOps.hpp`, `ops/DynamicUnrestrictedOps.hpp` | inline functions | file scope | OPERATOR | Provide `alpha_factor`, density build, BSH operator build, and coupled-equation application for each typed response-vector model. | world, manager, ground state, typed response vectors, points | typed vectors/operators/densities | `ResponseVectorKernels`, `ResponseSolverUtils` | `iterate`, property helpers |

### `src/apps/molresponse_v2/excited_ops/*.hpp`

These files are inline operator kernels for excited-state bundle iterations.

| Symbol(s) | Kind | Line(s) | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `BundleKernels.hpp` functions (`transpose_bundle`, `inner`, `metric_gram_schmidt`, `build_rotation_matrices`, `compute_excited_potentials`, `diagonalize_and_rotate_bundle`, `snapshot_bundle`, `compute_response_densities`, `ExcitedStateStep`, `apply_step_project_normalize`, `print_potential_diagnostics`) | inline functions/structs | file scope | SOLVER / OPERATOR | Implement the bundle-level linear algebra and potential application used by excited-state iterations. | world, bundles, ground state, configs | bundles, tensors, diagnostics | response/vector kernel helpers | `ExcitedResponse`, bundle solver |
| `RestrictedBundleCommon.hpp`, `RestrictedFullBundleOps.hpp`, `RestrictedTDABundleOps.hpp` helpers like `metric_inner` | inline functions | file scope | OPERATOR | Provide restricted full/TDA metric and operator specializations for excited-state bundles. | world, bundles | tensors/scalars | bundle kernels | excited bundle iterations |

### `src/apps/molresponse_v2/broadcast_json.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `broadcast_json_file` | function | 11 | IO | Loads a JSON file on rank 0 and broadcasts it to all ranks. | `World&`, filepath | `json` | file I/O, `broadcast_serializable` | `PropertyManager` ctor |

### `src/apps/molresponse_v2/VBCMacrotask.hpp`

| Symbol | Kind | Line | Layer | Purpose | Inputs | Outputs | Calls | Called by |
|---|---|---:|---|---|---|---|---|---|
| `SimpleVBCComputer` | class | 14 | PROPERTY | Computes and saves VBC derived states and exposes helper access to the B and C component vectors. | world, ground state | object | response I/O and VBC kernels | derived stage, beta/Raman assembly |
| key methods (`compute_and_save`, `get_BC_vecs`, `make_zeta_bc`) | methods | header inline | PROPERTY / OPERATOR | Build the VBC state, save it, and provide the auxiliary vectors used by beta/Raman contractions. | derived descriptor / response vectors | response vectors / `void` | `ResponseIO`, response-vector helpers | `run_derived_request`, `compute_beta` |

## Coverage Notes

- As in `legacy_inventory.md`, inline-heavy headers are grouped when they primarily provide container or helper APIs rather than independent algorithmic units.
- The current codebase is more modular than the legacy tree, but a direct one-symbol-per-row dump would be dominated by inline accessors, variant helpers, and JSON boilerplate; this inventory instead groups coherent inline method families the same way the legacy inventory grouped container logic.
- Excluded tests and markdown/reference files were intentionally left out so the comparison stays focused on production current response code.
