# Current Excited-State Implementation Report

Statements below are facts unless explicitly labeled `Inference:` or `Ambiguity:`.

Comparison baseline: the legacy solver documented in `docs/legacy_excited_state_report.md`.

## 1. Scope

This report analyzes the current excited-state reintegration work in the new `molresponse_v2` architecture. The primary code paths examined were:

- the response workflow driver and stage orchestrator in `src/madness/chem/MolresponseLib.hpp`
- the excited-stage adapter in `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp` and `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- supporting state, metadata, archive, scheduler, and test files under `src/apps/molresponse_v2`, `src/madness/chem`, and `src/apps/madqc_v2`

This report covers:

- current excited-state entry points and execution flow
- the current class/data model
- protocol planning, iteration, restart, metadata, and archive handling
- what parts of the old excited-state design appear to have been reintegrated already
- what is still partial or missing

This report does not cover:

- a redesign of the current architecture
- implementation work
- a full analysis of the linear-response solver beyond the interfaces and scheduling machinery that the excited-stage code depends on
- runtime validation beyond inspection of the existing regression tests

## 2. High-Level Architecture

The new code does not yet reintroduce the legacy all-roots `X_space` model. Instead, excited-state work is inserted as a distinct Stage 2c adapter inside the `molresponse_lib` workflow. The orchestrator plans an `ExcitedStateBundlePlan` from response parameters, solves all linear response states first, then runs a protocol-indexed excited-state bundle stage, and finally runs derived-state requests and property assembly. Reference: `src/madness/chem/MolresponseLib.hpp :: molresponse_lib::plan_required_states` (lines 999-1049), `src/madness/chem/MolresponseLib.hpp :: molresponse_lib::solve_all_states` (lines 3520-3568).

The current excited implementation is centered on the abstract interface `ExcitedStateBundleSolver`, with a concrete scaffold implementation `RestartAwareExcitedScaffoldSolver`. The factory returns the scaffold when a ground-state archive path is available; otherwise it falls back to a no-op solver. Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp :: ExcitedStateBundleSolver` (lines 44-51), `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: make_excited_state_bundle_solver_adapter` (lines 2520-2525).

Architecturally, the new excited stage already has:

- parameter plumbing (`excited.enable`, `excited.num_states`, `excited.tda`, `excited.guess_max_iter`, `excited.maxiter`, `excited.maxsub`, `excited.owner_group`)
- protocol-wise planning metadata
- protocol restart snapshots
- protocol-local root naming and root manifests
- an executable iterative scaffold solver
- tests for metadata presence, restart reuse, and protocol reprojection

The main incompleteness is numerical and downstream integration, not driver plumbing. The code comments and local checklist both describe the current solver as a scaffold or partial reintegration rather than finished legacy parity. References: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` comment at lines 333-335, `src/apps/molresponse_v2/EXCITED_STATE_EXECUTION_CHECKLIST.md` lines 63-98.

## 3. Relevant Files and Components

| File | Role | Key symbols |
| --- | --- | --- |
| `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp` | Excited-stage solver interface and protocol I/O contracts. | `ExcitedBundleSolverConfig`, `ExcitedBundleProtocolInput`, `ExcitedBundleProtocolResult`, `ExcitedStateBundleSolver`, `ExcitedStateBundleNoopSolver` |
| `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` | Current scaffold excited-state implementation: protocol prep, guess generation, restart selection, iteration, bundle rotation, snapshot I/O. | `RestartSnapshot`, `ExcitedTrialSpace`, `RestartAwareExcitedScaffoldSolver`, `ExcitedProtocolWorkflow`, `assign_excited_state_names`, `initialize_protocol_guess`, `iterate` |
| `src/madness/chem/MolresponseLib.hpp` | Top-level response workflow orchestration. Adds Stage 2c excited execution and writes excited metadata into final task metadata. | `ExcitedStateBundlePlan`, `build_excited_state_bundle_plan`, `execute_excited_state_bundle_stage`, `build_state_stage_metadata`, `solve_all_states`, `run_response` |
| `src/madness/chem/ResponseParameters.hpp` | User parameter surface for new excited controls. | `ResponseParameters::excited_enable`, `excited_num_states`, `excited_tda`, `excited_guess_max_iter`, `excited_maxiter`, `excited_maxsub`, `excited_owner_group` |
| `src/apps/molresponse_v2/ResponseRecord.hpp` | JSON metadata persistence schema. Contains the `excited_states` subtree and helper methods. | `ResponseRecord2`, `initialize_excited_bundle`, `record_excited_protocol_status`, `record_excited_protocol_timing`, `record_excited_protocol_energies` |
| `src/apps/molresponse_v2/ResponseVector.hpp` | Variant representation reused by the excited solver for one root at a time. | `ResponseVector`, `StaticRestrictedResponse`, `DynamicRestrictedResponse`, `StaticUnrestrictedResponse`, `DynamicUnrestrictedResponse` |
| `src/apps/molresponse_v2/ResponseIO.hpp` | Linear/derived response-vector archive I/O. Important mainly because the excited stage currently does not use it. | `save_response_vector`, `load_response_vector` |
| `src/apps/molresponse_v2/ResponseState.hpp` / `.cpp` | Linear response descriptors and perturbation vector builders used by the surrounding Stage 2 flow. | `LinearResponseDescriptor`, `LinearResponsePoint`, `perturbation_vector(...)` |
| `src/apps/molresponse_v2/FrequencyLoop.hpp` / `.cpp` | Existing linear solver loop. Important for understanding how the excited stage is inserted after the linear protocol/frequency workflow. | `iterate`, `computeFrequencyLoop`, `computeFrequencyPoint`, `StateSolvePersistence` |
| `src/apps/molresponse_v2/ResponseSolver.hpp` / `.cpp` | Linear-response numerical kernel. The excited scaffold reuses some low-level ingredients but not the same KAIN driver. | `response_solver`, `compute_density`, `make_bsh_operators`, `CoupledResponseEquations` |
| `src/apps/molresponse_v2/ResponseSolverUtils.hpp` | Shared utility functions used directly by the excited scaffold, especially protocol reprojection helpers. | `inner`, `infer_state_bundle_k`, `align_state_bundle_protocol` |
| `src/apps/molresponse_v2/StateGenerator.hpp` | Stage-1 linear state generator. Excited states are not generated here; only the plan is. | `GeneratedStateData`, `StateGenerator` |
| `src/apps/molresponse_v2/StateParallelPlanner.hpp` | Builds the linear channel/point scheduling plan. Excited owner-group data is derived from it, but not yet enforced in the excited solver. | `StateParallelPlan`, `PointOwnershipScheduler`, `StateParallelPlanner::build` |
| `src/madness/chem/WorkflowBuilders.hpp` | Connects `madqc --wf=response` and `molresponse2` into `ResponseApplication<molresponse_lib>`. | `add_response_workflow_drivers` |
| `src/madness/chem/Applications.hpp` | Response task wrapper that calls `molresponse_lib::run_response` and places metadata into the response task entry in `*.calc_info.json`. | `ResponseApplication<Library>::run`, `ResponseApplication<Library>::results` |
| `src/apps/madqc_v2/test_molresponse_excited_metadata_smoke.py` | Regression test for metadata/planning fields. | metadata smoke expectations |
| `src/apps/madqc_v2/test_molresponse_excited_restart_reuse.py` | Regression test for restart snapshot reuse. | restart reuse expectations |
| `src/apps/madqc_v2/test_molresponse_excited_protocol_projection.py` | Regression test for cross-protocol trial-space reprojection. | `EXCITED_TRIAL_REPROJECT` expectation |

## 4. Current Execution Flow

### 4.1 Driver-level entry point

The main production entry point is `madqc --wf=response`. `madqc.cpp` builds a workflow, and `WorkflowBuilders.hpp` adds:

1. an SCF driver
2. a response driver wrapping `ResponseApplication<molresponse_lib>`

References:

- `src/apps/madqc_v2/madqc.cpp :: main` (lines 140-183)
- `src/madness/chem/WorkflowBuilders.hpp :: add_response_workflow_drivers` (lines 75-84)
- `src/madness/chem/Applications.hpp :: ResponseApplication::run` (lines 252-265)

The compatibility executable `molresponse2` reaches the same path through `add_response_workflow_drivers`. Reference: `src/apps/molresponse_v2/molresponse2.cpp :: main` (lines 29-84).

### 4.2 Response workflow entry

`ResponseApplication::run` creates `task_<n>/molresponse`, switches into that directory, and calls `molresponse_lib::run_response(...)`. That means excited restart files are written under the response task directory, not the repository root. References:

- `src/madness/chem/Applications.hpp :: ResponseApplication::run` (lines 252-265)
- `src/madness/chem/MolresponseLib.hpp :: molresponse_lib::label` (line 106)

`Workflow::run` then serializes the response task metadata into `<prefix>.calc_info.json`. Reference: `src/madness/chem/Drivers.hpp :: qcapp::Workflow::run` (lines 97-125).

### 4.3 Stage 1: context and planning

`molresponse_lib::run_response` performs:

1. parse response inputs into `ResponseParameters`
2. build a `GroundContext`
3. build `PlannedStates`
4. solve all states
5. compute requested properties

Reference: `src/madness/chem/MolresponseLib.hpp :: molresponse_lib::run_response` (lines 4302-4354).

The excited-specific planning work happens in `plan_required_states(...)`:

1. generate linear states with `StateGenerator`
2. build derived-state plan with `DerivedStatePlanner`
3. build state-parallel plan for linear states with `StateParallelPlanner`
4. translate response knobs into `ExcitedStateBundlePlan`

References:

- `src/madness/chem/MolresponseLib.hpp :: plan_required_states` (lines 999-1049)
- `src/madness/chem/MolresponseLib.hpp :: build_excited_state_bundle_plan` (lines 977-997)

Important fact: Stage 1 does not create explicit excited-state root objects. It creates only an excited bundle plan. Root generation is deferred to Stage 2c inside the excited solver. Reference: `src/madness/chem/MolresponseLib.hpp :: plan_required_states` (lines 1004-1048) versus `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_fresh_guess` (lines 2062-2094).

### 4.4 Stage 2a/2b: linear states first

Before the excited stage runs, the existing linear-response machinery executes. `solve_all_states(...)`:

1. builds the runtime state/point ownership schedule
2. runs linear response in subgroup mode or serial mode
3. validates final-protocol linear readiness
4. runs the excited bundle stage
5. runs derived-state requests
6. assembles final stage metadata

Reference: `src/madness/chem/MolresponseLib.hpp :: solve_all_states` (lines 3520-3568).

The linear solver path uses `computeFrequencyLoop(...)` or `computeFrequencyPoint(...)`, which operate on `(state, protocol, frequency)` points and seed guesses from same-protocol archives, coarser-protocol archives, previous frequencies, or initializer guesses. References:

- `src/apps/molresponse_v2/FrequencyLoop.cpp :: computeFrequencyLoop` (lines 566-618)
- `src/apps/molresponse_v2/FrequencyLoop.cpp :: computeFrequencyPoint` (lines 620-650)
- `src/apps/molresponse_v2/FrequencyLoop.cpp :: seed_response_guess_for_point` (lines 187-260)

This matters because the excited stage is not integrated into that frequency loop. It is a separate protocol-only stage that starts after linear-state solving.

### 4.5 Stage 2c: excited bundle execution

The top-level excited entry point is `molresponse_lib::execute_excited_state_bundle_stage(...)`. It:

1. reads `PlannedStates::excited_state_bundle_plan`
2. constructs an `ExcitedBundleSolverConfig`
3. creates an excited solver adapter with `make_excited_state_bundle_solver_adapter(...)`
4. ensures the in-memory `state_metadata_json["excited_states"]` subtree exists
5. loops over every protocol threshold in the plan
6. creates or extends the per-protocol metadata node
7. decides whether the adapter needs to run
8. if needed, builds `ExcitedBundleProtocolInput` and calls `solver_adapter->solve_protocol(world, protocol_input)`
9. writes the returned status/energies/names/residuals back into `state_metadata_json`
10. accumulates an `excited_state_planner.execution` summary

Reference: `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage` (lines 2784-3092).

Per-protocol dispatch is controlled by:

- `saved` and `converged` booleans already present in the metadata node
- `plan.enabled`

If the node is already `saved && converged`, the stage can skip adapter execution and mark the protocol as `restart_ready_skip`. Reference: `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage` (lines 2908-2918, 2922-2956).

### 4.6 Adapter execution flow

`solve_protocol(...)` on the concrete scaffold implementation routes to `ExcitedProtocolWorkflow::solve_protocol(...)`. The protocol-level flow is:

1. read protocol restart snapshot file
2. if metadata says `restart_saved && restart_converged`, optionally load restart data and skip
3. call `prepare_protocol(...)` to set `FunctionDefaults<3>::k` and threshold
4. call `initialize_protocol_guess(...)`
5. call `iterate(...)`
6. prefix the final `stage_status` with the initialization strategy

Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: ExcitedProtocolWorkflow::solve_protocol` (lines 356-399).

### 4.7 Initialization flow inside the scaffold

`initialize_protocol_guess(...)` is the new-stage equivalent of legacy initialization. It:

1. ensures ground-state/preliminary data for the current protocol
2. chooses a restart seed using `select_restart_seed(...)`
3. if a restart seed exists, loads energies/state names/residuals/trial states with `load_restart_guess(...)`
4. calls `ensure_trial_space_matches_guess(...)`
5. otherwise, if an in-memory guess already exists, carries it over
6. otherwise, builds a fresh guess with `build_fresh_guess(...)`
7. for fresh guesses, writes a generic guess archive immediately

References:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: initialize_protocol_guess` (lines 2256-2307)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: select_restart_seed` (lines 2175-2254)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_fresh_guess` (lines 2062-2094)

Restart seed precedence is:

1. current protocol snapshot, if present and not considered stalled
2. nearest lower protocol snapshot
3. generic guess archive
4. in-memory carryover guess
5. fresh guess generation

Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: select_restart_seed` (lines 2212-2254), `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: initialize_protocol_guess` (lines 2260-2307).

### 4.8 Iteration flow inside the scaffold

`iterate(...)` is the current main excited bundle solve loop. It:

1. ensures ground data and trial-space/guess consistency
2. converts trial `x` states into a bundle of `ResponseVector` objects
3. loops up to `input.maxiter`
4. for restricted variants:
   - computes bundle potentials
   - diagonalizes and rotates the whole bundle
   - updates each root from rotated `v0`/`gamma`
5. for unrestricted variants:
   - skips bundle rotation
   - iterates each state independently
6. if no bundle rotation occurred, updates `omega_` from `estimate_state_energies(...)`
7. records per-root residual norms and per-iteration maximum residuals
8. checks convergence against `max_residual <= max(10 * threshold, 1e-7)`
9. writes a protocol restart snapshot and updates the generic guess archive
10. returns `ExcitedBundleProtocolResult`

Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate` (lines 2310-2503).

### 4.9 Finalization and output

After Stage 2c, the workflow:

1. runs derived-state execution
2. assembles metadata blocks including `excited_state_planner`
3. returns metadata to the response application
4. serializes the response task into `*.calc_info.json`

References:

- `src/madness/chem/MolresponseLib.hpp :: build_state_stage_metadata` (lines 3455-3509)
- `src/madness/chem/MolresponseLib.hpp :: solve_all_states` (lines 3550-3567)
- `src/madness/chem/Applications.hpp :: ResponseApplication::results` (lines 271-275)

Important fact: the property stage currently uses solved linear and derived states, but there is no corresponding consumer of `excited_states` metadata or excited restart files in property assembly. References:

- `src/madness/chem/MolresponseLib.hpp :: compute_polarizability` (lines 3603-3689)
- `src/madness/chem/MolresponseLib.hpp :: compute_hyperpolarizability` (lines 3691-3705)
- `rg` over `src/` shows `roots` is only created/merged in excited metadata, not consumed by property code

## 5. Class and Data Model

### 5.1 There is no standalone excited-state object

The new code currently has no dedicated `ExcitedState` or `ExcitedStateDescriptor` class in the runtime solver path. Instead, one excited root is represented in two different forms at different phases:

1. as one `vector_real_function_3d` inside `ExcitedTrialSpace::x_states`
2. as one `ResponseVector` inside a bundle `std::vector<ResponseVector>`

Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: ExcitedTrialSpace` (lines 326-331), `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_response_bundle_from_trial_space` (lines 1947-1955).

### 5.2 Bundle-level containers

The current bundle-level state is spread across:

- `trial_space_`
- `omega_`
- `state_names_`
- `residual_norms_`
- `iteration_max_residuals_`
- temporary `response_states` inside `iterate(...)`

Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: ExcitedProtocolWorkflow` member fields (lines 416-433).

This means the operative bundle is not owned by one container analogous to legacy `X_space`. The active representation is split across several parallel arrays plus one temporary response-vector bundle.

### 5.3 `ExcitedTrialSpace`

`ExcitedTrialSpace` contains:

- `num_states`
- `num_orbitals`
- `x_states`
- `omega`

Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: ExcitedTrialSpace` (lines 326-331).

Role:

- stores the protocol-independent or restartable guess state
- stores only `x`-channel functions, not full `x/y` response bundles
- is the only function-valued data persisted in excited restart snapshots

Important consequence: dynamic `y` channels are not persisted across protocol boundaries or restarts. The scaffold regenerates full `ResponseVector` objects from `x_states` each time, with `y` channels zero-initialized when `tda == false`. References:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: make_state_response_vector` (lines 1898-1935)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: sync_trial_space_from_response_bundle` (lines 1958-1964)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: write_restart_snapshot` (lines 310-324)

### 5.4 `ResponseVector` reuse

The excited stage reuses the linear `ResponseVector` variant types:

- `StaticRestrictedResponse`
- `DynamicRestrictedResponse`
- `StaticUnrestrictedResponse`
- `DynamicUnrestrictedResponse`

Reference: `src/apps/molresponse_v2/ResponseVector.hpp` (lines 16-18, 45-152).

The excited solver chooses the variant shape from two inputs:

- `tda`
- `ground_data_->isSpinRestricted()`

Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: make_state_response_vector` (lines 1898-1935).

### 5.5 Root identity

Excited roots are currently identified by:

1. their position in `omega_` / `trial_space_.x_states` / `response_states`
2. a protocol-local `state_names_` array
3. a per-protocol metadata manifest `roots = [{root_index, name, energy}, ...]`

References:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: assign_excited_state_names` (lines 110-140)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: ensure_state_names_for_protocol` (lines 2164-2173)
- `src/madness/chem/MolresponseLib.hpp :: build_excited_root_manifest` (lines 2768-2781)

Important fact: root names are assigned from energy ordering, with grouping tolerance `max(1e-12, 10 * protocol_threshold)`. Near-degenerate clusters get suffixes `a`, `b`, `c`, ... . This makes naming protocol-dependent. Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: assign_excited_state_names` (lines 118-139).

### 5.6 Metadata model

There are two excited metadata layers:

1. `metadata["excited_states"]`
2. `metadata["excited_state_planner"]`

`excited_states` stores the actual per-protocol plan and protocol results. `excited_state_planner` stores a summary of execution counters and protocol events. Reference: `src/madness/chem/MolresponseLib.hpp :: build_state_stage_metadata` (lines 3495-3508).

The protocol subtree can contain:

- `saved`
- `converged`
- `timings`
- `energies`
- `state_names`
- `roots`
- `iterations`
- `residual_norms`
- `iteration_max_residuals`
- `stage_status`
- `owner_group`

Reference: `src/madness/chem/MolresponseLib.hpp :: ensure_excited_protocol_placeholder_node` (lines 2716-2766), `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage` (lines 2997-3022).

Important fact: `ResponseRecord2` initializes only a smaller schema (`saved`, `converged`, `timings`, `energies`) when seeding `response_metadata.json`. The richer schema is added later by `ensure_excited_protocol_placeholder_node(...)` inside Stage 2c. References:

- `src/apps/molresponse_v2/ResponseRecord.hpp :: initialize_excited_bundle` (lines 84-108)
- `src/apps/molresponse_v2/ResponseRecord.hpp :: ensure_excited_protocol` (lines 563-582)
- `src/madness/chem/MolresponseLib.hpp :: ensure_excited_protocol_placeholder_node` (lines 2716-2766)

## 6. Implemented Functions

### `molresponse_lib::build_excited_state_bundle_plan`

- Location: `src/madness/chem/MolresponseLib.hpp` (lines 977-997)
- Purpose: convert response parameters into a protocol-indexed excited bundle plan.
- Inputs: `CalculationParameters`, `ResponseParameters`, `StateParallelPlan`.
- Outputs: `ExcitedStateBundlePlan`.
- Side effects: none.
- Called by: `plan_required_states(...)`.
- Calls: no helper beyond simple field extraction.
- Algorithmic role: binds user-facing excited controls to the Stage 2c orchestration layer.
- Notes / assumptions:
  - `owner_group` is clamped to `mapping_groups - 1`.
  - `protocols` comes directly from `calc_params.protocol()`.

### `molresponse_lib::execute_excited_state_bundle_stage`

- Location: `src/madness/chem/MolresponseLib.hpp` (lines 2784-3092)
- Purpose: run the excited bundle solver once per protocol and populate metadata.
- Inputs: `World`, `PlannedStates`, `GroundContext`, `ResponseParameters`, mutable `state_metadata_json`.
- Outputs: `ExcitedExecutionResult`.
- Side effects:
  - mutates `state_metadata_json["excited_states"]`
  - prints summary lines
  - calls the adapter, which may read/write restart snapshot files
- Called by: `solve_all_states(...)`.
- Calls:
  - `make_excited_state_bundle_solver_adapter(...)`
  - `ensure_excited_protocol_placeholder_node(...)`
  - `solver_adapter->solve_protocol(...)`
  - `build_excited_root_manifest(...)`
  - `broadcast_json_object(...)`
- Algorithmic role:
  - stage orchestration, not numerical solving
  - protocol-level dispatch and metadata aggregation
- Notes / assumptions:
  - execution is protocol-only; there is no excited frequency loop
  - `owner_group` is recorded in metadata and passed to the solver input, but the current code does not dispatch to a subgroup here
  - `saved` means a bundle snapshot exists, not necessarily convergence

### `ExcitedProtocolWorkflow::solve_protocol`

- Location: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 356-399)
- Purpose: one protocol solve entry point for the scaffold adapter.
- Inputs: `World`, `ExcitedBundleProtocolInput`.
- Outputs: `ExcitedBundleProtocolResult`.
- Side effects:
  - reads per-protocol restart snapshot
  - may update internal workflow state
  - may run `prepare_protocol`, initialization, iteration, and snapshot writes
- Called by: `RestartAwareExcitedScaffoldSolver::solve_protocol`.
- Calls:
  - `read_restart_snapshot(...)`
  - `load_restart_guess(...)`
  - `ensure_trial_space_matches_guess(...)`
  - `prepare_protocol(...)`
  - `initialize_protocol_guess(...)`
  - `iterate(...)`
- Algorithmic role:
  - top-level protocol controller
  - composes initialization strategy with iteration outcome
- Notes / assumptions:
  - if the caller marks the protocol already saved+converged, this function can skip solving
  - stage-status strings are concatenated with `__`

### `ExcitedProtocolWorkflow::prepare_protocol`

- Location: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 470-481)
- Purpose: set the MADNESS function defaults for the current excited protocol.
- Inputs: `World`, threshold.
- Outputs: none.
- Side effects:
  - `FunctionDefaults<3>::set_k(...)`
  - `FunctionDefaults<3>::set_thresh(...)`
- Called by: `solve_protocol(...)`.
- Calls: `protocol_k(...)`.
- Algorithmic role: protocol-to-wavelet-resolution bridge.
- Notes / assumptions:
  - hardcoded threshold-to-`k` mapping: `1e-2 -> 4`, `1e-4 -> 6`, `1e-6 -> 8`, `1e-8 -> 10`, else `12`
  - comment explicitly says full operator/ground setup is still to be moved here later

### `ExcitedProtocolWorkflow::select_restart_seed`

- Location: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 2175-2254)
- Purpose: choose the best available restart seed for a protocol.
- Inputs: `World`, `ExcitedBundleProtocolInput`, current protocol `RestartSnapshot`.
- Outputs: `RestartSeedChoice`.
- Side effects: prints skipped-stalled diagnostics.
- Called by: `initialize_protocol_guess(...)`.
- Calls:
  - `read_restart_snapshot(...)`
  - local stalled-snapshot filter
- Algorithmic role:
  - restart policy
  - same-protocol > lower-protocol > generic guess archive
- Notes / assumptions:
  - stalled snapshots are rejected if iteration residual history stays near `1.0`
  - this logic uses `iteration_max_residuals`, not metadata status

### `ExcitedProtocolWorkflow::initialize_protocol_guess`

- Location: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 2256-2307)
- Purpose: establish the initial bundle guess for one protocol.
- Inputs: `World`, `ExcitedBundleProtocolInput`, current protocol `RestartSnapshot`.
- Outputs: `ProtocolInitResult`.
- Side effects:
  - mutates `trial_space_`, `omega_`, `state_names_`, residual arrays
  - may write the generic guess archive
- Called by: `solve_protocol(...)`.
- Calls:
  - `ensure_ground_data(...)`
  - `select_restart_seed(...)`
  - `load_restart_guess(...)`
  - `ensure_trial_space_matches_guess(...)`
  - `build_fresh_guess(...)`
  - `ensure_state_names_for_protocol(...)`
  - `write_guess_archive(...)`
- Algorithmic role: scaffold equivalent of legacy `initialize`.
- Notes / assumptions:
  - restart reuse is tracked independently from convergence
  - fresh guesses are archived before the main iteration runs

### `ExcitedProtocolWorkflow::build_fresh_guess`

- Location: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 2062-2094)
- Purpose: create a new trial bundle when no restart/carryover guess exists.
- Inputs: `World`, `ExcitedBundleProtocolInput`.
- Outputs: none; mutates workflow members.
- Side effects:
  - populates `trial_space_`, `omega_`, `residual_norms_`, `last_iteration_count_`
- Called by: `initialize_protocol_guess(...)`, `iterate(...)`.
- Calls:
  - `create_trial_functions2(...)`
  - `make_nwchem_trial(...)`
  - `create_trial_functions(...)`
  - `make_random_trial(...)`
  - `iterate_trial(...)`
- Algorithmic role: guess generation and preconditioning.
- Notes / assumptions:
  - selection is heuristic and depends on `tda` and `num_states`
  - `make_nwchem_trial(...)` is only a naming holdover; it currently returns a derivative-based trial, not an NWChem-backed guess. Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` lines 1849-1853.

### `ExcitedProtocolWorkflow::iterate_trial`

- Location: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 1867-1896)
- Purpose: improve the trial `x`-space before the main solve.
- Inputs: `World`, `ExcitedBundleProtocolInput`, mutable `ExcitedTrialSpace`.
- Outputs: none.
- Side effects: mutates `trial.x_states` and `trial.omega`.
- Called by: `build_fresh_guess(...)`, `ensure_trial_space_matches_guess(...)`.
- Calls:
  - `project_and_orthonormalize(...)`
  - `estimate_state_energies(...)`
  - `sort_by_energy(...)`
- Algorithmic role:
  - cheap guess refinement
  - not the main solve loop
- Notes / assumptions:
  - energy updates are mixed `0.6 old + 0.4 estimate`
  - this is heuristic, not a generalized eigen solve

### `ExcitedProtocolWorkflow::align_trial_space_protocol`

- Location: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 2096-2118)
- Purpose: reproject stored trial states when the protocol changes `k` or threshold.
- Inputs: `World`.
- Outputs: none.
- Side effects: mutates `trial_space_.x_states`.
- Called by: `ensure_trial_space_matches_guess(...)`.
- Calls:
  - `ResponseSolverUtils::infer_state_bundle_k(...)`
  - `ResponseSolverUtils::align_state_bundle_protocol(...)`
  - `project_and_orthonormalize(...)`
- Algorithmic role: restart continuity across protocol ladders.
- Notes / assumptions:
  - tested indirectly by `test_molresponse_excited_protocol_projection.py`

### `ExcitedProtocolWorkflow::make_state_response_vector`

- Location: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 1898-1935)
- Purpose: turn one `x` trial state into a full `ResponseVector`.
- Inputs: `World`, one `vector_real_function_3d`, `tda`.
- Outputs: one `ResponseVector`.
- Side effects: none beyond function allocation.
- Called by: `build_response_bundle_from_trial_space(...)`.
- Calls: `align_component_to_ground(...)`, `zero_functions_compressed(...)`.
- Algorithmic role: bridge from trial-space-only representation to response-vector update representation.
- Notes / assumptions:
  - non-TDA dynamic channels start from zero
  - unrestricted channels also start from zero except `x_alpha`

### `ExcitedProtocolWorkflow::diagonalize_excited_bundle`

- Location: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 988-1140)
- Purpose: diagonalize the restricted bundle in the current basis.
- Inputs: `World`, overlap-like matrix `S`, operator-like matrix `A`, mutable `omega`, output rotation `U_out`.
- Outputs: `bool` success flag.
- Side effects: none outside output mutation.
- Called by: `rotate_excited_bundle_states(...)`.
- Calls:
  - `sygvp(world, A, S, 1, eigenvectors, eigenvalues)`
  - `svd(...)` for near-degenerate clusters
- Algorithmic role:
  - solves a generalized eigenproblem
  - orders, phase-fixes, and de-rotates degenerate clusters
- Notes / assumptions:
  - comment explicitly references legacy `excited_eig`
  - only used for restricted bundle rotation

### `ExcitedProtocolWorkflow::rotate_excited_bundle_states`

- Location: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 1143-1175)
- Purpose: rotate the restricted bundle and its accompanying potentials using the generalized-eigenvector matrix.
- Inputs: `World`, mutable restricted state bundle, mutable potentials, mutable `omega`, iteration index.
- Outputs: `bool` indicating whether rotation succeeded.
- Side effects: mutates state order and `omega`.
- Called by: `iterate_typed_bundle_legacy_sequence(...)`, `rotate_response_bundle(...)`.
- Calls:
  - `build_rotation_matrices(...)`
  - `diagonalize_excited_bundle(...)`
  - `rotate_bundle_states(...)`
- Algorithmic role: bundle-coupled root update and ordering step.
- Notes / assumptions:
  - unrestricted bundle rotation is explicitly skipped as work-in-progress. Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: rotate_response_bundle` (lines 1195-1225).

### `ExcitedProtocolWorkflow::iterate`

- Location: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 2310-2503)
- Purpose: current main excited solver loop.
- Inputs: `World`, `ExcitedBundleProtocolInput`, `restart_seed_reused`.
- Outputs: `ExcitedBundleProtocolResult`.
- Side effects:
  - mutates all workflow bundle state
  - writes protocol restart snapshot
  - writes guess archive
- Called by: `solve_protocol(...)`.
- Calls:
  - `ensure_ground_data(...)`
  - `ensure_trial_space_matches_guess(...)`
  - `build_response_bundle_from_trial_space(...)`
  - `iterate_typed_bundle_legacy_sequence(...)`
  - `iterate_response_state(...)`
  - `estimate_state_energies(...)`
  - `sync_trial_space_from_response_bundle(...)`
  - `write_restart_snapshot(...)`
  - `write_guess_archive(...)`
- Algorithmic role: scaffold equivalent of legacy `iterate`.
- Notes / assumptions:
  - convergence target is `max(10 * threshold, 1e-7)`
  - first iteration damping is `0.85`, later iterations use `0.65`
  - `result.failed = !converged`, so unconverged protocols are classified as failures

### `read_restart_snapshot` / `write_restart_snapshot`

- Location: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 267-324)
- Purpose: persist and restore the protocol-level excited restart state.
- Inputs:
  - read: `World`, path
  - write: `World`, path, `RestartSnapshot`
- Outputs:
  - read: `RestartSnapshot`
  - write: none
- Side effects:
  - reads/writes JSON snapshot file
  - reads/writes companion trial-state archive `<snapshot>.trial_states`
- Called by:
  - `solve_protocol(...)`
  - `select_restart_seed(...)`
  - `write_guess_archive(...)`
  - `iterate(...)`
- Calls:
  - `read_restart_trial_states(...)`
  - `write_restart_trial_states(...)`
- Algorithmic role: excited restart persistence.
- Notes / assumptions:
  - snapshot JSON stores energies, names, residuals, and iteration history
  - function-valued trial data lives in the companion archive, not in JSON

### `ResponseRecord2::initialize_excited_bundle`

- Location: `src/apps/molresponse_v2/ResponseRecord.hpp` (lines 84-108)
- Purpose: seed the metadata tree for the excited stage before execution.
- Inputs: plan fields and protocol list.
- Outputs: none.
- Side effects:
  - mutates `data_["excited_states"]`
  - writes the JSON file immediately
- Called by:
  - `JsonStateSolvePersistence::initialize_excited_bundle(...)`
  - then serial/subgroup linear stage setup
- Calls:
  - `ensure_excited_root(...)`
  - `ensure_excited_protocol(...)`
  - `write()`
- Algorithmic role: metadata schema initialization.
- Notes / assumptions:
  - this stores only plan data and minimal per-protocol placeholders
  - the richer stage-2c protocol fields are added elsewhere

## 7. Save/Load and Archive Structure

### 7.1 Files written by the current excited stage

The current excited implementation writes three distinct artifact types:

1. plan and status metadata inside `*.calc_info.json` response task metadata
2. optional placeholder metadata in `response_metadata.json`
3. excited restart snapshot files under the response task work directory

Driver placement:

- `ResponseApplication` runs in `task_<n>/molresponse`, so excited files live there. Reference: `src/madness/chem/Applications.hpp :: ResponseApplication::run` (lines 252-265).

Snapshot naming:

- per protocol: `<prefix>.excited_bundle.<protocol_key>.restartdata`
- companion trial-state archive: `<prefix>.excited_bundle.<protocol_key>.restartdata.trial_states`
- generic carryover guess: `<prefix>.excited_bundle.guess.restartdata`

References:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: protocol_restart_file` (lines 2506-2510)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: guess_archive_file` (lines 435-439)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: trial_space_archive_path` (lines 199-200)

The restart test explicitly expects these snapshot files under `task_1/molresponse`. Reference: `src/apps/madqc_v2/test_molresponse_excited_restart_reuse.py` (lines 72-79).

### 7.2 What the restart snapshot stores

`RestartSnapshot` stores:

- `converged`
- `iterations`
- `energies`
- `state_names`
- `residual_norms`
- `iteration_max_residuals`
- `trial_space_k`
- `trial_space_num_orbitals`
- `trial_states`

Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: RestartSnapshot` (lines 185-197).

Important fact: the stored function data is only `trial_states`, which correspond to `x`-space bundle functions. There is no archive of the full `ResponseVector` bundle, and no persisted `y` channels. References:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: write_restart_trial_states` (lines 247-265)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: sync_trial_space_from_response_bundle` (lines 1958-1964)

### 7.3 Metadata schema

The canonical final metadata shape exposed to users is:

```json
{
  "excited_states": {
    "plan": {
      "enabled": true,
      "num_states": ...,
      "tda": ...,
      "guess_max_iter": ...,
      "maxiter": ...,
      "maxsub": ...,
      "owner_group": ...,
      "protocols": ["1e-02", ...]
    },
    "protocols": {
      "1e-02": {
        "saved": true,
        "converged": false,
        "timings": {...},
        "stage_status": "...",
        "owner_group": 0,
        "iterations": ...,
        "energies": [...],
        "state_names": [...],
        "roots": [{"root_index": 0, "name": "es1", "energy": ...}, ...],
        "residual_norms": [...],
        "iteration_max_residuals": [...]
      }
    }
  },
  "excited_state_planner": {
    "plan": {...},
    "execution": {
      "solver_adapter": "...",
      "protocol_events": [...]
    }
  }
}
```

References:

- `src/madness/chem/MolresponseLib.hpp :: ensure_excited_protocol_placeholder_node` (lines 2716-2766)
- `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage` (lines 2798-3092)
- `src/madness/chem/MolresponseLib.hpp :: build_state_stage_metadata` (lines 3502-3508)

### 7.4 Restart behavior

The excited restart mechanism is snapshot-driven, not metadata-driven.

Fact:

- `execute_excited_state_bundle_stage(...)` can skip a protocol if metadata already says `saved && converged`. Reference: `src/madness/chem/MolresponseLib.hpp` (lines 2908-2915).
- the scaffold itself reuses restart snapshots by reading the protocol restart files and generic guess archive. Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: select_restart_seed` (lines 2212-2254).

Important current limitation:

- after the excited stage mutates `state_metadata_json`, the code does not write that updated excited subtree back to `response_metadata.json`; the final excited metadata is returned in `calc_info.json`, but the on-disk metadata file written during the linear stage is not refreshed by Stage 2c. References:
  - serial path captures a copy at `state_metadata_json = persistence.metadata_json();` in `src/madness/chem/MolresponseLib.hpp` line 2138
  - subgroup path writes merged `response_metadata.json` before excited execution in `src/madness/chem/MolresponseLib.hpp` lines 2626-2628
  - Stage 2c mutates only the in-memory `state_metadata_json` in `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage` (lines 2872-3090)

`Inference:` because of that write pattern, ordinary reruns likely reuse excited restart snapshots through `select_restart_seed(...)` rather than through the metadata-only `restart_ready_skip` path.

### 7.5 Naming model

Current naming is metadata-only:

- `state_names` is stored in snapshots and metadata
- `roots` is derived from `state_names` and `energies`
- no archive filename is keyed by state name
- no helper loads an excited root by state name

References:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: assign_excited_state_names` (lines 110-140)
- `src/madness/chem/MolresponseLib.hpp :: build_excited_root_manifest` (lines 2768-2781)
- `rg` over `src/` shows `roots` is only written/merged, not consumed elsewhere

## 8. Implemented vs Missing Features

### Implemented

- Parameter plumbing from input to plan metadata is implemented. References: `src/madness/chem/ResponseParameters.hpp` (lines 82-96, 184-201), `src/madness/chem/MolresponseLib.hpp :: build_excited_state_bundle_plan` (lines 977-997).
- Stage-2 workflow insertion is implemented. The excited bundle stage runs between linear solves and derived-state execution. Reference: `src/madness/chem/MolresponseLib.hpp :: solve_all_states` (lines 3547-3559).
- Protocol-level restart snapshots are implemented. References: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: read_restart_snapshot` / `write_restart_snapshot` (lines 267-324).
- Trial-space reprojection across protocol changes is implemented and regression-tested. References: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: align_trial_space_protocol` (lines 2096-2118), `src/apps/madqc_v2/test_molresponse_excited_protocol_projection.py` (lines 69-118).
- Protocol-local root naming and root manifests are implemented. References: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: assign_excited_state_names` (lines 110-140), `src/madness/chem/MolresponseLib.hpp :: build_excited_root_manifest` (lines 2768-2781).
- Restricted bundle rotation via generalized diagonalization is implemented. References: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: diagonalize_excited_bundle` (lines 988-1140), `rotate_excited_bundle_states` (lines 1143-1175).
- Metadata and summary tests are implemented. References: `src/apps/madqc_v2/test_molresponse_excited_metadata_smoke.py` (lines 84-173), `src/apps/madqc_v2/test_molresponse_excited_restart_reuse.py` (lines 97-133).

### Partial

- The numerical solver is executable but still scaffold-like. It performs bundle potential formation, restricted rotation, BSH-based updates, damping, and convergence checks, but comments and local checklist still frame it as incomplete legacy reintegration. References: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` comment at lines 333-335, `src/apps/molresponse_v2/EXCITED_STATE_EXECUTION_CHECKLIST.md` lines 63-98.
- Root solving is bundle-coupled only for restricted variants. Unrestricted variants fall back to per-state iteration and explicitly skip bundle rotation as `unrestricted_wip`. References: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: rotate_response_bundle` (lines 1195-1225), `iterate` (lines 2360-2400).
- Dynamic full TDHF/TDDFT channels are only partially restart-persisted. The live `ResponseVector` may contain `y` channels, but restart and cross-protocol state transfer keep only `x_alpha`. References: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: sync_trial_space_from_response_bundle` (lines 1958-1964), `make_state_response_vector` (lines 1898-1935).
- Excited metadata support exists in `ResponseRecord2`, but the current stage writes most excited results directly into `state_metadata_json` instead of using the `record_excited_protocol_*` helpers. References: `src/apps/molresponse_v2/ResponseRecord.hpp` (lines 110-137), `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage` (lines 2997-3022).
- Protocol result semantics are only partly refined. `saved` means "snapshot written", while unconverged results are marked `failed`. Reference: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate` (lines 2486-2503).

### Missing

- No dedicated excited-state descriptor/runtime object exists; state identity is still mostly index-based. Reference: there is no excited runtime state class beyond `ExcitedTrialSpace` and metadata arrays in `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 326-331, 416-433).
- `owner_group` is not enforced by execution. It is stored in metadata and passed through the input struct, but there is no use of `owner_group` inside `ExcitedStateBundleSolver.cpp`, and Stage 2c runs on the full `World`. References: `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage` (lines 2957-2967), `rg` over `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` finds no `owner_group` use beyond the input struct definition.
- `maxsub` is not used anywhere in the excited numerical kernel. References: `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp` line 26, `src/madness/chem/MolresponseLib.hpp` lines 988, 2967, and no use sites in `ExcitedStateBundleSolver.cpp`.
- There is no root-by-root archive save/load path analogous to `ResponseIO.hpp`. Excited roots are persisted only through bundle snapshots. References: `src/apps/molresponse_v2/ResponseIO.hpp` (lines 13-130), and no `save_response_vector` call sites inside `ExcitedStateBundleSolver.cpp`.
- No property stage consumes excited roots. References: `src/madness/chem/MolresponseLib.hpp :: compute_polarizability` / `compute_hyperpolarizability` / `compute_raman` (lines 3603-3725) and no `roots` consumer in `src/`.
- The build option `MADNESS_ENABLE_LEGACY_EXCITED_BUNDLE_ADAPTER` exists, but there are no source-level `#ifdef` or macro use sites that switch the implementation. References: `src/apps/molresponse_v2/CMakeLists.txt` (lines 42-51), and `rg` shows no other macro references.

## 9. Progress Assessment

State representation - partial  
Reason: bundle state exists, but there is no dedicated excited-state object; identity is still index- and metadata-based.

Protocol planning - implemented  
Reason: `ExcitedStateBundlePlan` is fully wired from `ResponseParameters` into Stage 2c metadata and execution.

Solver entry path - implemented  
Reason: `madqc --wf=response` and `molresponse2` both reach `execute_excited_state_bundle_stage(...)`.

Initialization / guess generation - partial  
Reason: restart selection, fresh-guess creation, and protocol reprojection exist, but the guess model is heuristic and stores only `x` channels.

Restricted bundle iteration - partial  
Reason: bundle potentials, generalized rotation, and per-state updates run, but the code still identifies itself as scaffold and does not claim full legacy numerical parity.

Unrestricted excited-state iteration - partial  
Reason: unrestricted roots can be represented and iterated, but bundle rotation is explicitly skipped as work in progress.

Convergence logic - partial  
Reason: a concrete residual-based stop criterion exists, but it is not the same structured density/residual/KAIN workflow as the linear solver and is still called out as incomplete in the local checklist.

Restart snapshots - implemented  
Reason: per-protocol and generic guess snapshots are read/written and exercised by tests.

Metadata / manifest system - partial  
Reason: final metadata contains plan, status, timings, root names, and root manifests, but schema ownership is split between `ResponseRecord2` and direct JSON mutation, and Stage 2c does not flush final excited metadata back to `response_metadata.json`.

Naming system - partial  
Reason: protocol-local names and root manifests exist, but names are energy-order-based and not backed by a stable runtime identity or archive naming scheme.

Owner-group execution - missing  
Reason: `excited.owner_group` is planned and recorded, but not enforced.

Subspace control (`maxsub`) - missing  
Reason: the parameter is threaded through plan/input/metadata but unused in solver code.

Property integration - missing  
Reason: no property routine consumes excited roots or excited metadata.

Build-time adapter switching - missing  
Reason: the CMake option exists but currently does not select a different compiled implementation.

## 10. Migration-Relevant Observations

The new architecture has already separated excited-state execution from the old monolithic response base. The reintegration seam is now explicit: `execute_excited_state_bundle_stage(...)` only requires a protocol-wise adapter contract (`ExcitedStateBundleSolver`), so future work can replace the scaffold without rewriting the whole Stage 2 driver. References: `src/madness/chem/MolresponseLib.hpp` (lines 2784-3092), `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp` (lines 44-77).

The current bundle model is intentionally different from legacy `X_space`, but it is still bundle-centric. All roots share `trial_space_`, `omega_`, `state_names_`, residual arrays, and restricted bundle rotation logic. This means reintegration work can still assume bundle-level coupling, but not a single legacy-style state container. References: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 416-433, 1981-2035, 2310-2503).

Cross-protocol restart continuity currently depends on `x`-space trial states, not full response bundles. That is a major design constraint if the final reintegration needs true TDHF/TDDFT `Y`-channel continuity across restarts or protocol ladders. References: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 1898-1935, 1958-1964, 2478-2481).

Root identity is still protocol-local and energy-order-dependent. Because `state_names` are reassigned from sorted energies with protocol-dependent degeneracy grouping, names are not yet a stable cross-protocol state identifier in the architectural sense. References: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: assign_excited_state_names` (lines 110-140), `iterate` (line 2458), `src/madness/chem/MolresponseLib.hpp :: build_excited_root_manifest` (lines 2768-2781).

The excited stage is better integrated into final task metadata than into the persistent `response_metadata.json` file. That split matters for restart, tooling, and downstream consumers. Today:

- user-facing final metadata in `*.calc_info.json` contains the rich excited subtree
- restart snapshots contain the actual reusable excited bundle data
- `response_metadata.json` is only partially aware of final excited results

References: `src/madness/chem/Applications.hpp` (lines 271-275), `src/madness/chem/MolresponseLib.hpp` (lines 2138, 2626-2628, 2872-3090), `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp` (lines 267-324, 2482-2484).

Relative to the legacy implementation, the strongest pieces already reintegrated are orchestration, restart plumbing, protocol metadata, and restricted bundle rotation. The weakest pieces are stable state identity, owner-group enforcement, full-bundle persistence, and downstream property consumption.

## 11. Code Reference Index

- `src/apps/madqc_v2/madqc.cpp :: main` - response workflow CLI entry (lines 140-183)
- `src/apps/molresponse_v2/molresponse2.cpp :: main` - compatibility entry (lines 29-84)
- `src/madness/chem/WorkflowBuilders.hpp :: add_response_workflow_drivers` - workflow registration (lines 75-84)
- `src/madness/chem/Applications.hpp :: ResponseApplication::run` - response task execution directory and `run_response` call (lines 252-265)
- `src/madness/chem/MolresponseLib.hpp :: run_response` - top-level response workflow (lines 4302-4354)
- `src/madness/chem/MolresponseLib.hpp :: build_excited_state_bundle_plan` - plan construction (lines 977-997)
- `src/madness/chem/MolresponseLib.hpp :: plan_required_states` - stage-1 planning (lines 999-1049)
- `src/madness/chem/MolresponseLib.hpp :: solve_all_states` - stage ordering (lines 3520-3568)
- `src/madness/chem/MolresponseLib.hpp :: execute_excited_state_bundle_stage` - stage-2c orchestration (lines 2784-3092)
- `src/madness/chem/MolresponseLib.hpp :: build_excited_root_manifest` - root manifest generation (lines 2768-2781)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp :: ExcitedBundleProtocolInput` - protocol solve input contract (lines 16-27)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp :: ExcitedBundleProtocolResult` - protocol solve result contract (lines 29-42)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: RestartSnapshot` - restart snapshot content (lines 185-197)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: ExcitedTrialSpace` - bundle trial-space container (lines 326-331)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: ExcitedProtocolWorkflow::solve_protocol` - protocol entry (lines 356-399)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: prepare_protocol` - protocol `k`/threshold setup (lines 470-481)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: ensure_ground_data` - ground/preliminary cache management (lines 491-519)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: estimate_state_energies` - heuristic energy estimator (lines 578-597)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: diagonalize_excited_bundle` - generalized bundle diagonalization (lines 988-1140)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: rotate_excited_bundle_states` - restricted bundle rotation (lines 1143-1175)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_trial` - guess refinement loop (lines 1867-1896)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: make_state_response_vector` - trial-to-response conversion (lines 1898-1935)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: sync_trial_space_from_response_bundle` - response-to-trial backprojection (lines 1958-1964)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: initialize_protocol_guess` - main initialization path (lines 2256-2307)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate` - main excited iteration loop (lines 2310-2503)
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: read_restart_snapshot` / `write_restart_snapshot` - snapshot I/O (lines 267-324)
- `src/apps/molresponse_v2/ResponseRecord.hpp :: initialize_excited_bundle` - excited metadata seeding (lines 84-108)
- `src/apps/molresponse_v2/ResponseVector.hpp` - runtime response-vector variants reused by excited roots (lines 16-230)
- `src/apps/molresponse_v2/ResponseSolverUtils.hpp :: align_state_bundle_protocol` - trial-space reprojection utility (lines 120-129)
- `src/apps/molresponse_v2/FrequencyLoop.cpp :: computeFrequencyLoop` / `computeFrequencyPoint` - surrounding linear protocol/frequency workflow (lines 566-650)
- `src/apps/madqc_v2/test_molresponse_excited_metadata_smoke.py` - metadata regression coverage (lines 55-173)
- `src/apps/madqc_v2/test_molresponse_excited_restart_reuse.py` - restart reuse regression coverage (lines 56-133)
- `src/apps/madqc_v2/test_molresponse_excited_protocol_projection.py` - protocol reprojection regression coverage (lines 55-118)
