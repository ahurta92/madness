# Molresponse v2 Function Reference

This document is a code-oriented map of the `molresponse_v2` stack.
It answers four questions for each major piece:

1. What does this function/type do?
2. Why does it exist?
3. Where does it fit in the end-to-end response pipeline?
4. What other functions are similar or adjacent?

Scope notes:
- This covers the new response pipeline centered on `src/apps/molresponse_v2/*` and `src/madness/chem/MolresponseLib.hpp`.
- "Major" functions are listed explicitly; trivial getters/setters and formatting helpers are grouped where appropriate.

## Complete File Map (New Molresponse Pieces)

This is the full "what lives where" index for the new stack.
Use it before diving into per-function details.

| File | Primary role in the pipeline | Closest sibling(s) |
|---|---|---|
| `molresponse2.cpp` | Compatibility CLI entrypoint, delegates to workflow driver. | `WorkflowBuilders.hpp` response driver wiring. |
| `StateGenerator.hpp` | Expands requested properties into linear state descriptors. | `DerivedStatePlanner.hpp` (derived-state expansion). |
| `StateParallelPlanner.hpp` | Computes static ownership plan for channels/points/subgroups. | Runtime schedule helpers in `MolresponseLib.hpp`. |
| `DerivedStatePlanner.hpp` | Builds and gates VBC-driven derived requests. | `PropertyManager.hpp` beta/raman triplet expansion. |
| `ResponseState.hpp/.cpp` | Core perturbation/state types and perturbation-vector assembly. | `Perturbation.hpp`, `ResponseVector.hpp`. |
| `ResponseVector.hpp` | Variant container for static/dynamic and restricted/unrestricted response vectors. | `FrequencyLoop.cpp` variant solve dispatch. |
| `GroundStateData.hpp/.cpp` | Ground-state checkpoint hydration, protocol preparation, preliminaries. | `ResponseManager.hpp/.cpp` protocol operators. |
| `ResponseManager.hpp/.cpp` | Protocol-level operator and parameter management. | `GroundStateData::prepareOrbitals(...)`. |
| `ResponseInitializer.hpp/.cpp` | Fallback initial guess construction for response points. | `seed_response_guess_for_point(...)` in `FrequencyLoop.cpp`. |
| `ResponseIO.hpp` | Archive read/write for linear and derived response vectors. | `ResponseRecord.hpp` metadata persistence. |
| `FrequencyLoop.hpp/.cpp` | Frequency-point solve loop, restart provenance, diagnostics recording. | `ResponseSolver.cpp` nonlinear update kernels. |
| `ResponseSolver.hpp/.cpp` | Coupled-response equations and variant-specific update kernels. | `ResponseSolverUtils.hpp`. |
| `ResponseSolverUtils.hpp` | Shared low-level solver helpers (inner products, BSH builders, logging helpers). | `ResponseSolver.cpp`. |
| `ResponseRecord.hpp` | Run metadata/status/timing/diagnostic persistence. | `ResponseDebugLogger.hpp` iteration logs. |
| `ResponseDebugLogger.hpp` | Iteration-level diagnostics and timing logs. | `ResponseDebugLoggerMacros.hpp`. |
| `ResponseDebugLoggerMacros.hpp` | Macro wrappers for low-overhead timed value/block logging. | direct `ResponseDebugLogger` calls. |
| `VBCMacrotask.hpp` | VBC state construction/persistence used by derived and property stages. | `DerivedStatePlanner::make_vbc_state(...)`. |
| `ExcitedStateBundleSolver.hpp/.cpp` | Excited-state stage adapter and restart scaffold. | Excited-stage orchestration in `MolresponseLib.hpp`. |
| `PropertyManager.hpp` | Alpha/beta/raman assembly and component-level dedupe/claims. | Stage-3 wrappers in `MolresponseLib.hpp`. |
| `InnerContributions.hpp` | Optional contraction contribution diagnostics sink. | `compute_response_inner_product_tensor(...)`. |
| `Perturbation.hpp` | Typed perturbation variant and canonical descriptors. | `ResponseState.cpp` perturbation operator builders. |
| `MolecularProperty.hpp` | Legacy/auxiliary property request descriptors/printers. | `ResponseParameters` + `StateGenerator` pipeline. |
| `broadcast_json.hpp` | MPI broadcast utility for JSON file payloads. | local JSON read/write helpers in `MolresponseLib.hpp`. |
| `MolresponseLib.hpp` | End-to-end workflow orchestrator for response calculation. | `run_response(...)` consumers in workflow drivers. |
| `ResponseParameters.hpp` | Response input schema, derived defaults, and validation. | `StateGenerator`, `DerivedStatePlanner`, `PropertyManager`. |
| `test_preliminaries.cpp` | Preliminaries regression check across protocols. | `GroundStateData`, `ResponseManager`. |
| `test_parameter_manager.cpp` | Parameter-manager/response parsing smoke test. | `ResponseParameters`. |
| `validate_property_component_claims.py` | Offline validator for subgroup component-claim behavior. | `try_claim_property_component_task(...)`. |

## Pipeline Map

1. **Entry and setup**
- `molresponse2` (compat wrapper) or `madqc --wf=response` enters workflow dispatch.
- `molresponse_lib::run_response(...)` orchestrates all stages.

2. **Stage 1: plan required work**
- Build ground context.
- Generate linear states from requested properties.
- Build state-parallel ownership plan.
- Build derived (VBC-driven) request plan.
- Build excited-state bundle plan metadata.

3. **Stage 2: solve states**
- Solve linear states over protocols/frequencies (serial or subgroup).
- Record metadata, solver diagnostics, restart provenance.
- Execute excited-state bundle stage adapter.
- Execute dependency-gated derived requests.

4. **Stage 3: assemble properties**
- Compute alpha/beta/raman from solved archives.
- Optional subgroup property-component precompute and merge.
- Return final metadata/properties/debug payload.

## Entry Points

### `src/apps/molresponse_v2/molresponse2.cpp`

| Function | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `print_help()` | Prints CLI usage for `molresponse2`. | Keeps compatibility executable self-describing. | Entry | `madqc --help`, workflow docs. |
| `main(int argc, char** argv)` | Initializes world, enforces `workflow=response`, constructs `Params`, dispatches response workflow builders, runs workflow. | Backward-compatible entry path while central logic lives elsewhere. | Entry | `molresponse_lib::run_response(...)` (actual response pipeline body). |

## Top-Level Orchestrator

### `src/madness/chem/MolresponseLib.hpp` (`struct molresponse_lib`)

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `run_response(...)` | Main response workflow: build context, plan, solve, compute properties, return `Results`. | Single canonical response execution API used by workflow driver. | All | `main(...)` in wrapper only dispatches here indirectly. |
| `GroundContext` | Carries molecule, ground data, response manager, archive/fock paths. | Centralizes shared state for all stages. | Stage 1+ | `PropertyContext` (stage-3 view). |
| `PlannedStates` | Bundles generated linear states + derived plan + excited plan + state-parallel plan. | Stable handoff between planning and solving. | Stage 1->2 | `SolvedStates` (post-solve handoff). |
| `SolvedStates` | Planned state context + merged metadata + merged debug log. | Feeds property stage and output serialization. | Stage 2->3 | `PropertyStageOutput`. |
| `PropertyStageOutput` | Properties JSON + vibrational + Raman artifacts. | Structured stage-3 output payload. | Stage 3 | Final `Results`. |
| `make_ground_context(...)` | Reads checkpoint/molecule context, initializes `GroundStateData` + `ResponseManager`. | Standardized ground bootstrap with consistent paths and metadata. | Stage 1 | `GroundStateData` constructors/load. |
| `build_excited_state_bundle_plan(...)` | Converts response knobs into excited-state per-protocol plan metadata. | Keeps excited stage parameterization explicit and restart-friendly. | Stage 1 | `execute_excited_state_bundle_stage(...)`. |
| `plan_required_states(...)` | Generates linear states, derived requests, state-parallel plan, excited plan; prints diagnostics. | Single planning convergence point used by downstream scheduler. | Stage 1 | `StateGenerator::generateStates`, `StateParallelPlanner::build`, `DerivedStatePlanner::build_vbc_driven_quadratic_plan`. |
| `build_state_solve_schedule_context(...)` | Combines planner output + restart metadata into runtime solve policy. | Decouples static planning from runtime restart-aware execution policy. | Stage 2 | `compute_runtime_point_ownership_policy`, `build_protocol_execution_policy`. |
| `execute_serial_state_solve(...)` | Runs protocol loops in one world communicator; solves pending linear points/channels; writes metadata/debug. | Deterministic fallback and baseline execution mode. | Stage 2 | `execute_subgroup_state_solve(...)`. |
| `execute_subgroup_state_solve(...)` | Runs state/point ownership in subgroup worlds, writes shard metadata/logs, merges shards. | Enables state-level and point-level parallelism. | Stage 2 | `execute_serial_state_solve(...)`. |
| `execute_excited_state_bundle_stage(...)` | Runs protocol-indexed excited-state adapter and records protocol status/timing/energies. | Reintegration seam for excited-state workflow with restart tracking now. | Stage 2c | `make_excited_state_bundle_solver_adapter(...)`. |
| `execute_derived_state_requests(...)` | Evaluates dependency gate for derived requests and executes ready requests (subgroup or serial deterministic lanes). | Prevents invalid derived solves and preserves progress reporting. | Stage 2d | `DerivedStatePlanner::evaluate_dependency_gate`, `run_derived_request(...)`. |
| `prepare_and_validate_final_protocol_state(...)` | Sets final protocol context and assesses full linear convergence status. | Ensures property stage runs at final protocol and records convergence state. | Stage 2->3 | `point_ready_in_metadata(...)`. |
| `build_state_stage_metadata(...)` | Assembles stage metadata (`state_parallel_runtime`, `derived_state_planner`, `excited_state_planner`). | Makes scheduler/derived/excited decisions inspectable and restartable. | Stage 2 finalization | Metadata merge helpers. |
| `compute_requested_properties(...)` | Serial property assembly loop for requested properties. | Baseline property execution path. | Stage 3 | `compute_requested_properties_with_property_group(...)`. |
| `compute_requested_properties_with_property_group(...)` | Subgroup-aware property component precompute + fallback + merge + property-group assembly. | Scales property component generation without giving up serial fallback safety. | Stage 3 | `compute_requested_properties(...)`. |
| `compute_polarizability(...)` | Stage-3 wrapper invoking `compute_alpha`. | Keeps property dispatch in orchestrator thin and explicit. | Stage 3 | `compute_hyperpolarizability(...)`, `compute_raman(...)`. |
| `compute_hyperpolarizability(...)` | Stage-3 wrapper invoking beta computation with triplet-mode params. | Routes response params (`beta.shg/beta.or/beta.all_triplets`) to property code. | Stage 3 | `::compute_hyperpolarizability(...)` in `PropertyManager.hpp`. |
| `compute_raman(...)` | Stage-3 wrapper invoking Hessian + Raman assembly and reporting. | Keeps Raman-specific flow in one path with shared context. | Stage 3 | `compute_hessian(...)`, `compute_Raman(...)`. |
| `parse_property_name(...)` | Converts string property names to enum dispatch. | Avoids stringly-typed control flow in property stage. | Stage 3 | `ResponseParameters::requested_properties()`. |
| `JsonStateSolvePersistence` | Adapter implementing `StateSolvePersistence` over `ResponseRecord2` + `ResponseDebugLogger`. | Decouples frequency solver from concrete metadata/log format. | Stage 2 infra | `StateSolvePersistence` interface in `FrequencyLoop.hpp`. |

### Key orchestration helpers in `MolresponseLib.hpp`

| Helper | What it does | Why it exists | Similar / related |
|---|---|---|---|
| `group_shard_file`, `group_console_file`, `group_derived_timing_file` | Derive deterministic shard filenames per subgroup id. | Stable merge and debugging artifacts. | Property/metadata shard merges. |
| `write_json_file`, `read_json_file_or_object` | I/O wrappers for JSON files. | Uniform error behavior for orchestration layer. | `broadcast_json_file(...)`. |
| `merge_state_metadata_json`, `merge_debug_log_json` | Merge subgroup shard metadata/debug logs into canonical JSON. | Needed for subgroup solve mode correctness. | Shard-based property merge. |
| `point_ready_in_metadata`, `point_needs_solving_from_metadata` | Readiness predicates against metadata. | Restart-aware scheduling decisions. | `point_needs_solving(...)` path using persistence. |
| `run_protocol_threshold_loop(...)` | Generic protocol-loop skeleton with callbacks (skip/setup/solve/finalize). | Removes duplicated protocol control flow across serial/subgroup paths. | `dispatch_owned_work_for_protocol(...)`. |
| `build_pending_work_manifest(...)` | Build pending owned channels/points for a protocol and lane set. | Central worklist logic for both serial preview and subgroup execution. | `build_pending_manifest_from_metadata(...)`. |
| `prepare_protocol_context(...)` | Applies protocol threshold to response manager + ground data and recomputes preliminaries. | Keeps protocol setup logic identical across serial/subgroup stage-2 paths. | protocol-setup lambdas in both execution paths. |
| `log_pending_manifest(...)` | Prints a standardized summary for pending channel/point work manifests. | Avoids duplicated logging format logic for serial and subgroup execution. | protocol policy diagnostics, subgroup pending-work logging. |
| `execute_manifest_work(...)` | Dispatches a `PendingProtocolManifest` in channel-series or channel-point mode. | Removes duplicated channel/point branching loops. | `run_frequency_loop_with_flush(...)`, `computeFrequencyPoint(...)`. |

## Response Parameters

### `src/madness/chem/ResponseParameters.hpp`

| Function / Group | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `ResponseParameters()` | Declares all response inputs with defaults and constraints. | Canonical parameter schema for response workflow. | Stage 1 input | `CalculationParameters` for general MRA/SCF settings. |
| `ResponseParameters(World&, parser)` | Reads input + CLI, derives property list, validates consistency. | One-step parse/normalize/validate entry. | Stage 1 input | `set_derived_properties()`, `validate_user_specified_properties()`. |
| `set_derived_properties()` | Fills `requested_properties` from legacy knobs when not explicitly set. | Backward compatibility for older input style. | Stage 1 input | Explicit `requested_properties` path bypasses this. |
| `validate_user_specified_properties()` | Ensures required dipole/nuclear blocks are present for explicit properties. | Fail-fast input validation. | Stage 1 input | Property dispatch in `MolresponseLib`. |
| `beta_shg()`, `beta_or()`, `beta_all_triplets()` | Beta triplet mode controls. | Limits beta frequency-pair expansion to desired physical process. | Stage 1->3 | `select_beta_frequency_triplets(...)`, `DerivedStatePlanner::select_beta_frequency_pairs(...)`. |
| Accessor family (`dipole_*`, `nuclear_*`, `state_parallel_*`, `excited_*`, etc.) | Typed reads of configured knobs. | Keeps business logic free of string keys. | All | Used heavily in planner and property dispatch. |

## State Generation and Planning

### `src/apps/molresponse_v2/StateGenerator.hpp`

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `GeneratedStateData` | Holds `states` vector + `state_map` by perturbation key. | Dual access pattern (ordered iteration + keyed lookup). | Stage 1 output | Used by solver and property loaders. |
| `GeneratedStateData::print_generated_state_map(...)` | Prints generated linear-state inventory. | Quick operator-facing sanity check. | Stage 1 diagnostics | Metadata JSON in `response_metadata.json`. |
| `StateGenerator::generateStates()` | Expands requested properties into merged, deduped linear perturbation descriptors and frequencies. | Converts high-level property intent into concrete linear solve targets. | Stage 1 | `DerivedStatePlanner::build_vbc_driven_quadratic_plan(...)` handles derived second-order requests. |

### `src/apps/molresponse_v2/StateParallelPlanner.hpp`

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `PerturbationChannelAssignment` | `(channel_index, label, owner_group)` mapping row. | Deterministic owner mapping and inspectable metadata. | Stage 1/2 | `owner_by_channel_index` vector in orchestration. |
| `StateParallelPlan` | Full planner decision payload (mode, groups, policy, assignments). | One struct captures all scheduler inputs/outputs. | Stage 1/2 | Runtime policy augmentation in `MolresponseLib`. |
| `StateParallelPlan::point_parallel_start_protocol_index` | Configures first protocol index where channel-point ownership may begin on fresh runs. | Keeps planner intent explicit while runtime policy can switch to final-only restart mode. | Stage 2 | `compute_runtime_point_ownership_policy(...)` in `MolresponseLib.hpp`. |
| `PointOwnershipScheduler` | Computes point owner by global linear point index modulo owner groups. | Deterministic frequency-point ownership in point mode. | Stage 2 | Channel ownership mapping in `StateParallelPlan`. |
| `StateParallelPlanner::build(...)` | Computes effective mode/groups/assignments from user knobs and state count. | Separates static planning from runtime execution logic. | Stage 1 | `build_state_solve_schedule_context(...)` (runtime layer). |

### `src/apps/molresponse_v2/DerivedStatePlanner.hpp`

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `DerivedStateRequest` | Encodes one derived request: target A, VBC source, frequencies, dependencies. | Stable unit for gating, execution, and metadata. | Stage 1/2 | `VBCResponseState` reconstruction. |
| `DerivedStatePlan` | List of derived requests. | Explicit derived workload handoff. | Stage 1->2 | `GeneratedStateData` for linear states. |
| `DerivedStateGateReport` | Per-request readiness and missing dependency report. | Transparent dependency gating before derived execution. | Stage 2d | Metadata under `derived_state_planner.dependency_gate`. |
| `build_vbc_driven_quadratic_plan(...)` | Builds deduped VBC-driven requests for beta/raman from response settings. | Converts property requirements into concrete derived requests. | Stage 1 | `StateGenerator::generateStates()` (linear counterpart). |
| `evaluate_dependency_gate(...)` | Checks if all required linear dependencies are ready per request. | Prevents invalid derived solves and makes blockers explicit. | Stage 2d | `point_ready_in_metadata(...)` and persistence checks. |
| `make_vbc_state(...)` | Rebuilds `VBCResponseState` from request metadata. | Execution bridge from planned metadata to runtime objects. | Stage 2d execution | `VBCMacrotask::compute_and_save(...)`. |
| `select_beta_frequency_pairs(...)` | Applies beta mode precedence (`all_triplets` > `or` > `shg`) for derived beta request expansion. | Ensures derived VBC planning matches property-stage beta triplet policy. | Stage 1 | `select_beta_frequency_triplets(...)` in property code. |
| `add_request(...)` | Internal helper to canonicalize frequencies, dedupe by request id, append dependencies. | Ensures stable request identity and no duplicate execution. | Stage 1 | `make_request_id(...)`. |

## State and Perturbation Model

### `src/apps/molresponse_v2/ResponseVector.hpp`

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `StaticRestrictedResponse`, `DynamicRestrictedResponse`, `StaticUnrestrictedResponse`, `DynamicUnrestrictedResponse` | Concrete response container layouts (x/y, alpha/beta) with `flat` + `sync/flatten`. | Supports solver operations over unified flat vectors while retaining semantic components. | Stage 2/3 | `ResponseVector` variant dispatch. |
| `ResponseVector` | `std::variant` over the four response container forms. | Type-safe runtime polymorphism for static/dynamic and restricted/unrestricted cases. | Stage 2/3 | `solve_response_vector(...)` visitation. |
| `make_response_vector(...)` | Factory selecting variant arm by static/unrestricted flags. | Removes repetitive variant construction logic. | Stage 2 I/O/initialization | `initialize_guess_vector(...)`. |
| `get_flat(...)` | Returns flattened function vector from any variant arm. | Uniform downstream tensor-contraction API. | Stage 2/3 | `sync()` on concrete arm restores structured fields. |

### `src/apps/molresponse_v2/ResponseState.hpp` + `.cpp`

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `canonicalize_response_frequency(...)` | Rounds frequencies to 1e-3 grid and normalizes near-zero. | Prevents duplicate scheduling/filenames from floating-point noise. | Stage 1/2/3 | Used in state generation and request IDs. |
| `LinearResponseDescriptor` | Describes one perturbation channel across frequencies/protocol thresholds. | Core linear workload descriptor for planning and solves. | Stage 1/2/3 | `LinearResponsePoint` selects one point from this channel. |
| `LinearResponsePoint` | Lightweight `(descriptor, protocol_index, freq_index)` view. | Point-level solve and metadata keying. | Stage 2 | `ChannelPoint` alias in planner terminology. |
| `SecondOrderResponseDescriptor` + `VBCResponseState`/`XBCResponseState` | Encodes second-order source states and naming for derived responses. | Carries derived response identity and associated B/C dependencies. | Stage 2d/3 | `DerivedStateRequest` metadata equivalent. |
| `make_perturbation_operator(...)` overloads | Builds raw real-space perturbation operators for dipole/nuclear/magnetic. | Canonical perturbation operator generation. | Stage 2/3 | `raw_perturbation_operator(...)` dispatches to these. |
| `raw_perturbation_operator(...)` | Variant dispatch to operator builder by perturbation type. | One entrypoint for all perturbation kinds. | Stage 2/3 | `perturbation_vector(...)` builds orbital-projected vectors from it. |
| `project_perturbation_onto_orbitals(...)` | Applies raw operator to ground orbitals, projects with `Qhat`, truncates. | Shared projection path for linear/property code. | Stage 2/3 | Used by both `perturbation_vector` overloads. |
| `perturbation_vector(...)` overloads | Produces solver perturbation RHS vectors for point or descriptor contexts. | Avoids duplicating static/dynamic duplication logic across call sites. | Stage 2/3 | Point overload (solver), descriptor overload (property stage). |

## Ground and Operator Context

### `src/apps/molresponse_v2/GroundStateData.hpp` + `.cpp`

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `GroundStateData(World&, archive, Molecule)` | Constructs from restart archive + molecule. | Main path for response workflow seeded from moldft checkpoint. | Stage 1 | `GroundStateData(World&, shared_ptr<SCF>)` for SCF-coupled contexts. |
| `load(...)` | Reads orbitals, energies, occupancy, and settings from archive and normalizes wavelet state. | Archive-to-runtime state hydration. | Stage 1 / protocol reload | `prepareOrbitals(...)` may trigger reload. |
| `prepareOrbitals(...)` | Reload/project/truncate orbitals to requested protocol `(k, thresh)`; rebuilds `Qhat`. | Ensures orbital basis matches protocol requirements. | Stage 2 protocol setup | `ResponseManager::setProtocol(...)`. |
| `computePreliminaries(...)` | Builds density, potentials, Hamiltonian (load from fock JSON or compute), `Hamiltonian_no_diag`. | Precomputes quantities needed repeatedly by coupled-response equations. | Stage 2/3 setup | `computeDensity`, `computeXCPotential`, `computeHFExchangeEnergy`, `computeKineticEnergy`. |
| `tryLoadHamiltonianFromJson(...)` | Attempts fast Hamiltonian restore from cached fock JSON protocol key. | Avoids recomputation where possible. | Stage 2/3 setup | fallback compute branch in `computePreliminaries(...)`. |
| Potential/energy helpers (`computeDensity`, `computeNuclearPotential`, `computeCoulombPotential`, `computeXCPotential`, `computeHFExchangeEnergy`, `computeKineticEnergy`) | Piecewise builders for preliminaries. | Isolate expensive physics kernels and simplify `computePreliminaries`. | Stage 2/3 setup | called only by `computePreliminaries(...)`. |

### `src/apps/molresponse_v2/ResponseManager.hpp` + `.cpp`

| Function | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `ResponseManager(World&, CalculationParameters)` | Stores calc params and response operator state. | Shared operator/config holder for solver stages. | Stage 1+ | `GroundStateData` holds ground-side counterparts. |
| `setProtocol(World&, L, thresh, override_k)` | Sets function defaults (k/thresh/box/refine), builds Coulomb and gradient operators, updates `vtol`. | Central protocol reconfiguration hook used before each threshold solve block. | Stage 2/3 setup | `GroundStateData::prepareOrbitals(...)` (orbital side of same protocol switch). |
| Getters (`getVtol`, `getCoulombOp`, `getGradOp`, `params`) | Expose protocol-configured operators and parameters. | Keeps callsites decoupled from internals. | Stage 2/3 | used throughout solver/property paths. |

## Linear Solver Flow

### `src/apps/molresponse_v2/ResponseInitializer.hpp` + `.cpp`

| Function | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `initialize_guess_vector(...)` | Creates a zero-compressed initial response vector in correct variant shape for point static/dynamic and spin mode. | Standard fallback initialization when no restart/continuation guess is available. | Stage 2 point solve init | `seed_response_guess_for_point(...)` calls this when archive/memory seeds unavailable. |

### `src/apps/molresponse_v2/ResponseIO.hpp`

| Function | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `save_response_vector(...)` templates | Writes response vector archive for descriptor/point. | Persist restart states and feed downstream property/derived stages. | Stage 2 output | `ResponseRecord2::record_status` records metadata side. |
| `load_response_vector(...)` templates | Loads archive into appropriately shaped response variant; handles `k` mismatch reconstruction/projection. | Restart reuse and dependency loading across stages/protocols. | Stage 2/3 input | `seed_response_guess_for_point(...)`, property loaders, VBC loaders. |
| point adapter overloads | Convenience wrappers for `LinearResponsePoint`. | Reduces boilerplate where point objects already exist. | Stage 2 | descriptor overloads are generic form. |

### `src/apps/molresponse_v2/FrequencyLoop.hpp` + `.cpp`

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `ResponseSolveDiagnostics` | Per-point solver outcome (`converged`, iterations, final residual, final density change). | Captures explicit diagnostics for metadata and run analysis. | Stage 2 | stored via `record_solver_diagnostics`. |
| `iterate<ResponseType>(...)` | KAIN/nonlinear iteration loop for one response variant with debug/timing hooks and convergence checks. | Shared iterative solve implementation across variant types. | Stage 2 | variant-specific kernels in `ResponseSolver.cpp`. |
| `solve_response_vector(...)` | Variant visitor dispatch to `iterate(...)` for supported types. | Keeps variant branching in one place. | Stage 2 | `computeFrequencyLoop`, `computeFrequencyPoint`. |
| `promote_response_vector(...)` | Converts/copies static/dynamic and restricted/unrestricted variants for continuation seeding. | Supports reuse of previous solutions across frequency/protocol contexts. | Stage 2 init | used by guess seeding and VBC BC-vector loading. |
| `GuessSeedResult` | Holds seeded response guess plus restart provenance chosen by the seed policy. | Returns explicit seed-result data to callers instead of side-effect-only status. | Stage 2 init | `RestartGuessProvenance`, `seed_response_guess_for_point(...)`. |
| `seed_response_guess_for_point(...)` | Guess precedence logic: same-point archive, coarser protocol, previous frequency memory/archive, initializer. Returns `GuessSeedResult`. | Encapsulates restart/continuation policy in one helper. | Stage 2 init | replacement for ad hoc seed logic in old loop. |
| `should_solve_point(...)` | Shared predicate for point solve gating (`!saved || (final && !converged)`). | Keeps solve eligibility behavior consistent across channel-loop and single-point paths. | Stage 2 scheduling | `point_needs_solving(...)` in orchestrator stage-2 logic. |
| `solve_and_record_point(...)` | Shared point lifecycle: iterative solve, timing, diagnostics, provenance, archive save, and status record. | Eliminates duplicate point-solve persistence/logging code between loop and point entrypoints. | Stage 2 execution | `computeFrequencyLoop(...)`, `computeFrequencyPoint(...)`. |
| `PointSolveResult` | Carries solved response vector and diagnostics from `solve_and_record_point(...)`. | Makes post-solve continuation handoff explicit (e.g., previous-frequency reuse). | Stage 2 execution | `computeFrequencyLoop(...)` previous-frequency continuation. |
| `computeFrequencyLoop(...)` | Solves an entire channel frequency series at one protocol index with persistence/logging. | Channel-series execution mode for serial/protocol-0 workflows. | Stage 2 | `computeFrequencyPoint(...)` for independent point mode. |
| `computeFrequencyPoint(...)` | Solves a single `(state, protocol, frequency)` point with same persistence/logging semantics. | Enables channel-point ownership mode for finer-grained subgroup parallelism. | Stage 2 | `computeFrequencyLoop(...)`. |
| `StateSolvePersistence` (interface) | Abstracts save/converged checks and metadata/debug recording from solve logic. | Solver stays storage-backend agnostic. | Stage 2 infra | implemented by `JsonStateSolvePersistence`. |

### `src/apps/molresponse_v2/ResponseSolver.hpp` + `.cpp`

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `response_vector_allocator` + `response_solver` alias | Allocator and solver type alias for `XNonlinearSolver`. | Standardizes nonlinear solver instantiation for response vectors. | Stage 2 | used in `iterate(...)`. |
| `alpha_factor(...)` overloads | Returns tensor prefactor by response variant. | Keeps variant-specific prefactors constexpr and centralized. | Stage 2/3 diagnostics | alpha scaling also in property assembly. |
| `compute_density(...)` overloads | Computes response-induced density for each variant. | Needed for convergence criterion and coupled equations. | Stage 2 | currently restricted variants implemented; unrestricted throws. |
| `make_bsh_operators(...)` overloads | Builds BSH operators for static/dynamic restricted variants. | Variant-specific kernel setup for coupled-response update. | Stage 2 | `ResponseSolverUtils::make_bsh_operators_response(...)`. |
| `CoupledResponseEquations(...)` overloads | Evaluates coupled-response RHS updates (exchange/local terms, projectors, apply BSH). | Core physics kernel in iterative update. | Stage 2 | called inside `iterate(...)`. |

### `src/apps/molresponse_v2/ResponseSolverUtils.hpp`

| Function | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `print_iteration_line(...)` | Prints iteration residual/drho/inner-product line with targets. | Human-readable convergence trace. | Stage 2 | debug JSON in `ResponseDebugLogger`. |
| `make_bsh_operators_response(...)` | Generic BSH operator builder from shifts/energies/freq. | Shared primitive used by solver variants and excited scaffold. | Stage 2 | `ResponseSolver::make_bsh_operators(...)`. |
| `inner(...)` | Dot-product helper for response vectors. | Avoids repeated manual loops. | Stage 2/3 | used in solver and excited scaffold. |
| `do_step_restriction(...)` | Optional step damping with max rotation threshold. | Stabilization hook for solver updates. | Stage 2 | gate currently disabled in `iterate` branch. |

## Persistence and Debug Logging

### `src/apps/molresponse_v2/ResponseRecord.hpp` (`class ResponseRecord2`)

| Function / Group | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| constructor `ResponseRecord2(World&, filepath)` | Loads/broadcasts existing metadata JSON or initializes empty root. | MPI-consistent metadata state for all ranks. | Stage 2 infra | `ResponseDebugLogger` for per-iteration logs. |
| `initialize_states(...)` | Seeds metadata tree for all planned linear states/protocols/frequencies. | Guarantees keys exist before solve updates. | Stage 2 init | `initialize_excited_bundle(...)` for excited stage. |
| excited recorders (`initialize_excited_bundle`, `record_excited_protocol_status/timing/energies`) | Manage excited-state protocol metadata subtree. | Integrates excited adapter outputs into shared metadata schema. | Stage 2c | `execute_excited_state_bundle_stage(...)`. |
| `record_status(...)` | Marks point saved + converged bool. | Core per-point status persistence. | Stage 2 | `record_timing`, `record_solver_diagnostics`, `record_restart_provenance`. |
| `record_timing(...)` | Stores per-point wall/cpu timing. | Performance and progress diagnostics. | Stage 2 | debug logger timings are iteration-level detail. |
| `record_solver_diagnostics(...)` | Stores final convergence diagnostics per point. | Postmortem for failed/slow points and tuning. | Stage 2 | `ResponseSolveDiagnostics`. |
| `record_restart_provenance(...)` | Stores seed source metadata (same/coarser/prev-frequency/initializer). | Explains restart behavior and convergence differences. | Stage 2 | guess logic in `seed_response_guess_for_point(...)`. |
| gating/query helpers (`is_saved`, `is_converged`, `missing_at_final_protocol`, etc.) | Query metadata for readiness and missing points. | Restart-aware orchestration and final checks. | Stage 2/3 | `point_ready_in_metadata(...)` helpers in orchestrator. |

### `src/apps/molresponse_v2/ResponseDebugLogger.hpp`

| Function / Group | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `start_state(...)`, `begin_iteration(...)`, `end_iteration(...)` | Opens logging context for point and appends one iteration entry. | Structured per-iteration diagnostics keyed by state/protocol/freq. | Stage 2 | macros in `ResponseDebugLoggerMacros.hpp`. |
| `log_value(...)`, `log_timing(...)`, `log_value_and_time(...)` | Stores step-level numeric values and timing. | Lightweight instrumentation points inside solver kernels. | Stage 2 | console timing in `TimedValueLogger`. |
| `write_to_disk()` | Writes aggregated debug JSON. | Persist run diagnostics for analysis. | Stage 2 | `to_json()` for in-memory merge. |
| `print_timing_table(...)`, `print_values_table(...)` | Pretty console summaries from logged iteration arrays. | Human-readable diagnostics without JSON parsing. | Stage 2 | persistent JSON logs for tooling. |
| `TimedValueLogger` | RAII-style wall/cpu/value logger helper around blocks. | Reduces instrumentation boilerplate. | Stage 2 | explicit `log_timing`/`log_value` calls. |

### `src/apps/molresponse_v2/ResponseDebugLoggerMacros.hpp`

| Macro | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `DEBUG_LOG_VALUE(...)` | Wraps an expression in `TimedValueLogger`, records its value and timing under `key`. | Lightweight instrumentation with one call-site macro and minimal boilerplate. | Stage 2 | explicit `TimedValueLogger` object + `log(...)`. |
| `DEBUG_TIMED_BLOCK(...)` | Times an arbitrary block and logs timing (with optional value hooks inside block). | Convenient scoped timing for solver sub-steps. | Stage 2 | explicit `TimedValueLogger` + manual `log()`. |

## Derived and Excited-State Helpers

### `src/apps/molresponse_v2/VBCMacrotask.hpp` (`class SimpleVBCComputer`)

| Function | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| constructor | Captures world, ground state, orbital count, spin mode. | Reusable VBC helper object for beta/raman tasks. | Stage 2d/3 | created in derived execution and property beta compute. |
| `get_BC_vecs(...)` | Loads B and C linear states and promotes to dynamic form. | Normalizes BC input representation for VBC contractions. | Stage 2d/3 | `promote_response_vector(...)`. |
| `make_zeta_bc(...)` | Computes zeta coupling term. | Shared VBC intermediate. | Stage 2d/3 | used in beta contraction path. |
| `compute_g(...)`, `compute_vbc_i(...)` | Internal Coulomb/exchange contraction kernels for VBC response build. | Encapsulates VBC intermediate math. | Stage 2d/3 | property-side `compute_beta(...)` also builds related zeta terms. |
| `compute_vbc_response(...)` | Assembles final VBC dynamic response from BC vectors/operators. | Core VBC response constructor. | Stage 2d/3 | `compute_and_save(...)` wraps persistence logic. |
| `compute_and_save(...)` | Load existing VBC archive or compute and persist a new one. | Restart-aware, reusable entrypoint for VBC state production. | Stage 2d/3 | `DerivedStatePlanner::make_vbc_state(...)` creates descriptor input. |

### `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp` + `.cpp`

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `ExcitedBundleSolverConfig`, `ExcitedBundleProtocolInput`, `ExcitedBundleProtocolResult` | Data contracts for excited bundle adapter. | Keeps excited stage interface explicit and serializable. | Stage 2c | metadata rows in `ResponseRecord2` excited subtree. |
| `ExcitedStateBundleSolver` (abstract) | Interface for protocol-wise excited solves. | Adapter seam for reintegration without altering orchestrator contracts. | Stage 2c | `ExcitedStateBundleNoopSolver`, scaffold adapter. |
| `ExcitedStateBundleNoopSolver` | Placeholder implementation returning minimal status. | Behavior-preserving fallback when full solver not active. | Stage 2c | `RestartAwareExcitedScaffoldSolver` in `.cpp`. |
| `make_excited_state_bundle_solver_adapter(...)` | Factory returning active adapter implementation. | Centralized adapter selection. | Stage 2c | invoked from `execute_excited_state_bundle_stage`. |
| `RestartAwareExcitedScaffoldSolver` | Restart-aware scaffold implementation with protocol prep, guess seeding, iterative placeholder updates, and snapshot I/O. | Provides executable excited-stage behavior before full legacy reintegration. | Stage 2c | no-op solver and future full solver adapter. |
| `read_restart_snapshot(...)`, `write_restart_snapshot(...)` | JSON snapshot persistence for excited protocol state. | Restart continuity and protocol skipping logic. | Stage 2c | `file_exists_world(...)` utility. |

## Property Assembly

### `src/apps/molresponse_v2/PropertyManager.hpp`

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `PropertyManager` class | Stores flat property rows (`alpha`, `beta`, `raman`) with load/save/lookup/setters. | Unified property persistence and dedupe checks. | Stage 3 | merged shard component rows in orchestrator. |
| `compute_response_inner_product_tensor(...)` | Builds contraction tensor from response vector sets; optionally logs per-function contributions. | Shared tensor contraction primitive used across alpha/beta/raman workflows. | Stage 3 | `global_inner_contributions()` backing store. |
| `compute_alpha(...)` | Loads dipole linear states and computes frequency-dependent alpha tensor. Skips missing frequency states. | Polarizability assembly from solved linear states. | Stage 3 | `compute_hyperpolarizability(...)`, `compute_raman(...)`. |
| `compute_hessian(...)` | Computes vibrational Hessian and derived vibrational quantities from nuclear response states. | Raman pipeline prerequisite. | Stage 3 | called by `compute_raman(...)` wrapper in orchestrator. |
| `compute_beta(...)` (pair-list overload) | Main beta/raman component task loop over explicit `(omegaB,omegaC)` pairs and BC perturbation pairs, with claim-file based dedupe in subgroup mode; skips per-task failures. | Shared engine for beta and Raman-component generation. | Stage 3 | cartesian overload wraps this; `compute_Raman_components(...)` uses Raman BC pairs. |
| `compute_beta(...)` (pair-of-vectors overload) | Convenience wrapper converting two frequency vectors to cartesian pair list. | Backward compatible call surface. | Stage 3 | explicit pair-list overload. |
| `select_beta_frequency_triplets(...)` | Applies triplet mode selection (`all_triplets`, `or`, `shg`). | Avoids unnecessary beta pair expansion and supports process-specific runs. | Stage 3 | `DerivedStatePlanner::select_beta_frequency_pairs(...)` for derived planning parity. |
| `compute_hyperpolarizability(...)` | Builds A and BC dipole perturbation sets and dispatches beta engine with selected triplets. | Top-level beta assembly entry used by orchestrator. | Stage 3 | `compute_Raman_components(...)`. |
| `compute_Raman_components(...)` | Builds dipole+nuclear BC set with `(omegaB,0)` pattern and dispatches beta engine for Raman components. | Reuses beta contraction core for Raman derivatives. | Stage 3 | `compute_hyperpolarizability(...)` (dipole-dipole variant). |
| `compute_Raman(...)` | Ensures Raman components exist then assembles component tensors by mode/frequency. | Produces final Raman tensor artifacts from component table. | Stage 3 | orchestrator wrapper `compute_raman(...)`. |
| `try_claim_property_component_task(...)` | Creates lock-file claims per component task. | Prevents duplicate subgroup task execution. | Stage 3 subgroup mode | claim-prefix logic in orchestrator component precompute stage. |

### `src/apps/molresponse_v2/InnerContributions.hpp`

| Function | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `global_inner_contributions()` | Returns process-global JSON object storing contraction contribution entries. | Optional diagnostics sink for `compute_response_inner_product_tensor(...)`. | Stage 3 diagnostics | disabled unless save flag enabled in contraction calls. |

## Supporting Utilities and Types

### `src/apps/molresponse_v2/Perturbation.hpp`

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `DipolePerturbation`, `NuclearDisplacementPerturbation`, `MagneticPerturbation` | Small typed perturbation structs. | Type-safe perturbation representation. | Stage 1/2/3 | consumed through `Perturbation` variant. |
| `Perturbation` | Variant wrapper over perturbation structs. | Uniform API across planning/solver/property code. | Stage 1/2/3 | `raw_perturbation_operator(...)` visits this type. |
| `describe_perturbation(...)` | Canonical string id for perturbation. | Stable key for filenames, metadata, and maps. | Stage 1/2/3 | `perturbation_type_string(...)` gives coarse class name only. |
| `perturbation_type_string(...)` | Human-readable class label (`Dipole`, `Nuclear`, `Magnetic`). | Diagnostics/printing convenience. | Stage 1 diagnostics | `describe_perturbation(...)` for unique IDs. |

### `src/apps/molresponse_v2/broadcast_json.hpp`

| Function | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `broadcast_json_file(World&, filepath)` | Rank 0 reads JSON file then broadcasts serialized content to all ranks. | Ensures all ranks share identical JSON payload without per-rank filesystem assumptions. | Stage 1/2 infra | file readers in orchestrator (`read_json_file_or_object`) are local-rank helpers, not world broadcast helpers. |

### `src/apps/molresponse_v2/MolecularProperty.hpp`

| Function / Type | What it does | Why it exists | Stage | Similar / related |
|---|---|---|---|---|
| `MolecularPropertyType`, `MolecularProperty` | Legacy/auxiliary explicit property request container. | Simple typed property request representation for tooling and debug print paths. | Stage 1 (aux) | main workflow now derives requests from `ResponseParameters` + `StateGenerator`. |
| `property_type_to_string(...)`, `print_requested_properties(...)` | Utility printers for property request vectors. | Quick visibility during development/tests. | Stage 1 diagnostics | richer metadata emitted via orchestrator JSON outputs. |

## Tests and Validation Tooling

### `src/apps/molresponse_v2/test_preliminaries.cpp`

| Function / Area | What it does | Why it exists | Similar / related |
|---|---|---|---|
| `tensor_approx_equal(...)` and `main(...)` flow | Validates preliminaries/Hamiltonian against reference fock data across protocol thresholds. | Regression check for protocol setup + preliminaries consistency. | Ground/protocol helpers in `GroundStateData` and `ResponseManager`. |

### `src/apps/molresponse_v2/test_parameter_manager.cpp`

| Function / Area | What it does | Why it exists | Similar / related |
|---|---|---|---|
| `main(...)` | Smoke test for parameter manager merge/print behavior for response-related parameter sets. | Input parsing regression sanity for mixed parameter groups. | `ResponseParameters` constructor/validation paths. |

### `src/apps/molresponse_v2/validate_property_component_claims.py`

| Function / Area | What it does | Why it exists | Similar / related |
|---|---|---|---|
| claim/shard/log parsers and checks | Offline validator for subgroup property-component scheduling artifacts, duplicate rows, and duplicate VBC writes. | Debugging aid when queue jobs complete but scheduling behavior needs audit. | lock-claim path in `try_claim_property_component_task(...)` and property subgroup component stage in `MolresponseLib.hpp`. |

## Similar-Function Quick Index

This section highlights frequent "which one should I use" pairs.

| If you need... | Prefer | Instead of | Reason |
|---|---|---|---|
| Solve one channel at one protocol in frequency order | `computeFrequencyLoop(...)` | Manual frequency loop around `computeFrequencyPoint(...)` | Keeps continuation/restart behavior and diagnostics consistent. |
| Solve individual points for point ownership mode | `computeFrequencyPoint(...)` | `computeFrequencyLoop(...)` | Enables independent point scheduling per subgroup lane. |
| Determine static parallel plan from inputs | `StateParallelPlanner::build(...)` | Hand-building group assignments in orchestrator | Centralized policy and diagnostics. |
| Runtime restart-aware protocol behavior | `build_state_solve_schedule_context(...)` + runtime policy helpers | Using planner output verbatim | Runtime must account for already-saved protocol data. |
| Build linear states from properties | `StateGenerator::generateStates()` | Manual perturbation/frequency expansion | Dedupes and canonicalizes frequencies consistently. |
| Build derived VBC requests | `DerivedStatePlanner::build_vbc_driven_quadratic_plan(...)` | Expanding VBC requests inside property kernels | Keeps derived execution independently schedulable and gateable. |
| Beta pair policy selection | `select_beta_frequency_triplets(...)` (property) and `select_beta_frequency_pairs(...)` (derived planner) | Full cartesian loops everywhere | Avoids unnecessary work and keeps planning/runtime parity. |
| Property assembly in serial | `compute_requested_properties(...)` | Calling individual property kernels ad hoc | Preserves ordered dispatch and shared context. |
| Property assembly with subgroup support | `compute_requested_properties_with_property_group(...)` | Rolling custom subgroup property logic | Includes lock-claims, fallback, merge, broadcast, and error handling. |

## Maintenance Checklist For New Features

When adding a new response feature, update these locations together:

1. Input surface
- `ResponseParameters.hpp` (default + getter + validation/derivation if needed)

2. Planning surface
- `StateGenerator.hpp` for new linear-state requirements
- `DerivedStatePlanner.hpp` for new derived request classes/dependencies
- `StateParallelPlanner.hpp` if ownership/scheduling policy changes

3. Solve path
- `FrequencyLoop.hpp/.cpp` and `ResponseSolver*` for new iteration/diagnostic behavior
- `ResponseRecord.hpp` schema if new metadata fields are persisted

4. Property path
- `PropertyManager.hpp` for assembly changes and storage schema
- `MolresponseLib.hpp` property dispatch wrappers for stage integration

5. Documentation
- `MOLRESPONSE_TUTORIAL.md` for stage-flow changes
- `STATE_PARALLEL_DESIGN.md` for scheduling semantics
- This file for function-level intent and similarities
