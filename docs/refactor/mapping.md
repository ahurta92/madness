# Legacy-to-Current `molresponse` Mapping

## Main Executable Path

### Step 1: Start MADNESS, parse the input filename, and read response plus ground-state archives
- **Legacy:** `main`, `read_and_create_density`, `ResponseParameters::read_and_set_derived_values`, `GroundParameters::read` — bootstrap the executable, parse input, and load response plus ground-state context.
- **Current:** `add_response_workflow_drivers`, `molresponse_lib::run_response`, `ResponseParameters::ResponseParameters(World&, parser)`, `make_ground_context` — route the workflow, rebuild response parameters, and load the SCF checkpoint plus ground runtime context.
- **What current adds:** protocol, naming
- **Mixed concerns in current:** none in the direct startup path; `run_response` is orchestration-heavy but not tagged `MIXED` in the current inventory.
- **Refactor action:** KEEP

### Step 2: Construct a `density_vector` matching the requested calculation mode, then construct `TDDFT`
- **Legacy:** `set_density_type`, `TDDFT::TDDFT` — choose the response container/model and build the monolithic solver object.
- **Current:** `plan_required_states`, `StateGenerator::generateStates`, `build_excited_state_bundle_plan`, `ResponseManager`, `GroundStateData` — split planning, runtime context, and solver support objects into separate stage-1 components instead of one `TDDFT` god-object.
- **What current adds:** naming, parallel
- **Mixed concerns in current:** `plan_required_states` is `SOLVER`; no explicit mixed flag, but it bundles linear, derived, excited, and parallel planning in one stage function.
- **Refactor action:** KEEP

### Step 3: Build a coarse protocol/load-balance setup if running with more than one rank
- **Legacy:** `TDDFT::set_protocol`, `TDDFT::make_nuclear_potential`, `TDDFT::initial_load_bal` — set numerical defaults and do coarse whole-world load balancing.
- **Current:** `ResponseManager::setProtocol`, `GroundStateData::prepareOrbitals`, `GroundStateData::computePreliminaries`, `build_state_solve_schedule_context` — set protocol numerics, rebuild ground preliminaries, and construct a restart-aware state-parallel execution schedule.
- **What current adds:** protocol, parallel, restart
- **Mixed concerns in current:** `build_state_solve_schedule_context` and its helper family are `SOLVER (MIXED: restart logic, protocol handling, state-parallel logic)`; `GroundStateData::computePreliminaries` is `OPERATOR (MIXED: IO, model-specific branching)`.
- **Refactor action:** EXTRACT

### Step 4: Select the first protocol threshold and dispatch to either excited-state or first-order response
- **Legacy:** `main` — choose the high-level solve path and enter the protocol loop.
- **Current:** `solve_all_states`, `execute_serial_state_solve` / `execute_subgroup_state_solve`, `execute_excited_state_bundle_stage`, `execute_derived_state_requests` — stage-2 orchestrator dispatches linear, excited, and derived stages under one restart-aware schedule.
- **What current adds:** restart, protocol, parallel
- **Mixed concerns in current:** `solve_all_states` is `SOLVER (MIXED: protocol handling, restart logic, state-parallel logic)`; `execute_serial_state_solve` and `execute_subgroup_state_solve` are also mixed.
- **Refactor action:** EXTRACT

### Step 5: For each protocol threshold, update defaults, reproject state, restart or guess, and iterate to convergence
- **Legacy:** `TDDFT::solve_excited_states`, `TDDFT::iterate_excited`, `TDDFT::solve_response_states`, `TDDFT::iterate_freq2` — own the full per-threshold lifecycle in one driver.
- **Current:** `run_protocol_threshold_loop`, `prepare_protocol_context`, `computeFrequencyLoop`, `computeFrequencyPoint`, `seed_response_guess_for_point`, `iterate`, `solve_and_record_point` — split protocol progression, protocol setup, restart seeding, point solve, and persistence recording into separate pieces.
- **What current adds:** naming, restart, protocol, parallel
- **Mixed concerns in current:** `execute_serial_state_solve`, `execute_subgroup_state_solve`, `computeFrequencyLoop`, `computeFrequencyPoint`, `seed_response_guess_for_point`, `solve_and_record_point`, and `iterate` all carry `MIXED` flags in the current inventory.
- **Refactor action:** REWRITE

### Step 6: Optionally save restart files and plot orbitals/densities
- **Legacy:** `TDDFT::save`, `do_vtk_plots`, plotting helpers — save a restart and optionally emit plots.
- **Current:** `save_response_vector`, `JsonStateSolvePersistence`, `ResponseRecord2`, `ResponseDebugLogger` — always write structured restart/status/debug artifacts; plotting is no longer a central part of the current flow.
- **What current adds:** naming, restart
- **Mixed concerns in current:** restart output is embedded in `solve_and_record_point` (`SOLVER (MIXED: restart logic, IO persistence)`), and the subgroup/serial stage-2 runners also mix persistence with orchestration.
- **Refactor action:** EXTRACT

### Step 7: For dipole response, compute and print the polarizability tensor
- **Legacy:** `TDDFT::polarizability`, `TDDFT::PrintPolarizabilityAnalysis` — assemble and print alpha directly from the in-memory solver state.
- **Current:** `compute_polarizability`, `compute_alpha`, `PropertyManager::save` — filter converged frequencies from metadata, reload final response archives, assemble alpha, and persist the flat property table.
- **What current adds:** restart, naming
- **Mixed concerns in current:** `compute_polarizability` is `PROPERTY (MIXED: IO archive loads, metadata/restart readiness filtering)`; `compute_alpha` is `PROPERTY (MIXED: IO archive loads, naming/state lookup)`.
- **Refactor action:** EXTRACT

## Excited-State Iteration Flow

### Step 1: Build projector, XC operator, KAIN solvers, residual workspaces, and convergence thresholds
- **Legacy:** `TDDFT::iterate_excited` — initializes the entire excited-state iteration environment inline.
- **Current:** `execute_excited_state_bundle_stage`, `ExcitedResponse::solve_protocol`, `RestartAwareExcitedScaffoldSolver::prepare_protocol`, excited bundle kernels in `excited_ops/*.hpp` — split excited protocol setup from the top-level stage dispatcher and the low-level bundle kernels.
- **What current adds:** restart, protocol, naming
- **Mixed concerns in current:** `execute_excited_state_bundle_stage` is `SOLVER (MIXED: restart logic, naming logic, protocol handling, model-specific branching)`; `ExcitedResponse::solve_protocol` is also mixed.
- **Refactor action:** EXTRACT

### Step 2: At each iteration, recompute transition densities and optionally redistribute work
- **Legacy:** `TDDFT::make_density`, `TDDFT::loadbal` — rebuild excited-state densities and rebalance whole-world data.
- **Current:** bundle kernels such as `compute_response_densities`, `compute_excited_potentials`, and subgroup ownership/runtime policy in `MolresponseLib.hpp` — density recomputation stays in kernel space, while redistribution is replaced by explicit subgroup scheduling rather than per-iteration global load balancing.
- **What current adds:** parallel
- **Mixed concerns in current:** subgroup scheduling sits in the stage-2 scheduling helpers and `execute_subgroup_state_solve`, both mixed with restart/protocol concerns.
- **Refactor action:** REWRITE

### Step 3: Form residual diagnostics from density changes and BSH residuals
- **Legacy:** `TDDFT::iterate_excited` — computes and prints residual diagnostics inline.
- **Current:** `ExcitedResponse::solve_protocol`, `ExcitedIterDiagnostics`, `ResponseDebugLogger`, and bundle-kernel diagnostics like `print_potential_diagnostics` — residual tracking is more structured and metadata-oriented.
- **What current adds:** restart
- **Mixed concerns in current:** `ExcitedResponse::solve_protocol` is `SOLVER (MIXED: restart logic, protocol handling, model-specific TDA/full branching)`.
- **Refactor action:** KEEP

### Step 4: Update the excited-state subspace and energies
- **Legacy:** `TDDFT::update_x_space_excited` — perform the nonlinear excited-state update and refresh energies.
- **Current:** bundle kernels such as `ExcitedStateStep`, `diagonalize_and_rotate_bundle`, `apply_step_project_normalize`, orchestrated by `ExcitedResponse::solve_protocol` and the scaffold solver.
- **What current adds:** restart, naming, protocol
- **Mixed concerns in current:** orchestration-side mixed flags remain in `execute_excited_state_bundle_stage` and `ExcitedResponse::solve_protocol`; the bundle kernels themselves are comparatively clean.
- **Refactor action:** KEEP

### Step 5: Sort final roots, print diagnostics, and run transition-property analysis
- **Legacy:** `TDDFT::sort`, `TDDFT::analysis`, `TDDFT::analyze_vectors` — finalize root ordering and print excited-state properties.
- **Current:** `ResponseRecord2::assign_excited_state_names`, `build_excited_root_manifest`, `execute_excited_state_bundle_stage` — root naming, stable root ids, and per-protocol excited metadata are explicit; there is no one direct current equivalent of the old inline analysis printout.
- **What current adds:** naming, restart
- **Mixed concerns in current:** `execute_excited_state_bundle_stage` is mixed with naming/restart/protocol/model branching.
- **Refactor action:** REWRITE

## Frequency-Response Iteration Flow

### Step 1: Build the property RHS (`PQ`) for dipole or nuclear perturbations
- **Legacy:** `TDDFT::PropertyRHS` — build the projected RHS for the active perturbation.
- **Current:** `make_perturbation_operator`, `raw_perturbation_operator`, `project_perturbation_onto_orbitals`, `perturbation_vector`, `initialize_guess_vector` — separate perturbation naming, operator construction, and projected-vector assembly from the solver loop.
- **What current adds:** naming
- **Mixed concerns in current:** none called out as mixed in this perturbation layer.
- **Refactor action:** KEEP

### Step 2: Build response BSH operators for `+omega` and `-omega`, including energy shifts when needed
- **Legacy:** `TDDFT::make_bsh_operators_response` — build the response BSH operator list inline in the legacy solver.
- **Current:** `ResponseSolverUtils::make_bsh_operators_response`, per-model `make_bsh_operators` in `ops/*.hpp`, and `ResponseKernels.hpp` — separate common BSH setup from model-specific operator kernels.
- **What current adds:** model-specific branching
- **Mixed concerns in current:** the model-specific branch is concentrated in `solve_response_vector` and `iterate`, both mixed on the solver side, not in the operator helper itself.
- **Refactor action:** KEEP

### Step 3: At each iteration, recompute transition density, load balance, and test density plus BSH residuals
- **Legacy:** `TDDFT::iterate_freq2` — does density rebuild, balancing, residual testing, and progress reporting in one loop.
- **Current:** `iterate`, per-model `compute_density`, `CoupledResponseEquations`, `ResponseDebugLogger`, plus subgroup schedule/runtime policy outside the iteration kernel — the numerical loop is cleaner, but runtime ownership and restart policy live around it.
- **What current adds:** restart, parallel, model-specific branching
- **Mixed concerns in current:** `iterate` is `SOLVER (MIXED: model-specific type dispatch via templates, convergence policy)`; surrounding stage-2 runners are mixed with restart/protocol/parallel concerns.
- **Refactor action:** EXTRACT

### Step 4: Apply the response update step with optional KAIN acceleration
- **Legacy:** `TDDFT::update_x_space_response` — do one response update and optional KAIN acceleration.
- **Current:** `iterate`, response-solver wrappers in `ResponseSolver.hpp/.cpp`, per-model kernels in `ops/*.hpp` — separate solver control from model-specific equation application.
- **What current adds:** model-specific branching
- **Mixed concerns in current:** `iterate` and `solve_response_vector` carry the model-specific mixed flags.
- **Refactor action:** KEEP

### Step 5: Save and/or plot the converged response orbitals and densities
- **Legacy:** `TDDFT::save`, `do_vtk_plots` — emit a restart archive and optional plots after convergence.
- **Current:** `save_response_vector`, `solve_and_record_point`, `JsonStateSolvePersistence`, `ResponseRecord2` — save canonical point archives and structured metadata; no central plotting analogue remains in the main current flow.
- **What current adds:** naming, restart
- **Mixed concerns in current:** `solve_and_record_point` is `SOLVER (MIXED: restart logic, IO persistence, failure policy)`.
- **Refactor action:** EXTRACT

## Current Features With No Legacy Equivalent

### Stable naming and metadata keys
- **Layer:** naming layer
- **Current code:** `LinearResponseDescriptor::make_key`, `LinearResponsePoint::response_filename`, `ResponseRecord2::freq_key`, `ResponseRecord2::protocol_key`, excited root-id/state-name helpers
- **Why new:** legacy relied on user-provided filenames and had no protocol-aware or frequency-aware canonical naming scheme.

### JSON-backed state persistence and restart provenance
- **Layer:** naming layer
- **Current code:** `JsonStateSolvePersistence`, `ResponseRecord2`, `record_restart_provenance`, solver diagnostics metadata
- **Why new:** legacy restart was just exact-file save/load with no point-level status tree or provenance recording.

### Restart-aware runtime scheduling
- **Layer:** solver layer
- **Current code:** `compute_runtime_point_ownership_policy`, `build_protocol_execution_policy`, `point_needs_solving_from_metadata`
- **Why new:** legacy had no metadata-driven skip/restart policy by protocol and point.

### Layered restart seed precedence
- **Layer:** solver layer
- **Current code:** `seed_response_guess_for_point`, `seed_retry_guess_from_previous_frequency_archive`, `promote_response_vector`
- **Why new:** legacy either loaded the exact restart file or started fresh.

### Protocol-indexed excited-state bundle stage
- **Layer:** solver layer
- **Current code:** `execute_excited_state_bundle_stage`, `ExcitedResponse::solve_protocol`, `ExcitedStateBundleSolver.cpp`
- **Why new:** legacy had excited-state iteration, but not a protocol-aware metadata-recorded bundle stage with stable root identity and restart-ready skip logic.

### Derived-state dependency gate and VBC execution stage
- **Layer:** solver layer
- **Current code:** `DerivedStatePlanner`, `execute_derived_state_requests`, `run_derived_request`, `SimpleVBCComputer`
- **Why new:** the current code adds a distinct derived-state stage driven by dependency gating rather than mixing all higher-order work into one procedural flow.

### State-parallel subgroup scheduling
- **Layer:** parallel layer
- **Current code:** `StateParallelPlanner`, `PointOwnershipScheduler`, `execute_subgroup_state_solve`, `MacroTaskQ::create_worlds`
- **Why new:** legacy only had whole-world load balancing, not subgroup/subworld scheduling or owner-lane work partitioning.

### Point-parallel execution with protocol-dependent mode switching
- **Layer:** parallel layer
- **Current code:** `use_channel_series_ownership_for_protocol_runtime`, `build_pending_work_manifest`, `computeFrequencyPoint`
- **Why new:** legacy never fanned out one channel into independent channel-frequency point tasks.

### Property-group execution and distributed component precompute
- **Layer:** parallel layer
- **Current code:** `compute_requested_properties_with_property_group`, shard merge of `properties_components.json`, property task claim files
- **Why new:** legacy property assembly always ran in the main world with no subgroup property stage.

### Structured Raman workflow with vibrational post-processing objects
- **Layer:** property layer
- **Current code:** `compute_hessian`, `compute_Raman_components`, `compute_Raman`, `print_raman_table`, `VibrationalResults`, `RamanResults`
- **Why new:** legacy inventory only highlighted direct dipole polarizability output; the current code has a formal stage-3 Raman pipeline with saved intermediate property rows and spectra objects.

### Preflight memory estimate before subgroup allocation
- **Layer:** parallel layer
- **Current code:** `preflight_memory_estimate`
- **Why new:** legacy had no explicit preflight memory guardrail before distributed execution.
