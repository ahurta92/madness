# Molresponse Tutorial (Developer-Focused)

This document explains how the current `molresponse` path is structured, how data flows through the stages, and how the new state-parallel refactor pieces fit together.

## 1) High-Level Pipeline

`molresponse_lib::run_response(...)` in `src/madness/chem/MolresponseLib.hpp` runs in three major stages:

1. Build runtime context and parse requested work.
2. Solve required states (linear + derived), including restart-aware state-parallel scheduling.
3. Assemble requested properties (polarizability, hyperpolarizability, Raman) from solved states.

The refactor goal has been to keep this orchestration readable while moving details into focused helpers.

## 2) Stage 1: Inputs and Planning

Planning primarily happens in:

- `src/apps/molresponse_v2/StateGenerator.hpp`
- `src/apps/molresponse_v2/StateParallelPlanner.hpp`
- `src/apps/molresponse_v2/DerivedStatePlanner.hpp`

Responsibilities:

- Convert response knobs (frequencies, directions, requested properties) into a concrete linear-state list.
- Canonicalize frequencies to avoid floating-point near-duplicate scheduling.
- Build state ownership and point ownership policies for state-parallel execution.
- Build derived-state request plans (currently VBC-driven).
- Provide dependency gating metadata for derived requests.

Output of Stage 1 is a `PlannedStates` object consumed by Stage 2.

## 3) Stage 2: Solving States

Stage 2 orchestration lives in `solve_all_states(...)` inside `src/madness/chem/MolresponseLib.hpp`.

### 3.1 Runtime Scheduling Context

`build_state_solve_schedule_context(...)` decides:

- serial vs subgroup path,
- state-ownership vs point-ownership mode by protocol,
- restart promotion behavior (if protocol-0 points are already saved),
- effective group counts for state and point modes.

### 3.2 Linear State Solves

Two execution paths:

- `execute_serial_state_solve(...)`
- `execute_subgroup_state_solve(...)`

Both use `computeFrequencyLoop(...)` for channel-series lanes and
`computeFrequencyPoint(...)` for channel-point lanes.

Inside the frequency solve path (`src/apps/molresponse_v2/FrequencyLoop.cpp`):

- `should_solve_point(...)` centralizes final-protocol re-solve gating.
- `seed_response_guess_for_point(...)` chooses best available initial vector:
  - exact point file,
  - coarser protocol point,
  - previous frequency continuation,
  - fresh initializer.
- `solve_and_record_point(...)` runs the iterative solve
  (`solve_response_vector(...)`) and records timing, diagnostics, restart
  provenance, save-status, and archives in one shared lifecycle helper.

Inside stage-2 orchestration (`src/madness/chem/MolresponseLib.hpp`):

- `prepare_protocol_context(...)` standardizes protocol setup work.
- `log_pending_manifest(...)` standardizes pending-work logging.
- `execute_manifest_work(...)` centralizes channel-vs-point manifest dispatch.

### 3.3 Persistence and Metadata

`JsonStateSolvePersistence` wraps:

- `ResponseRecord2` (`src/apps/molresponse_v2/ResponseRecord.hpp`) for status/timing metadata,
- `ResponseDebugLogger` for iteration-level debug values and timings.

In subgroup mode, each group writes shard files:

- `response_metadata.group<gid>.json`
- `response_log.group<gid>.json`

Universe rank 0 merges shards into canonical files. The merge now includes:

- `saved`
- `converged`
- per-point `timings`.

### 3.4 Derived State Solves (Stage 2c)

`execute_derived_state_requests(...)` handles derived requests:

- evaluate dependency gate from final linear-state readiness,
- execute ready derived requests (subgroup mode when available),
- fall back to deterministic serial lane execution if needed.

Each derived request records:

- success/failure,
- wall/cpu timing.

These are aggregated into `metadata["derived_state_planner"]["execution"]`, including `request_timings`.

## 4) Stage 3: Property Assembly

After state solving is complete and validated at final protocol, properties are assembled in the property stage using solved state archives and metadata.

Outputs land in task property JSON (`mad.calc_info.json` task 1 properties), including:

- `response_properties` list (alpha/beta/raman records),
- `raman_spectra`,
- `vibrational_analysis`.

## 5) Subworld / Macrotask Mechanics

Subgroup execution uses `MacroTaskQ::create_worlds(...)`:

- it partitions the universe communicator into `N` subgroup communicators,
- each rank gets a subgroup id (`rank % N`),
- subgroup collectives are independent,
- each subgroup can run its owned lane(s) concurrently.

Important design point:

- archives and communicator shape must remain consistent when reading/writing response vectors in a given lane strategy.

## 6) Where to Read First

For a concrete trace, read in this order:

1. `src/madness/chem/MolresponseLib.hpp`
2. `src/apps/molresponse_v2/StateParallelPlanner.hpp`
3. `src/apps/molresponse_v2/FrequencyLoop.cpp`
4. `src/apps/molresponse_v2/DerivedStatePlanner.hpp`
5. `src/apps/molresponse_v2/ResponseRecord.hpp`
6. `src/apps/molresponse_v2/MOLRESPONSE_FUNCTION_REFERENCE.md` (function intent/similarity index)

## 7) Current Refactor Status

Completed:

- state/point ownership planner extraction,
- serial/subgroup linear solve extraction,
- derived dependency gate + execution extraction,
- stage metadata assembly extraction,
- subgroup console redirection,
- per-point linear timing persistence,
- per-request derived (VBC) timing persistence,
- frequency canonicalization for duplicate-scheduling robustness.

Active validation focus:

- serial vs state-parallel numeric parity for alpha/beta/raman on restart-driven workflows,
- reproducibility and reduction-order sensitivity in property assembly.
