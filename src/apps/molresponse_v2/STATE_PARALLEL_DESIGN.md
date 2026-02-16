# Molresponse State-Parallel Macrotask Design

## Goal

Reduce wall time for response-property workloads by solving independent linear-response states concurrently on independent MPI processor groups (subworlds), while preserving existing numerical behavior.

This targets workloads where many states are independent and memory pressure per subworld remains acceptable.

## Current Constraints (from code)

1. `solve_all_states` in `src/madness/chem/MolresponseLib.hpp` has grown large and now mixes planning, runtime policy, subgroup orchestration, and metadata merge in one function body.
2. `ResponseRecord2` writes a shared JSON file with communicator collectives, so one global writer from multiple groups is unsafe.
3. `ResponseDebugLogger` currently writes a single log file (`response_log.json`), also unsafe for multiple groups.
4. Response vectors are stored with `ParallelOutputArchive` in `src/apps/molresponse_v2/ResponseIO.hpp`; load/save behavior depends on communicator shape.
5. Frequency continuation in `computeFrequencyLoop(...)` is protocol-local. At `ti>0`, each `(state,frequency)` point can restart from its own `ti-1` data, so post-`ti0` point fanout is valid.

## Proposed Execution Model

### 1) Add a state-parallel mode

New response knobs (proposed):

- `response.state_parallel = off|auto|on`
- `response.state_parallel_groups = <int>`
- `response.state_parallel_min_states = <int>`
- `response.state_parallel_property_group = <int>` (default `0`)

Selection logic:

- `off`: current serial behavior.
- `auto`: enable only when `n_states >= min_states` and `groups > 1`.
- `on`: force grouped solving (with safety checks and explicit fallback if invalid).

### 2) Split universe into subworlds

Use `MacroTaskQ::create_worlds(universe, groups)` from `src/madness/mra/macrotaskq.h`.

- `subgroup_id = universe.rank() % groups`
- each subgroup has local communicator and independent collectives

### 3) Deterministic state ownership

Assign each state once:

- static round-robin: `owner = state_index % groups`

Each subgroup only solves states it owns.
During point mode (`protocol_index >= state_parallel_point_start_protocol`), ownership switches to deterministic point ownership:

- `point_owner = linear_point_index % effective_point_groups`
- `effective_point_groups = min(mapping_groups, num_points)`

This avoids idle point lanes when requested groups exceed available `(state,frequency)` points.

### 4) Per-group ground/response contexts

Each subgroup builds its own runtime context:

- local `GroundStateData`
- local `ResponseManager`
- subgroup-specific fock cache filename to avoid file races, e.g. `prefix.fock.group<gid>.json`

### 5) Per-group metadata and debug logs

Use sharded files:

- `response_metadata.group<gid>.json`
- `response_log.group<gid>.json`

Then merge to canonical outputs on universe rank 0:

- merged `response_metadata.json`
- merged `response_log.json`

### 6) Property phase compatibility

To avoid communicator-mismatch issues when loading archived response vectors, run property assembly on one designated subgroup (`state_parallel_property_group`, default 0) using that subgroup communicator.

Then serialize property results/metadata and broadcast to universe rank 0 for final reporting.

## Data Flow Sketch

1. Generate all states once (same as now).
2. Build state ownership and point ownership lane counts in planner metadata.
3. For each protocol: run state ownership before point-start threshold, then point ownership after threshold.
4. On restart, if all protocol-0 points are already saved, promote point ownership to protocol index 0.
5. Each subgroup writes shard metadata/log.
6. Universe rank 0 merges shards.
7. Property group computes alpha/beta/raman.
8. Final JSON outputs are emitted as today.

## Failure Handling / Fallback

Force fallback to serial mode when:

- requested groups <= 1
- requested groups > `world.size()`
- state count too small (`auto` mode)
- subgroup setup fails

When fallback occurs, print an explicit reason on rank 0.

## Implementation Plan (PR-sized)

### PR1: Infrastructure and safe scaffolding

- Add new `ResponseParameters` knobs.
- Add group planner utility (state ownership, validation).
- Add metadata/log shard merge helpers.
- Keep solve path serial by default (`off`).

### PR2: Parallel state solve path

- Add subgroup execution branch in `solve_all_states`.
- Add subgroup-local context setup and shard outputs.
- Keep property phase serial (existing behavior) initially.

### PR3: Property phase on designated subgroup

- Run `compute_alpha`, `compute_hyperpolarizability`, `compute_hessian` on property subgroup.
- Broadcast final property JSON to universe rank 0.
- Validate parity against serial mode.

## Validation Plan

### Correctness

1. Compare serial vs grouped `response_metadata.json` for saved/converged coverage.
2. Compare alpha/beta/raman values within tolerances on existing molresponse tests.
3. Verify restart behavior (partially solved states) still works.

### Performance

On Seawulf HBM node (1 node, 8 tasks, NUMA mapping):

```bash
salloc -p hbm-short-96core --nodes=1 --ntasks-per-node=8 --time=01:00:00

MAD_NUM_THREADS=10 \
mpirun --map-by numa numactl --preferred-many=8-15 \
  ./build/src/apps/madqc_v2/madqc \
  --wf=response \
  --input=src/apps/madqc_v2/test_molresponse_h2o_alpha_beta_z.in
```

Run A/B timing for:

- serial (`state_parallel=off`)
- grouped (`state_parallel=on`, `state_parallel_groups=8`)

Track wall time and peak memory per rank.

## Notes

- This design keeps the current solver internals (`computeFrequencyLoop` + response vector restart logic) unchanged.
- Major risk is archive communicator compatibility across solve and property phases; the property-group step is meant to isolate this safely.
- Ongoing readability refactor: continue splitting `solve_all_states` into focused helpers (runtime ownership policy, local workset construction, serial/subgroup execution blocks) without changing numerical behavior.
- Current checkpoint: serial and subgroup state-solve loops are extracted into dedicated helper routines so `solve_all_states` primarily orchestrates planning, execution mode selection, and downstream derived/property stages.
- Current checkpoint: derived-state dependency-gate + execution path is also extracted into a helper, reducing `solve_all_states` responsibility to stage orchestration and metadata assembly.
- Current checkpoint: final-protocol validation and stage-2 metadata assembly are extracted into helpers, further narrowing `solve_all_states` to high-level orchestration.
- Current checkpoint: stage-2 setup/logging and runtime scheduling policy construction are extracted into a dedicated schedule-context builder to keep ownership policy logic centralized.
