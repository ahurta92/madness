# Molresponse State-Parallel Macrotask Design

## Goal

Reduce wall time for response-property workloads by solving independent linear-response channels concurrently on independent MPI processor groups (subworlds), while preserving existing numerical behavior.

This targets workloads where many channels are independent and memory pressure per subworld remains acceptable.

## Terminology

- Perturbation channel: one perturbation descriptor (e.g. `Dipole_x`).
- Frequency series: ordered frequencies solved for one perturbation channel.
- Channel point: one `(channel, frequency, protocol)` solve target.

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

### 3) Deterministic channel ownership

Assign each channel once:

- static round-robin: `owner = channel_index % groups`
- protocol-0 warmup owner lanes: `channel_owner_groups = min(groups, num_channels)`
- static round-robin during warmup: `owner = channel_index % channel_owner_groups`

Each subgroup only solves channels it owns.
During point mode (`protocol_index >= state_parallel_point_start_protocol`), ownership switches to deterministic point ownership:

- `point_owner = linear_point_index % effective_point_groups`
- `effective_point_groups = min(mapping_groups, num_points)`

This avoids idle point lanes when requested groups exceed available `(channel,frequency)` points.
It also allows runs with `groups > channels` to keep extra groups idle at protocol-0 and then activate them for point ownership after protocol-0 restart data is available.

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

1. Generate all channels once (same as now).
2. Build channel ownership and point ownership lane counts in planner metadata.
3. For each protocol: run channel-series ownership before point-start threshold, then point ownership after threshold.
4. On restart, if all protocol-0 points are already saved, promote point ownership to protocol index 0.
5. Each subgroup writes shard metadata/log.
6. Universe rank 0 merges shards.
7. Property group computes alpha/beta/raman.
8. Final JSON outputs are emitted as today.

## Current Allocation Model (By Stage)

The three execution stages do not yet use the exact same allocation strategy:

1. Stage 2b (linear channel solves):
- Protocol-aware deterministic ownership.
- Protocol 0 favors channel-series ownership.
- Later protocols can switch to channel-point ownership.
- Restart metadata controls pending manifests and runtime mode gates.

2. Stage 2c (derived/VBC requests):
- Deterministic round-robin ownership over ready derived requests.
- Owner map is static per ready-request list (`owner = ready_index % groups`).
- No dynamic request stealing yet in this stage.

3. Stage 3 (property components + downstream assembly):
- Two-phase model:
- Phase A: distributed component precompute (beta/raman components) across all subgroups.
- Phase A uses dynamic task claiming via per-task lock files so subgroups can grab remaining component work.
- Phase B: property-group downstream assembly/final reporting from merged component tables.

So the short answer is: no, all three stages are not allocated in the same way yet.
Current direction is intentional: Stage 3 now has dynamic component work pickup, while Stage 2b/2c are still deterministic-owner policies.

## Offline Validation (Queue Down)

When queues are unavailable, validate component scheduling artifacts offline:

```bash
python3 src/apps/molresponse_v2/validate_property_component_claims.py \
  --run-dir <run_root> \
  --strict
```

What it checks:
- claim lock coverage: each claimed component task maps to at least one output row,
- partial outputs: claimed tasks missing some inferred A-components,
- duplicate component rows across subgroup shards,
- optional duplicate VBC writes across subgroup console logs.

## Failure Handling / Fallback

Force fallback to serial mode when:

- requested groups <= 1
- requested groups > `world.size()`
- channel count too small (`auto` mode)
- subgroup setup fails

When fallback occurs, print an explicit reason on rank 0.

## Implementation Plan (PR-sized)

### PR1: Infrastructure and safe scaffolding

- Add new `ResponseParameters` knobs.
- Add group planner utility (channel ownership, validation).
- Add metadata/log shard merge helpers.
- Keep solve path serial by default (`off`).

### PR2: Parallel channel solve path

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
3. Verify restart behavior (partially solved channels) still works.

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
- Current checkpoint: state-parallel subgroup console logs are isolated into per-group files (`response_console.group<gid>.log`) to reduce output contention and improve traceability.
- Current checkpoint: per-point linear solve timings (`wall_seconds`, `cpu_seconds`) are persisted in `response_metadata` and merged from subgroup shards into canonical metadata.
- Current checkpoint: per-point restart provenance (`kind`, source protocol/frequency, disk/promotion flags) is persisted in `response_metadata` and merged from subgroup shards.
- Current checkpoint: derived-request timing aggregation is captured in `derived_state_planner.execution.request_timings` for both subgroup and serial fallback lanes.
- Current checkpoint: serial vs parallel parity checks on H2O indicate alpha is close, while beta/Raman still show non-trivial drift; next work is reproducibility isolation (parallel-vs-parallel, then stage-by-stage divergence localization).
- Testing status (February 20, 2026): in progress. Queue-backed end-to-end
  validation for the newest property-component scheduling changes is still
  pending.

## Refactor TODO Backlog

- Extract a dedicated stage-2 orchestrator component from `MolresponseLib.hpp` to reduce `solve_all_states` coupling.
- Unify serial/subgroup scheduling loops around a shared work-item generator to remove duplicated protocol/frequency ownership logic.
- Split `computeFrequencyLoop(...)` into focused helpers (`should_solve`, `initial_guess`, `solve_point`, `persist_point`) and remove `goto`-based flow.
- Split `StateParallelPlan` into smaller types (`planner config`, `planner decision`, `runtime adjustment`) to separate static policy from restart/runtime overrides.
- Continue trimming direct logging from solver internals by routing through structured metadata/debug sinks where possible.
