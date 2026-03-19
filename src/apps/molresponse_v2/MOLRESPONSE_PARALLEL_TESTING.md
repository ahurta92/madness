# Molresponse State-Parallel Testing

This directory tracks correctness and performance checks for state-parallel molresponse runs.

## Testing Status

As of February 20, 2026, runtime validation is still in progress.
The current implementation changes are committed, but full queue-backed
regression/performance testing remains ongoing.

## High-Level Design Intent

This test area is intended to capture both:
- behavior: how state-parallel scheduling should behave as work scales, and
- science workflow context: why molresponse needs many independent states for `alpha`, `beta`, and `raman`.

At a high level, molresponse follows three stages:
1. Parse response inputs and generate all required states/frequencies.
2. Solve missing linear/derived states (or restart from archives when available).
3. Compute requested response properties from converged states.

The current two H2O inputs are primarily stage-2 (state solve + scheduling/timing) tests.

## Protocol Design Intent

This section describes the intended scheduling behavior across protocols.

Core rule:
- A subgroup should execute one protocol phase at a time.
- Within that phase, it should process all manifest items assigned to that subgroup.
- Ground-state data should be read/initialized once per subgroup per protocol, then reused for all items in that protocol pass.

Restart rule for a `(state, frequency, protocol)` target:
- First try restart from the same protocol archive for that exact `(state, frequency)`.
- If missing, try nearest available lower protocol for that same `(state, frequency)`.
- If still missing, fall back to phase-specific guess strategy (described below).

Phase behavior:
- `phase 0` (lowest protocol): objective is to create protocol-0 seed solutions for every required `(state, frequency)`.
- In phase 0, each subgroup should own complete frequency chains for its assigned states so it can walk frequencies in order and reuse nearby-frequency guesses.
- `phase > 0`: once protocol-0 seeds exist, ownership can be relaxed to arbitrary `(state, frequency)` points because each point can restart from an existing seed/protocol archive.

Why phase 0 is different:
- The best initial guess for `(state, frequency_i)` is usually a nearby frequency of the same state.
- This favors per-state frequency continuity and linear frequency traversal in phase 0.
- More advanced restart heuristics are possible later, but this is the baseline policy.

Group-count policy intent when `groups > states` in phase 0:
- Use only a subset of groups in phase 0 (enough to cover state chains well).
- After phase 0 completes, expand to full group count for phase > 0 where independent point parallelism is higher.
- This avoids idle or poorly utilized groups during phase 0 while still exploiting maximum parallelism later.

Incremental implementation plan:
1. Enforce/verify per-protocol subgroup loop semantics and one-time ground-state read per protocol pass.
2. Formalize restart lookup order: same protocol -> nearest lower protocol -> phase fallback guess.
3. Strengthen phase-0 ownership policy (state-centric frequency chains).
4. Enable dynamic group activation: reduced active groups in phase 0, full groups in later phases.
5. Add metadata fields to report which restart source was used per solved point.
6. Add regression tests for clean start, partial restart, and mixed-protocol restart cases.

Status tracker (current):

| Step | Status | Notes |
| --- | --- | --- |
| 1 | partial | Per-protocol subgroup execution is in place; ground-state reuse semantics should be explicitly validated per subgroup/protocol pass. |
| 2 | partial | Restart-aware scheduling exists, but restart-source policy should be made explicit and traceable in code/metadata. |
| 3 | partial | Phase-0 pending work now forces state-mode with fixed state-owner mapping; restart with complete protocol-0 still promotes to point mode. |
| 4 | todo | Dynamic active-group resizing between phase 0 and later phases is not implemented yet. |
| 5 | partial | Linear-state metadata now records per-point restart provenance (`kind`, `source_protocol`, `source_frequency`, disk/promotion flags); restart-policy summaries can now be traced point-by-point. |
| 6 | partial | Manual restart and scaling checks exist; dedicated automated regression matrix is still missing. |

Code map (where to implement each step):

| Step | Primary files | Main symbols / edit hotspots |
| --- | --- | --- |
| 1 | `src/madness/chem/MolresponseLib.hpp` | `run_protocol_threshold_loop`, `execute_serial_state_solve`, `execute_subgroup_state_solve`; protocol `prepare` lambdas (one-time `computePreliminaries` intent per subgroup/protocol). |
| 2 | `src/madness/chem/MolresponseLib.hpp`, `src/apps/molresponse_v2/ResponseRecord.hpp` | `build_pending_work_manifest`, `point_needs_solving`; restart lookup helpers such as `best_converged_protocol` and protocol-key utilities. |
| 3 | `src/madness/chem/MolresponseLib.hpp`, `src/apps/molresponse_v2/StateParallelPlanner.hpp` | `use_state_ownership_for_protocol_runtime`, `dispatch_owned_work_for_protocol`, planner fields `point_parallel_start_protocol_index` and ownership policy selection. |
| 4 | `src/apps/molresponse_v2/StateParallelPlanner.hpp`, `src/madness/chem/MolresponseLib.hpp` | `StateParallelPlanner::build` (`state_owner_groups`, `execution_groups`, `effective_point_groups`) plus runtime activation in `build_state_solve_schedule_context`. |
| 5 | `src/apps/molresponse_v2/ResponseRecord.hpp`, `src/madness/chem/MolresponseLib.hpp` | Extend record schema for restart provenance; write/merge in persistence + metadata assembly (`state_parallel_runtime` block). |
| 6 | `test_molresponse/`, `src/apps/madqc_v2/`, `src/apps/znemo/` | Add restart/scaling regression drivers and comparison checks (clean start, partial restart, mixed protocol availability). |

## Property Component Parallel Status

Current property-stage behavior now follows a two-phase model:
- Phase A (component precompute): beta/raman frequency-dependent components are computed across all subgroups.
- Phase A uses dynamic claim files so idle subgroups can pick up remaining component tasks.
- Phase A writes subgroup shards: `properties_components.group<gid>.json`.
- Phase B (assembly): merged component rows are assembled into final property outputs on `state_parallel_property_group`.

Why this split is useful:
- Component contractions are independent and expensive, so they scale well with subgroup fanout.
- Final assembly is lighter and can remain centralized for stable output formatting and JSON emission.

Current implementation notes:
- Dynamic claiming is wired through `task_claim_prefix` in `compute_beta(...)` / `compute_hyperpolarizability(...)` / `compute_Raman_components(...)`.
- Component scheduling is dynamic in Phase A.
- Final property assembly still runs in one subgroup (intentional for now).

Status tracker (property parallel):

| Step | Status | Notes |
| --- | --- | --- |
| 1 | partial | Implicit work-item model exists in `(freqB,freqC,B,C)` loops; no explicit `PropertyComponentWorkItem` type yet. |
| 2 | partial | Manifest is implicit from loops + claim-index ordering; no standalone manifest object yet. |
| 3 | implemented | Multi-subgroup fanout for beta/raman components is active via dynamic task claiming. |
| 4 | implemented | Per-subgroup component shards are written as `properties_components.group<gid>.json`. |
| 5 | implemented | Global merge of component shards into `properties.json`, then downstream assembly pass. |
| 6 | partial | Component-level skip exists via `pm.has_*`; richer restart provenance for components is still pending. |
| 7 | partial | Runtime path exists; dedicated serial-vs-component-parallel regression harness still needs to be finalized. |

Code map for this future stage:
- `src/madness/chem/MolresponseLib.hpp`
- `compute_requested_properties(...)`
- `compute_requested_properties_with_property_group(...)`
- `src/apps/molresponse_v2/PropertyManager.hpp`
- `src/apps/molresponse_v2/MolecularProperty.hpp`

## Scope

Current focus:
- Validate state-parallel scheduling behavior (8 and 16 subgroups).
- Validate restart-aware point scheduling (`state_parallel_point_start_protocol 0`).
- Collect state timing data from `response_metadata*.json`.

Not yet the focus of these two inputs:
- Final property validation in `calc_info.json` (`alpha`, `beta`, `raman`) because these test inputs currently do not set `property true` / requested property list.

## Input Cases

- `tmp_h2o_alpha_beta_xyz_statepar8_f6.in`
- `tmp_h2o_alpha_beta_xyz_statepar16_f6.in`

Both use:
- H2O geometry
- `dipole.directions xyz`
- `quadratic true`
- input frequencies `[0, 0.02, 0.04, 0.06, 0.08, 0.1]`
- protocols `[1e-4, 1e-6]`
- point-level scheduling enabled from protocol 0

Important note:
- With quadratic enabled, the linear-state frequency set expands to include sum frequencies. For this case, metadata shows solved frequencies from `0.000` to `0.200` (11 points) for each dipole direction.
- This is the core reason parallel scheduling matters: property requests can expand into many independent state-frequency solves.

## Why State Counts Grow

For this style of input, requesting dipole directions `xyz` with quadratic terms means we need:
- linear dipole response states for each direction/frequency point, and
- additional frequency combinations needed by downstream quadratic workflows (for VBC and beta-like contractions).

So even with 6 user-entered frequencies, the effective linear solve set is larger (11 here, from metadata), and total work is roughly:
- `N_directions * N_effective_frequencies * N_protocols`

That scaling is exactly what this directory is trying to stress with 8-group and 16-group cases.

## Example Input Pattern

```txt
response
    kain true
    dipole true
    dipole.directions xyz
    quadratic true
    dipole.frequencies [0, 0.02, 0.04, 0.06, 0.08, 0.1]
    state_parallel on
    state_parallel_groups 16
    state_parallel_point_start_protocol 0
    state_parallel_min_states 1
    state_parallel_property_group 0
end
```

## How To Run (40-core Intel MPI setup)

```bash
cd /gpfs/projects/rjh/adrian/development/madness-worktrees/raman_branch_no_mul_sparse/test_molresponse
source ~/load_40core.sh

mpirun -machinefile hosts -n 16 \
  /gpfs/projects/rjh/adrian/development/madness-worktrees/raman_branch_no_mul_sparse/build-40core-intelmpi/src/apps/madqc_v2/madqc \
  --wf=response --input=tmp_h2o_alpha_beta_xyz_statepar16_f6.in
```

## If Queues Are Unavailable (Do This Offline)

1. Keep refactoring momentum with compile-only checks:
- `source ~/load_40core.sh`
- `cmake --build build-renamecheck --target madqc -j 12`

2. Validate scheduling logic statically from metadata/log wiring:
- Confirm stage-2 linear uses protocol-aware channel-series/channel-point manifests.
- Confirm stage-2 derived uses deterministic ready-request ownership.
- Confirm stage-3 component precompute uses dynamic claim files and subgroup shards.

3. Add/maintain documentation and regression scripts:
- Update this file and `src/apps/molresponse_v2/STATE_PARALLEL_DESIGN.md` when allocation logic changes.
- Prepare queued run commands/scripts so they can launch immediately once resources are available.

## What To Check In Logs

From `task_1/molresponse/response_console.group*.log`, verify:
- `State-parallel subgroup ...`
- `pending states=... pending points=... mode=point`
- `solving independent state-frequency points at thresh ...`

This confirms the subgroup is running point-based work (not serial fallback).

## Timing Data Collection

Primary source:
- `task_1/molresponse/response_metadata.json`

Useful summaries:

```bash
# Per-state, per-protocol timing summary
jq -r '
  .states as $S |
  ["state","protocol","count","min_wall_s","max_wall_s","avg_wall_s"],
  ($S | to_entries[] | .key as $state | .value.protocols | to_entries[] | .key as $p |
   [ .value.timings | to_entries[] | .value.wall_seconds ] as $w |
   [ $state, $p, ($w|length), ($w|min), ($w|max), (($w|add)/($w|length)) ])
  | @tsv' task_1/molresponse/response_metadata.json
```

```bash
# Quick check that all frequencies converged/saved
jq '.states' task_1/molresponse/response_metadata.json
```

## Current 16-Subgroup Run Status

Run directory:
- `tmp_h2o_alpha_beta_xyz_statepar16_f6/`

Observed behavior:
- Point-level scheduling is active in all subgroup logs.
- Protocol 0 and protocol 1 both execute in `mode=point`.
- Work is distributed as expected for 33 linear points (3 directions x 11 frequencies):
  - most groups: `pending points=2`
  - one group: `pending points=3`
- `response_metadata.json` shows all Dipole_x/y/z points converged and saved at both protocols.

Timing snapshot from metadata (wall seconds, average per point):
- Dipole_x: 1e-4 ~180.35 s, 1e-6 ~350.83 s
- Dipole_y: 1e-4 ~217.95 s, 1e-6 ~293.18 s
- Dipole_z: 1e-4 ~208.60 s, 1e-6 ~364.37 s

## Where We Are Vs Desired Design

Implemented:
- Restart-aware pending-work manifests.
- Point-level independent state-frequency scheduling from configurable protocol index.
- Fresh-run protocol-0 state-chain ownership with fixed state-owner mapping.
- Per-point restart provenance in `response_metadata.json` for solved linear points.
- Subgroup-local logging and metadata suitable for debugging/diagnostics.

Still to do:
- Add dedicated property-validation input(s) with `property true` and explicit requested properties.
- Add direct serial-vs-parallel property comparison step in this directory.
- Add a compact timing comparison table for 8 vs 16 groups in this document.
