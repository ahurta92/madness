# MolresponseLib Developer Guide

`src/madness/chem/MolresponseLib.hpp` is the orchestration header for the
`molresponse_v2` pipeline. Everything is a `static` method or nested type inside
`struct molresponse_lib`. This guide walks through each logical section and explains
what it does, what data it owns, and how it connects to the sections before and after.

---

## Quick navigation

The file is divided into 12 labelled sections (search for `SECTION` to jump):

| # | Section banner | What's inside |
|---|----------------|---------------|
| 1 | Context types | `Results`, `GroundContext`, planning bundles |
| 2 | State persistence | `JsonStateSolvePersistence` |
| 3 | Shard filenames | `group_shard_file`, `group_console_file` |
| 4 | Manifest/claim helpers | derived-state done-file tracking |
| 5 | Log capture | `FilteredLineStreambuf`, `ScopedRankLogRedirect` |
| 6 | JSON utils | file I/O, `broadcast_json`, `merge_state_metadata_json`, `with_subworld` |
| 7 | Stage 1 — Planning | ground context, state generation, schedule planning |
| 8 | Stage 2 schedule types | `StateSolveScheduleContext`, manifests, ownership policy |
| 9 | Stage 2b — Linear solve | serial and subgroup execution paths |
| 10 | Stage 2d — Derived solve | VBC-driven derived-state execution |
| 11 | Stage 2c — Excited bundle | excited-state protocol adapter loop |
| 12 | Stage 3 — Property stage | polarizability, hyperpolarizability, Raman |
| 13 | Top-level orchestration | `solve_all_states`, `run_response` |

---

## The three-stage pipeline

```
run_response(world, calc_params, response_params, scf_calc)
    │
    ├── Stage 1: make_ground_context + plan_required_states
    │     └── Returns: PlannedStates (linear + derived + excited + parallel plan)
    │
    ├── Stage 2: solve_all_states
    │     ├── 2a: serial OR subgroup linear state solve
    │     │     └── calls execute_serial_state_solve / execute_subgroup_state_solve
    │     ├── 2c: execute_excited_state_bundle_stage (per-protocol excited solve)
    │     └── 2d: execute_derived_state_requests (VBC hyperpolarizability prereqs)
    │           └── Returns: SolvedStates (PlannedStates + merged metadata + debug log)
    │
    └── Stage 3: compute_requested_properties_with_property_group
          ├── compute_polarizability
          ├── compute_hyperpolarizability
          └── compute_raman → print_raman_table
                └── Returns: PropertyStageOutput → Results
```

---

## Section-by-section breakdown

### Section 1 — Context types

These are the **value-type bundles** that flow through the pipeline.

```
GroundContext       — molecule + ground-state orbitals + archive filename
ExcitedStateBundlePlan — excited-state protocol parameters (enabled?, num_states, TDA, …)
PlannedStates       — generated linear states + derived plan + excited plan + parallel plan
SolvedStates        — PlannedStates + merged metadata JSON + debug log JSON
PropertyStageOutput — assembled properties JSON + vibrational + Raman results
Results             — public output: metadata, properties, debug_log, vib, raman
```

`Results` is the only **public** type — everything else is private to `molresponse_lib`.

---

### Section 2 — State persistence (`JsonStateSolvePersistence`)

Implements the `StateSolvePersistence` interface consumed by `computeFrequencyLoop()`
inside `FrequencyLoop.cpp`. It wraps two backends:

- **`ResponseRecord2`** (from `ResponseRecord.hpp`) — the JSON metadata file
  (`response_metadata.json`). Tracks `saved`/`converged` per `(state, protocol, freq)`.
- **`ResponseDebugLogger`** — writes `response_log.json` with per-iteration diagnostics.

The subgroup path creates one `JsonStateSolvePersistence` **per subgroup** pointing at
shard files (`response_metadata.group<N>.json`). Shards are merged back to the primary
file at protocol sync points.

Key methods you'll use:
```cpp
persistence.is_saved(pt)       // restart check
persistence.is_converged(pt)   // convergence gate
persistence.record_status(pt, converged)
persistence.record_timing(pt, wall, cpu)
persistence.metadata_json()    // get a snapshot for merging
```

---

### Section 3 — Shard filename generators

Simple helpers used to compute shard file paths for per-subgroup metadata, logs, and
derived-timing records. Nothing surprising here.

```cpp
group_shard_file("response_metadata.json", 2) → "response_metadata.group2.json"
group_console_file(1)                          → "response_console.group1.log"
```

---

### Section 4 — Derived-state manifest/claim helpers

The derived-state execution (Section 10) uses **file-based task claiming** so that
multiple subgroups don't duplicate work. Key files:

```
derived_request_manifest/<scope>/req<N>_<id>.done.json   ← written when request complete
<claim_prefix>.req<N>.lock                                ← claimed before execution
```

`derived_request_done_record_exists(...)` reads the `.done.json` and checks `"success": true`.
`write_derived_request_done_record(...)` atomically writes the completion record.

---

### Section 5 — Log capture

`ScopedRankLogRedirect` redirects `std::cout` and `std::cerr` to a per-subgroup file
for the duration of the subgroup solve. This prevents interleaved output from different
MPI groups on the same console.

`FilteredLineStreambuf` drops lines containing `"!!MADNESS: Hung queue?"` to avoid
log spam.

---

### Section 6 — JSON utilities

Small but heavily used. Functions you'll call frequently:

```cpp
// File I/O (rank-0 only, no MPI sync)
write_json_file(filename, json_data);
auto j = read_json_file_or_object(filename);   // returns {} on missing/corrupt

// Broadcast: modifies obj in-place on all ranks
broadcast_json(world, obj);                    // NEW in refactor

// Legacy: takes by value, returns broadcast result
auto j = broadcast_json_object(world, std::move(payload));

// Metadata merge (used after subgroup shards are collected)
merge_state_metadata_json(merged, shard);
merge_debug_log_json(merged_log, shard_log);

// Point readiness check (queries metadata JSON directly)
bool ok = point_ready_in_metadata(metadata, pt, require_saved, require_converged);

// Subworld RAII (added in refactor — not yet used at call sites)
with_subworld(world, ngroups, [&](World &sw) { ... });
```

**`broadcast_json` vs `broadcast_json_object`:** prefer `broadcast_json` for in-place
mutation; use `broadcast_json_object` only when you need a fresh copy from a specific
rank.

---

### Section 7 — Stage 1: Planning

Entry: `plan_required_states(world, calc_params, response_params)` → `PlannedStates`.

Steps inside Stage 1:
1. `make_ground_context(world, calc_params)` — loads archive, builds `GroundStateData`
2. `StateGenerator` — converts response input knobs to `GeneratedStateData`
   (one `LinearResponseDescriptor` per perturbation channel)
3. `DerivedStatePlanner` — inspects requested beta triplets, emits `DerivedStatePlan`
4. `build_excited_state_bundle_plan(...)` — reads `response.excited.*` params
5. `StateParallelPlanner::build_plan(...)` — decides serial vs subgroup, assigns channel
   ownership, computes `point_parallel_start_protocol_index`

The resulting `PlannedStates` is immutable from this point — later stages only read it.

---

### Section 8 — Stage 2 schedule context types

`StateSolveScheduleContext` aggregates everything the linear-solve loop needs:
- `linear_states` — the full state list
- `owner_by_channel_index` — which subgroup owns each channel
- `point_scheduler` — per-channel frequency-point ownership at fine-grained protocols
- `force_retry_removed_frequencies` — retry policy for failed points

`PendingProtocolManifest` / `PendingPointWorkItem` — the work list for one protocol step.
Built by `build_pending_manifest_from_metadata(...)`.

`build_state_solve_schedule_context(...)` is the factory; it reads the metadata file,
broadcasts it, and constructs the context used throughout Stage 2.

---

### Section 9 — Stage 2b: Linear state solve

Two execution paths selected by `StateParallelPlan::effective_mode`:

#### Serial path (`execute_serial_state_solve`)
One world, one `JsonStateSolvePersistence` per protocol step, all channels solved in
sequence. Simple but doesn't exploit parallelism across states.

#### Subgroup path (`execute_subgroup_state_solve`)
`MacroTaskQ::create_worlds(world, N)` splits `world` into N subgroups. Each subgroup:
1. Redirects its console to `response_console.group<N>.log`
2. Creates its own `JsonStateSolvePersistence` writing to shard files
3. Solves its assigned channels via `execute_manifest_work`
4. After the final protocol, participates in the **tail derived poll** loop
   (currently marked `// TODO(step4)` for extraction)

After each subgroup finishes a protocol, rank 0 **merges shards** into the primary
`response_metadata.json` and **broadcasts** the merged result before the next protocol.

**Inside `execute_manifest_work`:** calls `computeFrequencyLoop()` for each pending
`(channel, freq, protocol)` work item. Restart precedence per point:
1. Exact checkpoint for this protocol
2. Coarser-protocol snapshot (static→dynamic promotion)
3. Frequency continuation from nearest solved frequency
4. Fresh init

---

### Section 10 — Stage 2d: Derived state

`execute_derived_state_requests(world, plan, states, final_protocol, final_threshold,
                                 ground, response_manager, vbc_computer)`

Evaluates `DerivedStatePlanner::evaluate_dependency_gate(...)` to find requests whose
linear prerequisites are all converged, then executes ready requests via
`try_claim_property_component_task(...)` (file-lock based, safe for multi-subgroup).

`run_derived_request(...)` wraps a `SimpleVBCComputer::compute_vbc(...)` call and writes
the `.done.json` completion record.

In the subgroup path, this runs as the **tail poll loop** inside
`execute_subgroup_state_solve` while other subgroups may still be finishing their last
linear points.

---

### Section 11 — Stage 2c: Excited-state bundle

`execute_excited_state_bundle_stage(world, planned_states, ground_ctx, response_params,
                                     state_metadata_json)` → `ExcitedExecutionResult`

Iterates over `ExcitedStateBundlePlan::protocols` in order. For each protocol:
1. Reads existing metadata node, calls `ensure_excited_protocol_placeholder_node` to fill
   missing fields (now a 12-line merge loop using `ExcitedProtocolResult::to_json`)
2. Checks `restart_ready`: if saved + converged + bundle_state_present + restart_capable,
   sets `solver_needed = false` and skips the solver call
3. Otherwise calls `ExcitedResponse::solve_protocol(world, input)` → `ExcitedProtocolResult`
4. Broadcasts result via `broadcast_json`, normalises via
   `ResponseRecord2::normalize_excited_protocol_result`, records via
   `excited_metadata_record.record_excited_protocol_result`

`ExcitedExecutionResult` carries an `execution` JSON summary written to
`response_metadata.json["excited_state_planner"]`.

**Phase status:** Phases 1–3 complete (descriptors, archive/restart, solver init).
Phase 4 (iteration/convergence loop) is in progress in `ExcitedResponse`.

---

### Section 12 — Stage 3: Property assembly

`compute_requested_properties_with_property_group(...)` → `PropertyStageOutput`

Iterates over `PropertyType` (Polarizability, Hyperpolarizability, Raman) based on
what was requested in `response_params`. Optionally runs in a dedicated property subgroup
(`property_group` parameter) using `MacroTaskQ::create_worlds`.

#### `compute_polarizability`
Assembles the α tensor from converged linear response functions. Writes to
`properties.json["alpha"]`.

#### `compute_hyperpolarizability`
Assembles the β tensor from derived (VBC) states. Writes to `properties.json["beta"]`.
Handles all triplet modes: SHG, OR, all_triplets.

#### `compute_raman`
1. Calls `compute_hessian(...)` → `VibrationalResults` (normal modes, frequencies)
2. Calls `compute_Raman(...)` → α derivatives along nuclear displacement modes
3. Computes α² and β² invariants per mode
4. **`print_raman_table(world, raman, print_level)`** (extracted in refactor) —
   formats and prints the per-frequency Raman intensity table to stdout on rank 0

---

### Section 13 — Top-level orchestration

`run_response(world, calc_params, response_params, scf_calc)` → `Results`

The single public entry point called by `WorkflowBuilders.hpp`. Sequence:
```
Stage 1: plan_required_states(...)
Stage 2: solve_all_states(...)
Stage 3: compute_requested_properties_with_property_group(...)
         → aggregate into Results and return
```

`prepare_and_validate_final_protocol_state(...)` is called inside `solve_all_states`
after linear solves complete; it enforces the property gate
(`ResponseRecord2::enforce_ready_for_properties`).

`build_state_stage_metadata(...)` merges subgroup shards into the canonical
`response_metadata.json` at stage boundaries.

---

## JSON schema quick reference

All metadata flows through `response_metadata.json`.

```json
{
  "states": {
    "<state_id>": {
      "protocols": {
        "<1e-04>": {
          "saved":    { "<0.000>": true },
          "converged":{ "<0.000>": true },
          "timings":  { "<0.000>": { "wall_seconds": 1.2, "cpu_seconds": 1.1 } },
          "solver_diagnostics": { ... },
          "restart_provenance":  { ... }
        }
      }
    }
  },
  "excited_states": {
    "plan": { "enabled": true, "num_states": 3, ... },
    "protocols": {
      "<1e-04>": {
        "saved": true, "converged": false,
        "energies": [...], "roots": [...], "stage_status": "...",
        ...
      }
    }
  },
  "excited_state_planner": { "plan": {...}, "execution": {...} }
}
```

Key string formats (set by `ResponseRecord2` static helpers):
- State key: `pt.perturbationDescription()` e.g. `"Dipole_x"`
- Protocol key: `"%.0e"` e.g. `"1e-04"`, `"1e-06"`
- Frequency key: `"%.3f"` e.g. `"0.000"`, `"0.500"`

---

## Common gotchas

1. **All metadata writes go through `ResponseRecord2` or `JsonStateSolvePersistence`.**
   Never write `response_metadata.json` directly with `write_json_file` inside a live
   solve — use `persistence.record_status(pt, c)`.

2. **Subgroup shards are not merged until protocol sync.** If you read
   `response_metadata.json` during a subgroup solve you will see stale data.
   Use the shard merge pattern from `execute_subgroup_state_solve` to get live state.

3. **`broadcast_json` vs `broadcast_json_object`.**  `broadcast_json` modifies in-place
   (preferred); `broadcast_json_object` takes by value and returns a copy (legacy).

4. **`with_subworld` is available but call sites haven't been migrated yet.**
   New code using subworlds should use it; existing call sites will be migrated
   when the file is eventually split in Step 7 of the refactor.

5. **Protocol ordering is ascending threshold** (coarse → fine, e.g. 1e-4 → 1e-6).
   Functions that iterate protocols use `std::min_element` to find the final (tightest)
   protocol.

6. **Excited-state execution (Stage 2c) runs after linear solves (Stage 2a/b) but
   before derived states (Stage 2d).** The ordering in `solve_all_states` is 2a/b → 2c → 2d.

---

## Key types from other files used here

| Type | File | Role |
|------|------|------|
| `ResponseRecord2` | `ResponseRecord.hpp` | JSON metadata backend |
| `ExcitedProtocolResult` | `ResponseRecord.hpp` | Per-protocol excited-state result |
| `ExcitedRootDescriptor` | `ResponseRecord.hpp` | Stable root identity |
| `LinearResponseDescriptor` | `ResponseState.hpp` | One perturbation channel |
| `LinearResponsePoint` | `ResponseState.hpp` | One (channel, freq, protocol) point |
| `GeneratedStateData` | `StateGenerator.hpp` | Output of Stage 1 state gen |
| `DerivedStatePlan` | `DerivedStatePlanner.hpp` | VBC-driven requests |
| `StateParallelPlan` | `StateParallelPlanner.hpp` | Subgroup ownership plan |
| `ResponseManager` | `ResponseManager.hpp` | Protocol-level solver setup |
| `GroundStateData` | `GroundStateData.hpp` | Ground-state orbitals |
| `SimpleVBCComputer` | `VBCMacrotask.hpp` | VBC hyperpolarizability compute |
| `ExcitedResponse` | `ExcitedResponse.hpp` | Phase 4 excited solver |
| `PropertyManager` | `PropertyManager.hpp` | Property result storage |
