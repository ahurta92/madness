# madqc Interface and Timing Design

## Current Interface: molresponse ↔ madqc

### Entry Points

There are two ways to run response:

1. **`madqc --wf=response`** — the canonical entry point. madqc
   orchestrates the full workflow (SCF → response → properties).

2. **`molresponse2`** — compatibility wrapper that dispatches to the
   same code path as madqc.

Both converge to `molresponse_lib::run_response(...)` which runs the
three-stage pipeline and returns a `Results` struct containing:
- `metadata` — response state solve metadata (convergence, timings, restart)
- `properties` — computed property values (alpha, beta, raman)
- `vibrational_analysis` — hessian / normal modes
- `raman_spectra` — per-mode Raman intensities
- `debug_log` — iteration-level diagnostics

This gets written to `<prefix>.calc_info.json` by the madqc workflow
driver.

### What madqc Provides to molresponse

- `CalculationParameters` — protocol thresholds, prefix, DFT settings
- `ResponseParameters` — perturbation types, frequencies, directions,
  requested properties, state-parallel settings, excited-state knobs
- `shared_ptr<SCF>` — SCF handle for ground-state checkpoint access
- `outdir` — working directory for response outputs

### What molresponse Returns to madqc

The `Results` struct, which madqc serializes into the calc_info JSON.
This is the contract between the two.

---

## Current Timing: What's Recorded and Where

### Per-Point Timing (response_metadata.json)

Each linear response point `(state, protocol, frequency)` records:
- `wall_seconds` — wall time for the full point solve
- `cpu_seconds` — CPU time for the full point solve

Stored in:
```
states/<perturbation>/protocols/<protocol_key>/timings/<freq_key>
```

This includes guess loading, iteration, and save — the full cost
of producing one converged response vector.

### Per-Iteration Timing (response_log.json / debug logger)

Within each point solve, the debug logger records per-iteration:
- Step-level wall/cpu timings for individual operations
  (density, coulomb, exchange, BSH, KAIN, etc.)
- Per-iteration residual, density change, alpha values

Controlled by `TimedValueLogger` and `ResponseDebugLogger`.

### Excited-State Bundle Timing

Per-protocol excited-state solve records:
- `wall_seconds`, `cpu_seconds` per protocol level
- Stored in `excited_states/protocols/<protocol_key>/timings`

### Derived-State Timing

Per derived request records:
- `wall_seconds`, `cpu_seconds` per VBC request
- Stored in `derived_state_planner/execution/request_timings`

### Stage-Level Timing

Currently NOT explicitly recorded. The overall wall time is implicit
from the process start/end, but there's no structured record of:
- How long Stage 1 (planning) took
- How long Stage 2 (solving) took total
- How long Stage 3 (property assembly) took
- How long each protocol step took as a whole

---

## What's Missing for Proper Comparison

To compare legacy code, current serial, and state-parallel modes,
you need timing at three levels that are currently incomplete:

### Level 1: Stage-Level Timing

**Need:** Wall/CPU time for each of the three stages, plus per-protocol
totals within Stage 2. This tells you the shape of the run — is most
time in solving or in property assembly? Which protocol is expensive?

**Currently:** Not recorded. You'd have to grep timestamps from stdout.

**Proposed:**
```json
{
  "timing": {
    "stage_1_planning": { "wall_s": 2.1, "cpu_s": 1.8 },
    "stage_2_solving": {
      "total": { "wall_s": 1200.0, "cpu_s": 1150.0 },
      "protocols": [
        { "threshold": 1e-4, "wall_s": 300.0, "cpu_s": 280.0 },
        { "threshold": 1e-6, "wall_s": 900.0, "cpu_s": 870.0 }
      ],
      "excited_bundle": { "wall_s": 50.0, "cpu_s": 48.0 },
      "derived_states": { "wall_s": 30.0, "cpu_s": 28.0 }
    },
    "stage_3_properties": { "wall_s": 15.0, "cpu_s": 14.0 },
    "total": { "wall_s": 1217.1, "cpu_s": 1165.8 }
  }
}
```

### Level 2: Per-Point Timing (already exists, needs standardization)

**Currently:** Exists in response_metadata.json per point. Good.

**Improvement needed:** Make sure legacy code and the new serial
orchestrator produce timing in the same format so comparison is
automatic. The legacy code probably doesn't output structured JSON,
so either add it to legacy or build a parser for legacy stdout.

### Level 3: Per-Operation Timing (debug logger, already exists)

**Currently:** Exists in response_log.json per iteration. Good for
diagnosing individual point performance.

**Improvement needed:** Consistent operation names across legacy and
current code. If legacy calls it "compute_gamma" and current calls
it "coulomb_potential", the comparison breaks. Standardize the names
in the building block catalog from the type system design.

---

## Comparison Framework Design

### What We Want to Compare

| Comparison                    | What it tells you                       |
|-------------------------------|-----------------------------------------|
| Legacy vs new serial          | Did the refactor preserve correctness and performance? |
| New serial vs state-parallel  | What's the parallel speedup? Any numerical drift? |
| Protocol A vs protocol B      | How does protocol affect accuracy and cost? |
| Molecule A vs molecule B      | How does system size affect scaling? |

### Proposed Timing Report Structure

After each run, produce a standardized timing summary (in calc_info.json
or a separate timing.json) that captures:

```json
{
  "run_id": "h2o_alpha_serial_2026-03-25",
  "code_version": "current_serial | legacy | state_parallel",
  "molecule": "H2O",
  "properties": ["polarizability"],
  "parallelism": {
    "mode": "serial | state_parallel",
    "num_groups": 1,
    "total_ranks": 8,
    "threads_per_rank": 10
  },
  "protocol_schedule": [1e-4, 1e-6],
  "timing": {
    "stage_1_planning_wall_s": 2.1,
    "stage_2_solving_wall_s": 1200.0,
    "stage_3_properties_wall_s": 15.0,
    "total_wall_s": 1217.1,
    "per_protocol": [
      {
        "threshold": 1e-4,
        "wall_s": 300.0,
        "points_solved": 3,
        "avg_point_wall_s": 100.0
      },
      {
        "threshold": 1e-6,
        "wall_s": 900.0,
        "points_solved": 3,
        "avg_point_wall_s": 300.0
      }
    ]
  },
  "results": {
    "alpha_xx_static": -5.123456,
    "converged_points": 6,
    "total_points": 6
  }
}
```

### Gecko Integration

Gecko already compares MADNESS results against Gaussian basis set
calculations. The timing report should be designed so Gecko can
also consume it for performance comparison — not just accuracy.

This means Gecko needs to:
1. Read the timing summary from calc_info.json
2. Compare timing across runs (same molecule, different modes)
3. Compare properties across runs (same molecule, different modes)
4. Produce a combined report: "serial took X seconds and got alpha=Y;
   parallel took X/4 seconds and got alpha=Y±epsilon"

### Comparison Workflow (Seawulf Skill)

This is a natural addition to the skills database:

```
Skill: benchmark_serial_vs_parallel
Trigger: Want to validate state-parallel for a given molecule/property
Inputs: molecule, property type, protocol schedule, num_groups
Steps:
  1. Run madqc --wf=response in serial mode → timing_serial.json
  2. Run madqc --wf=response in state_parallel mode → timing_parallel.json
  3. Compare properties (must match within tolerance)
  4. Compare timing (report speedup)
  5. Optionally run Gecko for basis-set comparison
Output: Comparison report with pass/fail on correctness, speedup metrics
```

```
Skill: benchmark_legacy_vs_current
Trigger: Want to verify refactored code matches legacy results
Inputs: molecule, property type, protocol schedule
Steps:
  1. Run legacy code → legacy_results
  2. Run current serial code → current_results
  3. Compare properties (must match within tolerance)
  4. Compare per-point timing (identify regressions)
Output: Parity report
```

---

## Simplifying the madqc Interface

### Current Issue

The `run_response` function currently receives `CalculationParameters`,
`ResponseParameters`, `shared_ptr<SCF>`, and `outdir` as separate
arguments. The function itself is ~50 lines of setup (parsing response
input, writing response.in, re-reading it) before it even starts the
pipeline.

### Proposed Simplification

For the v3 serial orchestrator (Increment 5), the interface should be
cleaner:

```cpp
struct ResponseWorkflowInput {
    // Everything the response pipeline needs, pre-validated
    Molecule molecule;
    std::string archive_file;       // path to SCF checkpoint
    std::vector<double> protocols;  // truncation threshold schedule
    ResponseParameters response_params;  // perturbations, frequencies, etc.
};

struct ResponseWorkflowOutput {
    // Everything the response pipeline produces
    nlohmann::json properties;      // alpha, beta, raman values
    nlohmann::json metadata;        // state solve status, convergence
    nlohmann::json timing;          // structured timing at all levels
    nlohmann::json debug_log;       // iteration-level diagnostics
};

// Clean entry point
ResponseWorkflowOutput run_response(World &world,
                                     const ResponseWorkflowInput &input);
```

The madqc workflow driver is responsible for constructing the input
from its own parameters and SCF results. The response pipeline doesn't
know about madqc, SCF objects, or workflow builders — it receives a
self-contained input and returns a self-contained output.

This also makes it trivial to run the same pipeline from a test
harness, from legacy comparison scripts, or from a future Python
binding — you just construct the input differently.

### Timing Integration

The timing JSON is part of the output, produced by the orchestrator.
The orchestrator wraps each stage in wall/cpu timers:

```cpp
ResponseWorkflowOutput run_response(World &world,
                                     const ResponseWorkflowInput &input) {
    ResponseWorkflowOutput output;
    Timer total_timer;

    // Stage 1
    Timer stage1_timer;
    auto ground = load_ground(world, input);
    auto states = plan_states(world, input, ground);
    output.timing["stage_1_wall_s"] = stage1_timer.elapsed();

    // Stage 2
    Timer stage2_timer;
    solve_all(world, ground, states, input.protocols, output);
    output.timing["stage_2_wall_s"] = stage2_timer.elapsed();

    // Stage 3
    Timer stage3_timer;
    output.properties = compute_properties(world, ground, states, input);
    output.timing["stage_3_wall_s"] = stage3_timer.elapsed();

    output.timing["total_wall_s"] = total_timer.elapsed();
    return output;
}
```

Clean, readable, and every timing you need for comparison falls out
naturally.
