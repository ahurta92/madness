# Review the complete contract: SCF and static response

**Behavior and record layout accepted — 2026-11-12.** The contract is in
[design sections 1–7](design_draft.md#technical-plan). This walk-through uses the
[maintained LiH deck](src/examples/qc/response_lih_beta/response_lih_beta.in):
SCF, static alpha_zz and beta_zzz. The deck's existing scientific selectors are
retained; the new convergence/restart contract is not yet implemented.
All numbers below are invented to explain decisions, not reference results,
recommended tolerances or evidence that LiH converges with these settings.

## 1. What the user requests

Keep the deck's molecule block. Use the following proposed application blocks;
`madqc --wf=response --prefix=lih example.in` retains the existing CLI form.
There are no user-named stages. Fixed k=8 applies to both protocol thresholds
in this example; response inherits the workflow's protocol settings.

```text
dft
    xc hf
    k 8
    protocol [1e-4, 1e-6]
    l 50.0
    localize new
    maxiter 30
    convergence.criteria [bsh_residual, density_change, energy_change]
    convergence.bsh_residual 1e-6
    convergence.density_change 1e-6
    convergence.energy_change 1e-7
end

response
    dipole true
    dipole.directions z
    dipole.frequencies [0.0]
    quadratic true
    kain true
    maxiter 60
    convergence.criteria [bsh_residual, density_change, energy_change]
    convergence.bsh_residual 1e-6
    convergence.density_change 1e-6
    convergence.energy_change 1e-7
end
```

Omitted relative tolerances are zero. `restart auto` and endpoint-only saving
are the defaults. Add `checkpoint.every_iterations 5` to an application group
only if interior checkpoints are wanted. There is no `save`, keep-all, warmup,
stall or separate excited-state budget setting.

## 2. What executes and who owns it

| Work | Owner and completion rule |
| --- | --- |
| SCF `ground_state/scf` | Moldft owns orbitals/history. It evaluates the current prepared orbitals, then the shared controller decides whether another update is needed. Finish both requested protocols before releasing the reference. |
| Response `response/dipole-z` | Molresponse owns amplitudes/history and its internal DAG. The z state supplies A, B and C in this static one-direction example. It follows the same measurement/decision/checkpoint contract. |
| VBC `response/vbc-zz` | Direct task after B and C meet final criteria. Record both input roles even though they use the same state generation here. Save through the same checkpoint interface. |
| Properties | Molresponse computes alpha_zz from eligible inputs and beta_zzz after VBC and A are eligible. Both appear under calculation response. |

At 1e-4, use the documented model/metric intermediate absolute adjustments,
with the update cap as a fallback. At 1e-6, use exactly the requested limits.
Neither a small coarse residual nor reaching maxiter establishes final convergence.

Suppose the final evaluations produce:

| Metric | SCF | Response z | Required final limit |
| --- | --- | --- | --- |
| BSH residual | 3e-7 | 5e-7 | 1e-6 |
| Density change | 4e-7 | 7e-7 | 1e-6 |
| Energy change | 2e-8 | Inapplicable: driving frequency is fixed | 1e-7 where applicable |
| Kinetic change | 6e-8 | 2e-7 | Diagnostic; not selected |

The limits have each metric's declared model units. SCF uses maximum per-function
BSH, maximum physical spin-density L2 change and total-energy change. Response
uses the same BSH reduction and its model's induced density; its kinetic quantity
is the positive unit-weighted gradient sum. Definition/normalization review is
still required before the corresponding adapter is implemented.

Both solver tasks return `complete / final_criteria_met`. The application saves
those evaluated functions without a further KAIN update. SCF iteration 7 follows
seven final-protocol updates; response iteration 5 follows five. These finish
before their caps. KAIN history is not saved. If convergence is first established by a terminal evaluation after
the cap, that evaluation consumes no extra update.

## 3. What appears in results.json

This excerpt shows the whole requested property set and one metric snapshot
per state. To keep it readable, each calculation's required `metric_definitions`
and `metric_set` are omitted here. The real writer includes them: the former
holds full definitions/units; the latter maps metric names to those definition
IDs. The main contract specifies their fields. The root `attempt_id` names the
run writing this snapshot. Both calculations below omit their own `attempt_id`,
so each inherits `a1` as its immutable source-map owner. Results cannot be used as
standalone restart evidence.

```json
{
  "version": 2,
  "attempt_id": "a1",
  "calculations": [
    {
      "calculation_id": "ground_state",
      "application": "moldft",
      "model": "HF",
      "protocol": {"threshold": 1e-06, "k": 8},
      "properties": {"total_energy": {"value": -7.98, "units": "hartree", "value_type": "real", "status": "complete"}},
      "state_metrics": [
        {
          "states": {
            "ground": {"values": {"bsh_residual": 3e-07, "density_change": 4e-07, "energy_change": 2e-08, "kinetic_change": 6e-08}, "availability": {}}
          }
        }
      ]
    },
    {
      "calculation_id": "response",
      "application": "molresponse",
      "model": "TDHF",
      "protocol": {"threshold": 1e-06, "k": 8},
      "frequency_units": "atomic_angular_frequency",
      "properties": {
        "polarizability": {
          "units": "bohr^3",
          "value_type": "real",
          "values": [{"frequency": 0.0, "status": "complete", "components": {"zz": 20.0}, "availability": {}}]
        },
        "hyperpolarizability": {
          "units": "e^3*bohr^3/hartree^2",
          "value_type": "real",
          "values": [
            {
              "process": "static",
              "frequencies": {"B": 0.0, "C": 0.0},
              "output_frequency": 0.0,
              "status": "complete",
              "components": {"zzz": 300.0},
              "availability": {}
            }
          ]
        }
      },
      "state_metrics": [
        {
          "frequency": 0.0,
          "states": {
            "dipole-z": {
              "values": {"bsh_residual": 5e-07, "density_change": 7e-07, "kinetic_change": 2e-07},
              "availability": {"energy_change": {"status": "not_applicable", "reason": "fixed_driving_frequency"}}
            }
          }
        }
      ]
    }
  ]
}
```

The z state's metrics appear once, shared by alpha and beta. Full tensor requests
use this same shape with nine/twenty-seven named components and all contributing
states. An unavailable beta_zzz has no numeric value: its `components` map is
empty and `availability.zzz` contains `status: unavailable` and a reason. Alpha
may remain available if its own inputs are eligible; the response stage remains
incomplete until all requested properties are complete.

If SCF instead exhausts its final update allowance, the requested energy is
reported without a numeric result:

```json
{
  "units": "hartree", "value_type": "real",
  "status": "unavailable", "reason": "required_state_unconverged"
}
```

The last finite energy, for example -7.98 hartree, remains in analysis and the
human report as "Last SCF energy (unconverged)", tied to its measured sample.
Its metrics and usable checkpoint remain available for restart. Response waits
for the reference to qualify; it does not use this observation as a converged
reference. No new property status or alternate numeric tensor format is needed.

Each calculation's per-attempt `result_sources.json` binds these public state
IDs to their original iterations and generations. Here ground means calculation
`ground_state`, task `scf`, attempt `a1`, iteration 7, generation `g8`;
dipole-z means calculation `response`, task `dipole-z`, attempt `a1`, iteration 5,
generation `g6`. Both use protocol `{threshold: 1e-6, k: 8}`. Beta also uses
`response/vbc-zz/g1`. That map
remains tied to this published snapshot when later checkpoints are replaced.

### When electronic attempts differ from the run

This identity-only excerpt omits the required scientific fields to show the
lookup. Two optimizer evaluations share run `a1` but publish through different
electronic attempts:

```json
{
  "version": 2,
  "attempt_id": "a1",
  "calculations": [
    {"calculation_id": "opt_a1_eval_1_ground_state", "attempt_id": "a1_e1"},
    {"calculation_id": "opt_a1_eval_2_ground_state", "attempt_id": "a1_e2"}
  ]
}
```

The first entry resolves
`checkpoints/opt_a1_eval_1_ground_state/attempts/a1_e1/result_sources.json`.
If the optimizer later selects `a1_e3` for that same evaluation, this saved entry
still resolves `a1_e1`, whether `a1_e3` is unfinished or has published a checkpoint.
A new result entry can explicitly name `a1_e3`; it does not change the retained
older entry. Republishing the older entry in a later run also preserves its
explicit `a1_e1` binding. The map's state records may reference numerical samples
from still earlier attempts; map ownership does not relabel those samples.
A missing or mismatched map is an error, not a reason to follow the current
optimizer selector or latest checkpoint. Historical maps preserve this binding
but do not retain archives on their own. While this entry is in current
`results.json`, its exact sources are protected; after replacement, its payloads
may be reclaimed unless active work still needs them. Inspecting reclaimed
functions reports unavailability without rebinding. RESULT-SOURCES tests these cases.

## 4. What the checkpoint records

The following manifest excerpt expands one response metric. The full manifest
also has the stage/settings snapshot, metric definitions and archive inventory,
plus the other reported metrics. They are omitted here, not optional. `g6`
identifies this published generation; a generation number is not an update count.

```json
{
  "version": 2,
  "calculation_id": "response",
  "attempt_id": "a1",
  "generation": "g6",
  "task": {
    "task_id": "dipole-z",
    "kind": "iterative",
    "definition": {
      "subject": "dipole-z",
      "context_id": "response-context-1",
      "operator": "electric_dipole_z",
      "frequency": 0.0,
      "frequency_units": "atomic_angular_frequency",
      "protocol": {"threshold": 1e-06, "k": 8}
    },
    "inputs": [{"role": "ground_reference", "calculation_id": "ground_state", "task_id": "scf", "generation": "g8"}],
    "outcome": {"complete": true, "reason": "final_criteria_met"},
    "iteration": {
      "updates_this_attempt": 5,
      "updates_at_protocol": 5,
      "protocol_complete": true,
      "metrics": [
        {
          "metric": "density_change",
          "definition": "response.density.v1",
          "status": "measured",
          "value": 7e-07,
          "comparison": {"task_id": "dipole-z", "attempt_id": "a1", "protocol": {"threshold": 1e-06, "k": 8}, "iteration": 4, "subject": "dipole-z"},
          "scale": {"status": "measured", "value": 0.12},
          "evidence": {"current_norm": 0.12, "previous_norm": 0.1199998}
        }
      ]
    }
  }
}
```

The value is the measured density difference. Its scale is the larger model
norm of current and previous densities. The two norms alone cannot reconstruct
the difference. This metric entry contains no tolerance or pass/fail test.
The checkpoint retains original settings for provenance; checks against the
current request are written separately to this attempt's analysis.

For example, the density check is `7e-7 <= 1e-6`, so it passes. If a later
request adds relative tolerance 1e-6, the saved scale gives the limit
`1e-6 + 1e-6 * 0.12 = 1.12e-6`. No density recomputation is needed when the
saved definition, sample and scientific context remain compatible.

VBC has the same manifest envelope, with this direct task record:

```json
{
  "task_id": "vbc-zz",
  "kind": "direct",
  "definition": {
    "operators": {"B": "electric_dipole_z", "C": "electric_dipole_z"},
    "frequencies": {"B": 0.0, "C": 0.0},
    "frequency_units": "atomic_angular_frequency",
    "protocol": {"threshold": 1e-06, "k": 8}
  },
  "inputs": [
    {"role": "ground_reference", "calculation_id": "ground_state", "task_id": "scf", "generation": "g8"},
    {"role": "B", "calculation_id": "response", "task_id": "dipole-z", "generation": "g6"},
    {"role": "C", "calculation_id": "response", "task_id": "dipole-z", "generation": "g6"}
  ],
  "outcome": {"complete": true, "reason": "computed"}
}
```

There is no VBC iteration count or residual. Its own generation, archive
inventory and task/definition settings are validated and published together.
Both applications use checked publication and latest-usable retention. An
interrupted write cannot relabel an older archive with these new measurements.

## 5. Internal task timings

This is a timing excerpt: expand one response update, its final evaluation,
and the direct VBC task. Other iterations, SCF details and the coarse protocol
are omitted. Durations are illustrative. Each full record measures its parents
independently; omitted children are not implicit zero-cost work.

```json
{
  "version": 2,
  "attempts": [
    {
      "attempt_id": "a1",
      "execution": {"mpi_ranks": 2, "threads_per_rank": 2},
      "timing": {
        "name": "run",
        "clock_scope": "driver_rank",
        "seconds": 20,
        "children": [
          {"name": "moldft", "calculation_id": "ground_state", "seconds": 12},
          {
            "name": "molresponse",
            "calculation_id": "response",
            "seconds": 7,
            "children": [
              {"name": "load_reference", "seconds": 0.3},
              {
                "name": "protocol",
                "protocol": {"threshold": 1e-06, "k": 8},
                "seconds": 6.3,
                "children": [
                  {"name": "prepare_protocol", "seconds": 0.5},
                  {
                    "name": "solve",
                    "task_id": "dipole-z",
                    "seconds": 5.8,
                    "children": [
                      {
                        "name": "iteration",
                        "iteration": 4,
                        "seconds": 1.0,
                        "children": [
                          {"name": "prepare", "seconds": 0.1},
                          {
                            "name": "evaluate",
                            "seconds": 0.6,
                            "children": [{"name": "bsh_and_intermediates", "seconds": 0.5}, {"name": "metrics", "seconds": 0.1}]
                          },
                          {"name": "decide", "seconds": 0.01},
                          {
                            "name": "advance",
                            "seconds": 0.29,
                            "children": [{"name": "kain", "seconds": 0.24}, {"name": "restrict_step", "seconds": 0.05}]
                          }
                        ]
                      },
                      {
                        "name": "final_evaluation",
                                          "seconds": 0.7,
                        "children": [{"name": "prepare", "seconds": 0.1}, {"name": "evaluate", "seconds": 0.59}, {"name": "decide", "seconds": 0.01}]
                      },
                      {"name": "checkpoint", "iteration": 5, "seconds": 0.3}
                    ]
                  }
                ]
              },
              {
                "name": "vbc",
                "task_id": "vbc-zz",
                "protocol": {"threshold": 1e-06, "k": 8},
                "seconds": 0.2,
                "children": [{"name": "build_vbc", "seconds": 0.15}, {"name": "checkpoint", "seconds": 0.05}]
              },
              {"name": "properties", "task_id": "properties", "protocol": {"threshold": 1e-06, "k": 8}, "seconds": 0.2}
            ]
          }
        ]
      }
    }
  ]
}
```

There is one second of run-level overhead beyond the two application intervals.
The response task's 5.8 seconds include omitted iterations and the endpoint
checkpoint. Its evaluated state is published before VBC and properties consume
that generation. VBC and property nodes state their protocol explicitly because
they sit outside the protocol node.
BSH and metrics are split here only to illustrate separately measurable work;
keep a combined interval when asynchronous work prevents that distinction.
The final evaluation has no advance or update index. VBC has build/save
steps and no iteration. Shared bundle work would be timed once for the bundle.
Parallel state times may overlap; never sum child work to invent a parent's
elapsed time. Reuse in attempt a2 records a2 checking/loading cost without
copying a1's solver time into it.

## 6. What happens on restart

| New request or interruption | Action |
| --- | --- |
| Same request | Validate records/context, reassess saved measurements, reuse SCF and response without `app.evaluate()` or updates. Reuse VBC when exact sources/settings match. Return stored requested properties when their source snapshot and request remain valid. |
| Tighten only response BSH to 6e-7 | Saved 5e-7 still passes; other checks still pass. Reuse without numerical evaluation, update or checkpoint rewrite. |
| Tighten only response density to 5e-7 | Saved 7e-7 fails. SCF remains reusable. Resume the z state at 1e-6/k=8 with a fresh allowance and empty KAIN. A new z generation invalidates dependent VBC/property evidence. |
| A required change has no valid saved measurement/comparison evidence and its history cannot be restored | Obtain fresh evidence and establish a new baseline; never compare the loaded state with itself. Missing previous functions alone do not require reevaluation when valid saved scalar measurements already satisfy the request. |
| Add a frequency | Preserve existing independent states. Start the new state at the first requested protocol. It has its own update allowance. |
| Refining an excited root changes an existing 2PA photon frequency | For the same tracked root and perturbation, retain the saved response functions and protocol if still requested and the other scientific settings match. Rebuild equation work at the new frequency, reset comparisons/KAIN and use a fresh allowance. Old metrics cannot certify the new equation; affected 2PA sources/tensors are recomputed. Unmatched states and removed protocols follow the ordinary first-protocol rules. |
| One state's latest usable checkpoint is interrupted at 1e-4 | Resume that state at 1e-4; later saved independent states keep their progress. A completed intermediate checkpoint advances to its next requested protocol. |
| Crash mid-stage with endpoint-only saves | Resume from the last actually published checkpoint. Analysis history does not turn later unsaved iterations into resumable state. |
| Change operator/smoothing settings | Use supported continuation with fresh scientific checks; affected old metric values cannot certify the new problem. |
| Append a tighter final protocol | Execute that protocol for the affected solves even if the old measured values pass. Reference changes invalidate dependent work. |
| Expected checkpoint is corrupt or missing | Stop affected work and preserve files. No silent fresh start or older-generation fallback. |
| Request inspection/plotting | Read saved data in its original context, including unfinished work. No new convergence demand, solver work or checkpoint rewrite. |

`madqc --wf=response --prefix=lih --continue example.in` requires the existing
run. Automatic restart can also recognize and reuse it without `--continue`.
Explicit fresh-start or inspection mode conflicts with `--continue`.

If the response task were an excited-state bundle, the same response maxiter
and applicable metric tolerances would apply. Measurements would identify each
root, but the controller would return one bundle action and save all roots
together. Existing application checks must establish a valid bundle; every root
then passes its selected applicable metrics, without an extra aggregate test.
Ambiguous matching would make affected density/energy/kinetic changes unavailable;
it blocks completion only when those changes are selected. It would not freeze
roots, invalidate an otherwise usable bundle, or manufacture zero changes.

## 7. What this review can settle

This example connects user controls, app-owned state/history, shared decisions,
direct-task dependencies, property grouping, checkpoint evidence and timing.
Review the behavior and record meanings together. Exact scientific definitions
for the remaining engines and numerical-default validation retain the explicit
implementation prerequisites in the design. No solver was run to create this
example, and the consolidated revision has not received independent design review.


## Direct example: a VBC checkpoint

VBC uses the same native checkpoint publication interface but has no
`IterationMetrics`, iteration counter or solver residual. Molresponse first
qualifies B and C against the requested convergence criteria, then builds VBC.
The store validates their exact published generations and retains them while
VBC or current results still need them. Scientific eligibility stays in
molresponse's dependency manager.

This is the task portion of an illustrative manifest; the enclosing version-2
record also contains `calculation_id: response`, its producing `attempt_id`,
its own `generation`, the stage/settings snapshot, `metric_definitions: {}`
and the actual native archive inventory.

```json
{
  "task_id": "vbc_xy",
  "kind": "direct",
  "definition": {
    "B": {"operator": "electric_dipole_x", "frequency": 0.0},
    "C": {"operator": "electric_dipole_y", "frequency": 0.1},
    "frequency_units": "atomic_angular_frequency",
    "protocol": {"threshold": 1e-06, "k": 8}
  },
  "inputs": [
    {"role": "B", "calculation_id": "response", "task_id": "B_x", "generation": "g8"},
    {"role": "C", "calculation_id": "response", "task_id": "C_y", "generation": "g5"}
  ],
  "outcome": {"complete": true, "reason": "computed"}
}
```

B and C's manifests retain their original attempt/protocol/iteration and metrics;
VBC does not copy those metric histories. Its own attempt records when this
construction was performed. Molresponse may reuse VBC while the exact inputs
and relevant operator/settings context remain valid. A changed dependency or
invalidated construction setting requires a rebuild. A enters the dependency
check when the final quadratic property is assembled from A and VBC.


## Replacing a response after 2PA frequency refinement

Keep the current task-ID format. For a y-dipole response continued from
`0.1500` to `0.1501`, the replacement task's abbreviated record contains:

```json
{
  "task_id": "fd:dipole_y@f0.15010",
  "replaces": [
    {
      "calculation_id": "response",
      "task_id": "fd:dipole_y@f0.15000",
      "generation": "g4"
    }
  ]
}
```

The registration records this relationship before continuation; the new manifest
keeps it as provenance. Required ordinary task fields are omitted here. The
exact source also remains a pending input while its functions are needed for
restart. Once the new task has usable saved work, that temporary need can end.
The relationship itself does not keep old archives alive forever.

After replacement results are published and no active consumer or request needs
the old task, mark only that task `retired: true` in the run inventory. Keep its
small metadata and original result-source associations. Cleanup may then remove
its unreferenced payloads; the response calculation and its other tasks stay
active. Failure before those conditions are met leaves the old task protected.

The current frequency key rounds to five decimal places. A smaller shift may
retain the same task ID: record the actual new frequency in the new generation,
run fresh checks, and keep that task active. Older generations are reclaimed
through the ordinary reference rules. Matching task-ID text never proves that
two numerical frequencies are identical. These records and task-level retirement
are planned integration work, not newly implemented behavior.
