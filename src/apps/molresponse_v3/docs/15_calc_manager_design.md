# Calc Manager Design

Status: **design, agreed 2026-05-29.** The executor that turns a
`ResponsePlan` (`ResponsePropertyPlanner.hpp`, doc 14 Tier A) into actual
solves, using the persistence + restart layer (doc 13) and the per-state
instrumentation (`state_metrics.hpp`). Sits between the planner and
persistence in the doc-14 data flow.

---

## What the calc manager is

A **memory-aware scheduler over processor chunks.** It takes the response
states a plan requires and decides:

1. **What order** to run them in (sequencing), and
2. **How many processors** each calc gets — from *all processes for one calc*
   (large system, memory-bound) down to *N independent subgroups, one calc
   each* (small system, throughput-bound).

The governing resource is **memory per task for `calc(type, protocol)`**. A
single calculation either takes all processes or is sub-divided into chunks
(subworlds); the right granularity is a function of how much memory one state
needs at the active protocol versus what a node has. This is the v2
`state_parallel` subgroup idea, made *memory-driven* and *protocol-aware*
rather than a fixed group count.

It is NOT a property assembler (that contracts the solved states into
tensors — a separate Tier-A step) and NOT a chemical-property reporter
(doc 14 Tier B). It only schedules and runs state solves, and records results
through the existing `save_fd_state` / `save_es_roots`.

---

## The plan is a dependency DAG

A `ResponsePlan` is not a flat list — it is a DAG of state computations with
two distinct edge types.

### Hard edges (correctness)

A node cannot run until these prerequisites are **converged**. The defining
case is the quadratic / third-order path:

```
FD(B, ω_B) ─┐
            ├──► VBC(B, C, ω_B, ω_C) ──► [β assembly: ⟨x⁽¹⁾_A | VBC⟩ + ζ]
FD(C, ω_C) ─┘
```

VBC (the quadratic source built from converged linear states) is a
**first-class response state** — named, saved, loaded, restartable — not a
hidden intermediate. For β, v3 uses the **2n+1 / VBC-contraction**
formulation: β contracts the linear A states against VBC; there is *no*
explicit second-order solve. The explicit second-order state x⁽²⁾ (solve the
coupled equations with VBC as RHS — structurally an FD solve with `source` =
VBC, reusing `FDSolver`/`iterate`/save-load) is designed-for but **deferred to
γ / higher order**.

The two existing symbolic requests are also hard edges:
- `derived_fd[].es_root_id == "*"` — fires when its ES bundle converges, then
  expands to one FD per converged root (frequency = root energy).
- `nuclear_fd[].pert.atom < 0` — expands to one FD per atom against the
  concrete molecule (an *expansion*, not a convergence dependency).

### Soft edges (seeding / efficiency)

Improve cost but are not required for correctness; the scheduler treats them
as preferences and falls back to a fresh guess if a seed is unavailable:
- `dynamic(ω) ⟵seed⟵ static(0)` — frequency continuation.
- `P_{n+1} ⟵restart⟵ P_n` — coarser-protocol projection.

Both are already implemented as restart precedence in the persistence layer
(`try_load_fd_state` / `try_load_es_bundle`): exact key → lower protocol →
nearby frequency → fresh. The scheduler's job is to *order* work so those
seeds exist when a node runs.

### Dependency resolution is implicit-by-identity

Nodes do not store edge lists. A VBC node reconstructs its prerequisites'
identity keys — `FD(B, ω_B, P)`, `FD(C, ω_C, P)` at the same `protocol_key` —
and the **readiness gate** looks them up in `response_metadata.json`. The
aggregate metadata file *is* the dependency-resolution substrate. (This is the
v3 analog of v2's `DerivedStateGate`, minus the parallel edge bookkeeping.)

---

## The two-phase schedule

Sequencing is not obviously optimal, but this ordering both maximizes the
parallel region and exploits every restart path we built:

```
Phase A — SEED  (at the lowest protocol P0)
  A1.  static (ω = 0)                         one calc, robust, cheapest   [serial]
       (skipped if seeding from an ES guess / external initial state)
  A2.  all other frequencies @ P0,
       each seeded from static@P0                                         [parallel]

Phase B — RAMP  (P0 → P1 → … → P_max)
  every (state, freq) climbs its own ladder,
  each step restarting from its own lower-protocol solution    [embarrassingly parallel]
```

Rationale:
- Phase A front-loads the cheap, serializing work. `static@P0` is one state
  (hard to parallelize *across*), but it is the cheapest solve and the best
  seed for every dynamic frequency. It is also the cheapest **memory probe**
  for everything downstream (see below).
- The only strict serialization is `static@P0` first; once it exists, all
  other frequencies at P0 run in parallel, each seeded from it.
- Phase B is embarrassingly parallel: once every `(state, freq)` has a P0
  solution, each rung of each ladder is independent (restart from its own
  lower protocol). This is where parallelism maxes out — and where the
  expensive high-protocol work lives.

Policy knob — **seeding topology** at low protocol:
- *star-seed* (default): every dynamic ω seeded directly from `static@P0`
  (max parallel),
- *chain*: `static → ω₁ → ω₂ …` (closer guesses, more serial).
Default to star-seed where solves are cheap (low protocol).

---

## Reconcile against what's already on disk

Before scheduling, the manager diffs the plan against
`response_metadata.json`. Per requested state it decides:

| Disk state | Action |
|---|---|
| converged at the target `protocol_key` | **skip** |
| converged at a lower protocol only | **restart** (project up via the restart precedence) |
| present but not converged | **resume** from the partial archive |
| absent | **fresh** (seed per the soft edges, else perturbation/guess) |

This makes the manager idempotent and restart-safe: re-running a plan picks up
exactly where a killed run left off, at state granularity.

---

## Processor-group allocation (memory-driven, protocol-dependent)

Memory per state grows steeply with protocol —
`mem ≈ n_occ × n_leaves(k) × k³ × 8`, and CLAUDE.md measures `n_leaves` at
k=10 ≈ 4.6× the k=6 count. So the machine that fits 16 states in parallel at
k=6 may fit only 2 at k=10. The scheduler therefore **re-decides subdivision
at each protocol level**: wide fan-out at low protocol, collapsing toward
all-processes-per-calc at the expensive high protocol. (CLAUDE.md's "delay
full subworld expansion," made concrete.)

Group size at protocol P:

```
groups(P) = clamp( floor( node_mem_budget / mem_per_task(state, P) ), 1, n_natural )
```

where `mem_per_task(state, P)` comes from the instrumentation feedback loop
below, and `n_natural` is the number of independent ready nodes.

**This is the parallelism layer and lands AFTER the single-group skeleton.**
The first scheduler runs all-processes-per-calc (groups = 1): correct
sequencing, gates, reconcile, and logging, with no subdivision. Get the DAG /
schedule / restart loop right serially; add subworld sizing as a separate
increment once the loop is trusted and the memory data exists.

---

## Instrumentation feedback loop

The scheduler is both producer and consumer of the per-state metrics
(`state_metrics.hpp`, already landed):

- **Producer:** every solve records `{coeffs, bytes, rss_gb, iters}` (and,
  once the scheduler owns the timing loop, `wall_s`) at each protocol
  boundary into `response_metadata.json`.
- **Consumer:** `mem_per_task(state, P_high)` is needed *before* allocating
  the expensive protocol. The **memory-estimation study** predicts it from a
  measured low-protocol solve: because the low-protocol solve gives the
  *actual* adaptive leaf structure for this molecule (not a worst-case bound),
  scaling by `k³ × leaf-growth(P_low → P_high)` yields a real per-task
  estimate → a clean **pre-flight abort** instead of exit-137.

So `static@P0` does triple duty: cheapest solve, best dynamic seed, and the
memory probe that sizes every later group.

---

## Where it sits / what it reuses

```
ResponsePlan ─► CalcManager ─► (per node) ─► FDSolver / ESSolver via iterate_protocol
                  │                              │
                  │ reconcile vs metadata        │ restart: try_load_fd_state /
                  │ gate (hard edges ready?)     │          try_load_es_bundle
                  │ order (two-phase)            │ save:    save_fd_state /
                  │ size groups (later)          │          save_es_roots  (+ metrics)
                  ▼                              ▼
            response_metadata.json  ◄────────────┘   (status, restart, metrics, deps)
```

Reused as-is: `iterate_protocol` (+ its `post_step` hook for per-protocol
save/metrics), the restart precedence, the aggregate metadata, `protocol_key`,
`measure_state`. New: the DAG/gate, the two-phase ordering, the reconcile
step, and (later) subworld sizing.

---

## Increment plan

- **15a — single-group skeleton.** DAG build from a `ResponsePlan`, the
  readiness gate (implicit-by-identity against the aggregate), reconcile
  (skip/restart/resume/fresh), the two-phase ordering, all-processes-per-calc.
  Pure scheduling logic is unit-testable with a stubbed metadata/solve;
  end-to-end on the allocation drives real FD/ES solves. β still needs VBC
  (15b) to be fully runnable; α and ∂α/∂Q-resonant are runnable from existing
  states.
- **15b — VBC planner + persistence.** `VBCRequest` in the plan; the
  `derived_states/` subtree + `vbc__<B>__<C>__fB_fC__<protocol_key>` naming +
  save/load with restart precedence. Closes the β path (2n+1 contraction).
- **15c — processor-group allocation.** Memory-driven, protocol-dependent
  subworld sizing (v2 MacroTaskQ subgroups). Consumes the metrics.
- **15d — memory-estimation study + pre-flight abort.** Offline model from the
  logged data, fed into 15c's sizing and a clean OOM-avoiding abort.

x⁽²⁾ / second-order solve (γ) is beyond this sequence — it reuses `FDSolver`
with a VBC source and the 15b persistence, when γ lands.

---

## Why this shape

1. **Restart-safety is free.** Reconcile-by-identity + the persistence layer
   make every run idempotent at state granularity — kill and resume anywhere.
2. **Parallelism is maximized where it's expensive.** The two-phase schedule
   serializes only the cheap seed and parallelizes the costly ramp.
3. **Memory is a first-class input, not an afterthought.** Group sizing is
   driven by measured footprints; the cheapest solve doubles as the probe that
   prevents OOM on the expensive one.
4. **The hard/soft edge split keeps correctness and efficiency separable.**
   Gates enforce physics; seeds optimize cost and can always fall back.
