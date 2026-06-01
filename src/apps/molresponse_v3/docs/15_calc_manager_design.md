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

#### Two edge types, two bars — by purpose

The hard/soft split is a difference of *purpose*, which sets a different
acceptance bar and a different mechanism:

| edge | purpose | bar | mechanism |
|---|---|---|---|
| **prerequisite** (hard) | correctness — the dependent physically cannot be built from a partial input (e.g. VBC from a half-solved linear state) | dep must be **converged** | `prerequisites_converged` (or expansion) |
| **seed** (soft) | efficiency — just a good starting guess | best **non-diverged** snapshot, else fresh | `seed_from` / `try_load_*` |

Naming follows purpose: `CalcNode::prerequisites` + `prerequisites_converged()`
on the hard side; `seed_from` on the soft side. A prerequisite is *satisfied*
only when converged; a seed is *available* whenever a non-diverged snapshot
exists (and is optional — its absence just means a fresh guess).

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
  A1.  per perturbation, its SEED STATE solve @ P0          one solve per perturbation
       (ω = 0 if the perturbation has a static request, else its lowest |ω|).
       Anchors of DIFFERENT perturbations are independent.              [parallel across perturbations]
  A2.  all remaining (perturbation, ω) @ P0,
       each seeded from the NEAREST converged state               [parallel across everything]

Phase B — RAMP  (P0 → P1 → … → P_max)
  every (perturbation, ω) climbs its own ladder,
  each step restarting from its own lower-protocol solution    [embarrassingly parallel]
```

The unit of parallelism is the **perturbation perturbation** (dipole_x, dipole_y,
dipole_z, Nuc0_x, …) — these are physically independent solves and run
concurrently throughout. There is no single global "static" solve that
everything waits on.

Rationale:
- Phase A front-loads the cheap work. The only ordering is *within* a perturbation:
  one seed-state solve before the states that seed from it. Across perturbations, A1
  seed states run fully in parallel, and once each perturbation has an seed state, A2 fans
  out across every `(perturbation, ω)`.
- The seed state doubles as the cheapest **memory probe** for everything downstream
  (see below). Because `mem ≈ n_occ × n_leaves × k³` and `n_occ`/molecule are
  the same for every perturbation, a *single* low-protocol seed state already predicts
  the high-protocol footprint for all of them.
- Phase B is embarrassingly parallel: once every `(perturbation, ω)` has a P0
  solution, each protocol step of each ladder is independent (restart from its own
  lower protocol). This is where parallelism maxes out — and where the
  expensive high-protocol work lives.

**Seeding = nearest converged, not a fixed topology.** A fresh `(perturbation, ω)`
seeds from whatever converged state is *closest* in (perturbation, frequency,
protocol) space — typically the same perturbation's nearest-frequency or
coarser-protocol solution. This is exactly what the restart precedence in
`try_load_fd_state` already encodes (exact → coarser protocol → nearby
frequency → fresh); the scheduler's job is only to *order* work so a near
neighbor exists when a node runs. (Earlier drafts named fixed *star* / *chain*
topologies — retired: static-as-hub is just one guess and usually not the
closest. Seeding an FD from the *nearest excited state* is a plausible future
extension but is out of scope for now.)

Policy knob — **seed-selection strategy** (the distance metric / preferences
restart precedence uses), not the topology. Default: nearest converged.

### Future: cross-type seeding (ES ↔ FD)

`CalcNode::seed_from` is a free node-id string, deliberately *unconstrained* to
type — it can reference an ES node from an FD node or vice versa. The data
model already accommodates cross-type seeding; it is **not yet wired** (15a
`build_dag` sets it only for dynamic-FD ← same-perturbation-static, and neither
`schedule` nor the executor consume it). Two physically valuable cases motivate
finishing it, ideally alongside 15b:

- **FD near a resonance ← converged ES root.** As ω → ωₙ the linear response
  is dominated by root n's eigenvector (small denominator ωₙ−ω), so the ES
  vector is the best initial guess — and near-resonant FD is exactly where
  convergence is hardest. This *is* the DerivedFD path (resonant Raman / 2PA at
  ω = ωₙ): the natural guess for `DerivedFD(root n)` is root n's converged
  vector. The storage shapes line up — a TDA root is a `ResponseStateX`, which
  drops straight into FD Static (x=y=root for FD Full) — so it is mechanically
  a loader + projection, not a rewrite.
- **ES targeting a Cartesian direction ← FD response.** Seeding the ES
  iterative subspace with an x-polarized FD response biases toward bright
  states with oscillator strength in x. This additionally needs the ES solver's
  guess interface to accept an external seed (ES-side change, 15b+).

Wiring it requires: (1) the executor honors `seed_from` (load that node's
converged vector ahead of the same-pert FD fallback); (2) a small cross-type
seed loader (ES root → FD-shaped guess); (3) `schedule` treats `seed_from` as a
**soft ordering edge** — run the source first, fall back to fresh if absent, so
no correctness risk; (4) protocol re-projection (the prepare hook already does
this). Tracked as a 15b item because DerivedFD already wants ES-root seeding —
the same piece of work.

Related seed edges this unlocks (all node-id references — same mechanism):

- **TDA → full-RPA warmup (ES → ES).** Mixed ES+FD properties (2PA, resonant
  Raman) need the full-RPA bundle (X *and* Y orbitals), but TDA converges first
  and is a strong RPA guess. TDA and RPA bundles therefore legitimately
  coexist, with `ES(RPA).seed_from = ES(TDA)`. Consequence: the DerivedFD → ES
  hard edge must target the RPA bundle **by identity** (`tda=false, n_roots`),
  not the "first" bundle — 15a's `first_es_id` is a placeholder to replace, and
  there is no "≤1 ES bundle" constraint.
- **Nuclear FD → ES seed.** A nuclear-displacement state can seed an ES solve;
  nuclear nodes are already individually addressable (`nuc_<atom>_<axis>`), so
  targetability needs no extra machinery — only the cross-type loader above.
  (This is why `build_dag` taking `n_atoms` rather than the full `Molecule` is
  sufficient: the per-atom node id already carries the identity.)

### Derived states are FD states (promotion)

A DerivedFD at ω = ωₙ *is* an ordinary FD solve at that frequency. So when
`CalcManager::run` expands a symbolic `"*"` node, each concrete per-root node
should be **promoted to `CalcKind::FD`** (it already gets an `fd_node_id` in the
`fd_states` subtree). Promotion keeps every downstream check uniform —
existence/skip, nearest-frequency restart, and seeding all then treat a
promoted derived state exactly like any other FD point. Intended pipeline:

  1. FD at a range of frequencies (coarse scan);
  2. ES to obtain the excitation energies ωₙ;
  3. FD at the promoted derived frequencies ωₙ — each seeded from the nearest
     already-computed FD frequency *and/or* the converged ES root.

Two gaps to close when this lands (15b):

- 15a's `expand_converged_es` keeps `kind = DerivedFD`, which the FD executor
  does **not** solve (it stubs ES/DerivedFD). Promoting to `CalcKind::FD` on
  expansion is the fix — without it, promoted derived states never run.
- `try_load_fd_state` seeds only from the *same* frequency at a coarser
  protocol; nearest-*frequency* seeding (step 3) is not yet implemented.

### Resonant β at ES-derived frequencies (15b foresight)

Resonant SHG / hyper-Rayleigh is enhanced when **2ω hits an excited state** ωᵢ,
i.e. the driver is ω = ωᵢ/2. So the relevant linear states are `FD(B, ωᵢ/2)`
and `FD(C, ωᵢ/2)`, and VBC (the β source) is built from *those*. This composes
the two dependency mechanisms and needs two additions:

- **A frequency transform on the ES → derived expansion.** Today
  `expand_converged_es` sets the derived FD frequency to the root energy
  (`freq = w`). Resonant β wants `w/2`; other processes want `2ωᵢ` or
  `ωᵢ ± ωⱼ`. So the symbolic derived node must carry *how* to map a root energy
  to its FD frequency, not assume identity.
- **Two-stage expansion.** ES converges → expand `FD(ωᵢ/2)` (promoted, §above)
  → those converge → VBC gated on them via `prerequisites_converged`. The VBC's
  prerequisite FDs have no ids until the first expansion, so the VBC is itself
  symbolic until then. `CalcManager::run` already re-schedules every pass, so
  the second stage drops in: `expand` also emits VBC nodes once their derived
  inputs exist, and the strict-converged gate handles the rest.

---

## Reconcile against what's already on disk

Before scheduling, the manager diffs the plan against
`response_metadata.json`. Per requested state it decides:

| Disk state | Action |
|---|---|
| converged at the target `protocol_key` | **skip** |
| converged at a lower protocol only | **restart** (project up via the restart precedence) |
| present, not converged, **not diverged** | **resume** from the partial archive |
| present but **diverged** | **fresh** — discard the diverged archive and re-seed from the nearest *converged* neighbor (never resume a blown-up state) |
| absent | **fresh** (seed from the nearest converged state, else perturbation/guess) |

The diverged-vs-stalled split matters: a state that merely ran out of iterations
is a good warm start (resume), but a state flagged `diverged` (explosive growth /
NaN) would poison the next solve — so it is thrown away and re-seeded from a
known-good neighbor, not continued.

This makes the manager idempotent and restart-safe: re-running a plan picks up
exactly where a killed run left off, at state granularity.

### Reconcile / seeding refinements (review 2026-06-01)

- **Restart from a non-converged coarse seed (the key one).** A coarse-protocol
  partial — e.g. ~10 iterations at 1e-4 — is a perfectly good seed for the
  finer protocol step; it need NOT be converged. Today `reconcile_protocol` calls
  `has_coarser_converged_fd`, so a coarser *partial* (present, not converged)
  falls through to **Fresh** → the executor skips the load → the partial is
  discarded. Fix: unify the decision and the load on *"best coarser-or-equal
  snapshot that is not `diverged`"* (converged or not). `reconcile_protocol` returns
  `Restart` whenever such a snapshot exists; `try_load_fd_state` selects the
  same one (it currently loads any coarser but must additionally exclude
  `diverged`). Share one "pick best usable source" helper so the verdict and the
  load can never disagree.
- **Resume preserves residual continuity.** `Restart` and `Resume` are
  behaviorally near-identical (both load + iterate); the intended distinction is
  that `Resume` (exact-protocol step partial) carries `last_bsh_residual` forward from
  metadata so the convergence log spans the restart, whereas a coarser `Restart`
  starts the residual history fresh. Keep both actions; lean on `Resume`.
- **`k` is derived only — `override_k` deprecated.** This application always
  derives `k` from `thresh` via `default_k_for_thresh`. Overriding `k` would
  desync the reconcile key from the saved key (everything reads as "absent" →
  silent re-solves), so it is deprecated: drop the `override_k` plumbing from
  `ExecutorContext` / `solve_fd_protocol` / the test `--k` flag.

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

So a perturbation's **seed-state solve @ P0** does triple duty: cheapest solve, seed for
the rest of that perturbation, and — because `mem` is perturbation-independent — the
memory probe that sizes every later group across all perturbations.

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

## Proposed interface (15a)

Header: `calc/calc_manager.hpp` (scheduling — World-free, pure) +
`calc/calc_executor.hpp` (the FD/ES dispatch — World-bound). The split is the
whole point: everything that decides *what runs and in what order* is testable
with a stub executor and a hand-written `response_metadata.json`; everything
that touches `World`, `iterate_protocol`, and the archive lives behind one
injected interface.

### Schedulable unit — `CalcNode`

A node is one response state with its own protocol ladder. Identity is a stable
string (`id`) that doubles as the dependency key — the implicit-by-identity
contract of the "plan is a DAG" section.

```cpp
namespace molresponse_v3 {

enum class CalcKind { FD, ES, DerivedFD, NuclearFD /*, VBC (15b) */ };

struct CalcNode {
  CalcKind            kind;
  std::string         id;          // "fd:dipole_x@f0.05700" / "es:tda_n5" / ...

  // FD / DerivedFD / NuclearFD payload
  Perturbation        pert;
  double              freq = 0.0;  // DerivedFD: filled when its ES converges

  // ES payload
  bool                tda = true;
  int                 n_roots = 0;

  std::vector<double> protocols;   // coarse -> fine ladder (from the plan)

  // symbolic resolution (kept symbolic until the molecule / ES is known)
  std::string         es_root_id;  // DerivedFD: "*" or "es_root_0003"

  // hard edges: must be converged before this node may run. Stored as
  // node ids; resolution is the readiness gate's metadata lookup, NOT a
  // pointer graph — ids survive a restart with no in-memory DAG.
  std::vector<std::string> prerequisites;

  // soft edge: a seed HINT, not a binding. The executor always picks the
  // NEAREST converged state via try_load_fd_state's restart precedence; this
  // field only nudges ordering (run the hinted neighbor first). Empty -> pure
  // nearest-converged / fresh guess. Never required for correctness.
  std::string         seed_from;
};
```

`id` is built from the same string primitives persistence already uses, so a
node maps one-to-one onto a metadata path:
- FD/NuclearFD → `"fd:" + pert.description() + "@" + freq_key(freq)` ↔
  `fd_states/<pert>/<protocol_key>/<freq_key>`.
- ES → `"es:" + (tda?"tda":"full") + "_n" + n_roots` ↔ `excited_states/<protocol_key>`.

### DAG build — `ResponsePlan → std::vector<CalcNode>`

```cpp
/// Pure. Expands nuclear_fd (pert.atom < 0) per molecule into one node per
/// atom×axis; leaves derived_fd "*" symbolic (one node with es_root_id="*",
/// prerequisites = {the ES node id}). Seeds soft edges: every dynamic FD gets
/// seed_from = the static (freq=0) node of the same perturbation if present.
std::vector<CalcNode> build_dag(const ResponsePlan &plan, const Molecule &mol);
```

VBC (15b) adds `CalcKind::VBC` nodes here, with `prerequisites = {FD(B), FD(C)}`;
nothing else in the interface changes.

### Reconcile against disk — pure

```cpp
enum class NodeAction { Skip, Restart, Resume, Fresh };

/// Decide one node's action at one protocol protocol step, by inspecting the aggregate
/// metadata json (already loaded — no filesystem here). target_key is
/// protocol_key(thresh, k). Implements the doc-15 reconcile table,
/// including the diverged split: present+not-converged+not-diverged -> Resume,
/// but present+diverged -> Fresh (discard the blown-up archive, re-seed from a
/// converged neighbor).
NodeAction reconcile_protocol(const CalcNode &node,
                          const nlohmann::json &metadata,
                          const std::string &target_key);

/// Readiness gate: are all of node.prerequisites converged at this protocol step's key?
/// Implicit-by-identity — looks each dep id up in `metadata`. (v3 analog of
/// v2 DerivedStateGate.)
bool prerequisites_converged(const CalcNode &node,
                const nlohmann::json &metadata,
                const std::string &target_key);
```

Both take the metadata as a `const nlohmann::json&` so a test passes a literal
blob; the only filesystem touch is the manager's single `load_or_create`.

### The two-phase schedule — pure

The unit of execution is a **(node, protocol) step**, because the two-phase
schedule interleaves nodes *within* a protocol level (all of P0 before any P1):

```cpp
struct WorkItem {
  const CalcNode *node;
  double          thresh;   // this protocol step's protocol
  NodeAction      action;   // decided by reconcile_protocol at emit time
};

/// Emit work in dependency- and seed-respecting waves. A wave is a set of
/// items with no edges between them — safe to run concurrently (the 15c
/// group layer fans a wave across subworlds; the 15a skeleton runs a wave
/// sequentially). The unit of concurrency is the PERTURBATION CHANNEL:
/// dipole_x / dipole_y / Nuc0_x etc. are independent and share every wave.
/// Ordering realizes:
///   Phase A (protocol step = ramp.front()): each perturbation's seed state (ω=0 if present,
///     else lowest |ω|) — seed states of different perturbations run together — then
///     all remaining (perturbation, ω) at P0 in one wave.
///   Phase B (protocol step in ramp[1..]):   every (perturbation, ω) climbs; a protocol step-P wave
///     per level.
/// Seeding is nearest-converged (executor side via try_load_fd_state); a
/// node's seed_from only biases intra-perturbation ordering. Nodes whose ladder
/// doesn't include a protocol step are skipped for that protocol step. Gated nodes (deps not yet
/// ready — e.g. DerivedFD before its ES bundle converges) are held to a later
/// wave.
std::vector<std::vector<WorkItem>>
schedule(const std::vector<CalcNode> &dag,
         const std::vector<double> &global_ramp,
         const nlohmann::json &metadata);
```

`global_ramp` is the coarse→fine union of every node's `protocols` (reusing
`detail_planner::union_protocols`).

### Injected executor — the only World-bound surface

```cpp
struct NodeResult {
  bool                converged = false;
  std::string         reached_protocol_key;
  std::vector<double> es_root_freqs;  // ES only: drives DerivedFD expansion
};

/// One protocol step of one node. The real impl: set_response_protocol(thresh) ->
/// gs.prepare(...) -> try_load_fd_state/try_load_es_bundle for the seed ->
/// construct FDSolver/ESSolver, set_target -> iterate_protocol({thresh}) with
/// the save+measure_state post_step -> return convergence. Save+metadata are
/// the executor's job (via save_fd_state / save_es_roots), so the manager
/// never writes the aggregate itself.
struct ICalcExecutor {
  virtual NodeResult run_protocol(const WorkItem &item) = 0;
  virtual ~ICalcExecutor() = default;
};
```

Tests substitute a `StubExecutor` that records the call order and writes a
synthetic metadata entry — exercising gate/reconcile/ordering with zero MRA.

### Top level — `CalcManager`

```cpp
class CalcManager {
public:
  struct Policy {
    // Seed selection = nearest converged (restart precedence). The knob is
    // the strategy/metric, not a fixed topology.
    enum class SeedStrategy { NearestConverged } seed = SeedStrategy::NearestConverged;
    int max_iters_per_step = 25;
    // 15c: group sizing knobs land here (node_mem_budget, min/max groups).
  };

  CalcManager(ResponsePlan plan, std::string calc_dir, Policy policy = {});

  /// Build the DAG for this molecule (expands nuclear_fd, wires soft seeds).
  void build(const Molecule &mol);

  /// Drive every node to its target protocol against `exec`. Idempotent and
  /// restart-safe: re-loads metadata between waves, reconciles each protocol step, and
  /// resumes exactly where a killed run stopped. When an ES node's protocol step
  /// returns converged + es_root_freqs, expands its "*" DerivedFD node into
  /// concrete per-root nodes (freq = root energy at the final ES protocol),
  /// appends them to the DAG, and re-schedules.
  void run(madness::World &world, ICalcExecutor &exec);

private:
  ResponsePlan          plan_;
  std::string           calc_dir_;          // holds response_metadata.json
  Policy                policy_;
  std::vector<CalcNode> dag_;
};
```

`run()` is the only method that isn't pure, and even it is thin. Each pass:
reload the aggregate metadata → `schedule(...)` → run **only `waves.front()`**
(15a: sequentially; 15c: fan that one wave across processor groups) → handle ES
expansion → loop. Running just the front wave and re-scheduling each pass is
what makes every wave's reconcile actions final by the time it runs (the
front-wave-only contract above), and is also where the protocol-step barrier
lives. The expensive parts — MRA, MPI, archives — are entirely inside the
executor and the already-landed persistence/`iterate_protocol` layer.

### What this buys, mapped to 15a

- **Unit-testable scheduler:** `build_dag`, `reconcile_protocol`, `prerequisites_converged`,
  `schedule` are pure functions over `ResponsePlan` / `nlohmann::json`. The
  whole "skip/restart/resume/fresh" table and the two-phase ordering get gtest
  coverage with no `World`.
- **One World-bound seam:** `ICalcExecutor::run_protocol`. The real executor is a
  thin adapter over `iterate_protocol` + `save_fd_state`/`save_es_roots` +
  `try_load_*` — code paths that already exist and are tested.
- **Restart-safety falls out:** the manager holds no solved state in memory;
  every wave re-reads the aggregate and reconciles, so kill/resume works at
  protocol step granularity for free.
- **15b/15c slot in without interface churn:** VBC is a new `CalcKind` + a
  `build_dag` edge; group sizing is a `schedule`/`run` policy that turns a
  wave's sequential loop into a subworld fan-out. Neither touches `CalcNode`,
  the gate, or the executor seam.

### Decided

**Rung-level reconcile — the manager owns the protocol-ladder loop**
(agreed 2026-05-29). The manager drives every ready node at P0, then every
node at P1, …; one reconcile + one `exec.run_protocol({node, P})` per `(node,
protocol)`. This buys the all-of-P0-before-any-P1 barrier (so the cheap-level
memory probe sizes the expensive level), the widest same-level cross-perturbation
parallelism, and the finest restart granularity. `iterate_protocol` is invoked
once per protocol step (`thresholds = {P}`) rather than running its own full ramp — the
manager's `for P in ramp` loop replaces the ramp inside `iterate_protocol`.

The node-level alternative (each node carries its whole ladder; one call per
node; `iterate_protocol` climbs `[P0…Pmax]` internally) was rejected: it gives
no barrier between protocol levels, so the first node reaches Pmax before
others finish P0 and the memory probe can't do its job.

### Resolved (folded into the design above)

- **DerivedFD root-frequency provenance (2PA / resonance Raman).** Confirmed:
  these properties cannot be planned until the ES bundle is converged — the
  derived-FD frequencies *are* the converged root energies. The ES bundle is
  internally multi-state and is *not* parallelized across roots (one calc); the
  `"*"` DerivedFD nodes are gated on ES reaching `protocols.back()`, then
  expanded to one node per root at that final root energy.
- **`Molecule` in `build_dag`.** Confirmed needed: vibrational Raman requires
  ∂/∂Q for every atom in x/y/z, so `build_dag` must expand the `nuclear_all`
  sentinel into 3·N concrete nodes — the one place the otherwise
  molecule-independent planner has to see the `Molecule`. Accepted.

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
