# 19 — Excited-State Root-Parallel Design (R5-for-ES)

Status: **design, drafted 2026-06-12.** The first concrete target of the L2 /
R5 state-parallel layer (doc 16). It applies state-parallelism to the **ES
bundle's roots** — the cleanest place to prove the subworld fan-out/gather
pattern before tackling the FD wave scheduler. Read after `16_architecture.md`
(L2/R5) and `00_status.md`.

---

## 1. Motivation — the serial-over-roots bottleneck

The ES solver drives an M-root bundle, but every per-root step runs **serially
on the full communicator**. In `ESSolver::step_rotate_pieces`
(`solvers/es_solver.hpp:328`):

| Phase | Lines | Shape |
|---|---|---|
| per-root build (`density`, `gamma`, `V0x`, `T0x`, `E0x_full`, `E0x`) | ~350-364 | serial `for s`, each call collective on the universe |
| per-root Lambda assembly | ~369-372 | serial `for s` |
| subspace `A=⟨X|Λ⟩`, `S=⟨X|X⟩`, diagonalize, rotate U | ~389-405 | **couples all roots** |
| per-root Theta assembly | ~409-411 | serial `for s` |
| per-root BSH apply | ~425-433 | serial `for s`, each call collective |

For M roots on N ranks this is **M sequential collective passes**. The ranks are
not idle (all N work on root 0, then root 1, …); the waste is **strong-scaling
efficiency**: one response state is "small" relative to N ranks, so a single
root's exchange/Coulomb/BSH on all N nodes is communication-bound (the root
`CLAUDE.md` strong-scaling concern). Running G roots concurrently on G subworlds
of N/G ranks each keeps every subworld in its efficient ranks/work regime **and**
does G roots at once — the same rationale as v2's STATE_PARALLEL, applied to
roots instead of FD states.

The dominant cost per root is the exchange (now with the cached `K0`, doc commit
`e72c15349`) + Coulomb + BSH operator applies. Those are paid M times in series.

---

## 2. The decomposition — what parallelizes, what couples

- **Embarrassingly parallel across roots:** `compute_density/gamma/V0x/T0x/
  E0x(_full)`, Lambda assembly, Theta assembly, BSH apply. Root *s*'s pieces
  depend only on root *s*'s X and the (shared, read-only) ground state.
- **The single coupling point:** the subspace step (`A = ⟨X_i|Λ_j⟩`,
  `S = ⟨X_i|X_j⟩`, an M×M generalized eigensolve, rotate by U). It is **cheap** —
  M² scalar inner products + an M×M `sygvp`, **no operator applies**.

So the per-outer-iteration shape is a classic **fan-out → small gather →
fan-out**:

```
[parallel per-root BUILD of Λ + rotation pieces]   ← fan across subworlds
        ↓ gather X, Λ (and V0x/E0x/γ) to the universe
[universe subspace: A, S, diagonalize → U; rotate]  ← cheap, coupled
        ↓ scatter rotated pieces / assemble Θ
[parallel per-root BSH apply → new X]              ← fan across subworlds
        ↓ gather new X
[universe KAIN + residuals]                         ← cheap, coupled
```

KAIN and the subspace solve stay on the universe (cheap, and KAIN's history must
see the whole bundle). Only the two expensive fan-out phases move to subworlds.

---

## 3. Mechanism — MacroTaskQ (not raw subworlds)

Two precedents in-tree (see research, file:line below):

- **`MacroTaskExchangeRow`** (`src/madness/chem/exchangeoperator.h:357`,
  driver `exchangeoperator.cc:137-170`) — fans one exchange **row** (an
  independent collective `mul_sparse`/Poisson/`dot`) per subworld via a
  `MacroTaskOperationBase` functor with `max_batch_size=1`, shared orbital data
  passed as trailing full-vector args, results auto-reduced to the universe.
  **This is the direct template.** (It is also why we just gitignored
  `MacroTaskExchangeRow_task.*` — MADNESS already fans exchange *pairs* across
  subworlds; we are adding a fan over *roots* on top.)
- **v2 `execute_subgroup_state_solve`** (`src/madness/chem/MolresponseLib.hpp:2414`)
  — raw `MacroTaskQ::create_worlds` + per-subgroup archive reload + per-group
  disk shards + rank-0 JSON merge/broadcast. v2 used the *low-level* primitive
  because each FD state is a **long stateful iterative solve** that does not fit
  one functor call.

**Choice: MacroTaskQ.** The ES per-root work is **single-shot per outer
iteration** (build-one-root → return functions; BSH-one-root → return
functions), not a long stateful solve. That fits the `MacroTask` functor model
exactly, and we get cloud input-marshalling, dynamic load-balanced subworld
scheduling, and universe-side reduction for free. The outer iteration loop stays
on the universe; we fan *within* each iteration. (Contrast v2, which fanned the
*entire* per-state solve and merged via disk — heavier, and the wrong grain for
a coupled subspace iteration.)

### 3.1 The functor (Exchange-row pattern, root grain)

```cpp
// Build one root's Lambda (+ the V0x/E0x/gamma rotation pieces) in a subworld.
class MacroTaskESBuild : public MacroTaskOperationBase {
public:
  // leading vector = the partition axis: ALL roots' response orbitals,
  // flattened to M*n_occ Functions; trailing args are the shared ground state.
  typedef std::tuple<const std::vector<Function<double,3>>& /*X flat, M*n_occ*/,
                     const std::vector<Function<double,3>>& /*amo*/,
                     const Function<double,3>&               /*V_local*/,
                     /* focka, Q, c_xc, lo, n_occ ... */> argtupleT;
  using resultT = std::vector<Function<double,3>>;   // Lambda flat, M*n_occ

  // one ROOT per batch (n_occ functions), not one function — custom partitioner
  MacroTaskESBuild() { partitioner->set_dimension(n_occ); /* see §3.2 */ }

  resultT allocator(World& u, const argtupleT& a) const {
    return zero_functions_compressed<double,3>(u, std::get<0>(a).size());
  }
  resultT operator()(const std::vector<Function<double,3>>& Xroot /*n_occ*/,
                     const std::vector<Function<double,3>>& amo, ...) const {
    World& world = Xroot.front().world();          // == the subworld
    // build this root's density/gamma/V0x/T0x/E0x_full/E0x → assemble Lambda
    // ALL collectives here are on `world` (the subworld), never the universe.
    return lambda_root;                            // n_occ Functions
  }
};
```

A second functor `MacroTaskESBsh` does the BSH phase (input: rotated Θ pieces +
ε; output: new X per root). The driver mirrors `exchangeoperator.cc:146-162`:
build a shared `MacroTaskQ` with `set_nworld(G)`, submit, the universe result is
the gathered Λ (or new X).

### 3.2 Grain & partitioner

The partition axis is one `std::vector<Function>`; a root is `n_occ` functions.
Default `MacroTaskPartitioner` batches in `[min,max]`-sized chunks
(`macrotaskpartitioner.h:251`). We need **exactly `n_occ` per batch** (one root)
— subclass and fix the batch size to `n_occ`, so M batches map to G subworlds
round-robin. Result batches land disjointly (`batch.insert_result_batch`,
`macrotaskq.h:1160`), so root *s*'s Λ goes to its own slot — placed, not summed.

### 3.3 Subspace gather

After the build fan-out, the universe holds all X and all Λ as flat vectors →
reshape to `std::vector<State>` and call the existing
`rs::inner` / `rs::metric` / `rs::diagonalize` / `rs::transform` **unchanged**.
The rotation (`rs::transform`) and Θ assembly can stay universe-side initially
(cheap relative to the operator applies), or be folded into the BSH fan-out
later.

---

## 4. The hard constraint — ground-orbital replication

Each subworld needs the **full occupied-orbital set** (`amo`, `V_local`, `Q`,
`focka`) to build its roots. Passing them as trailing full-vector args means the
cloud **replicates** them to every subworld — the exact `n_occ`-replica memory
driver the root `CLAUDE.md` flags as the C6H6/naphthalene OOM wall. Root-parallel
**trades wall-time for memory replicas**, bounded by the same wall as FD
state-parallel. Therefore this layer is gated on:

1. **15d pre-flight memory estimate** — predict per-subworld footprint
   (`n_occ × n_leaves(k) × k³ × 8 B + intermediates`) from a measured
   low-protocol solve (`state_metrics`), pick `G = clamp(node_mem / per_task, 1,
   n_natural)`, and **abort cleanly** rather than OOM. G=1 ⇒ today's serial path
   (same code, no special case).
2. **(later) node-local shared ground replica** — one copy per node via MPI
   shared-memory windows instead of one per subworld, cutting the multiplier
   from `subworlds/node` to 1 for the ground contribution. Deferred L2 decision
   (doc 16 §"Open design decisions").

Per-subworld memory also *drops* the intermediates term: a subworld holding only
its roots' V0x/T0x/E0x/γ pays the ~8M-`Storage` peak **divided by G**, so
`stream_theta` (doc: `step_recompute_pieces`) becomes less necessary as G grows.

---

## 5. Parity contract

Serial (G=1) and parallel (G>1) results **must match within tolerance**; any
drift is a parallel-layer bug, not the solver (doc 16 L2 contract). The
reduction is a `gaxpy` accumulation into universe slots — order-independent up to
floating point, so per-root placement (disjoint batches) is bit-stable; the
subspace eigensolve runs identically on the universe regardless of G. Gate every
increment with a **serial-vs-parallel parity test** (`cm_*` with `G=1` vs `G=2`)
in addition to `cm_equiv` / the ES regression.

---

## 6. Gotchas (from the MacroTaskQ research)

1. **No universe collectives inside `operator()`.** `run_all` sets
   `forbid_fence` during execution; any `norm2(universe,…)` / `Function::trace()`
   / `inner` on the universe from inside a task deadlocks. Every collective in
   the functor must be on `Xroot.front().world()` (the subworld). (Same class as
   the `[[collective_in_rank_guard]]` rule.)
2. **Subworld fences + cloud cleanup at teardown.** A subworld container must be
   fenced before destruction or the run crashes at teardown (`cloud.h:170-177`)
   — this is the *same failure family* as the parked ES heap-OOB
   (`[[project_v3_es_analysis_parked]]`). Follow the Exchange driver's
   fence discipline exactly.
3. **G is fixed at queue creation** (`MacroTaskQFactory::set_nworld(G)`);
   `color = rank % G`, so co-locating a subworld on a node is the launcher's job
   (`--ntasks-per-node`). The 15d estimate picks G.
4. **Results must be valid task-result types** (`Function`, `vector<Function>`,
   `ScalarResult`); wrap energies as `ScalarResult<double>`, not raw doubles.
5. **Scheduling is dynamic / non-deterministic** (rank-0 first-come scheduler),
   fine for the commutative reduction but don't rely on which subworld runs which
   root.

---

## 7. Prototype increments (phased, each independently validated)

- **Inc 0 — parity harness.** Add a `--es-groups=G` knob (default 1) + a serial
  baseline capture. No behavior change at G=1. Establishes the parity gate.
- **Inc 1 — fan the per-root BUILD only.** Replace the serial build loop
  (`es_solver.hpp:350-372`) with a `MacroTaskESBuild` fan-out; subspace, rotate,
  Θ, BSH stay universe-side. Smallest validatable slice — proves cloud
  marshalling + the root-grain partitioner + parity (G=1 vs G=2 on h2o/c2h4).
- **Inc 2 — fan the BSH apply** (`es_solver.hpp:425-433`) via `MacroTaskESBsh`.
- **Inc 3 — 15d pre-flight + G sizing.** Wire the memory estimate; pick G;
  clean abort. Then sweep G on c6h6 (21 occ) for the strong-scaling curve.
- **Inc 4 (optional) — node-local shared ground replica** (MPI shm windows) to
  cut the replication cost; only if the replica memory is the binding constraint.

FD wave-scheduler parallelism (calc_executor.hpp:863) is the *same* layer applied
to the calc-manager waves and follows after the ES pattern is proven.

---

## 8. Open questions

1. **Rotation placement** — keep `rs::transform` universe-side (simple) or fold
   into the BSH fan-out (one less gather)? Start universe-side.
2. **Build result payload** — return only Λ (recompute V0x/E0x/γ post-rotation,
   the `step_recompute_pieces` trade), or return Λ + the rotation pieces (more
   data across the cloud, the `step_rotate_pieces` trade)? The cloud cost may
   flip the usual memory/compute trade — measure in Inc 1.
3. **Full vs TDA** — TDA roots are `n_occ`; Full carries x+y (`2·n_occ`). The
   grain doubles but the design is identical (the partitioner dimension is the
   per-root block size). Prototype on TDA first. (The symmetric-reduction
   Full-RPA solver was removed 2026-06 — only TDA and direct-Full remain.)
4. **Reuse for FD** — the same two functors generalize to FD states (no subspace
   coupling — even simpler). Keep the functor signatures FD-friendly.
