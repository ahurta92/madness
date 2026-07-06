# 21 — ES State-Parallel via Node-Aligned Subworlds (Stage 2, chosen path)

Status: **DRAFT** (design, no code). Chosen 2026-06-14 over the single-World
pmap design (doc 20) after the output-pmap finding (doc 20 §5.4). Companion to
doc 19 (the MacroTaskQ subworld track) and doc 20 (the pmap alternative).

The goal you asked for — **subworlds for free intermediate-locality + a
node-local ground state instead of per-task replication** — is achievable, but
NOT the way it was first framed. §1 is the hard infrastructure finding; §2 is
the design that gets the goal anyway; §3 costs; §4 reuse; §5 plan.

---

## 1. Infrastructure reality (verified in MADNESS source)

Two findings that shape everything:

**(a) There is no shared-memory Function storage.** MADNESS has no
`MPI_Win_allocate_shared` for coefficient data. `Split_type(MPI_COMM_TYPE_SHARED)`
exists (`safempi.h:674`) but only builds intra-node *communicators*, not shared
*memory*. `NodeReplicated` (`worlddc.h:81`, `replicate_on_hosts` at
`worlddc.h:656`) puts **one copy per node owned by the lowest rank on the host**,
and other ranks reach it by **intra-node active messages**, not a shared
pointer. So multiple subworlds on a node **cannot point at one physical φ**.
Building that = new MADNESS core (shm-windowed `FunctionImpl` storage) — out of
scope.

**(b) A subworld cannot reference universe data; it must be copied in.**
`copy_coeffs_different_world` (`funcimpl.h:1143`) fully serializes coefficients
across the world boundary. The Cloud (`cloud.h`) is exactly this:
`store`/`load` of Functions universe↔subworld via
`ContainerRecord{Output,Input}Archive` (`parallel_dc_archive.h`). Even the
`StoreFunctionPointer` policy ends in `arg = copy(subworld, arg)`
(`macrotaskq.h:1128`) — the copy is deferred, not avoided.

**Conclusion:** "share one φ across subworlds via shm" is impossible without new
core machinery. So we get the *goal* a different way (§2).

---

## 2. The design that reaches the goal: node-aligned subworlds

### 2.1 Why the v2 OOM happened, and how we avoid it
The v2 OOM (`CLAUDE.md`: "every MPI task holds a complete replica of all
occupied ground-state orbitals") came from φ being **RankReplicated** — a full
copy on *every rank*. That is a *choice*, not a requirement. If instead φ is
left **distributed within a subworld** (the subworld's normal pmap), then a
subworld holds **one copy of φ spread over its ranks**, and the per-rank φ cost
is `|φ| / (ranks in subworld)` — no replication penalty.

### 2.2 One subworld per node
Create the subworlds **node-aligned**: one subworld = one physical node, via
`universe.mpi.comm().Split_type(SHARED)` (`safempi.h:674`) instead of
MacroTaskQ's default `color = rank % nsubworld` (`macrotaskq.h:647`, which
*interleaves* subworlds across nodes). Then:

- Each node-subworld Cloud-copies φ in **once per protocol** → **one φ copy per
  node**, distributed over that node's ranks. Per-rank φ = `|φ| / ranks_per_node`
  — **identical to the current single-World cost. No memory penalty, no OOM.**
- This is the "node-local ground state" you wanted — realized by *placement*
  (one distributed copy per node-subworld), not by shm sharing.

### 2.3 Intermediate locality is free
Inside a node-subworld, `FunctionDefaults::get_pmap()` *is* that subworld's
pmap. So **every** intermediate (V0x, T0x, E0x, γ, θ, …) is automatically
subworld-local and mutually aligned — the assemble_lambda gaxpy just works. This
is the whole reason for going subworld instead of single-World pmaps (doc 20
§5.4: in one World those intermediates land on mismatched pmaps and crash).

### 2.4 Per-iteration flow (build local, couple global)
```
ONCE per protocol:
  • make node-aligned subworlds (Split_type SHARED)
  • Cloud φ, V_local, ε, focka into each node-subworld   (1 distributed copy/node)
  • assign roots round-robin to node-subworlds

EACH ES iteration:
  • FAN OUT: ship each root's current X to its node-subworld (Cloud)
  • LOCAL (per node-subworld, all roots it owns, batched as in Stage 1):
       build V0x/T0x/E0x/γ → assemble Λ → (Λ, and X) are the only things
       that must leave the subworld
  • GATHER: Cloud Λ + X for ALL roots back to the universe
  • COUPLE (universe): rs::inner(A), rs::metric(S), rs::diagonalize → ω,U;
       rs::transform(X, U)   ← rotation mixes roots, must be universe-side
  • SCATTER: ship rotated X back to node-subworlds
  • LOCAL: θ-assembly + BSH apply → new X per root        (stays in subworld)
```
The universe step is the same `rs::` code we already have — it runs on the
gathered, co-located bundle (0.5% of the work). Everything else is subworld-local.

---

## 3. Honest costs / risks

### 3.1 Per-iteration Cloud coupling is HEAVIER than the pmap shuffle
The gather/scatter is a **cross-world serialize/deserialize** (`copy_coeffs_
different_world`), not the cheaper intra-World `copy(f,pmap)` redistribution the
pmap design would use. We pay it every iteration (~40+). It only moves the
bundle (M·n_occ small functions, + Λ), which is tiny vs the build FLOPs — but it
must be **measured**; it is the main efficiency unknown (mirror of doc-20's
shuffle question, but more expensive per unit).

### 3.2 The Exchange operator inside a subworld
gamma's exchange `K[φ,x]` (36%, doc-19 §5.3) uses `madness::Exchange`, whose
own `MacroTaskExchangeRow` machinery is itself a subworld scheme — **nesting it
inside a subworld task is forbidden**. So inside a node-subworld the exchange
must use the **direct** `apply_exchange` path (no inner MacroTaskQ). That's the
K0=nullptr fallback already present in the kernels. Need to confirm the direct
path performs acceptably over a node-subworld's ranks.

### 3.3 Root-parallelism vs node count
One-subworld-per-node gives `nnodes` groups; each solves `M/nnodes` roots
(batched within the node). Root-parallelism is therefore capped at `nnodes`, not
`M`. For few nodes + many roots that under-parallelizes the *across-root* axis
(though each root still uses all the node's ranks data-parallel). More subworlds
→ more root-parallelism but more φ copies (memory). The knob is "subworlds per
node"; default 1 (memory-first). This is a tunable, not a fixed choice.

### 3.4 It's still a scaling bet
At small scale (h2o, 4 ranks, 1 node) there is exactly one subworld = the
universe, so this is a no-op + pure coupling overhead. The payoff is at
multi-node scale where confining each root's build+exchange to one node's
fabric beats spreading it over the whole machine. Must be shown with the
`--es-time` harness on a multi-node run before it's turned on.

---

## 4. What's reusable now (this flips doc 20 §6)

Going subworld means the **v2 / MacroTaskQ machinery becomes reusable**, which
it was not for the pmap design:

| Piece | Reusable here? |
|---|---|
| `MacroTaskQ` subworld pool + scheduler (`macrotaskq.h`) | **Yes** — but extend `create_worlds` to node-aligned `Split_type(SHARED)`. |
| `Cloud` store/load (`cloud.h`) | **Yes** — the fan-out/gather/scatter transport. |
| `ContainerRecord{Output,Input}Archive` | **Yes** — Function (de)serialization. |
| v2 `execute_subgroup_state_solve` / group ownership planning | **Concepts yes** — round-robin root→subworld assignment; v3 ES is one bundle so simpler. |
| v2 per-task φ **replication** | **No — deliberately drop it** (it's the OOM cause; we distribute φ per node-subworld instead). |
| v3 `ESSolver` / `kernels/` / `rs::` | **Yes, fully** — the subworld-local build is the existing kernels run in the subworld; the universe couple is the existing `rs::`. Stage-1 batching (`tda_batch`) runs *inside* each subworld. |

---

## 5. Phased plan

- **Inc 1 — node-aligned `create_worlds`.** A `Split_type(SHARED)` variant +
  a tiny test that prints the node→subworld map and confirms 1 subworld/node.
  No solver change. (If MPI ranks aren't node-contiguous, this also fixes the
  default interleaving.)
- **Inc 2 — φ into subworlds once, measure.** Cloud φ/V_local into each
  node-subworld; confirm **one distributed copy per node** (worldmem/MaxRSS) and
  that per-rank φ matches single-World. Pure memory validation.
- **Inc 3 — fan one phase (build) to subworlds, gather Λ, couple in universe.**
  Wire the §2.4 flow for the *build* only; BSH still universe-side. Validate
  bit-identity (converged roots, verify_es_batch.sh) at 1 node (subworld=universe)
  and 2 nodes; measure the Cloud coupling cost (the §3.1 unknown).
- **Inc 4 — fan BSH + exchange too; tune subworlds-per-node; scale study.**
  c2h4 / c6h6 across nodes; the decision point where this justifies itself.

**First step: Inc 1** — node-aligned subworld creation + a map-printing test.
Zero solver risk, and it's the prerequisite that makes the "one φ per node"
property real. Everything else builds on it.
```
