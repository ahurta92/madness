# 20 — ES State-Parallel via Locality Process Maps (Stage 2)

Status: **DRAFT** (design, not yet implemented). Companion to doc 19
(staged plan). Stage 1 (bundle batching, `--es-batch`) landed; this doc is
Stage 2 — the "two-pmap" locality design.

This doc is written to be readable without prior MADNESS-distribution
knowledge. §1–§3 build the mental model; §4 is the design; §5 the honest
trade-offs; §6 what is/isn't reusable from the v2 subworld state-parallel;
§7 the phased plan.

---

## 1. What a `Function` actually is, physically

A `madness::Function<double,3>` is **not** a grid or a contiguous array. It is
a **distributed octree**: a tree of nodes, one per `Key<3>` (a box at some
refinement `level` + translation), each node holding a small coefficient
tensor (`k³` doubles, k≈10). The set of nodes is stored in a
`WorldContainer<Key<3>, FunctionNode>` — a distributed hash map.

**Where does each node live?** A node for key `K` lives on exactly one MPI
rank, decided by a **process map (pmap)**: `pmap->owner(K) -> rank`
(`worlddc.h:135`). The default for MRA functions is `LevelPmap`
(`funcimpl.h:105`):

```cpp
ProcessID owner(const Key& key) const {
    if (key.level() == 0) return 0;
    hash = (level<=3 || odd level) ? key.hash() : key.parent().hash();
    return hash % nproc;          // <- scatter across ALL ranks
}
```

So **a single Function is spread across all ranks** — its tree nodes are
hashed out over the whole communicator. This is "data parallelism within one
function": even one root's orbital is distributed over every rank.

Key consequence to internalize: **the current v3 ES solver already uses every
rank for every operation.** Nothing about it is single-threaded-per-root in
the data sense. The "serial over roots" problem (doc 19 §1) was never about
data placement — it was about **fences** (next section).

---

## 2. Collectives, fences, and "owner computes"

Operations on Functions (`mul_sparse`, `apply`, `inner`, `truncate`, `gaxpy`)
are **collective over a `World`**: every rank in that World must call them, and
they end at a **fence** (`world.gop.fence()`) — a barrier that also flushes the
task queue. Between fences, work is asynchronous tasks.

**Owner computes:** a task touching node `K` is sent to `pmap->owner(K)` and
runs there (`worlddc.h:817`, `task(owner(key), ...)`). So *where a Function's
nodes live determines where its compute happens.* This is the lever Stage 2
pulls: move the data (via pmap) and the compute follows.

Why Stage 1 helped: the old loop did `M` kernel passes, each ending in a
fence — `M` barriers per phase, and each pass only offered the task queue
`n_occ` functions of work before stalling at its fence. Batching issues
`M·n_occ` functions before ONE fence → the queue is full, ranks stay busy, and
`M−1` barriers per phase vanish. Measured: −7% wall, all from the `build`
phase. **It did not change data placement at all.**

---

## 3. The one rule that shapes everything: binary ops assume aligned trees

`inner(f,g)`, `matrix_inner`, `gaxpy(f,g)`, `mul(f,g)` are computed
**node-local then reduced**. Look at `FunctionImpl::inner_local`
(`funcimpl.h:5670`), comment verbatim:

> `/// Returns the inner product ASSUMING same distribution`

Mechanically (`funcimpl.h:5637`): rank `r` walks `f`'s **local** nodes and, for
each key, *probes `g` locally* (`other->coeffs.probe(key)`). It sums
`f_node · g_node` only when **both** nodes are local to `r`. Then `gop.sum`
adds across ranks. `matrix_inner` (`vmra.h:985` → `inner_local` at
`funcimpl.h:6088`) builds per-rank local key→node maps for left and right and
multiplies where they coincide on a rank.

**Therefore:** if `f` and `g` are on *different* pmaps (a key on rank 2 for `f`
but rank 5 for `g`), the local probe misses, the term is dropped, and the
result is **silently wrong** — no error, just an undercount. So:

> **Binary ops require both operands on the SAME pmap (same key → same rank).**

This single rule is the entire reason Stage 2 has a "re-map for the subspace
step" cost (§4.4) and the reason the ground-state data must be *replicated*
(§4.2) — a localized root and a globally-distributed φ are NOT aligned, so
their product would be wrong unless φ is available everywhere.

Pmaps available out of the box (`worlddc.h:249–315`):
- `LevelPmap` / `WorldDCDefaultPmap` — hash scatter over all ranks (default).
- `WorldDCLocalPmap` — every rank owns every key (**RankReplicated**: a full
  copy on every rank).
- `WorldDCNodeReplicatedPmap` — owner = lowest rank on the host
  (**NodeReplicated**: one copy per physical node).

And the tools to move between them:
- `copy(f, pmap, fence)` (`mra.h:2111`) — redistribute `f`'s tree to a new pmap
  (this is communication: nodes physically move).
- `Function::replicate(DistributionType, fence)` (`mra.h:707`) — make a
  Rank/Node-replicated copy.
- `f.get_pmap()` (`mra.h:701`) — read a function's current pmap.

---

## 4. The design — locality pmaps ("two process maps")

### 4.0 The idea in response language

Give every excited-state root a **home**: root `i`'s response orbitals live
on a rank-group `Gᵢ` (a slice of the communicator), not scattered over all
ranks. The **ground state is shared everywhere** (replicated per node). Then:

- The per-root work (build, gamma, BSH) for root `i` runs **on its home group
  `Gᵢ`** — because owner-computes follows the data. Issue all roots' work
  together (Stage-1 batching) and the groups run **concurrently**, each
  confined to its own ranks → far less cross-rank chatter per operation.
- The only step that needs all roots together is the tiny `M×M` subspace
  solve (measured 0.5% of the iteration). For that one step we **briefly
  gather** roots+Λ onto a common pmap, diagonalize, rotate, and **scatter
  back** to the homes.

So: one *global* layout (all-ranks, for the collective subspace step) and one
*state-local* layout (per-group, for the embarrassingly-parallel build) — the
two pmaps you intuited.

### 4.1 The custom group pmap (lower level)

A new pmap class — confine a key to a contiguous rank range:

```cpp
class GroupPmap : public WorldDCPmapInterface<Key<3>> {
  int base_, size_;  Hash<Key<3>> h_;
public:
  GroupPmap(int base, int size): base_(base), size_(size) {}
  ProcessID owner(const Key<3>& k) const override {
    return base_ + (k.level()==0 ? 0 : h_(k) % size_);
  }
};
```

`G` groups over `P` ranks → group `g` owns ranks `[g·P/G, (g+1)·P/G)`. Root `i`
is assigned to group `i % G`. Place a root there with
`copy(x_i, group_pmap[i%G])` (once at seeding, and after each rotation).

### 4.2 Replicate the ground-state data (the alignment fix)

Every per-root kernel multiplies the root against ground data: `mul_sparse(V_local,
x_i)`, the Coulomb `J·φ`, the exchange `K[φ,x_i](φ)`. By §3, those operands
must be co-located with `x_i`. Since `x_i` is now confined to group `i`, the
ground data must be **available on every rank** — i.e. replicated:

```cpp
gs.V_local_alpha.replicate(NodeReplicated);
for (auto& phi : gs.amo) phi.replicate(NodeReplicated);
```

`NodeReplicated` = one physical copy per node (shared by all its ranks / all
groups on that node), not one per rank and not one per group. This is exactly
the "shared node-local orbital replica" improvement called out in the root
`CLAUDE.md` scaling section.

The kernels themselves **do not change** — they still call collective ops on
`world_`. The behavior changes only because the *operands' pmaps* changed:
replicated φ is local everywhere, localized `x_i` is local to `Gᵢ`, so each
root's products are computed on `Gᵢ`. (Two big caveats — the `Exchange`
operator, §5.3, and **intermediate output pmaps, §5.4, which turns out to be
the dominant difficulty**.)

### 4.3 Build/gamma/BSH: no code change, just placement

With roots localized + ground data replicated, the existing batched path
(`tda_batch::`) already issues all roots' ops in one pass. Owner-computes
routes root `i`'s tasks to `Gᵢ`. Groups run concurrently; the single fence
syncs all. This is Stage 1's batching *plus* placement — same kernels.

### 4.4 The subspace step: the irreducible coupling motion

`A = rs::inner(roots, Λ)` and `S = rs::metric(roots, roots)` need
`⟨root_i | Λ_j⟩` for all `i,j`. By §3 those operands must share a pmap — but
root `i` is on `Gᵢ`, Λ_j on `Gⱼ`. So, once per iteration:

1. `copy` roots and Λ onto a **common pmap** (the default all-ranks LevelPmap).
2. `rs::inner` / `rs::metric` / `rs::diagonalize` → ω, U (cheap: 0.5%).
3. `rs::transform(roots, U)` — the rotation **mixes roots**, so it must happen
   on the common pmap (it reads all roots to form each output).
4. `copy` the rotated roots back to their group pmaps for the next build.

Steps 1+4 are the real Stage-2 cost: a per-iteration **redistribution** of the
bundle (`M·n_occ` functions) between group-local and all-ranks layouts. The
*compute* in step 2–3 is negligible; the *data motion* is what we must measure.
(This is why doc 19 noted `step_recompute_pieces` becomes natural here: if we
recompute V0x/E0x/γ after rotation group-locally, only roots+Λ ever cross the
group boundary.)

---

## 5. Honest trade-offs (do NOT skip this)

### 5.1 It is a SCALING bet, not a free win
At the tested scale (h2o, 3 roots, 4 ranks) the default scatter pmap *already*
uses all ranks per op, and Stage 1 captured the fence win. Localization pays
off when **rank count is high enough that per-op all-to-all communication
dominates**, and when there are **enough roots to fill the groups**. At small
scale it can be *slower* (the per-iter re-map in §4.4 is pure overhead with no
locality benefit yet). The `--es-time` harness must show the per-op comm is
actually the bottleneck before this is worth turning on.

### 5.2 Memory: replication vs the current single distributed copy
Today (single World, default pmap) the ground state is **one** copy distributed
over all ranks — already memory-efficient. NodeReplication makes it **one copy
per node** — *more* memory than now, but far less than v2 subworlds (one per
task). Net: Stage 2 trades some ground-state memory for locality. Whether
that's a win depends on `n_occupied` vs node memory (the C6H6/naphthalene OOM
budget in `CLAUDE.md`).

### 5.3 The Exchange operator is the hard part (and it's 36% of the time)
`gamma`'s exchange `K[φ,x_i](φ)` goes through `madness::Exchange`, which has its
**own internal parallelization** (a MacroTaskQ row scheme, `exchangeoperator`).
It will not automatically respect our group pmaps, and nesting its task queue
under group placement is delicate. Localizing the exchange is both the biggest
prize (it's the dominant per-root cost) and the riskiest piece — likely its own
increment, possibly needing an exchange variant that honors a target pmap.

### 5.4 Output pmap discipline — the dominant difficulty (verified)
This is worse than "fiddly"; it's the crux. An op does **not** put its output on
"the per-root operand's" pmap — it inherits from a *specific* operand:
- `mul_sparse(world, a, v)` builds each result via `set_impl(left=a)`
  (`mra.h:1661` → `new implT(*a.get_impl(), a.get_pmap())`, `mra.h:set_impl`).
  So `V0x = mul_sparse(V_local, x_i)` lands on **V_local's** pmap, NOT `x_i`'s
  group. T0x (kinetic on `x_i`) lands on `x_i`'s group. So `V0x` and `T0x` end
  up on **different pmaps**.
- `assemble_lambda` then gaxpy-adds T0x + V0x − E0x_full + γ. Compressed gaxpy
  goes through `merge_trees`, which **asserts `get_pmap()==other.get_pmap()`**
  (`funcimpl.h:1233`). Mismatched pmaps → assert/crash, not a silent wrong
  number.

**Consequence:** you cannot just localize the roots — every per-root
intermediate (V0x, T0x, E0x, γ, θ) must be forced onto the *same* group pmap,
or assembly fails. MADNESS exposes no uniform "output on pmap P" knob, so this
means an explicit `copy(out, group_pmap)` after (or a pmap-aware variant of)
**every** kernel op. That is pervasive and fragile.

**This materially favors subworlds for *full* per-root locality** (§6 update):
inside a subworld, `FunctionDefaults::get_pmap()` *is* the subworld pmap, so
every intermediate is automatically subworld-local and mutually aligned —
locality is free, no per-op pmap fighting. The single-World pmap design gets
in-place global collectives + NodeReplicated φ (less memory), but pays by
fighting output-pmap conventions at every kernel step. The subworld design gets
locality for free, but pays with per-subworld φ replication (the `CLAUDE.md`
OOM driver) + Cloud coupling. **This trade-off is the real Stage-2 decision; it
must be settled before implementing.**

---

## 6. What can / cannot be reused from v2 state-parallel

v2 (`molresponse_v2`) has a mature state-parallel system, but it uses a
**different mechanism** — worth being precise about overlap.

**v2 mechanism = subworlds.** It splits the communicator into `G` independent
**sub-Worlds** (separate MPI sub-communicators), each a self-contained universe
that solves a subset of states, coordinated through the **Cloud** + MacroTaskQ,
with the ground state **replicated per subworld** and metadata **sharded** then
merged. (Files: `StateParallelPlanner.hpp`, `execute_subgroup_state_solve`,
MacroTaskQ.)

**This design = one World + pmaps.** No sub-communicators, no Cloud marshalling,
one ground replica per *node* (not per subworld).

| Piece | Reusable for the pmap design? |
|---|---|
| Subworld / sub-communicator machinery | **No** — different model. (It IS the doc-19 comparison track, Inc 0–4.) |
| Cloud marshalling of states between worlds | **No** — we never leave the one World; we re-map within it. |
| Metadata sharding + rank-0 merge | **No** — the ES bundle solve is one in-memory process. |
| `MacroTaskPartitioner` (root-grain batching) | Only if we take the subworld comparison route. |
| **Group-ownership *planning* concept** (how many groups, which states to which group, auto thresholds) | **Yes, conceptually** — but v3 ES is one bundle, far simpler than v2's frequency/channel planning; we need only `G` + round-robin. |
| **The memory lesson** (per-task φ replication is the killer) | **Yes** — it's exactly why we choose NodeReplicated, not per-group/per-task. |
| v3 ES solver itself (`ESSolver`, `kernels/`, `rs::`) | **Yes, fully** — Stage 2 changes *where functions live*, not the math. Kernels are pmap-agnostic (they call collectives on `world_`); only placement + ground replication + the §4.4 re-map are new. |

Bottom line: almost none of v2's *code* transfers (it's the subworld
mechanism, which is the alternative). What transfers is the *concept* of group
ownership and the *lesson* about ground-state memory — plus all of the v3 ES
kernel/solver stack, which is reusable as-is.

---

## 7. Phased implementation plan (each validatable with our harness)

Every increment is gated and checked with `verify_es_batch.sh` (bit-identity)
+ `time_es_phases.sh` (where the time goes). Default OFF until proven.

- **Inc A — `GroupPmap` + a `--es-groups=G` knob, placement only.** Add the
  pmap class; at seeding and after each rotation, `copy` root `i` onto
  `group_pmap[i%G]`; re-map roots+Λ to the common pmap for §4.4. **Ground data
  still default-distributed** (so build/gamma will be misaligned → keep G=1 a
  no-op and only exercise the re-map plumbing first). Validate: G=1 is
  bit-identical to today; the re-map round-trip (localize→gather→localize) is
  itself bit-identical.
- **Inc B — NodeReplicate the ground data.** `replicate(NodeReplicated)` φ +
  V_local in `set_gs`. Now build + V0x are genuinely group-local. Validate
  bit-identity at G=2; `--es-time` to measure build-phase comm change + the
  §4.4 re-map cost. Watch the per-node memory (worldmem / MaxRSS).
- **Inc C — localize the exchange (gamma).** The hard one (§5.3): make the
  ground-exchange apply honor the group pmap, or accept its internal scheme and
  measure. Biggest prize (36%).
- **Inc D — scale study.** Sweep `G` and rank count on c2h4 / c6h6 (8 / 21
  occ); plot strong-scaling. This is where Stage 2 either justifies itself over
  Stage-1-on-default-pmap or doesn't — decide based on data.

**Proposed first step: Inc A**, because it builds the GroupPmap + re-map
plumbing with zero numerical risk (G=1 no-op, round-trip bit-identity) and lets
us measure the §4.4 re-map cost *before* committing to the harder Inc B/C.
```
