# 25 — ES parallelism: the formulaic approach (runtime-grounded redesign)

Status: DESIGN THINKING (supersedes the *premise* of docs 21–24). Reconsiders ES
parallelism from the MADWorld runtime + the constraints: (a) do NOT replicate φ,
(b) parallelize the loop over states, (c) only the diagonalize couples states.

## 1. Problem structure

One ES iteration over M states sharing one ground state φ:
- **Per-state, independent (the ~99%):** build Λ_s (T0x,V0x,**γ=Coulomb+exchange**,E0x),
  θ_s, BSH_s, KAIN_s. γ (exchange) is the ~36% time-limiter.
- **Couples all states (the ~0.5%):** the M×M subspace `A=⟨X_i|Λ_j⟩`, `S`, diagonalize,
  rotation. A tiny reduction.

So the parallelism to exploit is the embarrassingly-parallel state loop; the only
genuine synchronization is an M×M scalar reduce.

## 2. What the MADWorld runtime actually gives us (code-grounded)

- A single Function's tree is distributed across **all** ranks by a pmap
  (`Key→owner`); work routes to the owner (`funcimpl.h` `task(coeffs.owner(key),…)`).
  So **φ in one World is one distributed copy: per-rank `|φ|/P`. It cannot OOM.**
- A vecfunc op over a flattened bundle is **one parallel wave + one fence**
  (`vmra.h` `apply`: submit all `apply_only(...,false)`, then a single `gop.fence()`).
  Flattening the M·n_occ bundle ⇒ the taskq spreads it over all ranks+threads. This is
  what `--es-batch` already does for V0x/T0x/BSH.
- `matrix_inner` is local-compute + `gop.sum` allreduce — the M×M coupling, native+cheap.
- **pmap menu:** `WorldDCDefaultPmap` (hash→all ranks: max parallelism, max comm),
  `WorldDCLocalPmap` (RankReplicated), `WorldDCNodeReplicatedPmap` (one/​node),
  `LBDeuxPmap` (load-balanced subtrees: locality). `LoadBalanceDeux` redistributes a
  tree weighted by a cost fn (SCF.cc `initial_load_bal`, mp2.cc) → balances work +
  subtree locality (less convolution comm).
- **Mixing pmaps in one op:** `matrix_inner`/`inner` tolerate different pmaps (local +
  allreduce). `gaxpy`/`transform`/`merge_trees` **require the same pmap** (asserted).
  So you cannot freely pin different states to different rank-subsets and still gaxpy
  φ·x_s across them.
- **`Group` (world/group.h):** subgroup bcast/reduce/allreduce over a rank-subset
  **without a separate World** — lightweight, no `MPI_Comm_create`. This is the right
  primitive for "processing groups," not subworlds.
- **MacroTaskQ ships φ by deep-copy into each subworld (`copy(subworld,φ)`), no
  cross-iter persistence** — i.e. it REPLICATES φ. Not what we want for an iterative
  solver that reuses φ every iteration.

## 3. The γ/exchange batching reality

γ_s = Q( J[ρ_s]·φ − c_xc·K[φ, x_s]·φ ).
- **Coulomb:** ρ_s is per-state (density `dot` can't batch numerically); but
  `J[ρ_s]·φ` is a per-function map → the **mul over the M·n_occ bundle batches** like
  V0x. Cheap.
- **Exchange:** `K[φ, x_s]` is a **distinct operator per state** (ket = x_s). A single
  `Exchange` over a flattened `[x_0..x_{M-1}]` ket would SUM across states (wrong). So
  exchange is **M independent operator-applies that share φ** — not one bundle apply.
  The lever is to run those M applies **concurrently** (one parallel wave, not M fenced
  passes) so the taskq overlaps them across ranks, with φ as one distributed copy.

## 4. The formulaic approach — ONE knob: G (number of processing groups)

Every option is a point on a single axis: **partition the M states into G groups; each
group does its states' builds; the M×M subspace is a global allreduce (the keystone).**
G sets the memory↔communication trade:

| | per-rank φ memory | exchange communication | state parallelism |
|---|---|---|---|
| **G = 1** (single World) | `|φ|/P` (one distributed copy, NO replication) | **global** (all-to-all across all P ranks) | flattened bundle over all ranks |
| **G = #nodes** | `|φ|/(P/G)` = one copy/node (G copies) | node-local | each node owns M/G states |
| **G = P** (v2 "state_parallel on") | `|φ|` per rank (P copies) → **OOM** | rank-local | one state per rank |

The unavoidable theorem (from §2's same-pmap rule): **localizing a state's exchange to a
group requires φ present in that group** — because `K[φ,x_s]` gaxpy/mul's φ against x_s,
and those operands must share a pmap. So **communication-locality ⟺ per-group φ copies**.
You cannot have both "no φ replication" and "localized exchange comm" — G is exactly the
dial between them. `per-rank φ ≈ |φ|·G/P`, and comm-locality rises with G.

This reframes docs 21–24: they are **the G = #nodes point** (one φ copy/node) — not a
separate paradigm, just one setting of G, and they used heavyweight subworlds where the
**`Group` API would do the per-group A/S allreduce without separate Worlds**.

**Corollary — what the user wants (no replication) = G = 1.** Then φ is one distributed
copy and we parallelize states by flattening the build (incl. γ) so the taskq spreads
M·n_occ work over all ranks. The open question is purely empirical: **does G=1 strong-
scale across nodes, or does the global exchange/convolution comm wall it?** That is the
strong-scaling experiment (`strong_scale_es.sh`). The answer picks G:
- G=1 scales → done, no replication, simplest. Extend `--es-batch` to γ; use
  `LoadBalanceDeux` on φ for locality; reduce fences.
- G=1 walls → raise G to the **smallest** value that scales within the memory budget,
  implemented with `Group` subgroup-collectives (not subworlds): partition states,
  redistribute φ to G group-pmaps (G copies), per-group build, **global allreduce of the
  M×M A/S** (keystone), global rotation, re-scatter. This is docs 21–24 generalized and
  made pmap-correct, with G a tunable (1 → #nodes), not hardwired to subworlds.

## 5. Recommended path

1. **Measure G=1 strong scaling NOW** (`strong_scale_es.sh`, c2h4, 1/2/4 nodes,
   single World, `--es-batch`, `--es-time`): does γ/build speed up with nodes? →
   `aggregate_strong_scale.py` reports γ-speedup + efficiency. This decides G.
2. **Batch γ in one World** (the time-limiter, G=1): a `compute_gamma_flat` mirroring
   `tda_batch.hpp` — batch the Coulomb `mul` over the M·n_occ bundle, and issue the M
   exchange applies as one un-fenced wave (overlap across ranks) instead of M fenced
   per-state passes. φ stays one distributed copy. A/B bit-identity vs the per-root path
   (the `--es-batch` discipline). This is the single highest-value, no-replication win.
3. **Load-balance φ** with `LoadBalanceDeux` (cost-weighted redistribute) to cut
   convolution/exchange internode traffic; batch compress/reconstruct fences.
4. **Only if (1) walls:** introduce G via the `Group` API (subgroup A/S allreduce),
   tuning G up from 1 to trade memory for comm-locality. Reuse the proven keystone
   allreduce; drop the subworld machinery.

## 6. What I'd have designed cold (the honest answer to "no context")

Single World, φ one distributed copy, parallelize states by flattening the **entire**
build (incl. γ) into bundle ops so the MADWorld taskq+pmap distribute M·n_occ work over
all ranks/nodes; the M×M diagonalize as a global `matrix_inner` allreduce; `LoadBalanceDeux`
for locality. Reach for processing groups (`Group`, not Worlds) **only** if measured
strong scaling forces it — and then as a tunable G, accepting `|φ|·G/P` memory as the
explicit price of comm-locality. No per-rank replication ever.
