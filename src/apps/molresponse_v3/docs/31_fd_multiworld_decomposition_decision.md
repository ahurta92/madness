# 31 — FD state-parallel: the multi-world decomposition decision

Status: **DECISION** (resolves the doc 24 vs doc 25 fork). Closed-shell FD first.
Supersedes the *framing* of docs 19–25; their proven primitives (Inc1/Inc2/
keystone/S1, the feasibility study) are reused unchanged. Drafted 2026-06-22.

> Doc-number note: 26–30 are claimed by other threads (26/27/28 exchange-tensor,
> 29 perf-model, 30 io-hdf5). This thread's docs continue at 31.

## 0. What this resolves

The parallel-runtime fork: **doc 24** (persistent node-aligned subworlds for the
ES bundle) vs **doc 25** (a single tunable `G` knob, with the `Group` API instead
of subworlds). Both asked "how do we parallelize *one* solve across nodes?" — and
the measurement says that is the wrong question. The right architecture is a
**two-axis decomposition with one solve per World**, and the cleanest first target
is **FD, not ES**.

## 1. The measurement that decides the fork

`strong_scale_es.sh` (c2h4, 6 roots, fixed maxiter, single rung 1e-4), aggregated
by `aggregate_strong_scale.py`, 2026-06-15:

| ranks | nodes | gamma (s) | total (s) | gamma-speedup | gamma-eff% |
|------:|------:|----------:|----------:|--------------:|-----------:|
| 8 | 1 | 8.03 | 25.15 | 1.00 | 100 |
| 32 | 4 | 10.86 | 31.45 | **0.74** | **18** |

(The 2-node point crashed, rc=134.) Spreading **one** ES bundle from 8→32 ranks
made the per-iteration cost *worse*. **G=1 strong scaling walls** — exactly doc 25
§4's failure mode. But this is *expected and correct*: 8 occupied orbitals over 4
nodes is 2 occ/node, far past the granularity floor. We do not want to push a
single state's spatial decomposition past ~**5 occupied orbitals per node**; comm
dominates below that. So the wall is a **boundary to respect**, not a bug to fix —
and it tells us the parallelism must come from a *different* axis.

## 2. The reframing — two axes, one solve per World

| Axis | What it partitions | Coupling | Has a strong-scaling ceiling? |
|---|---|---|---|
| **State `S`** | states across Worlds (one solve per World) | **FD: none.** ES: thin M×M subspace | no — embarrassingly parallel |
| **Spatial `R_state`** | ranks/nodes *within* one state's World (φ + response distributed by the pmap) | n/a | **yes** (the §1 wall; ~5 occ/node) |

`P ≈ S × R_state`. **The unit of work is one state per World; you size the World
to the state, then fill the machine with states.** This is v2's `state_parallel on`
done *right*: v2 fixed `R_state = 1 rank` → φ fully replicated per rank → OOM at
c6h6. The fix is `R_state ≥ 1 node` so φ is **distributed across the World's ranks**
(per-rank φ = |φ|/R_state — the feasibility study's exact 8× win, see doc 23 §0).
doc 24's subworld machinery is the *mechanism* for those Worlds; doc 25's `G` is the
*state axis*. Neither was wrong — they were each one axis of a 2-D grid.

## 3. The decision

1. **FD before ES.** FD states are *fully independent* (`calc_manager.hpp:180`:
   "distinct perturbations are independent and share every wave"). One World per FD
   state ⇒ **pure fan-out → solve → gather → fence, no per-iteration cross-world
   traffic at all.** The ES coupling (M×M allreduce + X re-broadcast, doc 24 §2) does
   not exist for FD. FD is the cleanest proof of the multi-world architecture and the
   highest-value win for the scaling goal.
2. **First cut = one FD state per node-aligned subworld** (the *medium-regime*
   decomposition, §4). Correct in every regime, optimal for medium, and built on the
   only proven primitive (`make_node_aligned_subworld`, Inc1; GS ship-in, S1).
3. **Simplest-first, measurement-driven.** The two refinements below are **deferred
   and gated on data** — we test into them, we do not build them up front:
   - **sub-node packing** (>1 state per node, `R_state` < 1 node) — a *small-regime*
     throughput optimization (φ is cheap, so several φ copies fit per node);
   - **multi-node-per-state** (`R_state` > 1 node) — a *large-regime* necessity (φ
     does not fit one node at the production protocol).

## 4. Regime map — pinned to the memory model, not atom count

The boundaries are where **per-state memory at the production protocol** crosses the
rank/node budgets (Xeon Max ≈ 576 GB/node; target ≤ 72 GB/rank at 8 ranks/node). We
have the *shape* (`rss ∝ n_occ`, ×k-growth, R4 study) and the *magnitude* (v2 MaxRSS,
root CLAUDE.md). The exact GB cutoffs are calibrated by the weak-scaling sweep (§7).

| Regime | Condition (per-state mem) | Decomposition | First-cut behavior | Today's molecules |
|---|---|---|---|---|
| **Small** | ≤ per-rank budget | `R_state`=1 rank; pack states/node | runs (1 state/node = over-decomposed but correct) | h2o (5), c2h4 (8) |
| **Medium** | per-rank < mem ≤ per-node | `R_state`=1 node; `S`=#nodes | **optimal** | ch3oh (9), c6h6 (21), naph (34)* |
| **Large** | > per-node budget | `R_state`=⌈n_occ/5⌉ nodes; `S`<#nodes | aborts pre-flight (§6) until R_state>1 lands | very large, or any mol at high k |

\*regime of c6h6/naphthalene at the real protocol is TBD by the §7 sweep.

**v2's OOM was the small-regime decomposition (1 state/rank, φ replicated) applied
to a medium molecule.** The fix is to *match the decomposition to the regime* — which
is the whole point of this design.

## 5. Mechanism — and why FD is the easy case

The code already anticipated this exact work:

- **The scheduler is already perturbation-parallel.** `CalcManager::build` emits a
  DAG; `run()` (`calc_executor.hpp:903`) schedules waves of independent states
  (`calc_manager.hpp:68`). FD states across distinct perturbations share every wave.
- **`run()` itself names the next step.** `calc_executor.hpp:859-867`: today it is
  "single-group — the whole communicator solves ONE state at a time"; "distributing a
  wave's states to rank SUBGROUPS sized to fit memory is the **15c design
  (STATE_PARALLEL / subworlds)**." The hard constraint stated there: *all ranks in a
  World must agree on which states they solve, or they deadlock on mismatched
  collectives.* Node-aligned subworlds satisfy it by construction — all ranks on a
  node solve that node's assigned states.
- **Binding a solve to a subworld is a constructor argument, not surgery.**
  `ExecutorContext` carries a borrowed `madness::World &world` (`calc_executor.hpp:131`);
  `solve_fd_protocol<Type,Shell>` (`:254`) and `FdResponseExecutor` (`:759`) take their
  World through it. Run a state in a subworld = build an `ExecutorContext` bound to
  that subworld with the GS shipped in (S1-proven).
- **Proven primitives, reused as-is:** `make_node_aligned_subworld`
  (`node_subworlds.hpp:52`, Inc1) for the World pool; GS ship-in under the
  `set_default_pmap(subworld)` discipline (S1, doc 23 §3); teardown-before-finalize
  (subworld reset before universe). **No allreduce, no re-broadcast, no Cloud** — FD
  is fan-out/gather only.

**The FD state-parallel loop (gated variant of `run()`):** build the node-subworld
pool once per protocol; for each wave, partition its WorkItems round-robin by node
index; each subworld solves its assigned FD states to convergence in its own World
(GS shipped in); gather per-state results + metadata to the universe; fence; next
wave. Metadata still flows through the response_metadata layer (sharded per group,
merged on rank 0 — the v2 subgroup pattern).

## 6. The flexible user solution

- **Knob:** `--fd-subworlds=G` (0 = current single-World path = the bit-identical
  reference; default). Mirrors how `--es-batch` / `--es-subworlds` thread through
  `ExecutorSettings → ExecutorContext`.
- **Auto-selector:** from `(n_occ, k, available nodes, mem budget)` pick `(S, R_state)`
  via the §4 regime model; manual override stays available.
- **Pre-flight memory estimate + clean abort** (already on the root CLAUDE.md priority
  list): before allocating the World pool, compute per-rank φ + response estimate from
  the mem model and abort with a message if the chosen `(S, R_state)` overflows the
  budget — far better than a silent SIGKILL. This is also what gates the "large" regime
  until `R_state > 1 node` lands.

## 7. The calibration experiment — weak scaling (new; the metric we actually want)

Strong scaling (fixed problem, more ranks) is the *wrong* metric here — §1 shows its
ceiling and we respect it. The right metric is **weak scaling**: hold work-per-World
≈ constant (1 state/node, fixed occ/node) and grow molecule + states + nodes
*together* — does wall-time stay flat? Design `weak_scale_fd.sh` as a sibling of
`strong_scale_es.sh`: sweep `(molecule, nodes)` keeping states/node fixed, capture
per-state `wall_s` + `MEMORY_HWM`, aggregate to a flat-line / efficiency plot. This
sweep also calibrates the §4 GB thresholds (where per-state mem crosses rank/node
budgets at the real protocol).

## 8. Next step + scope

**Next:** doc 32 — the FD state-parallel implementation plan (propose-diff-first),
staged like the ES work was:
- **F1** — standalone fan-out/gather proof (no solver change): partition a fixed set
  of independent FD states across the node-subworld pool, solve each in its World,
  A/B the converged α tensor vs the single-World path → bit-identical (converged FD
  states are path-stable; the `verify_es_batch.sh` discipline). The gate before any
  `run()` edit.
- **F2** — the gated `run()` variant + `--fd-subworlds` + pre-flight abort. A/B on
  converged states, 1 and 2 nodes.
- **F3** — weak-scaling sweep (§7) to calibrate regimes; then the auto-selector.

**Do-not-touch (this thread's standing list):** the reference kernels
(`compute_V0x`/`compute_gamma`/`compute_E0x` — exchange's gate-0 oracle); the
single-World `run()` path stays the `G=0` bit-identical reference; perf-model's
meter/profile schema; viz's `dump_mra_trees`/legacy export. Standing contracts:
collectives on all ranks (print-gate only); never write `response_metadata.json`
directly; **USER runs builds/solves on the alloc**; propose-diff-first on
solver/runtime edits; commit with `git -c core.hooksPath=/dev/null`.

**Deferred (test into, do not pre-build):** sub-node packing (small regime,
generalize to `make_subworld(R)`); multi-node-per-state (large regime); ES
state-parallel (reuses this World pool + the keystone allreduce, doc 24 §2).

**Cross-thread note:** the γ-batching idea from doc 25 §5 (batch the per-state
exchange in one World) lives in the **exchange** thread (doc 26 gamma_batching_design),
not here — it helps the G=1 path and composes with, but is independent of, this
state-axis decomposition.
