# 26 — γ batching (`compute_gamma_flat`): build spec

Status: DESIGN SPEC (propose-diff-first; diff to follow). Closed-shell. The
single highest-value, no-replication, per-iteration win — orthogonal to the
state-parallel scheduling (helps at G=1, G>1, and every save/load wave).

## 0. Why first
Strong scaling (c2h4, doc 25 data) showed the iteration is communication-bound —
adding nodes doesn't cut per-iter cost (8→32 ranks: 0.74× = slower). The per-iter
γ+build (the ~78%) is the cost. Batching attacks that cost directly, on a single
node, verifiable by bit-identity. The G>1 scheduler lands on top of a cheaper kernel.

## 1. The exchange contraction (verified, exchangeoperator.cc:182-187)
```
K[bra,ket](apply_to)_k = Σ_i ket_i · Poisson( bra_i · apply_to_k )
```
With `apply_to = φ`, the Poisson argument is `bra_i·φ_k`. So whether a term caches
hinges on whether its `bra` is φ (cacheable) or the response (not).

## 2. The γ term structure (per kernel — DO NOT over-claim caching)
- **TDA / ES** (`tda.hpp:148-157`): `pairs = {{φ, x}}` → ONE term `K[φ,x](φ)`,
  conv arg `φ_i·φ_k` → **CACHEABLE**.
- **Static FD** (`static.hpp:99-101`, comment 87-89 — "Y=X limit of Full, BOTH
  terms remain"): `{{φ,x},{x,φ}}` → `K[φ,x]` (direct, cacheable) **+** `K[x,φ]`
  (cross, conv arg `x_{s,i}·φ_k` → **NOT cacheable**).
- **Full FD** (`full.hpp:101-122`): γ_X = `K[φ,x]`+`K[y,φ]`; γ_Y = `K[φ,y]`+`K[x,φ]`.
  Direct (`K[φ,·]`) cacheable (same `g_ik` for X and Y); cross (`K[·,φ]`) not.

So: **TDA caches its whole exchange; FD/Full cache the direct half, batch the cross half.**

## 3. `compute_gamma_flat` (closed-shell), per iteration over the M-state bundle
γ_s = `Q( J[ρ_s]·φ − c_xc·(direct_s [+ cross_s]) )`.

**Coulomb** — batch the M Coulomb applies + the J·φ multiply:
```
J = apply(world, *coulop, rho_slab)        // M potentials, ONE wave
out[s] = mul_sparse(J[s], φ)               // J_s·φ_k  (per-state left operand)
```

**Direct exchange `K[φ,x](φ)` (all paths) — SHARED, cache-eligible:**
```
g_ik = Poisson(φ_i·φ_k)                     // φ-only, symmetric → n(n+1)/2 unique convs
exchange_dir_{s,k} = Σ_i x_{s,i} · g_ik     // per-state contraction, NO convolutions
```
`g_ik` depends only on φ ⇒ **cache per protocol** (φ fixed within a protocol; rebuild
on re-prep). TDA: this is the only exchange term ⇒ exchange ≈ free after iter 1.

**Cross exchange `K[x,φ](φ)` (FD/Full only) — NOT cacheable, batch per-iter:**
```
p_{s,i,k} = x_{s,i}·φ_k                      // M·n² products
P = apply(world, *coulop, flat(p))           // ONE Poisson wave (taskq-distributed)
cross_{s,k} = Σ_i φ_i · P_{s,i,k}            // contraction
```
Conv count stays M·n² (inherent — convolves the response), but as ONE wave instead
of M fenced per-state `apply_exchange` calls. (Full's `K[y,φ]` is the same with y.)

## 4. The no-mixing contract (the correctness invariant)
States are laid out state-major (`Xf[s*n+i]`); **`s` is NEVER a reduction index.**
- Direct: the convolution sees φ ONLY (no states) → cannot mix.
- Cross: products carry `(s,i,k)`; the contraction `Σ_i` is bounded to state s's
  own slice. `s` stays a free outer index.
- FORBIDDEN: handing the flattened `[x_0..x_{M-1}]` as a single Exchange ket — its
  internal `Σ_j` would run across states and sum them. `compute_gamma_flat` never
  does this; it builds the convolution from φ-only (direct) or per-(s,i,k) products
  (cross) and contracts within a state.

## 5. Memory knobs (gated)
- `g_ik` cache: `n²` distributed functions (`≈ n·|φ|` per-rank-ish, distributed).
  Two modes: (a) **per-iteration transient** (build/free each γ phase — M× fewer
  direct convs, low memory); (b) **per-protocol cache** (resident across iters —
  direct convs once per protocol, costs `~n²` resident). Tile `g_ik` in blocks to
  bound peak (mirrors MADNESS's exchange tiling). Default (a); (b) behind a flag.
- Cross products `M·n²` transient during the cross wave — tile by state if needed.

## 6. Verification
A/B bit-identity vs the per-state `apply_exchange` path (the `--es-batch`
discipline): same kernels, same inputs ⇒ reproduce to round-off. Single node, h2o
(TDA) and h2o-α (Static FD). Gate behind the existing `--es-batch` (extend it to γ)
or a new `--gamma-flat`; 0 = current per-state path = reference.

## 7. Where it slots
- `kernels/tda_batch.hpp`: add `compute_gamma_flat` (+ a `static_batch`/`full_batch`
  analog or a templated version covering the cross term).
- `solvers/es_solver.hpp` `step_*`: replace the per-state `compute_gamma` loop with
  `compute_gamma_flat` under the batch flag (the gamma pass at lines ~405-414).
- `solvers/fd_solver.hpp`: same for the FD γ (where Static/Full add the cross term).
- The per-protocol `g_ik` cache hangs off `ResponseGroundState` (built in
  `build_response_ground_state_*`, like the cached `K0`).

## 8. Expected impact
- **TDA/ES:** the exchange (a big chunk of the 36% γ) collapses to a per-protocol
  one-time `g_ik` build + cheap per-state contractions ⇒ γ ≈ Coulomb-only afterwards.
- **FD (static/dynamic):** direct half cached (M× fewer convs there) + cross half
  one wave instead of M fenced passes ⇒ fewer fences, better taskq overlap.
- Compounds with the G>1 state scheduler (doc 25): cheaper per-iter cost at any G /
  wave size; the memory-bounded `n_concurrent = min(M, ⌊budget/mem_per_state⌋)`
  scheduling (ES bundle + FD single-state save/load) sits above it.
