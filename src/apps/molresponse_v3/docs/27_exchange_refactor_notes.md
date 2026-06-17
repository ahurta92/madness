# 27 — Exchange refactor notes: the `K[a,b](c) = ⟨a | g(b,c)⟩` decomposition

Status: NOTES (toward the refactor; no code). Closed-shell. Working notes with
@user — formalize every two-electron/exchange op in one form, classify by whether
its convolution tensor caches, catalog the distinct `g(b,c)`, and find sharing.

## 1. The unified form
Verified contraction (exchangeoperator.cc:182-187):
`K[bra,ket](apply_to)_k = Σ_i ket_i · Poisson(bra_i · apply_to_k)`.

Write it as **`K[a,b](c)`** with `a = ket` (the contraction weight), `b = bra`,
`c = apply_to`, and the **convolution tensor** `g(b,c)_{ik} ≡ Poisson(b_i · c_k)`:
```
K[a,b](c)_k = Σ_i a_i · g(b,c)_{ik}          (an inner product of a with columns of g)
```
- The **convolution** `g(b,c)` (n² Poisson applies) is the expensive part — communication-heavy.
- The **contraction** `Σ_i a_i·g` is a cheap mul+reduce — local-ish.
Caching/sharing is entirely about the `g(b,c)` tensors.

## 2. Catalog — every closed-shell exchange term as K[a,b](c)
ρ_s = 2 Σ_p φ_p (x_p [+ y_p]); J[ρ_s] = Poisson(ρ_s) is the Coulomb/Hartree piece
(1 conv/state, separate). The exchange-like terms:

| # | term (where) | a (ket) | b (bra) | c (apply_to) | g(b,c) = Poisson(b·c) | result_k |
|---|---|---|---|---|---|---|
| 1 | V0x exch, X — `K[φ,φ](x)` (all paths) | φ | φ | x | `Poisson(φ_i·x_k)` = **T_x** | Σ_i φ_i T_x[i,k] |
| 2 | V0x exch, Y — `K[φ,φ](y)` (Full) | φ | φ | y | `Poisson(φ_i·y_k)` = **T_y** | Σ_i φ_i T_y[i,k] |
| 3 | γ direct, X — `K[φ,x](φ)` (all) | x | φ | φ | `Poisson(φ_i·φ_k)` = **g₀** | Σ_i x_i g₀[i,k] |
| 4 | γ direct, Y — `K[φ,y](φ)` (Full) | y | φ | φ | **g₀** (same!) | Σ_i y_i g₀[i,k] |
| 5 | γ cross (Static `K[x,φ](φ)`; Full γ_Y) | φ | x | φ | `Poisson(x_i·φ_k)` = **T_xᵀ** | Σ_i φ_i T_x[k,i] |
| 6 | γ cross, Full γ_X — `K[y,φ](φ)` | φ | y | φ | **T_yᵀ** | Σ_i φ_i T_y[k,i] |

(TDA = term 3 only. Static = 1,3,5. Full = 1,2,3,4,5,6.)

## 3. Two classes (by the convolution tensor g(b,c))
- **Class 1 — `g(b,c)` is φ-only ⇒ CACHEABLE.** Only `g₀ = Poisson(φ_i·φ_k)`
  (terms 3,4): both b and c are φ; the response enters only as the contraction
  weight `a`. φ is fixed per protocol ⇒ **cache g₀ once per protocol**; X and Y both
  reduce against the *same* g₀ columns. (Your "x and y do a local reduction with
  every column of g[φ,φ].") Symmetric ⇒ n(n+1)/2 unique convs.
- **Class 2 — `g(b,c)` contains a response orbital ⇒ NOT cacheable.** `T_x`,`T_y`
  (and their transposes): rebuilt every iteration (φ·x, φ·y change). n² convs each.

## 4. The sharing insight (why the K[a,b](c) view pays off)
`g(φ,x)` and `g(x,φ)` are **the same tensor transposed**:
`g(x,φ)_{ik} = Poisson(x_i·φ_k) = Poisson(φ_k·x_i) = T_x[k,i] = T_xᵀ_{ik}`.
So **one build of `T_x = Poisson(φ_i·x_k)` (n² convs) serves TWO exchange terms** —
the V0x ground exchange (term 1, contract `Σ_i φ_i T_x[i,k]`) AND the γ cross
(term 5, contract `Σ_i φ_i T_x[k,i]`). Likewise `T_y` serves terms 2 and 6.

Convolution-build count per state per iteration (Class 2), today vs fused:
| path | today (separate apply_exchange per term) | fused (build T_x[,T_y] once) |
|---|---|---|
| Static | term1 (n²) + term5 (n²) = **2n²** | **T_x: n²** |
| Full | terms 1,2,5,6 = **4n²** | **T_x + T_y: 2n²** |
Plus Class 1: g₀ once per protocol (vs every iter, every term, every state today).

So the refactor target: **build {g₀ (cached), T_x, T_y} once per state-iteration,
then express terms 1-6 as cheap contractions of those tensors.** Cuts Class-2
convs ~2×, and collapses Class-1 to a per-protocol one-time build.

## 5. Data movement of the two classes
- **Class 1 (g₀ cached):** per iter = contractions only (`Σ_i a_i g₀[i,k]`, a = x/y).
  Mul+reduce over the occupied index; g₀ resident & distributed. Cheap comm; the win
  is it removes convolutions from the hot loop entirely (TDA: all exchange; FD: half).
- **Class 2 (T_x, T_y):** the n² Poisson convolutions ARE the cost (internode
  convolution traffic). Fuse across terms (build once, contract twice). Batch over
  states into single waves (no per-term/per-state fences). This is the bulk of the
  remaining exchange comm and where the G>1 grouping (doc 25) localizes it.

## 6. What to measure (timing plan — instrument each type individually)
Add per-type timers (extend the `--es-time` `ES_TIMING` line) so reductions are
measurable, not assumed:
- `coulomb`  : Poisson(ρ_s) + J·φ (Hartree).
- `g0_build` : g₀ convolutions (amortized — report once per protocol).
- `g0_contract` : the Σ_i a_i g₀ reductions (Class 1, terms 3/4).
- `Tx_build`, `Ty_build` : the Class-2 response convolutions (terms 1,2 ≡ 5,6).
- `tx_contract` : the contractions for terms 1,2,5,6.
Baseline = current per-state path; then (a) cache g₀, (b) fuse T_x/T_y across
V0x+γcross, (c) batch over states — each measured as a delta.

## 7. Refactor shape (two problem classes)
Split the build into:
- a **cached/ground block** (g₀ per protocol; the V_local, focka, K0 already-cached
  pieces) — Class 1; built once per protocol, contracted per iter.
- a **per-iteration response block** (T_x, T_y; the Coulomb J[ρ_s]) — Class 2; built
  once per state-iteration, shared across the terms that need it, batched over states.
The iteration then = "contract cached-block columns with x/y" + "build response-block
tensors, contract". The K[a,b](c) catalog (§2) is the single source of which tensor
each term needs. Open: tracking ALL g(b,c) forms across FD/ES/VBC (β/Raman reuse the
same machinery) to find further shared tensors. Feeds: doc 26 (compute_gamma_flat is
the first concrete slice — g₀ cache + cross batching), doc 25 (G>1 localizes Class-2 comm).
