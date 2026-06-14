# 22 — Inc 3: Fan the ES build to subworlds (DRAFT for review)

Status: **DRAFT — design + skeleton, NO code applied.** Implements the
per-iteration fan/gather/couple/scatter of doc 21 §2.4 for the BUILD phase
only (BSH stays universe-side until Inc 4). Builds on Inc 1 (node-aligned
subworlds, commit 1cee76d36) + Inc 2 (one φ copy/node, ebb919f2e). This is the
first SOLVER-TOUCHING subworld change, so it lands only after review.

---

## 1. What changes, what doesn't

`step_rotate_pieces` today (per iter): top-of-iter Q+orthonormalize →
**per-root build loop (V0x/T0x/E0x/γ → Λ)** → subspace A,S,diagonalize → rotate
{X,V0x,E0x,γ} → θ → BSH → KAIN.

Inc 3 replaces ONLY the **bracketed build** with a MacroTask fan-out across
subworlds. The subspace / rotate / θ / BSH / KAIN stay **universe-side and
byte-identical** to today — that is deliberate: it confines the change to
"where the build ran" so bit-identity is easy to reason about (the universe
arithmetic is unchanged).

To keep the universe path identical, the fan-out returns, per root, the four
pieces the universe step consumes: **[Λ | V0x | E0x | γ]** (4·n_occ functions).
T0x and E0x_full are Λ-only and never leave the subworld. X is the *input*
(already universe-resident); we ship it out, build, ship the 4 pieces back.

---

## 2. The MacroTask build functor (skeleton)

New header `solvers/es_build_task.hpp` (closed-shell TDA first):

```cpp
class MacroTaskESBuildTDA : public MacroTaskOperationBase {
  long n_occ_;
  // one ROOT per batch: input stride n_occ, result stride 4*n_occ
  class RootPartitioner : public MacroTaskPartitioner { /* batch = n_occ funcs */ };
public:
  MacroTaskESBuildTDA(long n_occ) : n_occ_(n_occ) {
    partitioner.reset(new RootPartitioner(n_occ));   // see §4 (result strides)
    name = "MacroTaskESBuildTDA";
  }

  // batched first arg = ALL roots' X flattened (M*n_occ); trailing args are
  // Cloud-replicated whole to every subworld (the ground-state data).
  typedef std::tuple<
      const std::vector<real_function_3d>&,   // Xflat            [batched]
      const std::vector<real_function_3d>&,   // amo (φ)
      const real_function_3d&,                // V_local_alpha
      const Tensor<double>&,                  // focka
      const Tensor<double>&,                  // focka_no_diag
      const Tensor<double>&,                  // aeps
      const double&, const double&> argtupleT; // c_xc, lo

  using resultT = std::vector<real_function_3d>;        // [Λ|V0x|E0x|γ] per root
  resultT allocator(World& w, const argtupleT& a) const {
    std::size_t M = std::get<0>(a).size() / n_occ_;
    return zero_functions_compressed<double,3>(w, M * 4 * n_occ_);
  }

  // RUNS IN THE SUBWORLD (w = the subworld). args already Cloud-copied local.
  resultT operator()(const std::vector<real_function_3d>& Xbatch, /* one root */
                     const std::vector<real_function_3d>& amo,
                     const real_function_3d& Vloc,
                     const Tensor<double>& focka, const Tensor<double>& focka_nd,
                     const Tensor<double>& aeps, const double& c_xc, const double& lo) {
    World& w = amo.front().world();          // the subworld
    // Rebuild WORLD-BOUND objects subworld-side (they can't cross the Cloud):
    ResponseGroundState g0;
    g0.amo = amo; g0.V_local_alpha = Vloc; g0.focka = focka;
    g0.focka_no_diag = focka_nd; g0.aeps = aeps; g0.c_xc = c_xc; g0.lo = lo;
    g0.coulop = poperatorT(CoulombOperatorPtr(w, lo, 0.001*FunctionDefaults<3>::get_thresh()));
    g0.Qa = QProjector<double,3>(amo);
    g0.K0_alpha = nullptr;                   // => compute_V0x uses DIRECT apply_exchange
                                             //    (nesting Exchange's MacroTaskQ is forbidden)
    using K = Kernels<TDA,ClosedShell>;
    State x; x.x_alpha = Xbatch;             // this batch = one root
    auto rho = K::compute_density(w, g0, x);
    auto gamma = K::compute_gamma(w, g0, x, rho);
    auto V0x = K::compute_V0x(w, g0, x);
    auto T0x = K::compute_T0x(w, g0, x);
    auto E0f = K::compute_E0x_full(w, g0, x);
    auto E0x = K::compute_E0x(w, g0, x);
    auto lam = assemble_lambda(w, T0x, V0x, E0f, gamma);
    // pack [Λ | V0x | E0x | γ] -> 4*n_occ functions
    resultT out; out.reserve(4*n_occ_);
    append(out, lam.x_alpha); append(out, V0x.x_alpha);
    append(out, E0x.x_alpha); append(out, gamma.x_alpha);
    return out;                              // lands disjointly via result batch
  }
};
```

Key reuse: the kernels (`Kernels<TDA,ClosedShell>`, `assemble_lambda`) run
**unchanged** — `w` is the subworld, so every intermediate is subworld-local
and aligned (the whole reason for subworlds, doc 20 §5.4). This is the same
functor sketched pre-pivot; now it runs over the validated subworld machinery.

---

## 3. Solver integration (gated, reference preserved)

`es_solver.hpp`: add `int es_subworlds_ = 0` (+ setter, + `--es-subworlds=G`).
In `step_rotate_pieces`, when `es_subworlds_ > 0 && TDA/ClosedShell`:

```cpp
// replace the inline build loop with the fan-out
auto Xflat = flatten(out.roots);                       // M*n_occ
MacroTaskESBuildTDA op(n_occ);
MacroTask task(world_, op, MacroTaskQFactory(world_).set_nworld(es_subworlds_));
auto packed = task(Xflat, gs_.amo, gs_.V_local_alpha, gs_.focka,
                   gs_.focka_no_diag, gs_.aeps, gs_.c_xc, gs_.lo);  // gathered to universe
// unpack [Λ|V0x|E0x|γ] per root (stride 4*n_occ) into lambda[], V0x[], E0x[], gamma[]
```
Everything after (subspace, rotate, θ, BSH, KAIN) is the **current code,
untouched**. `if (es_subworlds_ == 0)` keeps the exact inline build (reference).

---

## 4. Open implementation points (the fiddly bits to nail)

1. **Custom partitioner result strides.** Default `MacroTaskPartitioner`
   batches the first vector by `max_batch_size` with input==result indexing.
   We need batch s → input `[s·n, (s+1)·n)`, result `[s·4n, (s+1)·4n)` (Batch
   has independent input/result `Batch_1D`; `insert_result_batch`). This is the
   one genuinely custom piece — mirror `MacroTaskPartitionerRow` but with the
   two strides. Verify against `macrotaskpartitioner.h`.
2. **World-bound objects.** `coulop`, `Qa`, `K0` cannot cross the Cloud —
   rebuilt subworld-side (above). `K0=nullptr` forces the DIRECT `apply_exchange`
   in `compute_V0x` (the exchange must not nest its own MacroTaskQ in a subworld
   task — doc 21 §3.2). The ground-exchange in compute_V0x then costs more per
   call (no cached operator); acceptable for Inc 3, revisit in Inc 4.
3. **Trailing-arg replication.** amo/V_local/focka ship to every subworld via
   the Cloud each call. With node-aligned subworlds (Inc 1) this is one
   distributed φ copy per node (Inc 2). Per-iter re-ship cost = the coupling
   unknown to measure.

---

## 5. Staging

- **Inc 3a — stock MacroTaskQ (interleaved subworlds), G = nnodes.** Validates
  the FLOW + bit-identity + Cloud coupling cost + φ copy COUNT (G copies, same
  as node-aligned). Reuses MacroTaskQ as-is — no core change. (Interleaved
  subworlds span nodes, so communication isn't yet node-local, but the
  arithmetic and memory count are already right.)
- **Inc 3b — node-aligned subworlds.** Add a `create_worlds` variant using
  `Split_type(SHARED)` (small MADNESS-core addition; Inc 1's
  `make_node_aligned_subworld` is the standalone proof) and feed it to the
  MacroTaskQ so each subworld = one node → communication locality + on-node φ.

---

## 6. Validation

- **Bit-identity:** `verify_es_batch.sh` (converged-root ω) with `--es-subworlds=2`
  vs reference (`0`), at 1 node and 2 nodes. Build-in-subworld must reproduce the
  universe build to round-off (it's the same kernels on the same inputs).
- **Cost:** `--es-time` — the build phase now includes fan/gather; compare its
  wall to the inline build. This quantifies the Cloud coupling (doc 21 §3.1) —
  the deciding measurement for whether subworld fan-out pays off, and at what
  node/root count.

---

## 7. Risks

- The result-stride partitioner is the main correctness risk (a wrong stride →
  pieces land on the wrong root → garbage; the bit-identity gate catches it).
- Per-iter Cloud coupling could dominate at small scale (expected; it's a
  multi-node bet — measure, don't assume).
- `MacroTask` lifecycle inside an iterative solver (created/destroyed each
  iter): confirm no leak / no residual fence state across iterations.
```
