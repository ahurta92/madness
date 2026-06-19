# Operator Contracts — molresponse_v3

The mathematical definition of every response operation, the **metric /
sign convention** it uses, and the **legacy source line** that is the
ground truth for it.

## Why this file exists

`molresponse_v3` is a refactor of a working legacy implementation
(`src/apps/molresponse/`). The legacy code is the base truth. A refactor
is only safe if every operation still matches that truth — and the
dangerous bugs are the ones in *conventions* (a sign, a metric, a factor)
that compile fine and converge to a plausible-but-wrong number.

We shipped exactly such a bug: the RPA subspace overlap and the
response-vector normalization were written with the **Euclidean** metric
`⟨X|X⟩ + ⟨Y|Y⟩` instead of the **indefinite** RPA metric
`⟨X|X⟩ − ⟨Y|Y⟩`. It compiled, it "converged," and it passed a loose
ω-tolerance check — because nothing tied the code to the legacy formula
tightly and automatically.

**How to use this file:**
- When you touch an operation, check the code against the formula here,
  **and** check the formula here against the cited legacy line.
- When you add an operation, add a row with its ground-truth line.
- If the code and the cited legacy line disagree, the code is wrong until
  proven otherwise — legacy is the reference.

Notation: `X` = excitation block, `Y` = de-excitation block. `φ` =
occupied ground-state orbitals. `⟨a|b⟩` over vecfuncs means
`Σ_p ⟨a_p|b_p⟩` (sum over orbitals). Closed-shell unless a β block is
named.

---

## Metrics and inner products

| Operation | Formula | Convention | Ground truth (legacy) | v3 location |
|---|---|---|---|---|
| Euclidean bundle inner | `⟨a\|b⟩ = Σ_p ⟨a_p\|b_p⟩` | positive-definite | `response_functions.h` `response_space_inner` | `rs::inner` (response_space_ops.hpp) |
| **RPA subspace overlap `S`** | `S_ij = ⟨Xi\|Xj⟩ − ⟨Yi\|Yj⟩` | **INDEFINITE** (Y block negative) | `ExcitedResponse.cpp:871` | `rs::metric` |
| **Response-vector norm** | `‖v‖ = √(⟨X\|X⟩ − ⟨Y\|Y⟩)` | **INDEFINITE** | `ResponseBase.cpp:2097-2099` | `rs::metric_inner` |
| TDA norm (x-only) | `‖v‖ = √(⟨X\|X⟩)` | standard | `ResponseBase.cpp:2077-2088` | `rs::metric_inner` (no Y block) |

The single most error-prone convention in the whole module: **the Y
(de-excitation) block enters the overlap and the norm with a minus
sign.** For x-only (TDA) states there is no Y block, so the indefinite
metric collapses to the standard `⟨X|X⟩` — `rs::metric` and
`rs::metric_inner` are written so this falls out automatically (`if
constexpr` on block presence).

---

## Subspace eigenproblem

| Operation | Formula | Convention | Ground truth | v3 location |
|---|---|---|---|---|
| subspace matrix `A` | `A_ij = ⟨Xi\|Λxj⟩ + ⟨Yi\|Λyj⟩` | Euclidean (+) | `ExcitedResponse.cpp:885` | `ESSolver::step_*` (`rs::inner(roots, lambda)`) |
| generalized eig | `A U = S U Ω` | `S` indefinite, solved with `sygvp` | `ExcitedResponse.cpp:1030` (`excited_eig`) | `rs::diagonalize` |
| eigenvalue ordering | dominance-swap → phase fix → cluster-unmix → **ascending sort** | stable ascending | `ExcitedResponse.cpp:1046-1122` | `rs::diagonalize` steps 3-5.5 |

Note the asymmetry: `A` uses the **Euclidean** (+) inner, only `S` uses
the **indefinite** metric. `sygvp` (symmetric-definite) is valid because
`S = ⟨X|X⟩ − ⟨Y|Y⟩` stays positive-definite for excitation-dominated
states (‖X‖ > ‖Y‖); it only loses definiteness at a triplet instability.
The **final ascending sort** is required for slot-stable KAIN history —
without it the greedy dominance swap flips slot identity between
iterations for near-degenerate roots (`ExcitedResponse.cpp:1122`
`sort_eigenvalues`).

---

## Per-root building blocks (kernels)

| Operation | Formula | Notes | v3 location |
|---|---|---|---|
| response density | `ρ = 2 · Σ_p φ_p (x_p + y_p)` | factor 2 = spin; `(x+y)` folds X,Y | `Kernels<Full,ClosedShell>::compute_density` (full.hpp:44) |
| Coulomb-exchange `γ` | `γ_X = J[ρ]φ − c_xc(K[φ,x]φ + K[y,φ]φ)` | **cross-channel** K[y,φ] couples X↔Y | `compute_gamma` (full.hpp); legacy `ResponseKernel.hpp:113` |
| `V0·x` | `V_local·x − c_xc K[φ,φ]·x` | acts on X and Y independently | `compute_V0x` |
| `E0·x` (off-diag) | `F_offdiag · x` | diag ε absorbed into BSH | `compute_E0x` |
| `E0·x` (full) | `F · x` | full Fock; for Λ assembly | `compute_E0x_full` |
| `T0·x` | `−½∇²·x` | kinetic | `compute_T0x` |
| **θ** (BSH driver) | `θ = V0x − E0x + γ` | FD adds `+ V_p` source | `assemble_theta` (assembly.hpp); legacy `ResponseKernel.hpp:271` |
| **Λ** (subspace) | `Λ = T0x + V0x − E0x_full + γ` | full Fock (not off-diag) | `assemble_lambda` (assembly.hpp); legacy `ResponseKernel.hpp:397` |
| BSH apply | `x_new = Q(BSH(ω)·(−2(θ + shift·x)))` | paired ±ω for Full (X:+ω, Y:−ω) | `bsh_apply` |

Density factor convention: `spin_factor × y_factor`. Restricted Static =
4, restricted Full = 2, unrestricted Static = 2, unrestricted Full = 1.
`alpha_factor = −(spin_factor × y_factor)`.

---

## RPA symmetric reduction — REMOVED (2026-06)

The `ESSolverFullRPA` symmetric-reduction solver (`(A−B)(A+B)u = ω²u`,
`u = X+Y`) and its `apply_AplusB` / `apply_AminusB` operators were
removed — not a direction we are pursuing. Full ES is the direct
paired-(X,Y) `ESSolver<Full, ClosedShell>` only.

---

## Top-of-iteration discipline

| Operation | Formula | Ground truth | v3 location |
|---|---|---|---|
| project | `Q·v` per spin (alpha→Qa, beta→Qb) on every block | `ExcitedResponse.cpp:2483-2517` | `rs::project` |
| orthonormalize | modified Gram-Schmidt in the **RPA metric** | `ResponseBase.cpp:2091` (normalize) | `rs::orthonormalize` |

`orthonormalize` MUST use `metric_inner` (indefinite), not the Euclidean
inner — so the bundle entering `diagonalize` is orthonormal in the same
metric the subspace overlap `S` uses (→ `S ≈ I`). This is the line that
had the bug.

---

## How this is enforced automatically

A contract table is a human aid; the gate is the test. See
`tests/test_v3_es_skeleton.cpp`:

1. **Convergence gate** — a run only PASSES if it actually converged
   (`!diverged` and `max_residual ≤ solver target`), not merely if ω
   landed near the reference. The metric bug stalled at a residual floor;
   this gate fails that.
2. **Legacy ω comparison** — converged ω must match the tabulated legacy
   value within `tol`. Legacy is an independent implementation, so a
   match is real signal.
3. **Two-solver cross-check** — retired with `ESSolverFullRPA` (removed
   2026-06). The remaining independent check is item 2 (legacy ω).

When adding an operation, add its row above **and** a check that ties it
to legacy — the table alone did not prevent the metric bug; the missing
piece was the automated tight comparison.

---

## Performance profile schema (v1)

Pinned by the **perf-model** thread (design: `docs/29_perf_model_design.md`).
The cross-thread contract: the machine-readable profile that
`WorldProfile::dump_json` emits, that `exchange` reports its Tx/tile counts +
phase timings into, and that `parallel-runtime` reads to settle the doc-24-vs-25
fork with measured numbers. **Stable key names — coordinate changes here.**

Emitted once per run by rank 0, **only** when built with `ENABLE_WORLD_PROFILE=ON`
**and** env `MADQC_PROFILE_JSON` is set (else absent — zero-effect contract).
Every per-phase statistic is a faithful copy of a reduced `WorldProfileEntry`
field: `{sum, min, max, pmin, pmax}` (`pmin`/`pmax` = the rank holding the min/max).

```jsonc
{
  "schema_version": 1,
  "world_size": 6,                 // P (MPI ranks)
  "total_cpu_s": 812.4,            // rank-0 cpu_time - WorldProfile::cpu_start
  "total_wall_s": 141.2,           // rank-0 wall_time - WorldProfile::wall_start
  "overhead_s_per_call": 3.0e-7,   // estimated profiling overhead/call
  "context": {                     // OPTIONAL; filled by the v3 caller, null in
                                   // core. The join key for the cost-model fit.
    "molecule": "h2o", "n_occ": 5, "k": 6, "thresh": 1e-4,
    "protocol": 0, "box_L": 30.0
  },
  "phases": [
    {
      "name": "FunctionImpl::apply",   // raw __FUNCTION__ / class::function key
      "phase": "apply",                // canonical taxonomy (PM-2); else "other"
      "count":      {"sum": 1.2e5, "min": 1.9e4, "max": 2.1e4, "pmin": 3, "pmax": 0},
      "cpu_excl_s": {"sum": 402.1, "min": 61.0, "max": 71.3, "pmin": 4, "pmax": 0},
      "cpu_incl_s": {"sum": 588.0, "min": 90.2, "max": 99.8, "pmin": 4, "pmax": 0},
      "nmsg_sent_excl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nmsg_sent_incl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nbyte_sent_excl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nbyte_sent_incl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nmsg_recv_excl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nmsg_recv_incl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nbyte_recv_excl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nbyte_recv_incl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0}
    }
    // … one object per registered profile entry
  ]
}
```

Canonical `phase` values (PM-2): `apply`, `compress`, `reconstruct`, `multiply`,
`inner`, `gaxpy`, `truncate`, `exchange` (γ build, block `rs_exchange_gamma`),
`projection` (`Q·v`, block `rs_projection`), `other`. The exchange thread's
Tx/tile counters aggregate under the `rs_exchange_gamma` block so they land in the
`exchange` phase. Cost-model consumers (§9 of the companion doc): use `cpu_excl_s`
+ `count` for `T_compute`, `nmsg_*`/`nbyte_*` for `T_comm`, and `max/sum` ratios
for the imbalance factor φ; recover wall by joining `context` with the coarse
`StateMetrics.wall_s` layer.
