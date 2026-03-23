# Response Solver Refactor Plan: Iteration Engine + Per-Type Ops Modules

**Status:** Approved design note — pending incremental implementation
**Scope:** `src/apps/molresponse_v2/` — response solver subsystem only

---

## 1. Motivation

The `molresponse_v2` response solver already has a shared nonlinear iteration
driver (`iterate<R>(...)` in `FrequencyLoop.hpp`) and per-type one-step updates
(`CoupledResponseEquations()` in `ResponseSolver.cpp`), but key mathematical
variation points are scattered across two templated helpers
(`compute_ground_exchange<R>` and `compute_gamma_response<R>` in
`ResponseVectorKernels.hpp`) that use `if constexpr` branches to handle all
response types in a single template body.

This means that to understand the complete mathematical definition of, say,
the static restricted response model, a developer must simultaneously track:

- `CoupledResponseEquations(..., StaticRestrictedResponse&, ...)` in `ResponseSolver.cpp`
- The `!response_has_y_channel_v<R> && !is_same_v<R, TDARestrictedResponse>` branch inside `compute_gamma_response<R>` in `ResponseVectorKernels.hpp`
- The `!response_has_y_channel_v<R>` branch inside `compute_ground_exchange<R>` in `ResponseVectorKernels.hpp`

The goal of this refactor is: **open one file, see all the math for one response model.**

---

## 2. Current Architecture

### 2.1 Call stack for one response point

```
MolresponseLib::run_response()
  └─ execute_serial_state_solve() / execute_subgroup_state_solve()
      └─ computeFrequencyLoop()   [FrequencyLoop.cpp:566]
          └─ solve_and_record_point()   [FrequencyLoop.cpp:373]
              └─ run_solver_attempt()   [FrequencyLoop.cpp:312]
                  └─ solve_response_vector()   [FrequencyLoop.cpp:485]
                      │  (std::visit dispatch)
                      └─ iterate<ResponseType>(...)   [FrequencyLoop.hpp:56]
                          ├─ make_bsh_operators(...)   [ResponseSolver.cpp]
                          └─ for each KAIN iteration:
                              ├─ CoupledResponseEquations(...)   [ResponseSolver.cpp]
                              │   ├─ compute_ground_exchange<R>()   [ResponseVectorKernels.hpp]
                              │   └─ compute_gamma_response<R>()   [ResponseVectorKernels.hpp]
                              ├─ solver.update(x, residuals)   [KAIN]
                              ├─ assign_flat_and_sync(rvec, kain_x)
                              ├─ compute_density(world, rvec, phi0)
                              └─ [convergence checks + logging]
```

### 2.2 Key files

| File | Role |
|------|------|
| `FrequencyLoop.hpp` | `iterate<R>(...)` — shared iteration driver (lines 56–295) |
| `FrequencyLoop.cpp` | `solve_response_vector()` variant dispatch; `computeFrequencyLoop()` |
| `ResponseSolver.hpp` | Free-function declarations: `alpha_factor`, `compute_density`, `make_bsh_operators`, `CoupledResponseEquations` |
| `ResponseSolver.cpp` | Implementations of above (one overload per response type) |
| `ResponseVectorKernels.hpp` | `compute_ground_exchange<R>` and `compute_gamma_response<R>` with `if constexpr` branches |
| `ResponseSolverUtils.hpp` | Shared utilities: BSH factory helper, KAIN allocator, norm helpers |
| `ResponseVector.hpp` | Type definitions + `assign_flat_and_sync()` |

### 2.3 Response types

| Type | Status | Flat channels |
|------|--------|---------------|
| `StaticRestrictedResponse` | Implemented | x (N) |
| `TDARestrictedResponse` | Implemented (phase 4 WIP) | x (N) |
| `DynamicRestrictedResponse` | Implemented | x+y (2N) |
| `StaticUnrestrictedResponse` | Not implemented (throws) | — |
| `DynamicUnrestrictedResponse` | Not implemented (throws) | — |

### 2.4 Variation points per response type

The following are already free-function overloads (compile-time dispatch by type):

| Function | Current file |
|----------|-------------|
| `alpha_factor(R&)` | `ResponseSolver.hpp` |
| `compute_density(world, R&, phi0)` | `ResponseSolver.hpp/cpp` |
| `make_bsh_operators(world, mgr, freq, eps, n, lgr, R&)` | `ResponseSolver.hpp/cpp` |
| `CoupledResponseEquations(world, gs, R&, vp, bsh, mgr, lgr)` | `ResponseSolver.hpp/cpp` |

The following are **templated helpers with `if constexpr` branches** — these are the problem:

| Function | Current file | Branches |
|----------|-------------|----------|
| `compute_ground_exchange<R>(world, R&, phi0)` | `ResponseVectorKernels.hpp` | `response_has_y_channel_v<R>`, unrestricted guard |
| `compute_gamma_response<R>(world, R&, phi0, Q, c_xc)` | `ResponseVectorKernels.hpp` | `response_has_y_channel_v<R>`, `is_same_v<R, TDARestrictedResponse>`, unrestricted guard |

### 2.5 What is genuinely shared (iteration engine)

The `iterate<R>()` template owns:
- KAIN nonlinear solver loop (`XNonlinearSolver`)
- Convergence checks (residual norm + density change vs. thresholds)
- Safety stops (non-finite detection, explosive growth, wall-time stall timeout)
- Iteration table logging (`ResponseDebugLogger`)
- `assign_flat_and_sync()` + density recompute
- `perturbation_vector()` call
- `ResponseSolveDiagnostics` accumulation and return

These **stay in `FrequencyLoop.hpp`**. No changes planned to the iteration engine.

### 2.6 Mathematical variation between types

| Aspect | Static Restricted | TDA Restricted | Dynamic Restricted |
|--------|------------------|----------------|-------------------|
| Channels | x only (N) | x only (N) | x + y (2N) |
| Response density | 2 Σ xᵢφᵢ | Σ xᵢφᵢ (no ×2) | Σ (xᵢ+yᵢ)φᵢ |
| Coulomb | 2J[ρ]φ | 2J[ρ]φ | 2J[ρ]φ (shared) |
| Exchange (x-ch) | K[φ₀,x] + K[x,φ₀] | K[φ₀,x] only | K[φ₀,x] + K[y,φ₀] |
| Exchange (y-ch) | — | — | K[φ₀,y] + K[x,φ₀] |
| BSH ops | N ops, μ=√(−2ε) | N ops, μ=√(−2(ε+ω)) | 2N ops: μˣ=√(−2(ε+ω)), μʸ=√(−2(ε−ω)) |
| alpha_factor | −4.0 | −4.0 | −2.0 |

The TDA and Static kernels share the same residual equation (`θ = −2(V₀x − εx + γx + vp)`);
they differ only in the density factor used inside `compute_gamma_response`.

---

## 3. Design Comparison: Templates vs. Virtual OO

### Template/overload approach (recommended)

The code already dispatches via free-function overloads. The "ops module"
concept becomes: one header per response type that groups all its overloads.

**Advantages:**
- Zero runtime overhead — no vtable, no heap allocation
- Matches the existing code style exactly
- Incremental migration: move code between files, no signature changes
- Compile-time type checking
- Future unrestricted types drop in by adding new overloads

**Disadvantages:**
- Full compilation of iterate<R>() for each R (already the case)
- Slightly harder to mock for unit tests (free functions vs. virtual methods)

### Virtual OO approach (not recommended)

Would require wrapping existing free functions in a class hierarchy, adding
heap-allocated ops objects, and threading them through `iterate()`. The
response type is always known at compile time, so runtime polymorphism buys
nothing here. This would increase complexity for no gain.

**Decision: keep and formalize the overload pattern.**

---

## 4. Target Design

### 4.1 File layout

```
src/apps/molresponse_v2/
├── FrequencyLoop.hpp/cpp          ← iteration engine (UNCHANGED)
├── ResponseSolverUtils.hpp        ← shared low-level helpers (UNCHANGED)
├── ResponseVector.hpp             ← type definitions (UNCHANGED)
├── ResponseVectorKernels.hpp      ← DEPRECATED after migration
│
├── ResponseSolver.hpp             ← facade: includes all ops headers
├── ResponseSolver.cpp             ← shared impls only (empty after migration)
│
└── ops/                           ← NEW
    ├── StaticRestrictedOps.hpp    ← static, closed-shell (ω = 0)
    ├── DynamicRestrictedOps.hpp   ← full TDDFT, closed-shell (ω ≠ 0)
    ├── TDARestrictedOps.hpp       ← TDA excited-state, closed-shell
    ├── StaticUnrestrictedOps.hpp  ← placeholder (stubs that throw)
    └── DynamicUnrestrictedOps.hpp ← placeholder (stubs that throw)
```

### 4.2 Ops interface (implicit contract per type R)

Each `ops/*.hpp` file defines these free functions for its concrete response type:

```cpp
// 1. Polarizability prefactor
constexpr double alpha_factor(const R& r);

// 2. Response density ρ¹ from current iterate
real_function_3d compute_density(World& world, const R& rvec,
                                  const vector_real_function_3d& phi0);

// 3. BSH operator set (N or 2N operators depending on R)
std::vector<poperatorT> make_bsh_operators(
    World& world, const ResponseManager& rm, double freq,
    const Tensor<double>& orbital_energies, int n,
    ResponseDebugLogger& lgr, const R& /* tag */);

// 4. One-step BSH iteration: θ → bsh·θ → Q·rsh
//    All response-type math lives here (exchange, gamma, residual, Q-proj)
vector_real_function_3d CoupledResponseEquations(
    World& world, const GroundStateData& g_s, const R& vecs,
    const vector_real_function_3d& v_p,
    const std::vector<poperatorT>& bsh,
    const ResponseManager& rm, ResponseDebugLogger& lgr);
```

No base class. No virtual functions. `iterate<R>()` calls these via ADL/overload
resolution, exactly as today.

### 4.3 Each ops file is self-contained

Each `ops/*.hpp`:
- `#include`s only: `../ResponseVector.hpp`, `../GroundStateData.hpp`, `../ResponseSolverUtils.hpp`
- Has **no `if constexpr` branches** (the type is fixed at the file level)
- Has **inline physics comments** matching the equations for that model
- Can be read in isolation to understand one complete response model

### 4.4 ResponseSolver.hpp becomes a facade

After migration:

```cpp
// ResponseSolver.hpp — facade
#pragma once
#include "ops/StaticRestrictedOps.hpp"
#include "ops/DynamicRestrictedOps.hpp"
#include "ops/TDARestrictedOps.hpp"
#include "ops/StaticUnrestrictedOps.hpp"
#include "ops/DynamicUnrestrictedOps.hpp"
// shared helpers still declared here
```

This preserves all existing `#include "ResponseSolver.hpp"` call sites.

---

## 5. Concern / Location / Target Table

| Concern | Current location | Target location | Shared or per-type |
|---------|----------------|-----------------|-------------------|
| `alpha_factor` | `ResponseSolver.hpp` | `ops/*.hpp` | per-type |
| `compute_density` | `ResponseSolver.hpp/cpp` | `ops/*.hpp` | per-type |
| `make_bsh_operators` | `ResponseSolver.hpp/cpp` | `ops/*.hpp` | per-type |
| `CoupledResponseEquations` | `ResponseSolver.hpp/cpp` | `ops/*.hpp` | per-type |
| `compute_ground_exchange` math | `ResponseVectorKernels.hpp` (if constexpr) | inlined in `CoupledResponseEquations` | per-type |
| `compute_gamma_response` math | `ResponseVectorKernels.hpp` (if constexpr) | inlined in `CoupledResponseEquations` | per-type |
| BSH factory helper | `ResponseSolverUtils.hpp` | `ResponseSolverUtils.hpp` (unchanged) | shared |
| KAIN loop | `FrequencyLoop.hpp:iterate()` | `FrequencyLoop.hpp` (unchanged) | shared |
| Convergence checks | `FrequencyLoop.hpp:iterate()` | `FrequencyLoop.hpp` (unchanged) | shared |
| Safety stops | `FrequencyLoop.hpp:iterate()` | `FrequencyLoop.hpp` (unchanged) | shared |
| Iteration logging | `FrequencyLoop.hpp:iterate()` | `FrequencyLoop.hpp` (unchanged) | shared |
| `assign_flat_and_sync` | `ResponseVector.hpp` | `ResponseVector.hpp` (unchanged) | shared |
| Restart seeding | `FrequencyLoop.cpp` | `FrequencyLoop.cpp` (unchanged) | shared |
| Variant dispatch | `FrequencyLoop.cpp:solve_response_vector()` | `FrequencyLoop.cpp` (unchanged) | shared |

---

## 6. Migration Plan

### Step 0 — Create directory structure (no behavior change)

- Create `src/apps/molresponse_v2/ops/`
- Add five empty placeholder `.hpp` files (one per type)
- Update `ResponseSolver.hpp` to `#include` them at the bottom (empty files → no-op)
- **Build must pass. No algorithmic change.**

### Step 1 — Migrate `StaticRestrictedResponse` (first slice)

1. Write `ops/StaticRestrictedOps.hpp` with:
   - `alpha_factor(StaticRestrictedResponse&)` — move from `ResponseSolver.hpp`
   - `compute_density(world, StaticRestrictedResponse&, phi0)` — move from `ResponseSolver.cpp`
   - `make_bsh_operators(..., StaticRestrictedResponse&)` — move from `ResponseSolver.cpp`
   - `CoupledResponseEquations(..., StaticRestrictedResponse&, ...)` — move from `ResponseSolver.cpp`
     - Inline the static branch of `compute_gamma_response<R>` (ρ = 2Σ xᵢφᵢ, both K terms)
     - Inline `compute_ground_exchange<R>` x-only logic
     - **No `if constexpr`. Type is fixed.**
2. Remove those declarations/definitions from `ResponseSolver.hpp/cpp`
3. Add `#include "ops/StaticRestrictedOps.hpp"` to `ResponseSolver.hpp`
4. Build and verify: same output from a known-good static calc

### Step 2 — Migrate `DynamicRestrictedResponse`

Same pattern as Step 1. Inline TDDFT (x+y) density and exchange.
Note: the shared Coulomb optimization (`J[ρ]` computed once for both x and y) lives here.

Verify with a known-good dynamic calc.

### Step 3 — Migrate `TDARestrictedResponse`

Same pattern. Inline TDA branch (ρ = Σ xᵢφᵢ, no ×2, single K term).

Note: TDA and Static kernels are mathematically nearly identical; modest duplication
(~5 lines) is acceptable here to avoid a shared helper that would re-introduce
cross-type coupling.

### Step 4 — Add unrestricted stubs

Create `ops/StaticUnrestrictedOps.hpp` and `ops/DynamicUnrestrictedOps.hpp`.
Define the 4 functions as `throw std::runtime_error("not yet implemented")`.
Remove the corresponding throw-stubs from `ResponseSolver.cpp`.

### Step 5 — Deprecate `ResponseVectorKernels.hpp`

After Steps 1–4, `compute_ground_exchange<R>` and `compute_gamma_response<R>`
have no callers. Mark the file as deprecated in a comment; remove in a subsequent
cleanup commit once confirmed zero remaining uses.

### Step 6 (optional) — Subdirectory reorganization

If per-type subdirectories are preferred over a flat `ops/` directory, this
can be done as a pure rename/move after the behavior is confirmed stable.

---

## 7. Key Risks

| Risk | Mitigation |
|------|-----------|
| Density factor difference TDA/Static missed during inlining | Explicitly note which branch is being inlined; review against `compute_gamma_response` before extraction |
| `if constexpr` semantics lost when inlined (silent wrong branch) | Inline one type at a time; compile with `-Wall`; run known-good numerical test |
| `flat`/`sync` contract broken by reordering calls | Keep `assign_flat_and_sync()` and `sync()` calls in the same positions; do not reorder |
| Q-projection placement changes | Keep `Qhat(rsh)` at the end of `CoupledResponseEquations`, same as today |
| BSH exponent sign for y-channel (`bsh_y` uses `−freq`) | Test Dynamic before and after migration with a known numerical result |
| `ExcitedResponse.cpp` also calls `CoupledResponseEquations` | Update its includes after declarations move; verify it compiles |
| `make_bsh_operators` for unrestricted types declared in `ResponseSolver.hpp` but not implemented | Stubs go in the unrestricted ops files; remove throw-only bodies from `ResponseSolver.cpp` |

---

## 8. Implementation Sketch

### Shared iteration engine (unchanged)

```cpp
// FrequencyLoop.hpp — no signature change
template <typename ResponseType>
ResponseSolveDiagnostics
iterate(World& world, const ResponseManager& response_manager,
        const GroundStateData& g_s, const LinearResponseDescriptor& state,
        const LinearResponsePoint& pt, ResponseType& response,
        ResponseDebugLogger& logger, size_t max_iter, double conv_thresh);
```

### Example: `ops/StaticRestrictedOps.hpp`

```cpp
// ops/StaticRestrictedOps.hpp
//
// Self-contained ops module for static (ω=0) restricted (closed-shell) response.
//
// Equation solved:  A·x = v_p
// where:
//   A·x_p = V₀·x_p − ε_p·x_p + Γˢᵗᵃᵗⁱᶜ[x]_p
//   Γˢᵗᵃᵗⁱᶜ[x] = 2·J[ρˢᵗᵃᵗⁱᶜ]·φ − c_xc·(K[φ₀,x]·φ + K[x,φ₀]·φ)
//   ρˢᵗᵃᵗⁱᶜ = 2·Σᵢ xᵢ·φᵢ           (y ≡ x in ω→0 limit, factor 2 from spin)
//   BSH: μ_p = √(−2·ε_p)             (ω = 0, no shift)
//
#pragma once
#include "../ResponseVector.hpp"
#include "../GroundStateData.hpp"
#include "../ResponseSolverUtils.hpp"

// ---- alpha_factor -----------------------------------------------------------
// Prefactor in α_{AB} = alpha_factor * ⟨x_A | V_B⟩
// Static/restricted: −4 = −2(spin) × 2(y≡x, B-matrix dropped)
constexpr double alpha_factor(const StaticRestrictedResponse&) { return -4.0; }

// ---- compute_density --------------------------------------------------------
// ρˢᵗᵃᵗⁱᶜ = 2·Σᵢ xᵢ·φᵢ  (y ≡ x in ω→0 limit → factor 2)
inline real_function_3d
compute_density(World& world, const StaticRestrictedResponse& rvec,
                const vector_real_function_3d& phi0) {
    auto xphi = mul(world, rvec.x_alpha, phi0, true);
    return 2.0 * sum(world, xphi, true);
}

// ---- make_bsh_operators -----------------------------------------------------
// N operators, one per occupied orbital, exponent μ_p = √(−2·ε_p)
std::vector<poperatorT>
make_bsh_operators(World& world, const ResponseManager& rm, double freq,
                   const Tensor<double>& orbital_energies, int n,
                   ResponseDebugLogger& lgr, const StaticRestrictedResponse&);

// ---- CoupledResponseEquations -----------------------------------------------
// One-step BSH iteration:
//   θ_p = −2·(V₀·x_p − ε_p·x_p + Γˢᵗᵃᵗⁱᶜ[x]_p + v_p^C)
//   x_new_p = bsh_x[p]·θ_p
//   result = Q̂·x_new   (project out occupied space)
vector_real_function_3d
CoupledResponseEquations(World& world, const GroundStateData& g_s,
                          const StaticRestrictedResponse& vecs,
                          const vector_real_function_3d& v_p,
                          const std::vector<poperatorT>& bsh_x,
                          const ResponseManager& rm, ResponseDebugLogger& lgr);
```

### Before/after pseudo-diff for `CoupledResponseEquations`

**BEFORE (ResponseSolver.cpp, lines 151–175):**
```cpp
vector_real_function_3d CoupledResponseEquations(World& world, const GroundStateData& g_s,
    const StaticRestrictedResponse& vecs, ...) {
    auto c_xc = g_s.xcf_.hf_exchange_coefficient();
    StaticRestrictedResponse k_0(vecs.num_orbitals());
    vector_real_function_3d g_x;
    DEBUG_TIMED_BLOCK(world, &logger, "g0_task",
        { k_0 = compute_ground_exchange(world, vecs, g_s.orbitals); });  // ← goes into RVK.hpp
    DEBUG_TIMED_BLOCK(world, &logger, "gx_task",
        { g_x = compute_gamma_response(world, vecs, g_s.orbitals, g_s.Qhat, c_xc); }); // ← goes into RVK.hpp
    // ... rest unchanged
}
```

**AFTER (ops/StaticRestrictedOps.cpp — no template, no if constexpr):**
```cpp
vector_real_function_3d CoupledResponseEquations(World& world, const GroundStateData& g_s,
    const StaticRestrictedResponse& vecs, ...) {

    const auto c_xc = g_s.xcf_.hf_exchange_coefficient();
    const auto& x   = vecs.x_alpha;
    const auto& phi0 = g_s.orbitals;
    const double lo  = 1.e-10;
    const double thresh = FunctionDefaults<3>::get_thresh();

    // Ground exchange:  K₀[x] = K[φ₀,φ₀](x)
    auto k0x = K(world, phi0, phi0)(x);   // same as compute_ground_exchange, static branch

    // Response coupling:  Γˢᵗᵃᵗⁱᶜ[x] = Q̂·( 2·J[ρ]·φ − c_xc·(K[φ₀,x]·φ + K[x,φ₀]·φ) )
    //   ρ = 2·Σᵢ xᵢ·φᵢ   (static: y≡x, factor 2)
    auto rho   = 2.0 * sum(world, mul(world, x, phi0, true), true);
    auto J_rho = apply(CoulombOperator(world, lo, thresh), rho);
    auto g_x   = g_s.Qhat(2.0 * (J_rho * phi0)
                  - c_xc * (K(world, phi0, x)(phi0) + K(world, x, phi0)(phi0)));

    auto v_local  = g_s.V_local * x;
    auto v0x      = v_local - c_xc * k0x;
    auto epsilonx = transform(world, x, g_s.Hamiltonian_no_diag, true);
    auto thetax   = -2.0 * (v0x - epsilonx + g_x + v_p);
    truncate(world, thetax);
    auto rsh = apply(world, bsh_x, thetax);
    return g_s.Qhat(rsh);
}
```

The key change: the two calls to `compute_ground_exchange` and `compute_gamma_response`
(which each contain `if constexpr` branches selecting the right math for this type)
are replaced by the literal math for this type, written out directly and commented.

---

## 9. Dependency / Control Flow Diagram

```
FrequencyLoop.hpp:iterate<R>()
│
├─ [SHARED — no change]
│   ├─ perturbation_vector()
│   ├─ make_bsh_operators()      ──────────────► ops/R_Ops.hpp
│   ├─ KAIN loop
│   │   ├─ CoupledResponseEquations()  ─────────► ops/R_Ops.hpp
│   │   │   ├─ ground exchange (inlined)
│   │   │   ├─ gamma/XC coupling (inlined)
│   │   │   ├─ V₀·x − ε·x
│   │   │   ├─ BSH apply
│   │   │   └─ Q̂ projection
│   │   ├─ KAIN update (XNonlinearSolver)
│   │   ├─ assign_flat_and_sync()
│   │   └─ compute_density()     ──────────────► ops/R_Ops.hpp
│   ├─ convergence checks
│   └─ safety stops / logging
│
└─ [SHARED — ResponseSolverUtils.hpp]
    └─ make_bsh_operators_response()   ◄── called from ops/R_Ops.hpp
```

---

## 10. First Refactor Slice Recommendation

Start with **`StaticRestrictedResponse` → `DynamicRestrictedResponse`** in that order:

1. Both are fully implemented and numerically validated
2. `Static` is simpler (x-only, 1 channel) — easiest to extract cleanly
3. `Dynamic` adds the x+y case and confirms the inlining pattern works for 2-channel math
4. Together they cover all currently active CI paths
5. After these two, `TDA` is easy (near-identical to Static), and unrestricted is stubs

After each extraction: `ninja molresponse2` must pass, and the Python smoke test
`test_molresponse_excited_metadata_smoke.py` should continue to pass.
