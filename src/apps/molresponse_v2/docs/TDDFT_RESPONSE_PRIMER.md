# TDDFT/TDHF Response Properties — A Developer Primer

This document explains the physics and code structure of the `molresponse_v2`
linear-response solver.  It is aimed at developers who want to understand how
the equations of TDDFT/TDHF map to the data structures and algorithm in the
codebase, or who want to extend the solver to a new property type.

---

## 1. Physical Background

### 1.1 Time-dependent density matrix

We work in the one-particle density-matrix (1PDM) formulation of TDDFT/TDHF.
The time-dependent 1PDM obeys the equation of motion:

```
i ∂γ/∂t = [F(t), γ(t)]
```

where `F(t) = h + g[γ(t)] + v(t)` is the time-dependent Fock operator:
- `h`    — one-electron (kinetic + nuclear attraction)
- `g`    — two-electron (Coulomb + exchange-correlation)
- `v(t)` — external perturbation (e.g. electric dipole along x/y/z)

The perturbation has the form:

```
v(t) = Σ_C λ_C [ v^{C†}(ω_C) e^{+iω_C t}  +  v^C(ω_C) e^{-iω_C t} ]
```

where `C` runs over perturbation channels (dipole components, nuclear displacements, …).

### 1.2 Perturbation expansion and first-order response equation

Expanding `γ` in powers of the coupling constants λ_C:

```
γ = γ^0  +  Σ_C λ_C γ^C  +  Σ_{BC} λ_B λ_C γ^{BC}  +  …
```

collecting the first-order terms in the frequency domain gives:

```
ω_C γ^C = [F^0, γ^C]  +  [F^C, γ^0]
```

This is the **coupled first-order response equation**.  Its solution gives `γ^C`,
from which all linear response properties (polarizability, optical rotatory
strength, Raman intensities, …) are assembled.

### 1.3 Orbital decomposition — the ResponseVector

The ground-state density matrix is `γ^0 = Σ_i φ_i φ_i†` (N occupied orbitals).
The first-order density matrix takes the form:

```
γ^C(r,r') = Σ_i [ x_i^C(r) φ_i†(r')  +  φ_i(r) y_i^C(r') ]
```

where:
- `x_i^C(r)` — **forward response orbital** (positive-frequency channel)
- `y_i^C(r)` — **backward response orbital** (negative-frequency channel)

These are stored in the `ResponseVector` data structure defined in
`ResponseVector.hpp`.  The four concrete types correspond to the four
combinations of spin restriction and TDA flag:

| Type | x channels | y channels | `alpha_factor` |
|------|-----------|-----------|----------------|
| `StaticRestrictedResponse`    | x_alpha | —      | −4.0 |
| `DynamicRestrictedResponse`   | x_alpha | y_alpha | −2.0 |
| `StaticUnrestrictedResponse`  | x_alpha, x_beta | — | −2.0 |
| `DynamicUnrestrictedResponse` | x_alpha, x_beta | y_alpha, y_beta | −2.0 |

`StaticRestricted` corresponds to the **Tamm-Dancoff Approximation (TDA)**:
setting `y ≡ 0` decouples the negative-frequency channel, giving

```
γ^C ≈ Σ_i x_i^C(r) φ_i†(r')
```

### 1.4 BSH integral-equation form

In MRA (Multiresolution Analysis) the differential equations are
ill-conditioned.  We solve instead the **Bound-State Helmholtz (BSH)**
integral form.  For each occupied orbital p:

```
x_p^C = -2 Ĝ(k_p^x) * [ V^0 x_p^C  -  Σ_{i≠p} ε_{ip} x_i^C  +  g_p'[γ^C] φ_p  +  V_p^C ]
y_p^C = -2 Ĝ(k_p^y) * [ V^0 y_p^C  -  Σ_{i≠p} ε_{ip} y_i^C  +  g_p'[γ^{C†}] φ_p  +  V_p^{C†} ]
```

where:
- `Ĝ(k)` is the BSH Green's function operator with parameter `k`
- `k_p^x = sqrt(-2(ε_p + ω_C))` — x-channel BSH exponent
- `k_p^y = sqrt(-2(ε_p - ω_C))` — y-channel BSH exponent
- `ε_{ip} = ∫ φ_i†(r) F^0 φ_p(r) dr` — off-diagonal Fock coupling (orbital localization)
- `V_p^C = Q̂ v^C φ_p`, `Q̂ = 1 - γ^0` — perturbation projected onto unoccupied space
- `g_p'[γ^C] φ_p` — response XC kernel (Γ_X in the paper)

Term-to-code mapping in `ResponseStatePotentials` (see `ResponseKernels.hpp`):

| Math term | Code symbol |
|-----------|------------|
| `V^0 x_p` | `v0` |
| `Σ_{i≠p} ε_{ip} x_i` | computed by `apply_hamiltonian_no_diag()` |
| `g_p'[γ^C] φ_p` | `gamma` |
| `v0 - ε + gamma` | `lambda` |
| `-2*(lambda + vp)` | `theta` (BSH source term) |

Convergence is driven by **KAIN** (Krylov-accelerated inexact Newton) until both
the orbital residual norm and the density change fall below protocol thresholds.
The `Q̂` projector is applied after every BSH step via
`project_response_channels()`.

### 1.5 Symplectic metric

For excited-state calculations the response states are rotated in a subspace
that respects the generalized eigenvalue structure.  The metric inner product is:

```
⟨Φ|Φ'⟩_metric = ⟨x_α|x'_α⟩ − ⟨y_α|y'_α⟩
```

The **minus sign on the y block is not arbitrary** — it arises from the
symplectic structure of the Casida (TDDFT) response matrix:

```
⎡ A   B  ⎤ ⎡X⎤       ⎡ 1   0 ⎤ ⎡X⎤
⎣ B*  A* ⎦ ⎣Y⎦  = ω  ⎣ 0  −1 ⎦ ⎣Y⎦
```

For TDA (`StaticRestricted`) `B = 0` and there is no y block; the metric
reduces to the plain overlap `⟨x|x'⟩`.

---

## 2. ResponseVector — Code Mapping

### 2.1 The four concrete types

```cpp
// ResponseVector.hpp

struct StaticRestrictedResponse {
    vector_real_function_3d x_alpha;   // x_i^C(r), i=0..N-1
    vector_real_function_3d flat;      // [ x_alpha[0..N-1] ]
};

struct DynamicRestrictedResponse {
    vector_real_function_3d x_alpha;   // forward  x-channel
    vector_real_function_3d y_alpha;   // backward y-channel
    vector_real_function_3d flat;      // [ x_alpha | y_alpha ]
};

struct StaticUnrestrictedResponse {
    vector_real_function_3d x_alpha;
    vector_real_function_3d x_beta;
    vector_real_function_3d flat;      // [ x_alpha | x_beta ]
};

struct DynamicUnrestrictedResponse {
    vector_real_function_3d x_alpha;
    vector_real_function_3d y_alpha;
    vector_real_function_3d x_beta;
    vector_real_function_3d y_beta;
    vector_real_function_3d flat;      // [ x_alpha | y_alpha | x_beta | y_beta ]
};

using ResponseVector = std::variant<
    StaticRestrictedResponse,
    DynamicRestrictedResponse,
    StaticUnrestrictedResponse,
    DynamicUnrestrictedResponse>;
```

### 2.2 Flat vector layout

Each struct maintains two synchronized views of the same data:
- **Typed channels**: `x_alpha`, `y_alpha`, `x_beta`, `y_beta`
- **Concatenated flat**: `flat`

The flat layout is:

```
StaticRestrictedResponse     (N slots):
  flat = [ x_alpha[0..N-1] ]
          channel 0 @ offset 0

DynamicRestrictedResponse    (2N slots):
  flat = [ x_alpha[0..N-1]  |  y_alpha[0..N-1] ]
          channel 0 @ 0        channel 1 @ N

StaticUnrestrictedResponse   (2N slots):
  flat = [ x_alpha[0..N-1]  |  x_beta[0..N-1] ]
          channel 0 @ 0        channel 1 @ N

DynamicUnrestrictedResponse  (4N slots):
  flat = [ x_alpha[0..N-1]  |  y_alpha[0..N-1]  |  x_beta[0..N-1]  |  y_beta[0..N-1] ]
          channel 0 @ 0        channel 1 @ N        channel 2 @ 2N     channel 3 @ 3N
```

**WARNING:** This layout is fixed by the binary archive format (`ResponseIO.hpp`).
Do NOT reorder channels without writing archive migration code.

### 2.3 Sync/flatten contract

```
After modifying flat directly (KAIN update, BSH apply, archive load):
    call sync()   — propagates flat → typed channels

After modifying typed channels directly (guess init, rotation result):
    call flatten() — propagates typed channels → flat

Use assign_flat_and_sync(r, new_flat) instead of the two-liner
    r.flat = new_flat; r.sync();
This makes the atomicity explicit and avoids the footgun of forgetting sync.
```

### 2.4 Helper methods

Every concrete struct provides:

| Method | Returns |
|--------|---------|
| `num_orbitals()` | N (number of occupied orbitals) |
| `num_channels()` | 1, 2, or 4 (static constexpr) |
| `num_flat_slots()` | N * num_channels() |
| `channel_offset(k)` | Start index of channel k in flat: k * N |
| `alpha_factor()` | Polarizability prefactor (−4.0 or −2.0) |

Free functions for use with the variant:

```cpp
response_alpha_factor(vec)    // runtime dispatch of alpha_factor()
response_num_orbitals(vec)    // N
response_num_channels(vec)    // num_channels()
response_num_flat_slots(vec)  // N * num_channels()
assign_flat_and_sync(r, f)    // atomic flat replacement + sync
get_flat(vec)                 // reference to flat member
```

### 2.5 alpha_factor derivation

```
α_{AB}(ω) = alpha_factor * ⟨x_A | V_B⟩

Physical derivation:
  Polarizability: α_{AB} = −⟨ρ¹_A | r_B⟩

  TDA/Restricted (StaticRestrictedResponse):
    ρ¹ = 2 Σ_i x_i φ_i*   (factor 2 from two spin channels)
    α  = −2·⟨Σ_i x_i φ_i*|r_B⟩·2 = −4·⟨x|v_B⟩   → alpha_factor = -4.0

  Dynamic/Restricted or Unrestricted:
    Both x and y channels contribute; the combined spin and channel factor
    gives a net −2.                                 → alpha_factor = -2.0
```

---

## 3. Minimal Working Example — Static Polarizability of Water

The following pseudocode sketches the key steps for computing the static
polarizability α(ω=0) of water using `molresponse_v2`.

```cpp
// 1. Load ground-state data from a completed moldft calculation.
GroundStateData gs = load_ground_state(world, "water.calc_info.json");

// 2. Choose the response type (TDA, closed-shell) and create a vector.
//    N = number of occupied orbitals; is_static = true (TDA); unrestricted = false.
size_t N = gs.getNumOrbitals();
ResponseVector rvec = make_response_vector(N, /*is_static=*/true, /*unrestricted=*/false);

// 3. Build the perturbation vector V_p for the dipole-x perturbation.
//    V_p = Q̂ (x * φ_p)  for each orbital p.
vector_real_function_3d vp = perturbation_vector(world, gs, dipole_x_state);

// 4. Construct the KAIN solver for the flat vector (size = N for TDA).
response_vector_allocator alloc(world, N);
response_solver kain(alloc);

// 5. Iterate: BSH → KAIN → sync.
for (int iter = 0; iter < max_iter; ++iter) {
    // Compute response potentials from the current state.
    auto pots = compute_response_potentials(world, rvec, gs);

    // Build BSH operators (x-channel only for TDA).
    auto bsh_ops = make_bsh_operators(world, manager, 0.0, gs.getEnergies(), N, logger, rvec);

    // theta = -2 * (lambda + vp)
    auto theta = pots.lambda.flat - 2.0 * vp;
    scale(world, theta, -2.0);

    // Apply BSH and project onto unoccupied space.
    auto new_flat = apply(world, bsh_ops, theta);
    project_response_channels(world, rvec, gs);   // Qhat

    // KAIN update.
    auto residual = sub(world, new_flat, get_flat(rvec));
    auto kain_flat = kain.update(get_flat(rvec), residual);

    // Update state atomically.
    assign_flat_and_sync(rvec, copy(world, kain_flat));

    // Check convergence.
    if (norm(residual) < threshold) break;
}

// 6. Assemble polarizability.
//    α_{xx} = alpha_factor * inner(x_alpha, vp_x)
double af = response_alpha_factor(rvec);  // -4.0 for StaticRestricted
double alpha_xx = af * inner(world, get_flat(rvec), vp).sum();
```

For real runs, steps 3–6 are handled automatically by the pipeline:
`FrequencyLoop.hpp` → `ResponseSolver.cpp` → `PropertyManager.hpp`.

---

## 4. Excited States — TDA Iteration Loop Step-by-Step

For excited states we solve for a **bundle of M states simultaneously** without
an external perturbation.  The excitation energies ω and transition vectors
(x_1, …, x_M) are the eigenvectors of the Casida generalized eigenvalue problem.

The entry point is `iterate_excited<R>(world, bundle, omega, gs, params)` in
`ExcitedResponse.cpp`.  It accepts any restricted closed-shell response type R:
- `TDARestrictedResponse` — x-only, TDA density ρ¹ = Σ xᵢφᵢ
- `DynamicRestrictedResponse` — x+y, full TDDFT density ρ¹ = Σ(xᵢ+yᵢ)φᵢ

Each iteration performs seven steps:

### Step 0 — Gram-Schmidt orthonormalization (`metric_gram_schmidt` in `BundleKernels.hpp`)

```cpp
metric_gram_schmidt(world, states, gs);
```

Called at the **start** of every iteration before computing potentials.  Mirrors
`gram_schmidt(Chi.X, Chi.Y) + normalize(Chi)` from the reference implementation.
Uses the symplectic metric ⟨x|x'⟩ − ⟨y|y'⟩ so S stays well-conditioned.

### Step 1 — Compute per-root potentials (`compute_excited_potentials` in `BundleKernels.hpp`)

For each bundle state i:
```
pots.v0[i]     = V^0 x_i  =  (V_local − c_xc·K₀) x_i
pots.gamma[i]  = Q̂·[response XC kernel]           (type-dispatched via ResponseKernels.hpp)
pots.lambda[i] = v0 − ε·x + gamma                  (no kinetic — T is implicit in BSH)
```

### Step 2 — Subspace rotation: S, A → ω, U (`diagonalize_and_rotate_bundle` in `BundleKernels.hpp`)

```
S_ij = metric_inner(states_i, states_j)           TDA: ⟨x|x'⟩;  Full: ⟨x|x'⟩−⟨y|y'⟩
A_ij = ½(⟨Φ_i | λ_j + T·x_j ⟩  +  ⟨Φ_j | λ_i + T·x_i ⟩)
```

**Important:** `T` (kinetic energy via `apply_kinetic_flat`) is applied to the
**state** `x_i`, not to `lambda_i`, and is only added for building A.  It must
NOT appear in the BSH source term (T is handled implicitly by the BSH Green's
function, which would double-count it).

`diagonalize_bundle` (`ResponseKernels.hpp`) solves `A c = ω S c` and applies:
1. Column swap to keep large coefficients near diagonal (root tracking)
2. Phase fixing so diagonal elements are positive
3. SVD rotation within near-degenerate clusters
4. Sort columns by ascending eigenvalue

After diagonalisation, states, lambda, v0, and gamma are all rotated by U, and
each state is renormalized under the metric.

### Step 3 — Snapshot bundle

```cpp
auto prev_bundle = snapshot_bundle(world, states);      // deep copy for KAIN
auto densities_before = compute_response_densities(...); // ρ¹ before BSH step
```

### Step 4 — Per-root BSH step (`ExcitedStateStep` in `BundleKernels.hpp`)

For each bundle state i:
```
θ_p = v0_p − ε_{ip}·x_i + γ_p               (= lambda_p; recomputed from rotated pots)
x_p^new = Ĝ(μ_p) · [−2·(θ_p + shift_p·x_p)]

μ_p = √(−2·(ε_p + ω_i + shift_p))
shift_p applied when ε_p + ω_i ≥ 0 to keep denominator negative (BSH stability)
```

TDA: N BSH operators at +ω.  Full: 2N operators (N at +ω, N at −ω).

### Step 5 — Bundle KAIN acceleration

```cpp
// Direct BSH step for the first two iterations; KAIN from iter ≥ 2.
auto kain_bundle = (iter >= 2)
    ? solver.update(prev_bundle, updated_bundle - prev_bundle)
    : std::move(updated_bundle);
```

The `bundle_solver<R>` is a `XNonlinearSolver` over `ResponseBundle<R>`, so KAIN
captures inter-state coupling from the subspace rotation step.

### Step 6 — Step restriction + Q-project + normalize (`apply_step_project_normalize` in `BundleKernels.hpp`)

For each root i:
1. If `||candidate[i] − prev[i]|| > max_rotation`, mix back toward prev
2. `project_response_channels(world, states[i], gs)` — Q̂ onto virtual subspace
3. `normalize_response_metric(world, states[i])` — metric norm = 1

### Step 7 — Convergence check

```cpp
const double max_res  = max over i of ||updated[i] − prev[i]||
const double max_drho = max over i of ||ρ¹_new − ρ¹_old||_2
converged = (max_res < params.dconv)
```

`ExcitedIterDiagnostics` carries both the per-root final residuals
(`residual_norms`) and the per-iteration max histories
(`iteration_max_residuals`, `iteration_max_density_changes`).

---

## 5. JSON Output — `excited_states.calc_info.json`

After convergence, `write_excited_state_calc_info()` (declared in
`ExcitedResponse.hpp`, implemented in `ExcitedResponse.cpp`) writes a JSON file
that follows the schema in `results/h2_tda_test/EXCITED_STATE_JSON_SCHEMA.md`:

```json
{
  "schema_version": "1.0",
  "metadata":    { "code": "molresponse_v2", "mpi_size": 4 },
  "ground_state":{ "model": "hf", "total_energy_au": -1.1336, "eigenvalues_au": [...],
                   "num_orbitals": 1, "wavelet_order_k": 8, "thresh": 1e-4, "basis": null },
  "calculation": { "type": "excited_state", "approximation": "tda",
                   "num_states_requested": 4, "dconv": 1e-4, ... },
  "convergence": { "converged": true, "converged_for_dconv": 1e-4,
                   "converged_for_thresh": 1e-4, "num_iterations": 12 },
  "excitations": [
    { "root": 1, "spin": "singlet", "irrep": null,
      "omega_au": 0.46789, "omega_ev": 12.73, "current_error": 8.3e-5,
      "oscillator_strength_length": null, "oscillator_strength_velocity": null,
      "converged": true },
    ...
  ]
}
```

`irrep` and `oscillator_strength_*` are `null` because MRA runs without symmetry
constraints and the current TDA implementation does not compute transition
dipoles.  Comparing against Dalton reference data should use only `omega_au`.

Both `test_tda_h2` and `ExcitedResponse::solve_protocol` write this file to the
ground-state directory (where `moldft.restartdata` lives).

---

## 6. Extending to Open-Shell (Unrestricted)

Unrestricted support (`StaticUnrestrictedResponse`, `DynamicUnrestrictedResponse`)
requires independent alpha and beta spin channels.  The data layout already
exists in `ResponseVector.hpp`; what is missing is the physics implementation
in the kernels.

Locations where `static_assert` guards block unrestricted use:

| Function | File |
|---------|------|
| `response_plain_inner` | `ResponseKernels.hpp` |
| `response_metric_inner` | `ResponseKernels.hpp` |
| `normalize_response_metric` | `ResponseKernels.hpp` |
| `project_response_channels` (free function) | `ResponseKernels.hpp` |
| `compute_response_density` | `ResponseKernels.hpp` |
| `make_excited_bsh_operators` | `ResponseKernels.hpp` |
| `apply_hamiltonian_no_diag` | `ResponseKernels.hpp` |
| `rotate_bundle` | `ResponseKernels.hpp` |

To implement unrestricted support:
1. Remove the `static_assert` for a function.
2. Add an `if constexpr (response_is_unrestricted_v<ResponseType>)` branch that
   handles alpha and beta channels independently.
3. For the metric inner product: `⟨Φ|Φ'⟩ = ⟨x_α|x'_α⟩ + ⟨x_β|x'_β⟩ − ⟨y_α|y'_α⟩ − ⟨y_β|y'_β⟩`
4. For `apply_hamiltonian_no_diag`: the Hamiltonian coupling between alpha and
   beta spins is not yet implemented — this requires the full UHF Fock matrix.
5. For `CoupledResponseEquations` in `ResponseSolver.cpp`: separate BSH operators
   per spin block.

Unrestricted restart is guess-only until Phase 6 of the excited-state
reintegration plan (see `EXCITED_STATE_REINTEGRATION_PLAN.md`).

---

## 7. Key File Index

| File | Role |
|------|------|
| `ResponseVector.hpp` | The four concrete response types + variant + free functions |
| `ResponseKernels.hpp` | Physics kernels: metric, density, potentials, BSH ops, rotation |
| `ResponseSolver.cpp/hpp` | Per-state iterative solver (`CoupledResponseEquations`, `make_bsh_operators`) |
| `FrequencyLoop.cpp/hpp` | Outer KAIN loop for frequency-dependent states |
| `ExcitedResponse.cpp/hpp` | TDA/TDDFT excited-state iteration (`iterate_excited`) |
| `ExcitedStateBundleSolver.cpp/hpp` | Full bundle solver with protocol loop, restart, metadata |
| `PropertyManager.hpp` | Property assembly (α, β, Raman) from solved states |
| `ResponseIO.hpp` | Archive I/O for response vectors |
| `GroundStateData.hpp` | Container for ground-state quantities (orbitals, energies, Qhat, V_local) |

---

## See Also

- `MOLRESPONSE_TUTORIAL.md` — pipeline walkthrough (input → solve → output)
- `EXCITED_STATE_REINTEGRATION_PLAN.md` — phase roadmap for excited-state support
- `STATE_PARALLEL_DESIGN.md` — subgroup scheduling design
- `docs/reintegration/` — detailed phase contracts and validation notes
