# GroundState Design — SCF Wrapper for Response

## Motivation

The v2 `GroundStateData` class duplicates almost everything from the core
`SCF` class: re-implements archive loading (same binary format), stores
local copies of orbitals, energies, occupations, molecule, k, L, xc.
All of this already lives in `SCF` with a well-tested implementation.

The v3 `GroundState` wraps `shared_ptr<SCF>` directly, eliminating all
duplication. It extends SCF with only the response-specific operators
that SCF doesn't provide.

## Architecture

```
GroundState
├── shared_ptr<SCF> scf_         ← owns or shares the SCF object
│   ├── amo, bmo                 ← orbitals (alpha, beta)
│   ├── aeps, beps               ← orbital energies
│   ├── aocc, bocc               ← occupations
│   ├── molecule                 ← geometry
│   ├── param                    ← CalculationParameters
│   ├── xc                       ← XCfunctional
│   └── potentialmanager         ← nuclear potentials
│
├── v_local_                     ← V_nuc + 2*V_coul + V_xc (cached)
├── hamiltonian_                 ← Fock matrix (cached)
├── hamiltonian_no_diag_         ← Fock with diagonal zeroed (cached)
├── q_projector_                 ← QProjector from orbitals (cached)
└── prepared_ flag               ← guards access to cached operators
```

**Key principle:** no local copies of anything SCF already stores. All
metadata accessors are one-line delegations to `scf_->` members.

## Construction Paths

### From archive (common case)

```cpp
auto gs = GroundState::from_archive(world, "moldft.restartdata", molecule);
```

The factory:
1. Reads archive header (version, L, k, xc, nmo_alpha, spin_restricted)
2. Builds `CalculationParameters` with L, xc, prefix, and nopen
   (inferred from `nmo_alpha` and nuclear charge)
3. Constructs `SCF(world, params, molecule)`
4. Calls `SCF::load_mos(world)` — handles all archive versioning,
   orbital loading, and k/thresh adjustment

**Why two-pass:** SCF::load_mos checks `L == param.L()`, so we must
read L from the archive before constructing SCF.

### From live SCF (madqc workflow)

```cpp
auto gs = GroundState(world, scf_ptr);
```

Direct shared ownership. Used when SCF was already run by madqc and
is available as a `shared_ptr<SCF>`.

## What SCF Provides (delegated, not copied)

| Data | SCF member | GroundState accessor |
|------|-----------|---------------------|
| Alpha orbitals | `scf->amo` / `scf->get_amo()` | `orbitals()`, `orbitals_alpha()` |
| Beta orbitals | `scf->bmo` / `scf->get_bmo()` | `orbitals_beta()` |
| Alpha energies | `scf->aeps` | `energies()`, `energies_alpha()` |
| Beta energies | `scf->beps` | `energies_beta()` |
| Alpha occupations | `scf->get_aocc()` | `occupations_alpha()` |
| Beta occupations | `scf->get_bocc()` | `occupations_beta()` |
| Molecule | `scf->molecule` | `molecule()` |
| Parameters | `scf->param` | `params()` |
| Spin restriction | `scf->is_spin_restricted()` | `is_spin_restricted()` |
| Box size | `scf->param.L()` | `L()` |
| HF exchange coeff | `scf->xc.hf_exchange_coefficient()` | `hf_exchange_coefficient()` |

## What GroundState Adds (response-specific)

These are quantities SCF computes during `solve()` but doesn't store
persistently. The response solver needs them pre-computed:

| Operator | Description | Built by |
|----------|-------------|----------|
| `V_local()` | V_nuc + 2*V_coul + V_xc | `prepare()` |
| `hamiltonian()` | Full Fock matrix | `prepare()` |
| `hamiltonian_no_diag()` | Off-diagonal Fock | `prepare()` |
| `Q()` | Virtual-space projector | `prepare()` |

## Protocol Stepping

The response solver iterates over protocol thresholds (e.g. 1e-4, 1e-6).
At each step:

1. **ResponseManager** sets `FunctionDefaults<3>` (k, thresh, truncate mode)
2. **GroundState::prepare()** is called:
   - Re-projects orbitals to current k/thresh (reloads from archive if k changed)
   - Rebuilds QProjector
   - Recomputes V_local (density → Coulomb → XC → sum)
   - Recomputes or loads Hamiltonian (checks fock.json first)
   - Builds off-diagonal Hamiltonian
   - Sets `prepared_ = true`

`prepare()` does NOT call `SCF::set_protocol()` — the response orchestrator
handles FunctionDefaults. This avoids creating redundant operators.

## Open-Shell Support

SCF natively stores both alpha and beta orbital sets. GroundState exposes
both through separate accessors:

- `orbitals_alpha()` / `orbitals_beta()`
- `energies_alpha()` / `energies_beta()`
- `occupations_alpha()` / `occupations_beta()`
- `num_alpha()` / `num_beta()`

The `orbitals()` convenience alias returns alpha orbitals (backward
compatible with v2 restricted-only code).

For `from_archive`: `nopen` is inferred from `nmo_alpha` in the archive
header: `nopen = 2 * nmo_alpha - total_nuclear_charge`.

## Response Solver Usage Pattern

After `prepare()`, the solver accesses:

```cpp
// Per-iteration kernel (from v2 StaticRestrictedOps.hpp)
auto v_local_x = gs.V_local() * x;                          // multiplicative
auto eps_x = transform(world, x, gs.hamiltonian_no_diag());  // Fock coupling
auto result = gs.Q()(exchange_term);                          // project to virtual

// Exchange computation
double c_xc = gs.hf_exchange_coefficient();
auto kx = K(world, gs.orbitals(), gs.orbitals())(x);
```

All access is read-only after `prepare()` completes.

## Files

| File | Role |
|------|------|
| `src/apps/molresponse_v3/GroundState.hpp` | Class declaration |
| `src/apps/molresponse_v3/GroundState.cpp` | Implementation |
| `src/madness/chem/SCF.h` | Wrapped class |
| `src/madness/chem/SCF.cc` | `load_mos()`, `make_nuclear_potential()`, etc. |

## What Was Eliminated

| v2 GroundStateData | v3 GroundState |
|--------------------|----------------|
| 85 lines of archive loading code | Delegates to `SCF::load_mos()` |
| Local copies of molecule, orbitals, energies, occ, k, L, xc | Accessors to `scf_->` members |
| Separate `XCfunctional xcf_` | Uses `scf_->xc` |
| Separate `PotentialManager` | Uses `scf_->potentialmanager` |
| `print_info()` duplicating SCF fields | Delegates to `scf_->` for values |

## Jira

BTS-62 — "v3 GroundState: thin SCF wrapper replacing v2 GroundStateData duplication" (Done)
