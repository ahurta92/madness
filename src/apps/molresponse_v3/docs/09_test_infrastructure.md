# Test Infrastructure & Fixture Management

## Overview

The molresponse_v3 test infrastructure separates **definitions** (checked into
git) from **data** (generated in scratch). A Python orchestration layer
(`fixture_manager.py`) bridges the two, generating inputs, submitting SLURM
jobs, and validating results.

```
Git (definitions)                         Scratch (data)
─────────────────                         ──────────────
tests/fixtures/                           $GECKO_SCRATCH/molresponse_fixtures/
├── numerical_settings.json               ├── systems/<name>/<tier>/
├── slurm_profiles.json                   │   ├── moldft.in
├── reference_db.json                     │   ├── job.sh
├── systems/<name>/                       │   ├── job_status.json
│   ├── geometry.xyz                      │   ├── slurm-<jobid>.out
│   ├── moldft.in (reference)             │   ├── moldft.restartdata.*
│   └── system.json                       │   └── moldft.calc_info.json
└── property/<category>/<name>.json       └── response/<fixture>/<tier>/
                                              ├── response.in
tests/fixture_manager.py                      └── *.calc_info.json
```

---

## Standard Numerical Tiers

Three accuracy tiers control all DFT and response calculations. The same
settings apply regardless of XC functional (HF, LDA, PBE, etc.).

Defined in `tests/fixtures/numerical_settings.json`.

### Tier definitions

| Parameter | Low | Medium | High |
|-----------|-----|--------|------|
| **dconv** | 1e-2 | 1e-4 | 1e-6 |
| **econv** | 1e-3 | 1e-5 | 1e-7 |
| **protocol** | [1e-4] | [1e-4, 1e-6] | [1e-4, 1e-6, 1e-8] |
| **eprec** | 1e-3 | 1e-4 | 1e-6 |
| **auto k (final)** | 6 | 8 | 10 |
| **Use case** | Smoke tests, dev | Validation, parity | Benchmarks, publication |

### How protocol interacts with convergence

The effective density convergence in the SCF solver (see `SCF.cc:2143`) is:

```
effective_dconv = max(FunctionDefaults<3>::get_thresh(), param.dconv())
```

This means the final protocol value must be <= dconv, otherwise the protocol
step is wasted — you can't converge density tighter than the MRA representation
allows. The auto k selection is:

| thresh | k |
|--------|---|
| >= 0.9e-2 | 4 |
| >= 0.9e-4 | 6 |
| >= 0.9e-6 | 8 |
| >= 0.9e-8 | 10 |

Note: `econv` is **not used** in the SCF convergence check — only `dconv`
matters. `econv` is recorded in metadata only.

### Box size (L)

L should be at least 10x the spatial extent of the molecule. L=200.0 bohr is
a safe default for most systems up to naphthalene. Small systems (H2, He, Li)
use smaller L values specified in their `system.json` `scf_settings`.

---

## System Fixtures

Nine molecular systems, each in `tests/fixtures/systems/<name>/`:

| System | Molecule | Spin | Occupied | Size Class |
|--------|----------|------|----------|------------|
| h2_hf | H2 | restricted | 1 | tiny |
| he_hf | He | restricted | 1 | tiny |
| li_uhf | Li | unrestricted | 2 (1α+2β) | tiny |
| lih_hf | LiH | restricted | 2 | small |
| oh_uhf | OH | unrestricted | 5 (4α+5β) | small |
| no_uhf | NO | unrestricted | 8 (7α+8β) | small |
| h2o_hf | H2O | restricted | 5 | small |
| c6h6_hf | C6H6 | restricted | 21 | medium |
| naphthalene_hf | C10H8 | restricted | 34 | large |

Each system directory contains:
- `geometry.xyz` — molecular geometry in angstrom
- `moldft.in` — reference input (for manual runs)
- `system.json` — metadata: atoms, electrons, size class, SCF settings

Geometry sources: H2O and LiH from v2 test inputs, C6H6 and naphthalene
from `GECKO_MOL_LIB`, open-shell systems from `08_test_fixtures.md`.

---

## Property Fixtures

Fifteen property test configurations across seven categories, each in
`tests/fixtures/property/<category>/<name>.json`.

Each fixture specifies:
- `system` — which system to use
- `tier` — which numerical tier (default: medium, protocol_study uses high)
- `property` — what to compute
- `frequencies` — optical frequencies
- `response_settings` — solver parameters
- `expected` — reference values for validation
- `reference_sources` — where reference values came from

### Categories

| Category | Fixtures | Purpose |
|----------|----------|---------|
| alpha_static | h2, he, h2o, oh_uhf | Static polarizability |
| alpha_dynamic | h2, h2o | Frequency-dependent alpha |
| beta_shg | h2o, lih | SHG hyperpolarizability |
| raman | h2o, lih | Raman spectra |
| excited_tda | h2, h2o | TDA excitation energies |
| excited_full | h2 | Full TDDFT excitations |
| protocol_study | h2o_alpha, lih_raman | Multi-protocol convergence (high tier) |

### Property composition

Response runs are naturally multi-property. Individual fixtures exist for
focused testing, but they can be composed into combined runs:
- Raman includes alpha (dipole) + nuclear responses
- Beta includes alpha (dipole) + quadratic
- A "production" run can compute alpha + beta + raman in a single invocation

The `merge_response_flags()` function in fixture_manager.py handles this
composition automatically.

---

## Reference Database

`tests/fixtures/reference_db.json` indexes existing MRA and Dalton reference
values from prior publication-quality calculations.

### Sources

| Source | Description | Accuracy | Systems |
|--------|-------------|----------|---------|
| Beta paper d06 | 72 molecules, HF alpha+beta | dconv=1e-6 | H2O, LiH, + 70 more |
| Raman paper mra-p07 | 7 molecules, HF raman | dconv=1e-4 to 1e-7 | H2O, CH4, NH3, SO2, DMSO, Anthracene, Pyrene |
| H2 TDA test | H2 excited states | legacy k=8 + Dalton d-aug-cc-pVQZ | H2 |

### Extracted values

**H2O** — Full alpha tensor at 9 frequencies (0.0 to 0.159 au):
- Static: xx=7.902, yy=9.192, zz=8.531 (isotropic=8.542)
- Beta tensor available via gecko `extract_beta`

**LiH** — Full alpha tensor at 9 frequencies (0.0 to 0.074 au):
- Static: xx=yy=25.272, zz=21.866 (isotropic=24.137)
- Beta tensor available via gecko `extract_beta`

**H2** — TDA excitation energies (4 roots):
- MADNESS legacy k=8: 0.46789, 0.47771, 0.48120, 0.48120 au
- Dalton d-aug-cc-pVQZ: 0.465376, 0.477069, 0.480889, 0.480889 au

### Known gap: Raman MRA values

The raman paper MRA runs have converged dipole and nuclear response functions
but the property assembly step did not produce a merged `calc_info.json` with
raman spectra. The data exists in the response archives — needs either a
re-run of PropertyManager or a gecko enhancement to reconstruct from archives.

### Using the reference database

Use gecko MCP to extract values from any of the 72 beta paper molecules:
```python
gecko.load_calculation("/path/to/beta_paper/d06/output/CH3OH")
gecko.extract_alpha("/path/to/beta_paper/d06/output/CH3OH")
gecko.extract_beta("/path/to/beta_paper/d06/output/CH3OH")
```

---

## SLURM Integration

### Cluster profiles

`tests/fixtures/slurm_profiles.json` defines three Seawulf cluster types:

| Cluster | Architecture | Cores | Memory | Build Dir |
|---------|-------------|-------|--------|-----------|
| xeonmax | Intel Sapphire Rapids | 96 | 128 GB HBM + 256 GB DDR5 | builds/.../debug |
| 40core | Intel Skylake | 40 | 192 GB DDR4 | builds/.../40core |
| 96core | AMD Milan | 96 | 256 GB DDR4 | builds/.../96core |

Each cluster specifies:
- `env_script` — environment loader (e.g. `~/load_xeonmax.sh`)
- `build_dir` — path to the compiled binaries
- `partitions` — SLURM partition names by time class
- `numa_mask` — for HBM-aware placement (xeonmax only)
- `account` — SLURM account (`pn_roha020218s`)

### Resource profiles

SLURM resources are determined by `(system_size_class, tier)`:

| Size Class | Systems | Medium: Nodes x Tasks | High: Nodes x Tasks |
|------------|---------|----------------------|---------------------|
| tiny | H2, He, Li | 1x1 (10 min) | 1x2 (30 min) |
| small | LiH, OH, NO, H2O | 1x2 (1 hr) | 1x4 (4 hr) |
| medium | C6H6 | 2x4 (4 hr) | 4x4 (12 hr) |
| large | naphthalene | 4x4 (8 hr) | 6x4 (24 hr) |

### Environment scripts

The `~/load_*.sh` scripts set up modules and Intel oneAPI. They guard
the Intel `setvars.sh` sourcing with `set +u`/`set -u` because Intel's
scripts have unbound variables that fail under bash strict mode.

### Generated job scripts

Each SLURM job script:
- Sources the cluster's environment script
- Sets `MAD_NUM_THREADS` and `OMP_NUM_THREADS`
- Changes to the work directory (all files co-located)
- Logs start/end timestamps and metadata
- Directs stdout/stderr to the work directory (`slurm-%j.out/err`)

---

## fixture_manager.py

### CLI commands

```bash
# Set GECKO_SCRATCH first (or source setenv.sh)
export GECKO_SCRATCH=/gpfs/scratch/ahurtado/gecko_calcs

# Show inventory: systems, fixtures, checkpoint status
python3 fixture_manager.py status

# Generate moldft.in for all systems at a tier
python3 fixture_manager.py generate-ground-states --tier medium

# Generate SLURM job.sh for all systems
python3 fixture_manager.py generate-job-scripts --tier medium --cluster xeonmax

# Generate + submit all pending ground states
python3 fixture_manager.py submit --tier medium --cluster xeonmax

# Check SLURM job status
python3 fixture_manager.py check --tier medium

# Validate results against expected values
python3 fixture_manager.py validate alpha_static/h2o --tier medium
```

### Options

| Flag | Description |
|------|-------------|
| `--tier` | low, medium, or high |
| `--system` | Single system name (default: all) |
| `--cluster` | xeonmax, 40core, or 96core |
| `--xc` | Exchange-correlation functional (default: hf) |
| `--force` | Re-submit even if checkpoint exists |

### Python API

```python
from fixture_manager import FixtureManager

fm = FixtureManager()

# Generate inputs
fm.generate_moldft_input("h2o_hf", "medium")
fm.generate_job_script("h2o_hf", "medium", cluster="xeonmax")

# Submit and track
job_id = fm.submit_ground_state("h2o_hf", "medium", cluster="xeonmax")

# Check
results = fm.check_jobs("medium")

# Validate
result = fm.validate("alpha_static/h2o", tier="medium")
```

---

## Raman Benchmark Database

Separate from the v3 test fixtures, a raman benchmarking database lives at
`$GECKO_SCRATCH/raman_benchmark/`. It mirrors the raman paper data structure
with updated numerical settings for publication-quality results.

### Structure

```
$GECKO_SCRATCH/raman_benchmark/
├── <Molecule>/
│   ├── <Molecule>.mol          # geometry file
│   ├── benchmark.json          # metadata
│   ├── mra/
│   │   ├── <mol>_raman.in      # madqc input
│   │   ├── job.sh              # SLURM script
│   │   └── (outputs after run)
│   └── dalton/                 # symlinks to raman paper Dalton results
│       ├── aug-cc-pVDZ -> ...
│       ├── aug-cc-pVTZ -> ...
│       └── aug-cc-pVQZ -> ...
```

### Molecules

| Molecule | Atoms | Occ | Dalton bases |
|----------|-------|-----|--------------|
| H2O | 3 | 5 | DZ, TZ, QZ, d-QZ |
| CH4 | 5 | 5 | DZ, TZ, QZ, d-QZ |
| NH3 | 4 | 5 | DZ, TZ, QZ |
| SO2 | 3 | 16 | DZ, TZ, QZ, d-DZ, d-QZ |
| DMSO | 10 | 21 | (none yet) |
| Anthracene | 24 | 47 | DZ, TZ, QZ |
| Pyrene | 26 | 55 | DZ, TZ, QZ |

### Numerical settings

Modified high-accuracy protocol for raman benchmarking:

| Parameter | Value |
|-----------|-------|
| dconv | 1e-5 |
| econv | 1e-6 |
| protocol | [1e-4, 1e-6, 1e-7] |
| eprec | 1e-6 |
| gopt | true |
| state_parallel | on |
| state_parallel_groups | 64 (8 nodes x 8 tasks/node) |
| frequencies | [0.0, 0.02, 0.04, 0.06, 0.08, 0.1] |

### Dalton results

Dalton basis set results from the original raman paper are **symlinked**, not
copied. These are not recomputed. They serve as finite-basis-set reference
points for the MRA basis-set-limit comparison.

---

## Jira Tracking

- **BTS-36** — "v3 Inc 0: Skeleton app + test fixture library + orchestration"
  (In Progress, under epic BTS-21 "molresponse_v3 — Ground-Up Rebuild")

---

## Files Created in This Increment

### Skeleton app (`src/apps/molresponse_v3/`)
- `CMakeLists.txt` — build config, links MADresponse2 + MADchem
- `main.cpp` — loads ground state, prints summary
- `README.md` — description and build instructions

### Test fixtures (`src/apps/molresponse_v3/tests/`)
- `fixture_manager.py` — orchestration (generate, submit, check, validate)
- `generate_systems.sh` — shell script for manual checkpoint generation
- `fixtures/numerical_settings.json` — three accuracy tiers
- `fixtures/slurm_profiles.json` — cluster configs and resource profiles
- `fixtures/reference_db.json` — MRA reference values index
- `fixtures/systems/*/` — 9 system directories (geometry, input, metadata)
- `fixtures/property/*/` — 15 property fixture JSONs across 7 categories

### GroundState (SCF wrapper)
- `GroundState.hpp` — thin wrapper around `shared_ptr<SCF>`, response-specific cached operators
- `GroundState.cpp` — `from_archive()` factory, `prepare()` for response preliminaries

### Build integration
- `src/apps/CMakeLists.txt` — added `add_subdirectory(molresponse_v3)`

---

## Known Issues

### Li (open-shell) fixture
The Li ground-state checkpoint appears to have been generated as spin-restricted
despite `nopen 1` in the input. The moldft input may need explicit
`spin_restricted false`. The GroundState wrapper handles this correctly
(infers nopen from archive header), but the fixture data itself may need
regeneration.

### Raman benchmark MPI errors
H2O and NH3 raman benchmarks computed all spectra successfully but crashed at
MPI finalization with `MPI_ERR_TRUNCATE`. CH4 hit a hung queue error from
state_parallel load imbalance. SO2 likely OOM at 8 ppn. Mitigations: reduce
`state_parallel_groups` and `maxiter`, use 1 ppn for SO2.
