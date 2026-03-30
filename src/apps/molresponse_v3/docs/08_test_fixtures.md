# molresponse_v3 Test Fixture Library

## Strategy

Organize all test fixtures into a structured library with two layers:
1. **System library** — geometries and ground-state checkpoints for each molecule
2. **Property fixtures** — specific test configurations per property per system

The system library is shared across all testing (v2, v3, legacy, gecko).
Property fixtures are specific test setups that reference a system and
define what to compute and what the expected results are.

Gecko can be used to collect geometries for new systems and to generate
Gaussian basis set reference values for comparison.

---

## System Library

### Directory Structure

On Seawulf, a central fixture library:
```
/gpfs/projects/rjh/adrian/test_fixtures/
├── systems/
│   ├── h2_hf/
│   │   ├── geometry.xyz
│   │   ├── moldft.in
│   │   ├── moldft.restartdata       (generated)
│   │   ├── moldft.calc_info.json    (generated)
│   │   └── system.json              (metadata)
│   ├── he_hf/
│   ├── h2o_hf/
│   ├── lih_hf/
│   ├── c6h6_hf/
│   ├── naphthalene_hf/
│   ├── oh_uhf/                      (open-shell)
│   ├── no_uhf/                      (open-shell)
│   └── li_uhf/                      (open-shell, atom)
├── property_fixtures/
│   └── (see below)
└── reference_values/
    └── (gecko-generated comparisons)
```

In the repo, we check in the inputs and metadata (not the archives):
```
src/apps/molresponse_v3/tests/
├── fixtures/
│   ├── systems/
│   │   ├── h2_hf/
│   │   │   ├── geometry.xyz
│   │   │   ├── moldft.in
│   │   │   └── system.json
│   │   ├── he_hf/
│   │   ├── h2o_hf/
│   │   ├── lih_hf/
│   │   ├── c6h6_hf/
│   │   ├── naphthalene_hf/
│   │   ├── oh_uhf/
│   │   ├── no_uhf/
│   │   └── li_uhf/
│   └── property/
│       ├── alpha_static/
│       ├── alpha_dynamic/
│       ├── beta_shg/
│       ├── raman/
│       ├── excited_tda/
│       ├── excited_full/
│       └── (future: tpa, resonant_raman, ...)
├── generate_systems.sh
└── generate_property_fixtures.sh
```

### System Metadata (system.json)

Each system gets a metadata file describing it:
```json
{
  "name": "h2o_hf",
  "description": "Water, Hartree-Fock, closed-shell",
  "molecule": "H2O",
  "method": "hf",
  "spin": "restricted",
  "n_atoms": 3,
  "n_electrons": 10,
  "n_occupied": 5,
  "geometry_source": "optimized at HF/aug-cc-pVQZ",
  "notes": "Standard small polyatomic test case",
  "scf_settings": {
    "xc": "hf",
    "k": 6,
    "L": 200.0,
    "protocol": [1e-4, 1e-6]
  }
}
```

---

## Closed-Shell Systems

### H2 — Minimal (1 occupied orbital)

```
H  0.0  0.0  -0.7
H  0.0  0.0   0.7
```

SCF: `xc=hf, k=5, L=20.0, protocol=[1e-4]`

Purpose: fastest possible iteration. Development workhorse. Every
increment gets tested here first.

### He — Single Atom (1 occupied orbital)

```
He  0.0  0.0  0.0
```

SCF: `xc=hf, k=6, L=20.0, protocol=[1e-2]`

Purpose: spherical symmetry edge case. No bond, no nuclear displacement
modes (for Raman). Tests that operators handle atoms correctly.
Already used in excited-state restart tests.

### H2O — Small Polyatomic (5 occupied orbitals)

```
O   0.0000   0.0000   0.2217
H   0.0000   1.4309  -0.8867
H   0.0000  -1.4309  -0.8867
```

SCF: `xc=hf, k=6, L=200.0, protocol=[1e-4, 1e-6]`

Purpose: standard validation target. Reference values exist for alpha,
beta, raman from v2. Multi-protocol for protocol-accuracy studies.

### LiH — Small Heteronuclear (2 occupied orbitals)

```
Li  0.0  0.0  0.0
H   0.0  0.0  3.015
```

SCF: `xc=hf, k=6, L=40.0, protocol=[1e-4]`

Purpose: different electronic structure from H2O. Existing Raman test.
Large permanent dipole — good for testing beta.

### C6H6 — Medium (21 occupied orbitals)

Geometry: use existing fixture (from your current test data).
Source via gecko if needed.

SCF: `xc=hf, k=6, L=200.0, protocol=[1e-4, 1e-6]`

Purpose: medium-sized system for scaling tests. Enough orbitals to
stress-test memory and parallelism. High symmetry.

### Naphthalene — Larger (34 occupied orbitals)

Geometry: use existing fixture.

SCF: `xc=hf, k=6, L=200.0, protocol=[1e-4, 1e-6]`

Purpose: larger system for state-parallel scaling. Realistic size for
production response calculations.

---

## Open-Shell Systems

### OH Radical — Minimal Open-Shell (4α + 3β)

```
O   0.0  0.0  0.0
H   0.0  0.0  1.832
```

SCF: `xc=hf, k=6, L=40.0, protocol=[1e-4], nopen=1`

Purpose: simplest open-shell diatomic. Tests that the type system
handles spin-unrestricted building blocks correctly. Small enough
for fast iteration.

### NO — Open-Shell Heteronuclear (7α + 6β)

```
N   0.0  0.0  0.0
O   0.0  0.0  2.175
```

SCF: `xc=hf, k=6, L=40.0, protocol=[1e-4], nopen=1`

Purpose: more electrons than OH, different symmetry. Standard benchmark
for open-shell response properties.

### Li Atom — Open-Shell Atom (1α + 0β in minimal, or 2α + 1β)

```
Li  0.0  0.0  0.0
```

SCF: `xc=hf, k=6, L=20.0, protocol=[1e-4], nopen=1`

Purpose: open-shell atom. Spherical symmetry + spin = double edge case.
Tests that open-shell operators work for atoms.

---

## Property Fixtures

Each property fixture defines: which system, what perturbation, what
frequencies, what protocol, and what the expected output looks like.

### Structure

```
src/apps/molresponse_v3/tests/fixtures/property/
├── alpha_static/
│   ├── h2.json        # H2 static polarizability test config
│   ├── h2o.json       # H2O static polarizability test config
│   └── oh_uhf.json    # OH open-shell static polarizability
├── alpha_dynamic/
│   ├── h2.json        # H2 dynamic polarizability at specific frequencies
│   ├── h2o.json
│   └── oh_uhf.json
├── beta_shg/
│   ├── h2o.json       # H2O SHG hyperpolarizability
│   └── lih.json
├── raman/
│   ├── h2o.json       # H2O Raman
│   └── lih.json
├── excited_tda/
│   ├── h2.json        # H2 TDA (matches test_tda_h2)
│   └── h2o.json
├── excited_full/
│   ├── h2.json
│   └── h2o.json
└── protocol_study/
    ├── h2o_alpha.json  # Multi-protocol alpha for T3 workflow
    └── lih_raman.json  # Multi-protocol Raman for T3 workflow
```

### Property Fixture Format

```json
{
  "name": "h2o_alpha_static",
  "system": "h2o_hf",
  "property": "polarizability",
  "perturbation": "dipole",
  "directions": "xyz",
  "frequencies": [0.0],
  "protocol": [1e-4, 1e-6],
  "response_settings": {
    "maxiter": 15,
    "kain": true,
    "dconv": 1e-4
  },
  "expected": {
    "alpha_tensor_1e-6": {
      "xx": null,
      "yy": null,
      "zz": null,
      "note": "fill from v2 reference or gecko"
    }
  },
  "reference_sources": {
    "v2_serial": "path to v2 calc_info.json",
    "gecko_gaussian": "path to gecko comparison",
    "legacy": null
  }
}
```

The `expected` values start as null and get filled in as we generate
reference data. The `reference_sources` track where each reference
came from.

---

## Fixture Generation

### Agent Task: Organize Existing Fixtures

```
I need you to organize existing test fixtures into the v3 fixture
library structure. Here is what currently exists:

Existing checkpoint directories (on Seawulf):
[FILL IN: paths to your existing h2, h2o, c6h6, naphthalene checkpoints]

Existing mol files (on Seawulf):
[FILL IN: path to your mol file library]

For each existing system:
1. Create the system directory under tests/fixtures/systems/
2. Extract or create geometry.xyz from the checkpoint
3. Create moldft.in that reproduces the checkpoint
4. Create system.json with the metadata
5. Note the Seawulf path to the existing checkpoint

For new systems (He, LiH, OH, NO, Li):
1. Create the system directory
2. Write geometry.xyz (use geometries from this document)
3. Write moldft.in
4. Create system.json
5. Generate the checkpoint on Seawulf using generate_systems.sh

Create generate_systems.sh that can regenerate all checkpoints.
Create generate_property_fixtures.sh as a placeholder.

Output: report what was created, what existing data was found,
and what needs to be generated on Seawulf.
```

### Gecko Integration for New Systems

For systems where you need optimized geometries or don't have mol files:

```bash
# Use gecko to look up geometry for a molecule
gecko geometry NO --basis=aug-cc-pVQZ --method=hf

# Use gecko to generate a MADNESS input from a known geometry
gecko madness-input NO --basis=aug-cc-pVQZ --method=hf --output=no_uhf/moldft.in
```

For generating Gaussian basis set reference values:

```bash
# Run Gaussian calculation via gecko for comparison
gecko run NO --basis=aug-cc-pVQZ --method=hf --property=polarizability
gecko compare NO --madness=no_uhf/results.json --gaussian=no_uhf/gaussian.json
```

---

## Which Fixtures for Which Increment

| Increment | Systems Needed | Properties Tested |
|-----------|---------------|-------------------|
| 0 | H2, He | ground-state loading only |
| 1 | H2 | building blocks (no solve) |
| 2 | H2, H2O | static alpha (FD static solver) |
| 3 | H2, H2O | dynamic alpha (FD full solver) |
| 4 | H2 | naming, save/load, restart |
| 5 | H2O | end-to-end alpha via orchestrator |
| 6 | H2O, LiH | alpha + beta + raman |
| 7 | H2 | TDA excited states |
| 8 | H2 | full excited states |
| 9 | H2 | ES naming and restart |
| 10 | H2O | madqc integration |
| 11 | H2O, C6H6 | state-parallel scaling |

Open-shell fixtures (OH, NO, Li) come into play when open-shell
support is added — likely as a sub-increment after Inc 3 or as a
dedicated increment.

---

## Reference Value Collection

For each fixture × property combination, collect references from
multiple sources and track them:

| System | Property | v2 Serial | Legacy | Gecko/Gaussian | Status |
|--------|----------|-----------|--------|----------------|--------|
| H2 | static alpha | TBD | TBD | TBD | generate at Inc 2 |
| H2 | TDA omega_1 | TBD | TBD | TBD | generate at Inc 7 |
| He | static alpha | TBD | — | TBD | generate at Inc 2 |
| H2O | static alpha | exists | — | TBD | extract from ref JSON |
| H2O | dynamic alpha | exists | — | TBD | extract from ref JSON |
| H2O | beta (SHG) | exists | — | TBD | extract from ref JSON |
| LiH | raman | exists | — | TBD | extract from ref JSON |
| C6H6 | static alpha | TBD | — | TBD | generate at Inc 11 |
| OH | static alpha | — | — | TBD | open-shell milestone |
| NO | static alpha | — | — | TBD | open-shell milestone |

Values marked "exists" can be extracted from the v2 reference JSONs
already checked into the repo (e.g.,
`mad_madqc_test_molresponse_h2o_alpha_beta_z.py.calc_info.ref.json`).

This table grows with each increment and becomes the living record
of the protocol-accuracy work (T3).
