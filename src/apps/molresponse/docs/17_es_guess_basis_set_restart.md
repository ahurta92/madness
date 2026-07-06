# 17 — ES guess: virtual-orbital ("NWChem") guess + basis-set restart

Status: **research, 2026-06-09.** Motivated by the ES-convergence problem: the
v3 ES solver *converges*, but to the wrong / missing states unless the initial
guess already overlaps the target subspace. The fix is a better guess. This doc
maps what exists and lays out the paths to (a) a virtual-orbital ("NWChem-style")
guess and (b) restarting an MRA ES calc from a basis-set (Dalton/NWChem) result.
Feeds the ES-convergence work (the deliverable that unblocks ES/2PA/resonant-Raman
benchmarks and R2's ES-density export).

## The problem

ES is an eigenproblem; iterative diagonalization (KAIN/Davidson over the response
subspace) finds the states with the largest overlap with the *guess*. v3's two
guesses are both "cold":
- **Random** (`ESSolverGuess.hpp::make_initial_guess_tda_rhf`) — atom-localized
  noise, Q-projected to virtual space. No energy/symmetry information.
- **SolidHarmonics** (default, `create_solid_harmonics_guess`) — single-excitation
  trials `Y_lm(r)·φ_occ`, Q-projected. For L=1 this is `{x,y,z}·φ_occ` (a dipole
  guess) — good for bright dipole-allowed states, blind to others.

Neither orders trials by excitation energy or targets specific states, so the
lowest roots can be missed or scrambled (exactly the NWChem caveat: "low-lying
roots can be missed when they lack overlap with the initial guess — request
several more roots than needed"). Hence: a guess built from **virtual orbitals**,
ordered by orbital-energy difference (the CIS diagonal), is the standard remedy.

## What exists in MADNESS (the infrastructure is already here)

| Capability | Where | Note |
|---|---|---|
| Gaussian AO basis → MRA functions | `SCF::project_ao_basis` (SCF.cc:603), `AtomicBasisSet` (molecularbasis.h:149) | XML basis defs (sto-3g … aug-cc-pVQZ); battle-tested |
| AO coeffs → MRA MOs | `transform(world, ao, MOs)` | column=MO, row=AO |
| **External-result seam** | `ES_Interface` (ESInterface.h:56) | data model `{basis_set, atoms, energies, MOs, occupancies}` α+β, slymer-backed |
| Ground-state restart from NWChem | `SCF::initial_guess_from_nwchem` (SCF.cc:766), `NWChem_Interface : ES_Interface` (NWChem.h) | reads `.out`/`.movecs`, aligns geometry, builds MRA MOs |
| **Legacy virtual-orbital ES guess** | `molresponse/ExcitedResponse.cpp`: `create_virtual_ao_guess` (3030), `make_nwchem_trial` (178) | project AO→virtual (Q), diagonalize Fock there, seed singles. **Not in v3.** |
| Dalton excitation energies + response tensors | gecko `plugins/dalton/parse.py` (`parse_electronic_excitations`, polarizability) | parses `DALTON.OUT` text; gives ω/symmetry, **not** MO/excitation vectors |
| v3 guess storage | `ResponseStateX<ClosedShell>.x_alpha[n_occ]` | occupied-indexed response orbitals; virtuals implicit via `Q` |

**No** Molden reader and **no** Dalton `ES_Interface` exist in MADNESS yet;
slymer is the Gaussian parser the NWChem path uses.

## The paths (increasing fidelity / effort)

### Path A — in-house virtual-orbital guess (the "NWChem guess") — RECOMMENDED FIRST
Self-contained; needs only an AO basis file (already supported). Port the legacy
`create_virtual_ao_guess` to v3:
1. `project_ao_basis` a chosen basis (e.g. aug-cc-pVDZ) → MRA AO functions.
2. `Q`-project to the virtual space; build the Fock/Hessian in that space and
   diagonalize → approximate virtual orbitals `{φ_a, ε_a}`.
3. Emit ES trial vectors as the **lowest (ε_a − ε_i) single excitations** — trial
   for `i→a` = `φ_a` in occupied slot `i` of `x_alpha` (the v3 X shape). Order by
   ε_a−ε_i; oversample (already have `run_oversampled_tda_warmup` to downselect).

This *is* NWChem's TDDFT Davidson initial guess (CIS diagonal). It directly fixes
"converges to the wrong states" because the trial subspace is the actual
low-energy excitation space, energy-ordered. New `ESGuessMode::VirtualAO`.

### Path B — restart from an external job's MOs (NWChem now; Dalton via new reader)
Reuse `ES_Interface` → occupied+virtual MOs projected to MRA (better virtuals
than an in-house minimal basis), then build the Path-A guess from the *external*
virtuals. NWChem works today (`initial_guess_from_nwchem`); Dalton needs a
`Dalton_Interface : ES_Interface` — **or** a Molden→`ES_Interface` adapter (Dalton
writes Molden: basis + MO coeffs + occupancies). The seam + AO→MRA pipeline
already exist; the gap is the reader.

### Path C — seed directly from the external EXCITED-STATE vectors (the answer)
Read the basis-set excitation amplitudes `c_ia^(n)` (Dalton response vectors) and
reconstruct each response orbital `x_i^(n) = Σ_a c_ia^(n) φ_a` on the MRA grid →
seed the v3 ES guess with the basis-set excited states directly. Highest fidelity
(it *is* the answer in a finite basis; MRA refines it) and the strongest form of
state targeting. Most plumbing: needs the MO coeffs **and** the excitation
vectors — gecko parses Dalton's ω's but not yet the vectors; would need Dalton's
Molden (MOs) + response-vector output (`RSPVEC`/SIRIUS) or the daltonproject API.

### State targeting (cheap, immediate, orthogonal)
gecko **already** extracts Dalton's excitation energies. Even with no vectors, the
Dalton ω-ordering/symmetries can (i) be the correctness reference, and (ii) set
`n_roots` / select-by-symmetry / order the roots. Path C is the strongest version.

## Recommended sequencing
1. **Path A — DONE (`18f853182`, 2026-06-09).** `ESGuessMode::VirtualAO` in
   `ESSolverGuess.hpp` (aug-cc-pvdz, energy-ordered singles). Validated: h2o
   recovers all four roots in order at the coarse rung (0.3196/0.3807/0.4101/0.4212
   vs target 0.317/0.378/0.399/0.409) — **incl. the 0.378 root SolidHarmonics
   missed** (it found a spurious 0.46). Upper two tighten at 1e-6/k8 (solver, not
   guess). Basis defaulted; `--es-guess-basis` knob = follow-up. Open-shell falls
   back to solid harmonics.
2. **Targeting** — feed Dalton ω's (gecko) as reference + to set/select roots.
3. **Path B/C (Dalton)** — `Dalton_Interface : ES_Interface` (or Molden adapter) →
   external virtuals (B), then extend gecko to extract the excitation vectors for a
   direct state seed (C).

## Status (2026-06-09): Path A landed; B + C → FUTURE FEATURES
Path A (`ESGuessMode::VirtualAO`) is the active, validated win. **Paths B (Dalton
`ES_Interface`/Molden reader) and C (seed from Dalton excitation vectors) are
DEFERRED to future features** — the seam (`ES_Interface` + AO→MRA projection)
exists, so they can be picked up when external-restart / state-targeting from a
basis-set job is needed. Not on the near-term path.

## Gaps to fill
- v3 has no virtual-orbital guess (Path A) — the legacy code is the template.
- No Molden reader / Dalton `ES_Interface` in MADNESS (Path B).
- gecko parses Dalton ω/tensors but not MO coeffs or excitation vectors (Path C).

## Sources
- NWChem TDDFT guess (Davidson, orbital-energy-difference singles, oversample roots):
  [NWChem Excited-State Calculations](https://nwchemgit.github.io/Excited-State-Calculations.html)
- Dalton / Molden interface (MO coeffs, response/CI vectors):
  [Dalton Project (JCP 2020)](https://pubs.aip.org/aip/jcp/article/152/21/214115/198790),
  [grid-data interop / Molden](https://arxiv.org/pdf/2401.17925)
