# Computing response properties with madqc

How to compute molecular **response properties** — polarizability, hyper­polar­iz­ability,
Raman, and excited states — through `madqc --wf=response`. This is the per-property
companion to the general driver guide ([`README.md`](README.md)); for the engine itself
and its current capability matrix see
[`../molresponse/README.md`](../molresponse/README.md).

> **Engine.** `madqc --wf=response` runs the unified **molresponse_v3** engine —
> it is the default and the only engine (the older `engine v2` / MolresponseLib
> was removed; a deck that still names it errors out). The recipes below target
> v3; no engine line is needed. See
> [`../molresponse/MIGRATION_FROM_V2.md`](../molresponse/MIGRATION_FROM_V2.md).

---

## How property selection works

The `response` block selects properties two ways (you can mix them):

1. **High-level flags** — `dipole`, `quadratic`, etc. madqc derives the
   `requested_properties` list from them.
2. **Explicit list** — `requested_properties [polarizability, hyperpolarizability, raman]`.

A response run always needs the **perturbation directions and frequencies**
(`dipole.directions`, `dipole.frequencies`) and a **protocol ladder** (the
truncation-threshold schedule, set in the `dft` block as `protocol`). The run
solves the ground-state SCF first, then the response.

```text
dft
    xc        hf
    protocol  [1e-4, 1e-6]      # threshold ladder (k climbs 6 → 8)
    dconv     1e-4
end
response
    engine              v3
    dipole              true
    dipole.directions   xyz
    dipole.frequencies  [0.0, 0.04]
end
molecule
    units angstrom
    O   0.0  0.0     0.1173
    H   0.0  0.7572 -0.4692
    H   0.0 -0.7572 -0.4692
end
```

---

## Polarizability α  ✅ supported (static + dynamic)

The linear dipole–dipole response α(ω). Static (ω=0) and frequency-dependent.

| Knob | Meaning |
|------|---------|
| `dipole true` | request linear dipole response (→ `polarizability`) |
| `dipole.directions xyz` | Cartesian components to solve (`x`, `y`, `z`, or any subset) |
| `dipole.frequencies [0.0, 0.04]` | driving frequencies ω (a.u.); `[0.0]` = static only |

```text
response
    engine              v3
    dipole              true
    dipole.directions   z          # α_zz only (fast); use xyz for the full tensor
    dipole.frequencies  [0.0]
end
```

**Result** (`<prefix>/results.json` → `response` entry): `polarizability`, one
item per `frequency` with `components` such as `zz`, its `status` and
`convergence`. Example (H₂O/HF): α_zz(0)=8.53.

> Requesting a single direction (`z`) solves one perturbation channel — much
> faster, and the natural choice for a quick check.

---

## First hyperpolarizability β  ✅ supported (static + dynamic)

The quadratic response β, static and dynamic (SHG/OR) — a benchmark-quality
capability; validate against published reference β values. **Operational
caveat:** for dynamic β, make sure all VBC contributions reach the top protocol
(check the convergence diagnostics) — otherwise the assembled β is incomplete.

| Knob | Meaning |
|------|---------|
| `quadratic true` | compute β from the defined perturbations (adds `hyperpolarizability`) |
| `beta.shg true` | second-harmonic generation β(−2ω; ω, ω) (default) |
| `beta.or true` | optical rectification β(0; ω, −ω) |
| `dipole.frequencies` | the ω used to form the β process |

```text
response
    engine              v3
    dipole              true
    dipole.directions   z
    dipole.frequencies  [0.0]      # static β
    quadratic           true
end
```

**Result:** `hyperpolarizability`, one item per pair of driving frequencies
(`frequency_b`, `frequency_c`) with `components` such as `zzz`. Example
(H₂O/HF): β_zzz(0,0)≈7.76.

---

## Vibrational Raman  🟡 preliminary (full spectrum); single component ✅

Polarizability gradients ∂α_AB/∂Q from β(dipole; dipole, nuclear) — the response
to a nuclear displacement — and, for the full set, the Hessian, normal modes and
per-mode Raman activities. Closed shell only.

| Knob | Meaning |
|------|---------|
| `requested_properties [polarizability, raman]` | request the Raman (polarizability-gradient) property |
| `dipole.directions` / `dipole.frequencies` | the optical perturbation (frequencies of the derivative) |
| `raman.nuc_atom N`, `raman.nuc_axis 0/1/2` | one nuclear coordinate: nine ∂α_AB/∂Q rows |
| `raman.nuc_atom -1` | all 3N coordinates + Hessian + normal-mode spectrum |

```text
response
    dipole               true
    dipole.directions    xyz
    dipole.frequencies   [0.0]
    requested_properties [polarizability, raman]
    raman.nuc_atom       -1
    dconv                1e-3
end
```

The nuclear-displacement responses are ~10³ larger than the dipole ones and the
FD gate is absolute (`bsh < 5*dconv`), so a response `dconv 1e-3` is needed for
those legs at `protocol [1e-4]`; the derivative source is cusped and sharpens
only by climbing (e.g. `protocol [1e-4, 1e-6]`). **Result:**
`polarizability_gradient`, one item per optical `frequency` with `components`
keyed by coordinate `nuc_<atom>_<axis>` and then dipole pair (∂α_AB/∂Q_C); with
the full set also `vibrational_analysis` (Hessian, frequencies in hartree,
normal modes, reduced masses) and `raman_spectrum` (per optical frequency, per
mode `freq_cm1`, `alpha2`, `beta2`, `pol_int` in Å⁴/amu, `dep_ratio`, plus the
same invariants in a.u.). The legacy
`nuclear.*` keys are not read by this engine. Working decks:
`src/examples/qc/response_lih_raman_full`, `response_h2o_raman_tpa`.

---

## Excited states (TDA / RPA)  🟡 experimental

Excitation energies and (with analysis) transition properties. Eigenpairs
compute; **convergence and climbing across protocols are being stabilized**, so
treat ES results as experimental and check residuals.

| Knob | Meaning |
|------|---------|
| `excited.enable true` | run the excited-state bundle |
| `excited.num_states N` | number of roots to target |
| `excited.tda true` | Tamm-Dancoff approximation; `false` = full RPA (X,Y) |
| `excited.maxiter` / `excited.maxsub` | iteration cap / subspace size |

```text
response
    engine          v3
    excited.enable  true
    excited.num_states 3
    excited.tda     true
end
```

**Scope:** closed-shell only. Open-shell ES / open-shell full-RPA are out of
scope. ES / diffuse-state references use **d-aug-cc-pVQZ**.

---

## Convergence & protocol knobs

| Where | Knob | Effect |
|-------|------|--------|
| `dft` | `protocol [1e-4, 1e-6]` | truncation-threshold ladder; each rung raises k and reprojects |
| `dft` | `dconv` | SCF density convergence |
| `response` | `maxiter` | response iteration cap |
| `response` | `dconv` | response convergence target |
| `response` | `print_level` | 0–3, response verbosity |

A coarse single rung (`protocol [1e-4]`) is fastest for a smoke test; climb
(`[1e-4, 1e-6]`) for converged values (k 6 → 8).

---

## The response result (`results.json`)

Every response value is in `<prefix>/results.json`, in the `response` entry:
`polarizability`, `hyperpolarizability`, `polarizability_gradient`,
`excited_states`, `two_photon_absorption`, `vibrational_analysis`,
`raman_spectrum`, each item with its `status` and `convergence`; the
scheduler's verdict in `run`; each saved state in `state_metrics`. The format
is described in
[`../molresponse/docs/guides/results_file.md`](../molresponse/docs/guides/results_file.md),
and `<prefix>/summary.txt` prints it for people.

The response task in `<prefix>.calc_info.json` holds only `type` and the run's
times (`metadata`: `elapsed_time`, `mpi_size`, `nthreads`, `clock_scope`,
`v3_timing`), which the wall-time lines at the end of the log print. For the full top-level schema (SCF/CIS/correlated
tasks, tensor-json encoding) see the **Output** section of [`README.md`](README.md).

---

## See also

- [`README.md`](README.md) — the general madqc driver: workflows, CLI, input
  format, full `calc_info.json` schema, visualization.
- [`../molresponse/README.md`](../molresponse/README.md) — the v3 response
  engine: capability matrix, status, v2 relationship.
- [`../molresponse/MIGRATION_FROM_V2.md`](../molresponse/MIGRATION_FROM_V2.md)
  — switching a v2 deck to `engine v3` (and what to keep on v2 for now).
- `cm.sh` (`/gpfs/scratch/ahurtado/madness_es_bench/README.md`) — the build/run/
  validate harness; `cm_mq <mol>` drives exactly this response path.
