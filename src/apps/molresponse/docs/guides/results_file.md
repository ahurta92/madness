# The results file: `results.json`

`<run>/results.json` holds every value a person wants from a run folder, and
next to each value how well converged the states it was computed from were.
Facts about one saved state are in that state's node folder under `states/`;
time and memory are in `timings.json`. Each application replaces its own entry
when it finishes and leaves the others alone. A file of another version is
replaced, not merged.

## Top level

```json
{
  "version": 4,
  "provenance": {
    "code_version": "<git commit the program was built from>",
    "launched_at": "2026-11-25T09:12:03.512Z",
    "mpi_processes": 1,
    "threads_per_process": 16,
    "archive_io": "native"
  },
  "calculations": [ <entry>, ... ]
}
```

`provenance` describes the launch that wrote the file last. `code_version` is
recorded when the program is built, so it names the commit of the last build
(`git rev-parse HEAD` then), with `-dirty` appended when tracked files had
uncommitted changes.
`threads_per_process` counts the pool threads and the main thread, as
`timings.json` and `calc_info.json` do. `archive_io` is `native` or `hdf5`.

## One entry

```json
{
  "calculation_id": "response",
  "application": "molresponse",
  "protocol": {"threshold": 0.0001, "k": 6},
  "targets": {
    "driven_states": {"bsh_residual": 0.0005, "density_change": 0.0005},
    "excited_states": {"bsh_residual": 0.0005, "density_change": 0.0005, "energy_change": 0.0005}
  },
  "properties": { "<name>": <property>, ... },
  "run": {"stop_reason": "complete", "dropped_work": [], "stalled": []},
  "state_metrics": {
    "dipole_z-f0.00000": {
      "protocol": {"threshold": 0.0001, "k": 6},
      "subjects": {"response": {"values": {"bsh_residual": 2.3e-06, "density_change": 3.4e-05}, "availability": {}}}
    }
  }
}
```

- `calculation_id` is `ground_state` (written by `moldft`) or `response`.
- `protocol` is the final level of the request.
- `targets` is the final-level target of each metric, as the convergence
  decision module resolved it for the application. When the response entry
  solves driven and excited states with different targets, each kind is named
  (`driven_states`, `excited_states`); otherwise `targets` is one object of
  metric to target. A driven state's residual and density targets use the
  scaled rule: the state met them when its value was at most target × (1 +
  scale). Targets are written to three significant figures.
- `convergence` appears on the entry instead of on its properties when every
  value of the entry comes from the same state: the `ground_state` entry has
  one `convergence`, the worst measured values of its saved state, and its
  properties have none.
- `run` appears only on the response entry. `dropped_work` and `stalled` are
  always lists, empty when nothing was dropped.
- `state_metrics` is keyed by node name: the per-state table behind every
  `convergence`. Each node gives the level it was saved at and the measured
  values of each of its subjects. A driven state has one subject (`response`),
  an excited-state set one per root (`es_root_0000`, …), the ground state one
  (`ground`). A second-order source is built rather than iterated, so its
  `subjects` is empty.

## One property

A property with one answer has `value` (on the `ground_state` entry the
`convergence` is the entry's, so the property has none):

```json
"total_energy": {
  "units": "hartree",
  "value_type": "real",
  "value": -1.1335935503913517,
  "status": "complete"
}
```

A property with several answers (one per frequency, frequency pair or root) has
a `values` list. Each item has the keys that say what it belongs to, then
`value` or `components`, then `status` and `convergence`:

```json
"polarizability": {
  "units": "bohr^3",
  "value_type": "real",
  "values": [
    {"frequency": 0.0,
     "components": {"zz": 6.4364},
     "status": "complete",
     "convergence": {"bsh_residual": 2.25e-06, "density_change": 3.42e-05}}
  ]
}
```

The excited states are one property, `excited_states`, with one item per
root. Each value of a root is its own key, so `units` maps each to its unit:

```json
"excited_states": {
  "value_type": "real",
  "units": {"excitation_energy": "hartree", "oscillator_strength": "dimensionless",
            "transition_dipole": "atomic_units", "transition_quadrupole": "atomic_units"},
  "definitions": {"transition_quadrupole": "nontraceless second moment"},
  "values": [
    {"root": "es_root_0000",
     "excitation_energy": 0.46583,
     "oscillator_strength": 0.28548,
     "transition_dipole": {"x": 0.0, "y": 0.0, "z": -0.95878},
     "transition_quadrupole": {"xx": 0.0, "xy": 0.0, ...},
     "status": "complete",
     "convergence": {"bsh_residual": 3.85e-05, "density_change": 0.000109, "energy_change": 1e-07}}
  ]
}
```

A root value that was not computed is absent from its item and named in
`availability`; the item is then `partial`, or `unavailable` when nothing was
computed.

- A property holds the value at the final protocol level only. Earlier levels
  are in the node histories and the log.
- `status` (`complete`, `partial` or `unavailable`) is on every property or
  item, decided when the file is written. An `unavailable` value has no
  `value`; its components, if any, are empty. For driven states `status` is
  the verdict, because the scaled rule lets a value slightly above its target
  still pass.
- `availability` appears only when something is missing: it lists each missing
  component (or `value`) with the reason, which names the state that fell short
  by its folder name.
- `convergence` is the worst measured value of every metric that has a
  target in the entry's `targets` (`bsh_residual`, `density_change`, and
  `energy_change` where the ground state or an excited root has an energy
  target), over every state the value was computed from, at the final level,
  taken from `state_metrics`, written to three significant figures
  (`state_metrics` and every property value keep full precision). For a value
  of an excited root, that root's own metrics count. **It is a measurement of
  the states, not an error bar on the property.** A metric no state measured is
  left out (a driven state has no `energy_change`), so a value whose states
  were never saved has an empty `convergence`.
- A property that was not requested is absent, not `unavailable`.

Which states a value used follows from what it is:

| Property | Item per | Identifying keys | States its convergence covers |
|---|---|---|---|
| `total_energy`, `orbital_energies`, `energy_components`, `dispersion_energy`, `dipole`, `gradient`, `populations` | (one value) | | the ground state (written once, on the entry) |
| `polarizability` | frequency | `frequency` | the dipole states at that frequency |
| `hyperpolarizability` | frequency pair | `frequency_b`, `frequency_c` | for each component ABC: the A state at ω_B + ω_C, the B and C states, and the second-order source |
| `polarizability_gradient` | laser frequency | `frequency` | as β; components are keyed by nuclear coordinate, then dipole pair |
| `excited_states` | root | `root` | the root; the item holds `excitation_energy`, `oscillator_strength`, `transition_dipole`, `transition_quadrupole` |
| `two_photon_absorption` | root | `root`, `energy`, `frequency` | the root and the three dipole states at `frequency` (half the root's energy) |
| `vibrational_analysis` | (one value) | | the static nuclear-displacement states; `value` holds `hessian`, `frequencies`, `normal_modes`, `reduced_mass` |
| `raman_spectrum` | laser frequency | `frequency` | the Hessian's states and that frequency's polarizability-gradient states; `value` is the list of mode rows |

The ground entry's `energy_components` holds the final `kinetic`,
`nuclear_attraction`, `coulomb`, `exchange_correlation` and
`nuclear_repulsion` terms, and `pcm` when PCM is on. `dispersion_energy` is
present when a dispersion correction is on, `populations` when a population
analysis was requested, `dipole` and `gradient` when computed.

## A state name is its folder name

Every key of `state_metrics`, and every name in an `availability` reason, is a
node folder name:

| Name | Folder |
|---|---|
| `ground_state` | `states/ground_state/` |
| `dipole_z-f0.00000` | `states/response/dipole_z-f0.00000/` |
| `excited-full` | `states/response/excited-full/` |
| `<vbc name>` | `states/response/<vbc name>/` |

## The text summary: `summary.txt`

`<run>/summary.txt` is `results.json` formatted for people. It is written from
the parsed `results.json` alone, every time `results.json` is written, so the
two never disagree; nothing reads it back. `madqc` prints the same text last in
the log when the workflow's `print_level` is 1 or more (the response group's for
`--wf=response`, the `dft` group's for `--wf=scf`); the file is written at every
`print_level`.

- Tables carry values and convergence only, no status column. A missing value
  prints as `-`; an item that has no value still shows each requested component.
- Convergence values and targets are in scientific notation to three
  significant figures. A value above its target is marked `*`; one line at the
  end of the section says `NOT COMPLETE` at the top is the verdict (driven
  states use the scaled rule).
- Values in hartree are printed to 8 decimals, response property values to 6
  decimals, other values to 6 significant figures, response frequencies to 5
  decimals as in the state folder names. A value below 1e-10 in magnitude
  prints as `0`.
- Energy components print in the order kinetic, nuclear attraction, Coulomb,
  exchange-correlation, nuclear repulsion.
- If any property or item is not `complete`, or the run stopped for another
  reason, dropped work or has stalled states, the file opens with a
  `NOT COMPLETE` block: one line per item with its reason from `availability`.
- The Hessian, normal modes, transition quadrupoles, polarizability-gradient
  components and anything else too large for a table are named on one
  `not printed (see results.json)` line.

The summary of the H2 deck `src/examples/qc/response_h2_es_rpa`:

```text
madqc summary
  code      900a79af8-dirty
  launched  2026-11-26 10:50 UTC
  run       1 process x 16 threads

Ground state (moldft, hf), level 1.0e-04 k 6
  reached   bsh_residual 2.16e-05   density_change 2.79e-05   energy_change 2.73e-08
  targets   bsh_residual 5.00e-04   density_change 1.00e-04   energy_change 1.00e-05

  property      value        units
  total_energy  -1.13359355  hartree

  energy_components (hartree)
  term                  value
  kinetic               1.12695922
  nuclear_attraction    -3.63456716
  coulomb               1.31782227
  exchange_correlation  -0.65891220
  nuclear_repulsion     0.71510434

  orbital_energies (hartree)
  spin   orbital  value
  alpha  0        -0.59489421

Response (molresponse), level 1.0e-04 k 6
  targets   driven:  bsh_residual 5.00e-04   density_change 5.00e-04
            excited: bsh_residual 5.00e-04   density_change 5.00e-04   energy_change 5.00e-04

  excited states (energy in hartree)
  root          energy      osc. strength  bsh_residual  density_change  energy_change
  es_root_0000  0.46583341  0.285480       3.85e-05      1.09e-04        1.00e-07

  transition dipole (atomic_units)
  root          x  y  z
  es_root_0000  0  0  -0.958778

  polarizability (bohr^3)
  frequency  component  value     bsh_residual  density_change
  0.00000    zz         6.436397  2.25e-06      3.42e-05

  not printed (see results.json): transition_quadrupole

  run: complete. Nothing dropped. Nothing stalled.
```
