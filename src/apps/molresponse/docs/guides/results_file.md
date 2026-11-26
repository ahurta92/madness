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
  scale).
- `run` appears only on the response entry. `dropped_work` and `stalled` are
  always lists, empty when nothing was dropped.
- `state_metrics` is keyed by node name: the per-state table behind every
  `convergence`. Each node gives the level it was saved at and the measured
  values of each of its subjects. A driven state has one subject (`response`),
  an excited-state set one per root (`es_root_0000`, …), the ground state one
  (`ground`). A second-order source is built rather than iterated, so its
  `subjects` is empty.

## One property

A property with one answer has `value`:

```json
"total_energy": {
  "units": "hartree",
  "value_type": "real",
  "value": -1.1335935503913517,
  "status": "complete",
  "convergence": {"bsh_residual": 2.2e-05, "density_change": 2.8e-05}
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
     "convergence": {"bsh_residual": 2.3e-06, "density_change": 3.4e-05}}
  ]
}
```

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
- `convergence` is the worst measured value of `bsh_residual` and
  `density_change` over every state the value was computed from, at the final
  level, taken from `state_metrics`. For a value of an excited root, that
  root's own metrics count. **It is a measurement of the states, not an error
  bar on the property.** A metric no state measured is left out, so a value
  whose states were never saved has an empty `convergence`.
- A property that was not requested is absent, not `unavailable`.

Which states a value used follows from what it is:

| Property | Item per | Identifying keys | States its convergence covers |
|---|---|---|---|
| `total_energy`, `orbital_energies`, `energy_components`, `dispersion_energy`, `dipole`, `gradient`, `populations` | (one value) | | the ground state |
| `polarizability` | frequency | `frequency` | the dipole states at that frequency |
| `hyperpolarizability` | frequency pair | `frequency_b`, `frequency_c` | for each component ABC: the A state at ω_B + ω_C, the B and C states, and the second-order source |
| `polarizability_gradient` | laser frequency | `frequency` | as β; components are keyed by nuclear coordinate, then dipole pair |
| `excitation_energy`, `oscillator_strength` | root | `root` | the root |
| `transition_dipole`, `transition_quadrupole` | root | `root` | the root |
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
