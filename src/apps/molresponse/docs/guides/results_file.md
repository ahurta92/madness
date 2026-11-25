# The results file: `results.json`

`<run>/results.json` holds every value a person wants from a run folder. Facts
about one saved state are in that state's node folder under `states/`; time and
memory are in `timings.json`. Each application replaces its own entry when it
finishes and leaves the others alone. A file of another version is replaced, not
merged.

## Top level

```json
{
  "version": 3,
  "provenance": {
    "code_version": "<git commit, as node descriptions store it>",
    "launched_at": "2026-11-25T09:12:03.512Z",
    "mpi_processes": 1,
    "threads_per_process": 16,
    "archive_io": "native"
  },
  "calculations": [ <entry>, ... ]
}
```

`provenance` describes the launch that wrote the file last. `threads_per_process`
counts the pool threads and the main thread, as `timings.json` and
`calc_info.json` do.
`archive_io` is `native` or `hdf5`.

## One entry

```json
{
  "calculation_id": "response",
  "application": "molresponse",
  "protocol": {"threshold": 0.0001, "k": 6},
  "properties": { "<name>": <property>, ... },
  "run": {"stop_reason": "complete", "dropped_work": [], "stalled": []},
  "state_metrics": {
    "dipole_z-f0.00000": {
      "protocol": {"threshold": 0.0001, "k": 6},
      "subjects": {"response": {"values": {"bsh_residual": 2.3e-06}, "availability": {}}}
    }
  }
}
```

- `calculation_id` is `ground_state` (written by `moldft`) or `response`.
- `protocol` is the final level of the request.
- `run` appears only on the response entry. `dropped_work` and `stalled` are
  always lists, empty when nothing was dropped.
- `state_metrics` is keyed by node name. Each node gives the level it was saved
  at and the measured values of each of its subjects. A driven state has one
  subject (`response`), an excited-state set one per root (`es_root_0000`, …),
  the ground state one (`ground`). A second-order source is built rather than
  iterated, so its `subjects` is empty.

## One property

```json
"polarizability": {
  "units": "bohr^3",
  "value_type": "real",
  "values": [
    {"frequency": 0.0, "status": "complete",
     "components": {"zz": 6.4364}, "availability": {},
     "states": ["dipole_z-f0.00000"]}
  ]
}
```

A property holds the value at the final protocol level only. Earlier levels are
in the node histories and the log. Its `status` (`complete`, `partial` or
`unavailable`) is decided when the file is written. A missing component is
listed under `availability` with the reason. A property that was not requested
is absent, not `unavailable`.

Where a property has one value per frequency or per root, each item in `values`
has its own `states`. Otherwise `states` is on the property.

## A state name is its folder name

Every name in a `states` list, and every key of `state_metrics`, is a node
folder name:

| Name | Folder |
|---|---|
| `ground_state` | `states/ground_state/` |
| `dipole_z-f0.00000` | `states/response/dipole_z-f0.00000/` |
| `excited-full` | `states/response/excited-full/` |
| `<vbc name>` | `states/response/<vbc name>/` |

A `states` list names the nodes the value needs. When the value is complete,
every name in it is a folder that exists. When it is `unavailable`, a state it
needs may not have been saved.
