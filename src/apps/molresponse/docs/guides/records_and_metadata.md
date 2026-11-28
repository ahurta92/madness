# Records: results, states, and restart

Where a response calculation's results and state live, and what the files
guarantee. If you consume MADNESS response output programmatically, this is the
contract.

## The files

| file | scope | audience |
|---|---|---|
| `<run>/results.json` | the whole run folder, one entry per application (`ground_state`, `response`) | every value a person or a script wants, each naming the saved states it came from — see [results_file.md](results_file.md) |
| `<run>/summary.txt` | the same | `results.json` formatted for people; written with it, never read back |
| `<run>/states/…/<node>/` | one saved state (ground state, driven state, excited-state set, second-order source) | the restart: `description.json` (settings, level, stored metrics, updates) and `history.jsonl` (one row per update) |
| `<run>/timings.json` | every launch in the run folder | where the time went, per state, level, iteration and step, including each process group's work |
| `<prefix>.calc_info.json` | the whole workflow — one entry per task | the driver's record: what ran, in what order, whether it failed, and each task's wall time |
| `<prefix>.out` | the whole workflow | humans |

For the response task, `calc_info.json` holds only `type` and the run's times
(`metadata`: `elapsed_time`, `mpi_size`, `nthreads`, `clock_scope`, and
`v3_timing` for plan_build, solve and assemble), which the wall-time lines at
the end of the log print. Every response value is in `results.json`.

## What `results.json` guarantees

- **One level per entry.** Each entry holds the values at the final level of the
  request (`protocol`). The earlier levels of a ladder are in each node's
  `history.jsonl`, not in `results.json`.
- **Values say how far to trust them.** Each value or item has a `status`
  (`complete`, `partial`, `unavailable`), and a `convergence` computed from the
  states it came from, beside the entry's `targets`. Each state in
  `state_metrics` has its own `status` and `updates`.
- **Items are identified by fields, not position.** A polarizability item by its
  `frequency`, a hyperpolarizability item by `frequency_b` and `frequency_c`, an
  excited state or a two-photon item by its `root`. Select by those fields.
- **A failed calculation is recorded.** Its entry has `"status": "failed"` and
  the error.

## Restart

Every state is saved as a node folder under `<run>/states/`. On a later launch
each requested state reads its node: it is reused when its stored metrics meet
the current request at the final level (the same test that decides `status` in
`results.json`), resumed from its saved functions when it does not, and started
from a nearby saved frequency of the same perturbation, a seed, or a guess when
it has no node. A state saved at a coarser level seeds the next level.

- **Rank-agnostic.** Native archives store their own layout and are
  redistributed on read, so a state written on one process count reloads on
  another.
- **One scheduling per level.** Within a launch each state is scheduled at most
  once per level, so no state is solved again and again; work the run cannot
  do (for example a quadratic source on an open-shell ground state) is listed in
  `run.dropped_work` with its reason.

## Atomicity

A node's archive and its description are replaced together in one save;
`results.json`, `summary.txt` and `timings.json` are written to a temporary file
and renamed over the target. On POSIX that rename is atomic within a filesystem,
so a crash during a write leaves either the previous file intact or the new one
complete. These writes happen on rank 0 and are collective: a failure on rank 0
is raised on every rank.

## Consuming it

```python
import json
results = json.load(open("h2o/results.json"))
response = next(c for c in results["calculations"] if c["calculation_id"] == "response")

# alpha_zz at omega = 0, if it was computed from converged states
item = next(v for v in response["properties"]["polarizability"]["values"]
            if abs(v["frequency"]) < 1e-9)
if item["status"] == "complete":
    print(item["components"]["zz"], item["convergence"])
```

Filter on `status` rather than assuming, and select items by their identifying
fields. See [example output](example_output.md) for a real run, and
[`madqc/HDF5_IO.md`](../../../madqc/HDF5_IO.md) for the optional HDF5 backend.
