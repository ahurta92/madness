# molresponse

A unified MADNESS solver for **molecular response properties** — linear and
non-linear — built on multiresolution analysis (MRA). It is the response engine
behind `madqc --wf=response` (the default and only engine), designed to scale to
larger systems and to export research-quality, MRA-fidelity data. Because it is
basis-set-free, its properties converge in resolution toward the complete-basis
limit rather than being extrapolated in a Gaussian basis.

## What it computes

| Property | Status | Guide |
|---|---|---|
| **Polarizability α** (static + dynamic) | ✅ supported (closed shell) | [guides/polarizability.md](docs/guides/polarizability.md) |
| **First hyperpolarizability β** (static + SHG) | ✅ supported (closed-shell HF) | [guides/hyperpolarizability.md](docs/guides/hyperpolarizability.md) |
| **Vibrational Raman** (components and full normal-mode analysis) | ✅ supported (closed-shell HF) | [guides/raman.md](docs/guides/raman.md) |
| **Excited states** (TDA / Full RPA) | ✅ supported (closed shell) | [guides/excited_states.md](docs/guides/excited_states.md) |
| **Two-photon absorption** | 🟡 preliminary (closed-shell singlets) | [guides/two_photon_absorption.md](docs/guides/two_photon_absorption.md) |
| **DALTON warm-start** (seed ES/SCF from DALTON) | ✅ supported | [guides/dalton_warm_start.md](docs/guides/dalton_warm_start.md) |
| Resonance Raman · open-shell property assembly and excited states | 🔬 future | — |

**Scope:** property assembly and excited-state execution are currently closed-shell.
Open-shell driven dipole states can be solved and saved, but their property
assembly is skipped. Quadratic DFT response is rejected because its second XC
kernel is not implemented. Full Raman analysis requires the complete nuclear
and derivative state set; partial selections return derivative components.

## Documentation

- **Feature guides** — one page per property, in [`docs/guides/`](docs/guides/):
  [polarizability](docs/guides/polarizability.md) ·
  [hyperpolarizability](docs/guides/hyperpolarizability.md) ·
  [raman](docs/guides/raman.md) ·
  [excited states](docs/guides/excited_states.md) ·
  [two-photon absorption](docs/guides/two_photon_absorption.md) ·
  [DALTON warm-start](docs/guides/dalton_warm_start.md).
- **[`docs/guides/example_output.md`](docs/guides/example_output.md)** — what a run
  prints and writes (real output).
- **[`docs/guides/records_and_metadata.md`](docs/guides/records_and_metadata.md)** —
  the results/metadata schema, restart, and provenance: the interface to script against.
- **[`docs/guides/formalism.md`](docs/guides/formalism.md)** — the shared response
  formalism the guides specialize.
- **[`docs/guides/parallelism.md`](docs/guides/parallelism.md)** — subworld /
  state-parallel design.
- **[`docs/FEATURE_STATUS.md`](docs/FEATURE_STATUS.md)** — the capability matrix.
- **[`../madqc/RESPONSE_PROPERTIES.md`](../madqc/RESPONSE_PROPERTIES.md)** — run
  recipes and the deck reference; **[`../madqc/HDF5_IO.md`](../madqc/HDF5_IO.md)** —
  optional HDF5 I/O.

## Using it

**Through madqc (recommended).** Put a `response` block in a madqc deck and run
the response workflow (v3 is the default — no engine line needed):

```text
response
    dipole              true
    dipole.directions   z
    dipole.frequencies  [0.0, 0.04]
    quadratic           true        # also compute β
end
```
```bash
madqc --wf=response h2o.in
```

Outputs land in `<prefix>.calc_info.json` (machine-readable) and `<prefix>.out`
(human-readable), with the property tensors in `response_metadata.json`. See
[`../madqc/RESPONSE_PROPERTIES.md`](../madqc/RESPONSE_PROPERTIES.md) for
per-property knobs and the result schema.

**Standalone binary** (development / low-level testing) — drives the response seam
directly on a ground-state archive:

```bash
molresponse --archive=<path_to_moldft.restartdata>
```

## Relationship to the other response codes

`molresponse_v2` (`MolresponseLib.hpp`) was this engine's predecessor and parity
reference during development; **it has since been removed** — this engine is now
the sole response engine and the default for `madqc --wf=response`, and a deck that
still carries `engine v2` errors out. Migrating a v2 deck is a one-line edit:
delete the `engine v2` line. See [`MIGRATION_FROM_V2.md`](MIGRATION_FROM_V2.md).

<!-- EGM architecture contract: begin -->

## Architecture and build

The public entry point is `madqc --wf=response`; the standalone `molresponse`
main and executable are retired. The response engine remains in this directory,
along with its private numerical tests and import/analysis tools.

[orchestrator/response_workflow.hpp](orchestrator/response_workflow.hpp) defines
the collective engine input/output seam. [madqc_adapter.hpp](madqc_adapter.hpp)
constructs requests and loads the exact qualified ground-state checkpoint.
Property/state dependencies stay in the response planner and calculation manager;
SCF-to-response dependencies stay in the madqc workflow.

Response runs use the shared lifecycle for driven states, coupled TDA/RPA roots,
VBC and the supported alpha, beta, Raman and 2PA properties. Applications own
functions and comparison history; shared controllers consume their measurements.
The existing DAG owns response dependencies. Native checkpoints record exact
sources, and properties load only inputs that meet the current request.

One response convergence block applies to driven and excited states. The default
criteria are maximum per-function BSH residual, density change and energy change
(where applicable). Provisional absolute defaults are `5e-4`, with relative
terms zero. Set an explicit tolerance when selecting kinetic change. Final
protocol requirements are mandatory; intermediate protocols retain the legacy
threshold floors. These defaults still need the later full accuracy calibration.
`dconv`, `excited.guess_max_iter`, separate excited-state budgets and experimental stall/seed-stage controls
are not public inputs. Restart uses saved functions and fresh KAIN/history; exact
saved evidence can satisfy a tighter request without another evaluation.

Results contain grouped properties and state metrics once; analysis and shared
steady-clock timing records remain separate. Open-shell HF driven states have
native metrics/restart, but their property contractions and excited-state theory
remain unavailable. Quadratic XC and meta-GGA response limitations are unchanged.
Full Raman and 2PA retain their documented preliminary scientific status.

Build `madqc` for application runs. CMake continues to build the response
private test targets and analysis/import tools, including `dump_mra_trees`,
`seed_from_dalton`, `seed_moldft_from_dalton` and `tpa_from_dalton`. Existing
optional HDF5 utilities remain supported; shared HDF5 restart is deferred.

Run `ctest -R madness/test/molresponse` for engine tests and the `response_*`
[qctest cases](../../examples/qc/README.md) for scientific workflow checks.
