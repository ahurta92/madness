# MADQC Convergence Evidence and Policy Design

## Goal

Use controlled ground-state and response trajectories, including HDF5 field comparisons, to determine whether response convergence is limited by the ground-state reference, the response solver, or stopping gates before changing convergence behavior across moldft and molresponse.

## Scope and constraints

This phase is analysis-first. It must not modify the Fable or Astra source worktrees, numerical reference files, or application branches. Analysis artifacts may be added to the checkpoint repository and the session folder. Forced runs are diagnostic experiments and must not be treated as calibrated convergence references.

Dalton integration remains in scope for the same release, after the core metric and checkpoint work stabilizes. YAML parser investigation is deferred until the metric and checkpoint schema are stable.

## Evidence matrix

Run a small controlled matrix using H₂, LiH, and water:

- ordinary protocol-driven calculations;
- forced ground-state calculations with a fixed iteration budget at each protocol rung;
- response calculations started from each exhausted ground-state endpoint;
- at least one static response and one frequency-dependent or excited-state response per molecule where supported.

Every run must identify molecule, application, workflow, protocol value, `k`, threshold, frequency/state, iteration phase, sampling phase, and whether it is ordinary, forced, or configured-progress.

## Metrics and machine-readable evidence

Record the following per iteration when available:

- total-density L2 change;
- spin-density L2 change for open-shell calculations;
- per-function BSH residuals and the aggregate gate value;
- energy change;
- kinetic trace change, with the scalar trace defined as `sum_i T_ii`;
- response density and BSH histories;
- response properties and iteration counts.

Application JSON must retain per-function contributions while keeping aggregate gate values explicit. Missing or incomparable measurements must be encoded as unavailable, never as zero.

## HDF5 field comparison

At exhausted ground-state endpoints, preserve native MADNESS HDF5 checkpoints and create structured `.mad.h5` exports with `dump_mra_trees --coeffs` for `pymra`. Compare only matched fields with compatible cell, representation, spin/state identity, protocol metadata, and sampling phase. Use common-grid L2 comparisons for the analysis artifact and label them as sampled comparisons.

The current HDF5 bridge does not make native response archives directly readable by `pymra`; response field comparisons remain optional until an explicit export path exists. Ground-state orbital-derived density comparisons must not be presented as response densities.

## Analysis and decision criteria

For each exhausted ground-state endpoint, compare response iteration count, response residual trajectories, and final properties.

- Ground-state-reference limited: improving the ground-state endpoint consistently reduces response work or changes the response property.
- Response-solver limited: response work and residual behavior remain similar across ground-state endpoints.
- Gate limited: requested metrics continue improving when the run stops, or the gate disagrees with the requested metric.

The analysis notebook must show protocol-rung shading and switch markers, response trajectories grouped by frequency and rung, endpoint metric tables, response work versus ground-state endpoint, and ordinary versus forced full histories.

## Policy questions to settle before implementation

1. Use total-density L2 change as the common charge metric and retain spin-density change separately for open-shell cases.
2. Use the unweighted kinetic trace `sum_i T_ii` as the scalar kinetic metric and retain per-function contributions.
3. Prefer maximum aggregation for gate evaluation unless the evidence supports another policy.
4. Remove or de-emphasize the unexplained quadrature gate.
5. Keep user-requested criteria separate from protocol progress detection.
6. Avoid universal defaults when applications cannot share meaningful values.
7. Treat a stalled requested metric as a reportable condition that can trigger a protocol transition or request tighter convergence.

## Deliverables

- Reproducible run manifest and portable JSON evidence.
- Durable native and structured HDF5 analysis artifacts.
- Updated matplotlib notebook with unavailable-data handling.
- Session analysis note recording measurements, interpretations, and limitations.
- A reviewed metric and gate contract before modifying solver behavior.

## Success criteria

The evidence phase is complete when each selected molecule has matched ground-state endpoints and response trajectories, the notebook can reproduce the endpoint and field comparisons, and the team can classify the dominant limitation for at least one representative closed-shell and one open-shell response case. Only then should the shared convergence implementation plan be executed.
