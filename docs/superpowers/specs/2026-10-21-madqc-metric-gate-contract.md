# MADQC Metric and Gate Contract — Evidence-Phase Proposal

## Status

Reviewable proposal based on the H₂ and LiH HDF5 runs, the existing forced trajectories, and the time-bounded H₂O attempt. This contract does not change solver behavior. It defines the evidence and interfaces required before that change.

## Metric definitions

### Charge and spin density

- `total_density_l2` is the application’s L2 change in total charge density.
- `spin_density_l2` is recorded separately when alpha and beta densities can vary independently.
- Open-shell calculations carry both metrics when available; the spin metric is not folded into the total metric.
- A missing metric is `{ "available": false, "value": null }` and is never encoded as numeric zero.

### BSH residual

- Each orbital or response function retains its raw L2 residual.
- The aggregate gate value is recorded separately from the per-function list.
- The proposed aggregate for gate evaluation is the maximum absolute per-function residual.

### Energy

- `energy_change` is the scalar change in the application energy between consecutive samples.
- Its units and sampling phase are recorded with the value.

### Kinetic energy

- The scalar orbital/response kinetic metric is the trace `sum_i T_ii`.
- Per-function `kinetic_energy_i` contributions remain in the application record.
- Occupation-weighted physical kinetic energy may be reported as a separate diagnostic, but it is not substituted for the trace metric.

## Gate and progress policy

1. User-requested metric criteria and protocol-progress detection are separate decisions.
2. The solver reports requested metric values and whether each is improving, stalled, satisfied, or unavailable.
3. A configured protocol may advance when requested metrics stall or fail to meet the current rung, subject to the application’s explicit continuation policy.
4. No universal numeric defaults are introduced in this phase. Tests and decks must provide criteria explicitly when numerical equivalence is required.
5. Maximum aggregation is the proposed default for collections of function metrics because it preserves the worst member. RMS/minimum alternatives remain reporting options, not implicit gate substitutions.
6. The quadrature gate is not part of the shared metric contract unless a future experiment demonstrates a documented, application-independent interpretation.
7. Forced runs are never used to calibrate universal stopping thresholds.

## HDF5 and field evidence

Native application HDF5 archives are checkpoint containers. Structured `.mad.h5` exports created with `dump_mra_trees --coeffs` are the current `pymra` input. A field comparison is valid only when cell, representation, spin/state identity, endpoint, and sampling phase are compatible. Ground-state orbital-derived density proxies are labeled separately from response-state densities.

## Current evidence classification

- H₂ ground-state field differences are measurable and concentrated in the molecular core.
- H₂ and LiH demonstrate that the HDF5 response backend can save and reload response state archives.
- The current evidence does not yet classify response convergence as reference-limited, response-solver-limited, or gate-limited because the same response observable has not yet been run from two fully exhausted ground-state endpoints for both a closed-shell and an open-shell case.
- H₂O remains a partial, time-bounded run and is not a convergence reference.

## Required acceptance tests before implementation

- Schema tests cover unavailable values, response-to-ground-state endpoint linkage, ordinary/forced/configured-progress classification, total versus spin density, per-function residuals, and kinetic trace contributions.
- Regression decks explicitly set criteria needed to reproduce numerical references.
- A matched endpoint experiment shows the response iteration count, residual trajectory, and final property for at least one closed-shell and one open-shell observable.
- The implementation preserves native HDF5 checkpoint compatibility and keeps structured export metadata explicit.

## Approval boundary

This proposal must be reviewed before changing moldft or molresponse stopping behavior. Review approval authorizes a separate implementation plan; it does not itself modify solver gates.
