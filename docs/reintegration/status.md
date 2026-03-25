# Reintegration Status Dashboard

Master roadmap:

- [`../excited_state_reintegration_gap_analysis.md`](../excited_state_reintegration_gap_analysis.md)

Last reconciled:

- 2026-03-16

## Overall Status

Current phase focus:

- Phase 4 — Iteration and Convergence Logic
- Phase 4 contract and prompt are instantiated under `docs/reintegration/`

Last completed phase:

- Phase 3 — Solver Initialization Path

Current state:

- Root identity and metadata persistence foundation are complete.
- Restart/archive semantics are reliable for restricted variants.
- Fresh-start initialization is now explicit, reproducible for the covered cases, and consistent with guess-archive reuse.
- Phase 4 implementation is in progress:
  restricted-shell dual-gate convergence metadata, per-root accelerator scaffolding, and
  explicit step restriction landed in code and now compile under the `xeonmax` sandbox build.
- A real MPI baseline now exists in `analysis/results/h2_excited/v4` for the dynamic restricted
  path. That run shows Stage 2c, restart sidecars, and metadata persistence are live, but
  convergence is still blocked by repeated post-iteration-1 bundle-rotation failure
  (`sygvp_exception`) plus a step-restriction policy that currently differs from legacy.
- Remaining work is primarily restricted-shell numerical parity, unrestricted completion,
  scheduling integration, and downstream consumers.

## Phase Dashboard

| Phase | Name | Status | Notes |
| --- | --- | --- | --- |
| 1 | Foundational Data Structures | Completed | Stable root descriptors, slot permutation tracking, unified metadata write path |
| 2 | Archive and Restart System | Completed | Typed restart schema, bundle archive round-trip, same/lower/guess restart flow |
| 3 | Solver Initialization Path | Completed | Explicit sort/truncate/seed flow, `guess_bundle` fresh-start contract, and guess-archive consistency coverage landed |
| 4 | Iteration and Convergence Logic | In Progress | Dual-gate convergence metadata and restricted `maxsub` scaffolding landed; `v4` shows live runtime path but unresolved rotation/convergence parity |
| 5 | Naming and Protocol-State Mapping | Planned | Downstream naming policy and protocol-edge-case hardening remain open |
| 6 | Property Integration | Planned | No excited-root property consumer yet |
| 7 | Verification and Cleanup | Planned | Numerical parity harness and cleanup remain open |

## Current Support Boundary

| Area | Status |
| --- | --- |
| Stage 2c metadata parity | In place |
| Stable root identity | In place |
| Restricted full restart | In place |
| Restricted convergence diagnostics | In place |
| Restricted `maxsub` scaffolding | In place, runtime validation pending |
| Unrestricted full restart | Not implemented |
| Owner-group execution semantics | Not implemented |
| Property-stage excited-root consumption | Not implemented |
| Numerical parity against legacy | Partially translated, not yet runtime-validated |

## Immediate Follow-Up

Before starting the next implementation phase:

1. Use the checked-in Phase 4 contract and prompt as the starting point.
2. Use `analysis/results/h2_excited/v4` as the current runtime baseline and fix the restricted
   dynamic numerical blockers, starting with repeated bundle-diagonalization failure after
   iteration 1 and the current step-restriction behavior.
3. Compare the corrected restricted-shell iteration/convergence path against legacy expectations in
   an MPI-capable environment.
4. Reconcile Phase 4 as completed only after runtime validation, then update notes, validation,
   roadmap reconciliation, and active contracts.
