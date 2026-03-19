# Phase 3 - Solver Initialization Path

Master roadmap:

- [`../../excited_state_reintegration_gap_analysis.md`](../../excited_state_reintegration_gap_analysis.md)

## Objective

Bring the Stage 2c fresh-start initialization path in `molresponse_v2` closer to the legacy excited-state solver by making the initial trial-bundle construction, refinement, truncation, and root-order selection behavior explicit and reproducible, while preserving the Phase 1 identity contract and the Phase 2 restart/archive semantics.

## Scope

In scope:

- fresh-start excited-state initialization behavior in Stage 2c
- trial-bundle sizing, projection, normalization, refinement, truncation, and initial ordering behavior
- consistency between fresh-start initialization and the existing guess/archive snapshot contract
- focused tests for initialization consistency and restart-from-guess behavior

Out of scope:

- full restricted-shell numerical parity beyond initialization
- convergence-gate redesign
- KAIN / `maxsub` runtime activation
- owner-group execution semantics
- property-stage consumers
- widening unrestricted full-restart support unless fully implemented and validated

## Primary Goals

1. Port the legacy initialization sequence into the current `ResponseVector`-based Stage 2c model as faithfully as practical.
2. Make the fresh-start path and the guess-archive path produce consistent initial root ordering and naming for the supported cases.
3. Keep the current restart schema and stable-root metadata contract intact while improving initialization behavior.

## Files To Inspect First

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/apps/molresponse_v2/ResponseVector.hpp`
- `src/apps/molresponse_v2/ResponseRecord.hpp`
- `src/madness/chem/MolresponseLib.hpp`
- `src/apps/molresponse/ExcitedResponse.cpp`
- `src/apps/molresponse/ResponseBase.cpp`
- `src/apps/madqc_v2/test_molresponse_excited_metadata_smoke.py`
- `src/apps/madqc_v2/test_molresponse_excited_restart_reuse.py`
- `src/apps/madqc_v2/test_molresponse_excited_guess_fallback.py`

## Relevant Current Symbols

- `ExcitedProtocolWorkflow::build_fresh_guess`
- `ExcitedProtocolWorkflow::iterate_trial`
- `ExcitedProtocolWorkflow::initialize_protocol_guess`
- `ExcitedProtocolWorkflow::ensure_trial_space_matches_guess`
- `ExcitedProtocolWorkflow::build_response_bundle_seed`
- `assign_excited_state_names`
- `ExcitedTrialSpace`
- `ExcitedRootDescriptor`
- `ExcitedBundleProtocolResult`

Legacy reference symbols:

- `ExcitedResponse::initialize`
- `ExcitedResponse::iterate_trial`
- `ResponseBase::orthonormalize_fock`
- `ResponseBase::sort`

## Standing Guardrails

- Treat the master roadmap as authoritative.
- Preserve Stage 2c as the excited-state entry point.
- Do not reintroduce legacy `X_space` as a public or persistent model.
- Keep stable root identity separate from runtime slot order.
- Keep restart behavior explicit and variant-aware.
- Do not silently widen unrestricted restart support.

## `ResponseVector` Variant Notes

This phase primarily targets the initialization path, not restart-support expansion.

- `static_restricted`
  - in scope for initialization refinement
  - must remain `full_bundle_resume`
- `static_unrestricted`
  - initialization logic may be improved if needed
  - must remain explicitly `guess_only_unrestricted_variant` unless fully implemented and validated
- `dynamic_restricted`
  - in scope for initialization refinement
  - must remain `full_bundle_resume`
- `dynamic_unrestricted`
  - initialization logic may be improved if needed
  - must remain explicitly `guess_only_unrestricted_variant` unless fully implemented and validated

## Implementation Requirements

- Reconstruct the legacy initialization sequence before editing code:
  oversized trial bundle, occupied-space projection, Gram-Schmidt/normalization, trial refinement, energy sorting, and lowest-root selection.
- Keep the fresh-start path consistent with the Phase 2 snapshot-kind model:
  `guess_trial_seed`, `guess_bundle`, and `protocol_bundle`.
- Preserve stable root descriptors and `slot_permutation`.
- If initialization changes affect guess archive contents, keep them restart-compatible with the existing typed snapshot schema.
- Document any intentional deviations from the legacy initialization flow in the phase notes.

## Non-Goals

- Do not redesign the Stage 2c orchestration model.
- Do not change the Phase 2 restart precedence.
- Do not implement property consumers.
- Do not claim full unrestricted restart support unless the code and tests actually support it.
- Do not jump ahead to full convergence/KAIN parity unless a small change is strictly required by initialization logic.

## Suggested Implementation Steps

1. Audit the current fresh-start path in `build_fresh_guess(...)` and `iterate_trial(...)` and compare it directly with legacy `ExcitedResponse::initialize(...)` / `iterate_trial(...)`.
2. Identify which legacy initialization steps are already present and which are currently approximated or missing.
3. Adjust the new initialization path so the trial bundle sizing, refinement, truncation, and initial root ordering are explicit and reproducible.
4. Reconcile fresh-start output with guess-archive reuse so tested no-op guess reuse does not silently reshuffle names/order.
5. Add or update focused tests for fresh-start consistency and guess-archive consistency.
6. Record the actual implementation decisions in Phase 3 notes.

## Required Deliverables

- code changes limited to the Phase 3 initialization scope
- tests covering fresh-start initialization behavior and restart-from-guess consistency
- a filled Phase 3 notes file
- a Phase 3 validation record
- a roadmap reconciliation after implementation

## Validation Requirements

- compare fresh-start ordering/selection behavior against the intended legacy flow on a small closed-shell case
- confirm a fresh run and a guess-archive-seeded run preserve stable root ids and deterministic names for the tested case
- confirm `response_metadata.json` and workflow metadata still agree
- record the execution environment used, including whether validation was done under `load_40core.sh` or `load_xeonmax.sh`

## Completion Criteria

- Phase 3 produces a documented, reproducible fresh-start initialization path aligned with the roadmap.
- The tested fresh-start and guess-archive paths produce consistent root identity and naming behavior.
- The implementation does not weaken the existing Phase 1 or Phase 2 contracts.
