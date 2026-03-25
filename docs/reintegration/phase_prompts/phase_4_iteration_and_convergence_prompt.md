# Phase 4 - Iteration and Convergence Logic

Master roadmap:

- [`../../excited_state_reintegration_gap_analysis.md`](../../excited_state_reintegration_gap_analysis.md)

## Objective

Restore the core legacy excited-state update and convergence behavior inside the Stage 2c `ResponseVector`-based solver for the restricted-shell paths, using the now-explicit Phase 3 fresh-start seed as the baseline, while preserving the Phase 1 identity contract and the Phase 2 restart/archive semantics.

## Scope

In scope:

- restricted-shell Stage 2c iteration and convergence behavior
- legacy-style bundle update flow: bundle potentials, bundle rotation, `theta` assembly, BSH update, residual update, normalization, and any step restrictions needed for parity
- convergence-gate translation and machine-readable per-iteration diagnostics
- `maxsub` / accelerator behavior for the excited-state solver if it is meant to remain a supported input
- focused tests for restricted-shell iteration/convergence behavior and metadata visibility

Out of scope:

- owner-group execution semantics
- downstream property-stage consumers
- widening unrestricted full-restart support unless fully implemented and validated
- unrelated linear-state solver changes
- rewriting the roadmap or redesigning Stage 2c

## Primary Goals

1. Port the restricted-shell legacy update kernel into the current `ResponseVector`-based Stage 2c model closely enough to support parity debugging.
2. Replace the current simplified excited-state convergence gate with an explicit legacy-informed convergence contract and expose the diagnostics needed to validate it.
3. Make `maxsub` behavior explicit: either operational for the excited solver or clearly diagnosed/disabled rather than silently ignored.

## Files To Inspect First

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/apps/molresponse_v2/ResponseVector.hpp`
- `src/apps/molresponse_v2/ResponseRecord.hpp`
- `src/madness/chem/MolresponseLib.hpp`
- `src/apps/molresponse_v2/ResponseSolver.hpp`
- `src/apps/molresponse_v2/ResponseSolver.cpp`
- `src/apps/molresponse/ExcitedResponse.cpp`
- `src/apps/molresponse/ResponseBase.cpp`
- `src/apps/madqc_v2/test_molresponse_excited_metadata_smoke.py`
- `src/apps/madqc_v2/test_molresponse_excited_restart_reuse.py`
- `src/apps/madqc_v2/test_molresponse_excited_protocol_projection.py`
- `src/apps/madqc_v2/test_molresponse_excited_guess_fallback.py`
- `src/apps/madqc_v2/test_molresponse_excited_restart_mode_compat.py`

## Relevant Current Symbols

- `ExcitedProtocolWorkflow::iterate`
- `ExcitedProtocolWorkflow::iterate_typed_bundle_legacy_sequence`
- `ExcitedProtocolWorkflow::compute_bundle_response_potentials`
- `ExcitedProtocolWorkflow::rotate_excited_bundle_states`
- `ExcitedProtocolWorkflow::diagonalize_excited_bundle`
- `ExcitedProtocolWorkflow::iterate_state_from_potentials`
- `ExcitedProtocolWorkflow::iterate_state`
- `ExcitedProtocolWorkflow::build_response_bundle_seed`
- `ExcitedProtocolWorkflow::build_fresh_guess`
- `ExcitedProtocolWorkflow::write_restart_snapshot`
- `ResponseRecord2::record_excited_protocol_result`

Legacy reference symbols:

- `ExcitedResponse::iterate`
- `ExcitedResponse::update_response`
- `ExcitedResponse::bsh_update_excited`
- `ExcitedResponse::rotate_excited_space`
- `ExcitedResponse::excited_eig`
- `ResponseBase::compute_response_potentials`
- `ResponseBase::update_residual`
- `ResponseBase::kain_x_space_update`

## Standing Guardrails

- Treat the master roadmap as authoritative.
- Preserve Stage 2c as the excited-state entry point.
- Do not reintroduce legacy `X_space` as a public or persistent model.
- Keep stable root identity separate from runtime slot order.
- Keep restart behavior explicit and variant-aware.
- Do not silently widen unrestricted restart support.
- Do not rewrite the roadmap during implementation.

## `ResponseVector` Variant Notes

This phase primarily targets restricted-shell iteration and convergence parity.

- `static_restricted`
  - in scope for kernel/convergence translation
  - must remain `full_bundle_resume`
- `static_unrestricted`
  - do not claim feature-complete iteration parity in this phase unless it is actually implemented and validated
  - must remain explicitly `guess_only_unrestricted_variant` otherwise
- `dynamic_restricted`
  - in scope where the restricted-shell kernel is shared or intentionally differentiated
  - must remain `full_bundle_resume`
- `dynamic_unrestricted`
  - do not claim feature-complete iteration parity in this phase unless it is actually implemented and validated
  - must remain explicitly `guess_only_unrestricted_variant` otherwise

## Implementation Requirements

- Reconstruct the legacy restricted-shell update sequence before editing code:
  bundle potentials, generalized subspace diagonalization, bundle rotation, `theta` assembly, BSH application, residual evaluation, normalization, and any accelerator hook.
- Replace the current simplified convergence gate with an explicit convergence contract grounded in the legacy code and the current implementation evidence.
- If `maxsub` remains an exposed excited-state input, make it operational or emit explicit metadata/diagnostics that it is unsupported.
- Preserve the Phase 3 fresh-start seed and the Phase 2 typed restart schema while changing the runtime kernel.
- Keep any new iteration/convergence metadata on the authoritative `ResponseRecord2` path.
- Document any intentional deviations from the legacy iteration sequence in the phase notes.

## Non-Goals

- Do not redesign the Stage 2c orchestration model.
- Do not change the current restart precedence.
- Do not implement property consumers.
- Do not claim full unrestricted restart or numerical support unless the code and tests actually support it.
- Do not jump ahead to owner-group execution semantics unless a tiny diagnostic hook is strictly required.

## Suggested Implementation Steps

1. Audit the current restricted-shell `iterate(...)` path against legacy `ExcitedResponse::iterate(...)` and `update_response(...)`, and write down the exact missing numerical/convergence pieces before editing code.
2. Translate the missing restricted-shell update-kernel steps into the typed bundle helpers in `ExcitedStateBundleSolver.cpp`.
3. Port or intentionally replace the legacy convergence gates, and make the active criteria visible in machine-readable metadata and logs.
4. Make `maxsub` behavior explicit in the excited solver: implement the intended accelerator path or add a clear unsupported diagnostic if that is the deliberate outcome for now.
5. Add or update focused tests that prove the restricted-shell iteration path, convergence metadata, and restart-visible metadata remain coherent.
6. Record the actual implementation choices, support boundaries, and any remaining gaps in Phase 4 notes and roadmap reconciliation.

## Required Deliverables

- code changes limited to the Phase 4 iteration/convergence scope
- tests covering restricted-shell iteration/convergence behavior and metadata visibility
- a Phase 4 contract-consistent notes file after implementation
- a Phase 4 validation record after implementation
- a roadmap reconciliation after implementation

## Validation Requirements

- compare restricted TDA per-iteration energy/residual behavior against the intended legacy flow on a small closed-shell case
- validate metadata parity between `response_metadata.json` and workflow metadata after the updated iteration path runs
- verify restart-visible metadata remains coherent after the kernel/convergence changes
- record the execution environment used, including the interactive job and whether validation used `load_40core.sh` or `load_xeonmax.sh`
- explicitly state which `ResponseVector` variants were validated and which remain outside the supported Phase 4 boundary

## Completion Criteria

- Phase 4 produces a documented restricted-shell iteration/convergence path aligned with the roadmap.
- The solver exposes enough machine-readable iteration data to support parity debugging.
- The tested restricted-shell paths preserve the existing Phase 1 and Phase 2 contracts.
- `maxsub` behavior is no longer ambiguous for the excited-state solver.
