# Reintegration Context

This file captures the standing architectural context for the excited-state reintegration work in `molresponse_v2`.

Master roadmap:

- [`../excited_state_reintegration_gap_analysis.md`](../excited_state_reintegration_gap_analysis.md)

Current implementation reports:

- [`../legacy_excited_state_report.md`](../legacy_excited_state_report.md)
- [`../current_excited_state_implementation_report.md`](../current_excited_state_implementation_report.md)

## Current State

As of the latest roadmap reconciliation:

- Phase 1 is complete.
- Phase 2 is complete.
- The next planned implementation phase is Phase 3.

What is already in place:

- Stage 2c is the excited-state entry point in `molresponse_v2`.
- Stable root identity is represented by `ExcitedRootDescriptor`.
- Root ordering changes are tracked with `slot_permutation`.
- Stage 2c metadata writes are unified through `ResponseRecord2::record_excited_protocol_result(...)`.
- Final excited metadata is written both to workflow metadata and to `response_metadata.json`.
- Restart snapshots are explicit and variant-aware.
- Restricted variants support full-bundle resume.
- Unrestricted variants remain explicitly guess-only for restart unless and until that support is implemented and validated.

## Architectural Guardrails

These guardrails apply to every remaining phase unless the roadmap is explicitly reconciled to change them.

1. Keep Stage 2c as the public excited-state workflow entry point.
   - Do not move excited solving back to a standalone legacy-style driver model.

2. Do not reintroduce legacy `X_space` as a public or persistent architecture.
   - Bundle-coupled math can stay local to `ExcitedStateBundleSolver`.
   - External interfaces must remain manifest- and descriptor-driven.

3. Treat stable root identity as separate from runtime slot order.
   - `root_id` and `stable_index` are identity.
   - `slot_index` is transient runtime placement.
   - `display_name` is presentation only.

4. Keep metadata authoritative through the persistence layer.
   - Stage 2c writes should flow through `ResponseRecord2`.
   - `response_metadata.json` and workflow metadata must agree on the excited subtree.

5. Keep restart semantics explicit and variant-aware.
   - Do not silently claim equivalent restart support for all variants.
   - Unsupported restart modes must be marked in code and metadata.

6. Preserve the restart precedence order.
   - `current protocol snapshot -> nearest lower protocol snapshot -> generic guess archive -> carryover -> fresh guess`

7. Do not expose solver-internal bundle slots as the public contract.
   - Property consumers, manifests, and future APIs must use stable root descriptors.

8. Preserve the current workflow layering.
   - linear states first
   - excited Stage 2c second
   - derived/property stages after that

## `ResponseVector` Variant Guidance

Treat these four variants as the true typed solver-state categories.

| Variant | Meaning | Practical interpretation | Current restart status |
| --- | --- | --- | --- |
| `static_restricted` | restricted-shell, static response | TDA-like restricted excited-state solving | `full_bundle_resume` |
| `static_unrestricted` | unrestricted-shell, static response | TDA-like unrestricted excited-state solving | `guess_only_unrestricted_variant` |
| `dynamic_restricted` | restricted-shell, full response | RPA-like restricted excited-state solving | `full_bundle_resume` |
| `dynamic_unrestricted` | unrestricted-shell, full response | RPA-like unrestricted excited-state solving | `guess_only_unrestricted_variant` |

Variant rules:

- Do not assume `x`-only persistence is sufficient for dynamic variants.
- Do not silently reconstruct missing channels as zero if that changes restart semantics.
- If a phase changes support boundaries, reconcile both the roadmap and `contracts/active_contracts.md`.

## Practical Repository Focus

The main files for remaining reintegration work are:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/apps/molresponse_v2/ResponseVector.hpp`
- `src/apps/molresponse_v2/ResponseRecord.hpp`
- `src/madness/chem/MolresponseLib.hpp`
- `src/apps/madqc_v2/test_molresponse_excited_*.py`

## Validation Environment Notes

Use repository-practical validation rather than generic assumptions.

- Use the interactive Slurm allocation for multi-rank / 40-core validation when required.
- In the current environment, `load_40core.sh` is the reliable path for MPI-enabled excited-state test execution.
- `load_xeonmax.sh` is still useful for serial smoke validation, but it should not be assumed to provide a usable MPI launcher in every allocation.

If the validation environment changes, record that change in phase notes and validation output.
