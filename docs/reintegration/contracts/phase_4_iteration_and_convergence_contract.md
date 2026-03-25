# Phase 4 Contract - Iteration and Convergence Logic

Master roadmap:

- [`../../excited_state_reintegration_gap_analysis.md`](../../excited_state_reintegration_gap_analysis.md)

## Scope Contract

This phase is allowed to change:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp` if a small solver-local contract adjustment is required for iteration or convergence metadata
- `src/apps/molresponse_v2/ResponseRecord.hpp` if per-iteration or convergence metadata needs to be normalized through the existing authoritative write path
- `src/madness/chem/MolresponseLib.hpp` only if Stage 2c needs a small wiring change for iteration diagnostics or convergence result propagation
- focused excited-state Stage 2c tests under `src/apps/madqc_v2/`
- documentation or notes needed to record Phase 4 outcomes

This phase must not change:

- the public Stage 2c orchestration model in a way that redefines later phases
- the stable root identity contract (`root_id`, `stable_index`, `slot_permutation`)
- the Phase 2 restart precedence or the current explicit restart support boundary
- owner-group execution semantics beyond diagnostics needed to observe the current behavior
- downstream property consumers or unrelated linear-state solver behavior

## Architectural Invariants

- Preserve Stage 2c as the excited-state entry point.
- Do not reintroduce public/persistent legacy `X_space`.
- Keep stable root identity separate from slot order.
- Keep metadata authoritative through `ResponseRecord2`.
- Keep restart semantics explicit by `ResponseVector` variant.
- Preserve the current restart precedence:
  `current protocol snapshot -> nearest lower protocol snapshot -> generic guess archive -> carryover -> fresh guess`

## Variant Contract

For Phase 4, the goal is restricted-shell numerical-kernel and convergence parity. This phase is not allowed to silently widen restart or feature support for unrestricted variants.

| Variant | Allowed change in this phase | Must remain true |
| --- | --- | --- |
| `static_restricted` | Legacy-style bundle update, convergence logic, and `maxsub` / accelerator behavior may be implemented or clarified | `full_bundle_resume` restart support remains intact |
| `static_unrestricted` | Diagnostics or explicit guardrails may be improved if touched indirectly | Restart support must remain explicitly `guess_only_unrestricted_variant` unless fully implemented and validated |
| `dynamic_restricted` | Legacy-style restricted full-response iteration and convergence logic may be improved | `full_bundle_resume` restart support remains intact |
| `dynamic_unrestricted` | Diagnostics or explicit guardrails may be improved if touched indirectly | Restart support must remain explicitly `guess_only_unrestricted_variant` unless fully implemented and validated |

## Persistence / Metadata Contract

- Phase 4 may extend per-iteration or convergence metadata, but it must preserve the current protocol-result contract for `roots`, `state_names`, `slot_permutation`, `response_variant`, `restart_support_mode`, `snapshot_kind`, `bundle_state_present`, and `restart_capable`.
- Any new iteration or convergence diagnostics must flow through the existing `ResponseRecord2`-based authoritative write path rather than ad hoc JSON mutation.
- If `maxsub` or an accelerator becomes live, restart artifacts and metadata must make the active behavior explicit rather than silently changing solver semantics.

## Validation Contract

Required validation before declaring Phase 4 complete:

- at least one restricted TDA case showing the translated bundle update and convergence path is active
- at least one restricted dynamic case if the shared restricted kernel changes affect that path
- metadata parity check between workflow metadata and `response_metadata.json`
- explicit statement of variant coverage and any variants intentionally left unsupported
- validation under the interactive 40-core job using the established environment workflow

## Exit Criteria

- Restricted-shell Stage 2c iteration and convergence behavior are documented and aligned with the reconciled roadmap.
- The solver exposes machine-readable convergence diagnostics sufficient for parity debugging.
- `maxsub` behavior is either operational in the excited solver or explicitly diagnosed as unsupported.
- No unrestricted restart/support boundary is widened implicitly.
