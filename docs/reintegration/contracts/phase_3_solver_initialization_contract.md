# Phase 3 Contract - Solver Initialization Path

Master roadmap:

- [`../../excited_state_reintegration_gap_analysis.md`](../../excited_state_reintegration_gap_analysis.md)

## Scope Contract

This phase is allowed to change:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp` if a small solver-local metadata or helper adjustment is required for initialization behavior
- focused excited-state Stage 2c tests under `src/apps/madqc_v2/`
- documentation or notes needed to record Phase 3 outcomes

This phase must not change:

- the public Stage 2c orchestration model in a way that redefines later phases
- the stable root identity contract (`root_id`, `stable_index`, `slot_permutation`)
- restart precedence or the Phase 2 restart support boundary
- property consumers or downstream property APIs
- unrelated linear-state solver behavior

## Architectural Invariants

- Preserve Stage 2c as the excited-state entry point.
- Do not reintroduce public/persistent legacy `X_space`.
- Keep stable root identity separate from slot order.
- Keep metadata authoritative through `ResponseRecord2`.
- Keep restart semantics explicit by `ResponseVector` variant.
- Preserve the current restart precedence:
  `current protocol snapshot -> nearest lower protocol snapshot -> generic guess archive -> carryover -> fresh guess`

## Variant Contract

For Phase 3, the goal is initialization-path parity and fresh-start consistency. This phase is not allowed to silently widen restart support.

| Variant | Allowed change in this phase | Must remain true |
| --- | --- | --- |
| `static_restricted` | Fresh-start initialization and guess refinement may be improved | `full_bundle_resume` restart support remains intact |
| `static_unrestricted` | Initialization heuristics may be clarified or improved if needed | Restart support must remain explicitly `guess_only_unrestricted_variant` unless fully implemented and validated |
| `dynamic_restricted` | Fresh-start initialization and bundle seeding may be improved | `full_bundle_resume` restart support remains intact |
| `dynamic_unrestricted` | Initialization heuristics may be clarified or improved if needed | Restart support must remain explicitly `guess_only_unrestricted_variant` unless fully implemented and validated |

## Persistence / Metadata Contract

- Phase 3 may change how fresh guesses are built or truncated, but it must preserve the current metadata contract for `roots`, `state_names`, `slot_permutation`, `response_variant`, `restart_support_mode`, `snapshot_kind`, `bundle_state_present`, and `restart_capable`.
- If initialization changes affect guess archives or protocol snapshots, they must remain compatible with the Phase 2 typed restart schema.
- Fresh-start logic must remain consistent with the snapshot-kind distinction:
  `guess_trial_seed`, `guess_bundle`, and `protocol_bundle`.

## Validation Contract

Required validation before declaring Phase 3 complete:

- at least one fresh-start restricted case that exercises the refined initialization path
- at least one restart-from-guess case showing name/order stability relative to the fresh-start result
- explicit statement of variant coverage
- metadata parity check between workflow metadata and `response_metadata.json`

## Exit Criteria

- The fresh-start path in `molresponse_v2` is documented and behaves consistently with the reconciled roadmap.
- Initial root ordering and naming are reproducible for the tested cases.
- No Phase 2 restart support boundary is widened implicitly.
