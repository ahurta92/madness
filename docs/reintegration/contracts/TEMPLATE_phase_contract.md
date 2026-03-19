# Phase N Contract - <Title>

Master roadmap:

- [`../../excited_state_reintegration_gap_analysis.md`](../../excited_state_reintegration_gap_analysis.md)

## Scope Contract

This phase is allowed to change:

- <allowed area>
- <allowed area>

This phase must not change:

- <forbidden area>
- <forbidden area>

## Architectural Invariants

- Preserve Stage 2c as the excited-state entry point.
- Do not reintroduce public/persistent legacy `X_space`.
- Keep stable root identity separate from slot order.
- Keep metadata authoritative through `ResponseRecord2`.
- Keep restart semantics explicit by `ResponseVector` variant.

## Variant Contract

For each variant, state what this phase may change.

| Variant | Allowed change in this phase | Must remain true |
| --- | --- | --- |
| `static_restricted` | <text> | <text> |
| `static_unrestricted` | <text> | <text> |
| `dynamic_restricted` | <text> | <text> |
| `dynamic_unrestricted` | <text> | <text> |

## Persistence / Metadata Contract

- <metadata rule>
- <restart rule>

## Validation Contract

Required validation before declaring the phase complete:

- <validation item>
- <validation item>

## Exit Criteria

- <criterion>
- <criterion>
