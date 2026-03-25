# Phase N - <Title>

Master roadmap:

- [`../../excited_state_reintegration_gap_analysis.md`](../../excited_state_reintegration_gap_analysis.md)

## Objective

<State the single phase objective in one paragraph.>

## Scope

In scope:

- <item>
- <item>

Out of scope:

- <item>
- <item>

## Primary Goals

1. <goal>
2. <goal>
3. <goal>

## Files To Inspect First

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/apps/molresponse_v2/ResponseVector.hpp`
- `src/apps/molresponse_v2/ResponseRecord.hpp`
- `src/madness/chem/MolresponseLib.hpp`
- <tests/docs specific to this phase>

## Relevant Current Symbols

- `<symbol>`
- `<symbol>`

## Standing Guardrails

- Treat the master roadmap as authoritative.
- Preserve Stage 2c as the excited-state entry point.
- Do not reintroduce legacy `X_space` as a public or persistent model.
- Keep stable root identity separate from runtime slot order.
- Keep restart behavior explicit and variant-aware.

## `ResponseVector` Variant Notes

State how this phase affects each of:

- `static_restricted`
- `static_unrestricted`
- `dynamic_restricted`
- `dynamic_unrestricted`

If support is not changed for a variant, say so explicitly.

## Implementation Requirements

- <requirement>
- <requirement>
- <requirement>

## Non-Goals

- <non-goal>
- <non-goal>

## Suggested Implementation Steps

1. <step>
2. <step>
3. <step>

## Required Deliverables

- <code change>
- <tests>
- <notes/reconciliation items>

## Validation Requirements

- <test or validation expectation>
- <test or validation expectation>

## Completion Criteria

- <criterion>
- <criterion>
