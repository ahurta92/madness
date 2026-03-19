# Generate Next Phase Prompt

Use this guide to generate the next implementation prompt from the reconciled state.

Master roadmap:

- [`../excited_state_reintegration_gap_analysis.md`](../excited_state_reintegration_gap_analysis.md)

Supporting inputs:

- [`status.md`](status.md)
- [`contracts/active_contracts.md`](contracts/active_contracts.md)
- latest phase notes
- latest validation record

## Prompt Generation Rules

1. Use the master roadmap as the source of truth.
2. Scope the prompt to one roadmap phase unless there is a deliberate, documented split.
3. Preserve existing architectural guardrails.
4. State support boundaries explicitly, especially for `ResponseVector` variants and restart behavior.
5. Separate primary goals from non-goals.
6. Require code references, tests, and reconciliation notes in the output.

## Minimum Prompt Inputs

Before writing the prompt, extract:

- the next planned phase name and purpose
- the unresolved gaps that phase is meant to address
- dependencies already completed
- active contracts and support boundaries
- validation expectations for that phase

## Prompt Structure

Use this structure:

1. Task title
2. Objective
3. Scope
4. Primary goals
5. Files to inspect first
6. Relevant current symbols
7. Implementation requirements
8. Non-goals
9. Design guidance
10. Suggested implementation steps
11. Required deliverables
12. Validation requirements
13. Completion criteria

## Phase Prompt Guardrails

Every next-phase prompt should say:

- the master roadmap remains authoritative
- do not rewrite the roadmap during implementation
- do not redesign unrelated parts of the system
- do not silently widen support boundaries

## After Prompt Generation

- save the new prompt in the phase work log or task record used for the implementation
- update `status.md` if the active phase changed
