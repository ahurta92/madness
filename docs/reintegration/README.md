# Excited-State Reintegration Workflow

This directory is the working scaffold for the remaining excited-state `molresponse_v2` reintegration work.

Master roadmap:

- [`../excited_state_reintegration_gap_analysis.md`](../excited_state_reintegration_gap_analysis.md)

That roadmap is the source of truth for:

- completed versus planned phases
- remaining gaps
- implementation order
- architectural risks
- acceptance criteria

This scaffold does not replace the roadmap. It organizes the phase-by-phase execution loop around it.

## Workflow Loop

Use this loop for each implementation phase:

1. Prompt
   - Generate the next phase prompt from the master roadmap and current status.
2. Implementation
   - Implement only the scoped phase work.
3. Notes
   - Record what actually changed, what differed from the plan, and what was discovered.
4. Reconciliation
   - Update the master roadmap so it reflects the repository as implemented, not as previously planned.
5. Contract Update
   - Update active contracts, invariants, and support boundaries.
6. Validation
   - Record what was tested, where, and what remains unvalidated.
7. Next Prompt
   - Generate the next implementation prompt from the reconciled state.

## Start Here

1.docs/excited_state_reintegration_gap_analysis.md
2.docs/reintegration/context.md
3.docs/reintegration/contracts/active_contracts.md
4.latest phase_notes
5.the current phase prompt

## Files

- [`context.md`](context.md)
  - architectural guardrails, current implementation state, and `ResponseVector` variant guidance
- [`status.md`](status.md)
  - lightweight dashboard for phase progress
- [`reconcile_phase_workflow.md`](reconcile_phase_workflow.md)
  - how to reconcile completed work back into the roadmap
- [`generate_next_phase_prompt.md`](generate_next_phase_prompt.md)
  - how to write the next phase prompt from the current state
- [`phase_prompts/TEMPLATE_phase_prompt.md`](phase_prompts/TEMPLATE_phase_prompt.md)
  - template for future implementation prompts
- [`phase_notes/TEMPLATE_phase_notes.md`](phase_notes/TEMPLATE_phase_notes.md)
  - template for recording actual implementation decisions and deviations
- [`contracts/TEMPLATE_phase_contract.md`](contracts/TEMPLATE_phase_contract.md)
  - template for a phase-specific execution contract
- [`contracts/active_contracts.md`](contracts/active_contracts.md)
  - standing contracts and currently active support boundaries
- [`validation/TEMPLATE_phase_validation.md`](validation/TEMPLATE_phase_validation.md)
  - template for per-phase validation records

## Usage Rule

Before starting any new phase:

- read the master roadmap
- read [`status.md`](status.md)
- read [`contracts/active_contracts.md`](contracts/active_contracts.md)
- create or update the phase contract
- generate the phase prompt from the reconciled state

After finishing any phase:

- write phase notes
- reconcile the roadmap
- update status
- update active contracts
- write validation
