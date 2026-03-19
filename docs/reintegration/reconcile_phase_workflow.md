# Reconcile Phase Workflow

Use this workflow after a phase implementation is complete.

Master roadmap:

- [`../excited_state_reintegration_gap_analysis.md`](../excited_state_reintegration_gap_analysis.md)

## Goal

Make the roadmap describe the repository as it now exists.

Do not rewrite the roadmap from scratch.
Update only the sections changed by the completed work.

## Required Inputs

Before reconciling, gather:

1. The master roadmap
2. The code changes that were actually implemented
3. Phase notes
4. Validation results
5. Any new support boundaries or constraints discovered

## Reconciliation Checklist

### 1. Identify the completed phase

- Confirm the roadmap phase number and title.
- Confirm whether the implementation fully completed the planned phase or only part of it.

### 2. Verify actual implementation

- Read the modified code, not just the prompt or summary.
- Confirm the real symbols, files, and metadata paths used.
- Record any deviations from the planned design.

### 3. Update the roadmap minimally

Update only the sections that changed:

- Executive Summary
- Current Solver Call Graph and Control Flow
- Comparison Matrix
- Gap Inventory
- Ordered Implementation Roadmap
- Risks and Design Tensions
- Codex Implementation Brief
- Code Reference Index

### 4. Record completed work precisely

For the completed phase:

- mark `Status: Completed`
- add an `Implementation Summary`
- record actual design decisions
- record deviations from the original proposal
- record new constraints discovered
- record follow-up work moved to later phases

### 5. Update remaining work

- Remove resolved gaps from the live gap inventory or convert them to residual follow-up items.
- Adjust later phases if the dependency order changed.
- Tighten support boundaries if the implementation revealed a narrower valid scope.

### 6. Synchronize the scaffold

After the roadmap update:

- update [`status.md`](status.md)
- update [`contracts/active_contracts.md`](contracts/active_contracts.md)
- ensure the next phase prompt will use the reconciled state

## Output Expectations

The reconciled roadmap should let the next agent answer:

- what is actually complete
- what remains incomplete
- what design decisions were already taken in code
- what constraints must be preserved in the next phase
