# Active Contracts

Master roadmap:

- [`../../excited_state_reintegration_gap_analysis.md`](../../excited_state_reintegration_gap_analysis.md)

This file records the standing contracts that are currently active for the excited-state reintegration work.

## Standing Architectural Contracts

1. Stage 2c remains the excited-state workflow entry point.
2. Legacy `X_space` must not return as a public or persistent architecture.
3. Stable root identity is defined by descriptor data, not by bundle slot.
4. Stage 2c metadata is authoritative through `ResponseRecord2`.
5. Workflow metadata and `response_metadata.json` must remain in parity for the excited subtree.
6. Restart precedence remains:
   `current protocol snapshot -> nearest lower protocol snapshot -> generic guess archive -> carryover -> fresh guess`

## Current Identity Contract

- `root_id` and `stable_index` are identity.
- `slot_index` is runtime placement only.
- `display_name` is presentation only.
- `slot_permutation` records the current bundle-slot to stable-root mapping.

## Current Variant / Restart Contract

| Variant | Current restart support mode |
| --- | --- |
| `static_restricted` | `full_bundle_resume` |
| `static_unrestricted` | `guess_only_unrestricted_variant` |
| `dynamic_restricted` | `full_bundle_resume` |
| `dynamic_unrestricted` | `guess_only_unrestricted_variant` |

Rules:

- Do not advertise broader restart support than the code actually implements.
- If a future phase changes any variant support boundary, update this file and reconcile the roadmap.

## Phase-Specific Contract State

Current active implementation phase:

- Phase 4 — Iteration and Convergence Logic

Current phase-specific contract:

- [`phase_4_iteration_and_convergence_contract.md`](phase_4_iteration_and_convergence_contract.md)

Most recently completed phase-specific contract:

- Phase 3 — Solver Initialization Path
- [`phase_3_solver_initialization_contract.md`](phase_3_solver_initialization_contract.md)
