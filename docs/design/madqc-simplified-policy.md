# madqc simplified iteration and restart policy

This amendment starts from `madqc-drop-6` (9179f5d32076a2ba555906dd146dc330ffc040d7). It records the November 6, 2026 design discussion and intentionally revises the older B5/B6/B8 behavior; it is not a behavior-preserving cleanup or a claim that the entire new policy is implemented.

## Accepted direction

Intermediate protocol stages perform their configured iteration count. The requested final representation is mandatory even if convergence criteria pass earlier. At the final protocol, all requested applicable criteria should govern convergence, with the iteration cap reported separately from convergence. Every executed stage receives a fresh allowance on each run; prior counts remain history. KAIN starts with an empty subspace.

Partial intermediate work resumes its saved stage. Completed intermediate work advances only when the numerical setup remains valid; changed numerical settings restart the saved stage. An otherwise unchanged appended tighter stage is executed next. Unconverged final checkpoints may resume. Valid unchanged converged ground states can be loaded without new iteration.

The ground state must satisfy the final request before dependent calculations. Independent frequency states retain completed valid work and share protocol boundaries. Coupled excited-state roots remain active together; larger bundles reuse old roots and add independent guesses at the saved stage. Refining the ground state invalidates dependent convergence while compatible functions remain useful guesses. Fresh comparable observations are required after numerical changes.

Retain only the latest usable numerical checkpoints, preserving the prior checkpoint until its replacement is successfully published. Keep lightweight iteration, timing and metric histories across stages and runs. Disabling restart means fresh initialization, with checkpoint writing still enabled.

## First implementation slice

This slice implements fixed intermediate iteration counts and removes the final-stall refusal to resume. The shared exit policy receives an explicit final-protocol flag. Intermediate stages ignore convergence and plateau exits until their configured count, unless the engine declares divergence. At completion they report the actual existing convergence verdict or iteration limit. The diagnostic count, when explicitly supplied, retains precedence.

Response execution must determine the role from the whole requested node ladder, because the scheduler invokes a single protocol at a time. Direct multi-stage response driver calls likewise treat every stage before the last as intermediate. Named SCF uses its supplied ladder. Standalone SCF behavior is unchanged.

Final convergence selection and optional live final plateau exits retain their existing behavior in this slice. A saved stalled final stage with unmet criteria receives a resume plan, including its historical count without deducting that count from the new allowance. A saved final state that now satisfies the applicable existing checks can still be loaded.

| Scenario | Required result in this slice |
| --- | --- |
| Intermediate checks pass after 2 of 10 iterations | Continue through iteration 10 |
| Intermediate plateau after 2 of 10 iterations | Continue through iteration 10 |
| Intermediate count reached, checks unmet | Record iteration limit; advance to the supplied next stage |
| Early intermediate convergence | Still execute the requested final stage |
| Divergence during intermediate work | Stop; do not advance |
| Single-threshold response wave with later requested stage | Treat as intermediate |
| Saved final plateau, criteria still unmet | Resume at final with a fresh allowance |
| Saved final checks satisfied | Load, under existing validity rules |

## Remaining alignment work

1. Make requested final criteria govern every live solver and reuse decision, with explicit metric definitions and defaults; remove optional final plateau termination under the simplified policy.
2. Add numerical-setup evidence and conservative validity rules for completed-stage advancement and whole-model restarts.
3. Keep all coupled excited-state roots active, validate the existing no-lock numerical cases, and make their live stopping agree with their published verdict.
4. Implement compatible bundle expansion and response-guess reuse after ground-state refinement, clearing stale convergence evidence.
5. Test mixed-progress restart, within-stage checkpoint cadence, safe replacement and retention, and durable per-state histories. The current response save boundary is solve completion.

The existing vault design and decision log remain historical records. B5's early intermediate exit and B6's stalled-final load-only rule are superseded by this slice. B8, D-025, D-092 and the remaining B6 changes are tracked above, not silently claimed complete.

## Validation

Pending implementation and validation. Policy tests alone do not establish molecular convergence, root identity, or restart crash safety.
