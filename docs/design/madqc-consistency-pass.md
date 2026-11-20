# Final design consistency pass

2026-11-13. Compared the reviewed [acceptance checklist](madqc-acceptance-checklist.md)
with the [design and implementation plan](../../design_draft.md),
[input mapping](../../input_option_mapping.md),
[functional walkthrough](../../functional_design_walkthrough.md),
[SCF/response records](../../output_record_examples.md) and
[grouped tensor example](../../tensor_result_examples.md).
This is a document consistency pass, not an independent Goldfish review,
code adoption audit or numerical validation.

## Corrections applied

| Accepted case | Before | After |
| --- | --- | --- |
| A01 | Input mapping and functional walkthrough still prescribed a custom migration message for `excited.maxiter`. | Both use generic unknown-input rejection, including mixed inputs. The separate `dconv`/`econv` migration guidance is unchanged. |
| A04, B01 | One restart example treated absent comparison history as a reason to establish a new baseline; the input mapping loosely said to reuse completed states. | Valid saved scalar evidence is assessed first even when previous functions are absent. Final qualifying states are reused; completed intermediate work advances; unfinished work resumes. |
| A06, C03 | The functional example used the undeclared reason `no_comparable_sample` and an outdated description of record approval. | It uses `no_previous_sample` and points to the accepted record contract. |
| B02, C01, C04 | Controller prose referred to a manifest update as publication. The timing example put the response checkpoint after VBC/properties. | Publication atomically selects the validated generation through `latest`. The example publishes the response sample before dependent work, under its task timing. Direct-task timings carry their final protocol explicitly. |
| C03, C04 | VBC prose put build timing in its checkpoint, while the record contract keeps timing in `timings.json`. | VBC build timing stays in `timings.json` with the task identity; checkpoint fields retain payload, completion and exact dependencies. |
| B03 | The response implementation row grouped ambiguous matching with invalid-state conditions. | Ambiguity makes affected comparisons unavailable; bundle validity remains a separate check. Unselected change diagnostics cannot block completion. |
| C02, D03 | The shared retention-test paragraph mentioned only optimizer snapshots. | It also covers finite-difference parent plans and their publication-before-retirement rule. |
| All cases | The plan did not link the reviewed acceptance IDs into test curation and still said to create the existing metric evidence document. | The design links the checklist; test curation maps its IDs to existing/adapted checks. The evidence document is extended. No passing coverage is implied. |

## Coverage checked

| Acceptance cases | Contract and planned implementation coverage |
| --- | --- |
| A01–A06 | Input routing, prepared-state metrics, intermediate/final criteria, saved-evidence reuse and fresh update allowances; steps 2–3 and adapter checks in steps 6–7. |
| B01–B03 | Independent progress, exact dependencies and coupled-root validity/comparisons; restart planning in step 2, workflow boundaries in step 5 and response adoption in step 7. |
| C01–C02 | Strict selected-checkpoint recovery, atomic publication, endpoint/optional saves and protected-generation retention; step 4 and repeated-evaluation checks in step 5. |
| C03–C04 | Common property/availability formats, immutable source associations and shared metrics; step 4. Timing hierarchy and overlap use controlled-clock checks in step 3; numerical-reference mapping remains in step 11. |
| D01–D03 | Saved geometry identities, optimizer save boundaries, history validity, fresh step allowance and finite-difference plan reuse; deterministic-provider checks in step 5, using step 4 publication/retention. |

The grouped tensor example required no changes. Scientific values, tolerance
examples, numerical references, source-code descriptions and the release
inventory are unchanged. Only illustrative timing-parent durations changed to
include the checkpoint within the response solve.

## Validation and remaining work

Validation checks cover embedded JSON syntax, property component accounting,
measurement reason codes, timing nesting/publication order, local links,
acceptance IDs and unchanged scientific prerequisites. A before/after diff and
machine-readable checks are retained with the review artifacts. All six JSON
blocks, 163 local links and 11 document anchors passed the document checks.

Florian's Nemo confirmations and other model-specific prerequisites remain open.
CC2/6D work stays deferred from the current review. The six previously applied
review repairs and these consistency corrections still need a fresh independent
review. No source implementation, numerical test or new review round was run.
