# madqc focused independent verification — round 5

**Author resolution — 2026-11-16:** the retention policy below was subsequently
settled and applied to the design and test plan. Protect latest non-retired task
checkpoints, current published result sources and active dependencies transitively;
historical metadata alone does not pin payloads. Cleanup follows successful
publication and complete reference checks. The user authorized shared-core
implementation at this review boundary. The round-5 findings below remain a
historical review of their stated hash; this amendment is not a new independent
review or a waiver of scientific/adoption prerequisites.

**Completed 2026-11-16.** Three fresh reviewers read the unchanged 3,348-line
draft, with the implementation plan at lines 2512–3348, and checked the focused
result-attempt clarification against its companion records, tests and source.

**The requested fix is applied and independently verified by all three
reviewers.** Each result calculation now identifies its immutable source-map
owner explicitly or inherits the enclosing run attempt. Numerical producing
attempts remain distinct. The examples and planned old-result tests agree.

Two reviewers found one separate Major question about which historical results
keep their numerical archives protected from cleanup. The binding itself is
ready to implement. Retention-dependent acceptance remains gated by that policy.

## Verified clarification

- Root `results.json.attempt_id` identifies the run writing the snapshot.
- Calculation `attempt_id` identifies the source-map owner. Omission inherits
  the root; a differing owner must be explicit.
- Calculation ID and effective attempt resolve one immutable source-map path.
  Invalid, missing or mismatched maps fail without current-selector fallback.
- Numerical states inside that map retain their own original attempts, samples
  and generations. A later optimizer selection or checkpoint cannot relabel them.
- RESULT-SOURCES/C03/C04 and the owning serializer/shared/native-reader plans
  cover inherited, explicit and mixed owners; an older result during newer
  unfinished and concluded work; older numerical producers; and strict errors.

See [the field contract](design_draft.md:1814),
[publication and source retention](design_draft.md:1903),
[record example](output_record_examples.md:181), and
[planned tests](docs/design/madqc-test-index.md:293).

The former readiness-4-1 finding is closed by this verification. The earlier
seed/protocol, regression curation, native-reader migration and finite-difference
plan revisions remain verified at their reviewed scopes. This does not claim
that their planned implementation or numerical validation has been performed.

## Remaining Major: define the lifetime of retained results

The design requires preserving archives needed by retained result snapshots.
However, it does not yet define which superseded snapshots count as retained,
where that membership is authoritative, or when their archive references end.
A and B may have the same calculation ID but different electronic attempts:
calculation-level retirement cannot distinguish their lifetimes.

After B replaces the current result, protecting only current `results.json` may
lose an intentionally retained A. Protecting every historical source map forever
would prevent the intended bounded cleanup. The existing retention helper accepts
a caller-supplied reference set; it does not choose this policy.

**Required before retention integration:** define that live reference set and
its publication/release order. Test cleanup and restart while A is still required,
then release A and verify reclamation only after its last dependency ends. The
source binding must remain an immutable record of origin even when an unneeded
numerical payload is later reclaimed.

A bounded policy based on current published results and live checkpoint/workflow
dependencies may be sufficient; explicit historical-payload pinning should only
be introduced if required. This choice and any adjusted old-snapshot loading
expectations remain proposed, not applied. No new registry or keep-all policy
has been added.

## Independent outcomes

| Reviewer | Blocker | Major | Nit | Binding clarification |
| --- | --- | --- | --- | --- |
| Implementation readiness | 0 | 1 | 0 | Verified |
| Codebase details | 0 | 1 | 0 | Verified |
| Clarity | 0 | 0 | 0 | Verified |

The two Major observations are one issue, readiness-5-1. No additional finding
was reported. The round is not clean because the retention decision remains
open. This is the fifth round, reaching the EGM loop's configured cap; no further
review round or retention policy change is authorized automatically.

## Programming boundary

**Startable now:** curated core tests, typed measurements/shared arithmetic,
pure convergence/restart policy, result-attempt serializers/validators and
synthetic source-map lookup tests, preserving existing applications.

**Still gated:** retention-dependent publication/cleanup and full native loading
of intentionally retained older snapshots require the lifetime decision. Owner-
confirmed scientific definitions and numerical calibration/accuracy/performance
evidence remain the existing gates for their affected adapters. Optimizer
replacement/history, shared HDF5 and deferred CC2/6D science are not brought into
this work by the review.

The approved changes were to the design, examples and planned test descriptions.
No C++, executable tests or numerical references changed. Structural checks
passed; no builds, MPI tests or numerical calculations ran. The draft and selected
context remained unchanged during this independent verification.

Independent reports: [readiness](/tmp/madqc-final-independent-review-yapnzvuv/stage9-round5/implementation_report.md),
[codebase](/tmp/madqc-final-independent-review-yapnzvuv/stage9-round5/codebase_report.md),
[clarity](/tmp/madqc-final-independent-review-yapnzvuv/stage9-round5/style_report.md).

Reviewed SHA256: `9d4c8cc8d561b103b345bfc2ef226f8047751862ad063a28a419cf3ab3341179`.
