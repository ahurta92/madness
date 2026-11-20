# Goldfish annotated verification report — round 5

**Author resolution — 2026-11-16:** the retention policy below was subsequently
settled and applied to the design and test plan. Protect latest non-retired task
checkpoints, current published result sources and active dependencies transitively;
historical metadata alone does not pin payloads. Cleanup follows successful
publication and complete reference checks. The user authorized shared-core
implementation at this review boundary. The round-5 findings below remain a
historical review of their stated hash; this amendment is not a new independent
review or a waiver of scientific/adoption prerequisites.

**Target:** [design_draft.md](design_draft.md), focused result/source-map binding
with full-contract consistency reading; plan lines 2512–3348.
**Status:** requested binding fix applied and verified by all three independent
reviewers. One separate retention-policy finding remains proposed and unapplied.
The fifth-round cap is reached; no further automatic review is scheduled.

## Closed: readiness-4-1 — calculation source-map attempt

The approved root/calculation attempt distinction, inheritance/override rule,
exact map path, strict validation and distinction from numerical producing
attempts are explicit at design_draft.md:1814–1839 and 1903–1914. Examples and
RESULT-SOURCES coverage agree. Old snapshots resolve their original maps despite
new electronic selections, including unfinished and concluded newer attempts.
All three reviewers verified this clarification; source implementation remains
planned work. The scoped serializer/reader lookup is ready to implement.

## readiness-5-1 — authoritative retained-result membership and release

**Severity:** Major. **Reviewers:** implementation and codebase.
**Category:** PROPOSE FIX. **Status:** policy decision proposed, not applied.

The contract requires protecting archives referenced by retained old snapshots
but does not specify what establishes or ends that retained status. A and B
can share a calculation ID, so calculation retirement alone cannot represent
their different lifetimes. Current results alone can drop A's protection;
retaining every historical map indefinitely prevents bounded cleanup.

**Evidence:** design_draft.md:1216–1233,1391–1405,1632–1662,1832–1839;
madqc-test-index.md:310–320; CheckpointRetention.h:43–51;
StageCheckpoint.h:444. Raw readiness-5-1 and codebase-5-1 are deduplicated here.

**Proposed resolution:** select the authoritative set of live result references
and the acquire/release ordering used by publication and cleanup. Prefer an
explicit bounded rule for existing records if only current published results and
active checkpoint/workflow dependencies require archives; do not make every
historical source map a permanent retention root. If old snapshots must remain
reloadable after replacement, define their explicit membership and release.
These are alternatives for the next decision, not new implemented features.

Add actual cleanup and interruption coverage: preserve A while it is required
after B publishes, reload the same authoritative reference state on restart,
then release A and reclaim only after no dependency requires it. Preserve its
original source identity even after a permitted payload reclamation; never
substitute B's state. Align inspection expectations with the selected policy.

**Affected phases:** step 4 retention/publication, step 5 optimizer integration,
steps 6–7 native old-snapshot acceptance. This does not block the measurement/
controller foundation or the now-verified source-map lookup and synthetic tests.

## Review boundary

The requested fix passes. The remaining issue concerns the lifetime of archive
protection, not ambiguity in the new attempt-reference fields. It was exposed
by checking the specified old-result retention interaction. No retention policy
was changed and the original binding finding is not reopened.

The carried scientific-definition prerequisite (readiness-1-1) and deferred
CC2/6D work (readiness-1-2) retain their existing scope. The previously recorded
editorial Nit about repeated finite-difference prose remains declined and does
not block. This round has one unique Major and no new Blocker or Nit.

The review loop has reached five rounds and is at its decision checkpoint:
accept the review boundary and proceed with unaffected programming, explicitly
authorize further review, or return the retention question to drafting. No risk
acceptance or full-release readiness is inferred. See
[the readiness report](implementation_readiness_report.md) for the programming boundary.
