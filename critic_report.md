# Design critique — Stage 7, round 3

**Target:** [design_draft.md](design_draft.md), feature specification lines 1–1894.\
**Status:** all three independent reviewers completed; two consecutive clean rounds reached.\
**Reviewed SHA-256:** `5e3b700cece542ac971ee68a0f5cf2262164dbb5f05b077ee5c53a9d180d52da`

| Reviewer | Status | Blocker | Major | Nit |
| --- | --- | --- | --- | --- |
| Design | Completed | 0 | 0 | 0 |
| Structural | Completed | 0 | 0 | 0 |
| Style | Completed | 0 | 0 | 2 |

Rounds 2 and 3 are clean under EGM's rule: no unsettled Blocker or Major findings
and no missing reviewers. The critique has converged; a fourth round is not
required. Both minor wording edits were approved and applied before the Stage 8
handoff. Their source sentences predate the approved cleanup, and neither fix
requires reverting an earlier edit.

The target and all 72 linked context files remained unchanged during this round.
Required read receipts are included below. The detailed implementation plan was
outside this review scope. No production code was edited or numerical tests run.

[Applied cleanup summary and combined diff](/tmp/madqc-final-independent-review-yapnzvuv/stage7-round3/applied_cleanup.md).\
[Annotated findings and exact proposed replacements](goldfish_review_report.md).

## Design reviewer — verbatim report

Status: **completed**. Findings: **none**—no Blocker, Major or Nit findings; no finding IDs assigned.

The active feature specification addresses the standard’s problem, scope, architecture, alternatives and risk requirements. Scientific definitions and numerical defaults still requiring validation are explicitly identified as adapter prerequisites, rather than presented as implemented or validated behavior.

Read receipt:

- Read the complete persona: `/home/ahurtado/.codex/plugins/cache/egm/egm/0.1.1/agents/design_critic.md`, lines 1–63.
- Read the complete standard: `/home/ahurtado/.codex/plugins/cache/egm/egm/0.1.1/resources/document_standards/feature-spec-standard.md`, lines 1–39.
- Read **all target lines 1–1894** of [design_draft.md](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/design_draft.md), in chunks. Reread the return-structure passage affected by output truncation.
- Read the supplied `stage7-round3/reading_list.txt`.
- Target SHA-256, identical before and after review: `5e3b700cece542ac971ee68a0f5cf2262164dbb5f05b077ee5c53a9d180d52da`.

Actual context inspected, relative to the repository root:

- `src/madness/chem/StageCheckpoint.h`: lines 1–80 and 341–404, plus targeted publication, validation, generation and retention searches.
- `src/madness/chem/CheckpointRetention.h`: complete file, lines 1–54.
- `src/apps/molresponse/calc/calc_manager.hpp`: lines 567–650, plus protocol/planning searches.
- `src/apps/molresponse/calc/calc_executor.hpp`: lines 1725–1815 and 1880–1940, plus protocol/subworld searches.
- `src/apps/molresponse/solvers/iterate_protocol.hpp`: targeted protocol/preparation search results.
- `src/apps/molresponse/solvers/es_kinetic_evidence.hpp`: complete file, lines 1–68.
- `src/apps/molresponse/solvers/es_root_identity.hpp`: lines 1–110.
- `output_record_examples.md`: lines 150–190 and 345–380, plus targeted source/snapshot/publication searches.

Excluded the Detailed Implementation Plan. Read no prior reviews or progress artifacts. Made no edits, ran no tests and spawned no agents.

## Structural reviewer — verbatim report

Status: **completed**.

**Findings: none.** No finding IDs assigned; 0 Blockers, 0 Majors, 0 Nits.

The active feature specification structurally covers the standard’s required areas:

- Problem: lines 29–62.
- Goals and non-goals: lines 64–100.
- Architecture, responsibilities and data flow: lines 102–1670.
- Alternatives and concrete tradeoffs: lines 1625–1628 and 1672–1685.
- Failure modes, scalability limits, security implications and mitigations: particularly lines 579–587 and 917–1027.
- Validation and release requirements: lines 1712–1893.

The standard does not prescribe exact section names or require a standalone risks heading. The distributed risk discussion therefore does not constitute a missing-header finding.

Heading levels are sequential, all five fenced blocks close, and table rows have consistent delimiters. All 75 Markdown links resolve to existing local files or identifiable section anchors. Planned contracts are explicitly distinguished from implemented behavior; outstanding adapter work was not treated as a structural defect or implementation-readiness finding.

Read receipt:

- Read the complete structural persona and complete feature-specification standard.
- Read every active target line, **1–1894**, in untruncated chunks.
- Excluded the Detailed Implementation Plan’s body from review. The initial heading inventory exposed its headings only.
- Read the complete supplied GoldfishContext reading list.
- Additional context inspected: the heading index of `functional_design_walkthrough.md`, including its referenced evaluation-return anchor; and `output_record_examples.md`, lines 1–110. Checked existence of the active specification’s linked files.
- Read the session `using-superpowers` skill, whose explicit subagent exemption applies.
- No prior reviews, progress files, source implementations or other scratch artifacts were inspected. No edits or tests performed.
- Target SHA-256, unchanged across inspections: `5e3b700cece542ac971ee68a0f5cf2262164dbb5f05b077ee5c53a9d180d52da`.

## Style reviewer — verbatim report

Status: **completed**. Mode: REVIEW; no separate style standard. Two Nit findings; no Blocker or Major findings.

Target: [design_draft.md](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/design_draft.md), lines 1–1894.

- **style-3-1 — Nit — lines 1439–1440.** “The user-facing criteria, metric definitions, tolerance policy and checkpoint cadence are unchanged” does not identify the baseline. The specification changes several existing conventions, so a fresh reader could interpret this as a claim of compatibility with current behavior. **Editorial action: CALIBRATE.** State the intended equivalence between the two submission mechanisms: “Direct submission and callbacks use the criteria, metric definitions, tolerance policy and checkpoint cadence specified above.”

- **style-3-2 — Nit — lines 1642–1645.** Adjacent sentences repeat the absence of old-checkpoint compatibility requirements. The second sentence adds useful distinctions about reusing archive primitives and preserving convergence claims, but repeats the format boundary. **Editorial action: consolidate.** State the compatibility boundary once while retaining those distinctions and the separate collaborator-approval requirement.

Read receipt:

- Read the complete `/home/ahurtado/.codex/plugins/cache/egm/egm/0.1.1/agents/style.md` persona.
- Read every active target line, 1–1894, in chunks. Reread lines 1441–1680 separately after aggregate output truncation. Excluded the Detailed Implementation Plan.
- Target SHA-256, identical before and after review: `5e3b700cece542ac971ee68a0f5cf2262164dbb5f05b077ee5c53a9d180d52da`.
- Read the complete supplied `reading_list.txt`.
- Context inspected, relative to the repository root:
  - `docs/design/madqc-simplified-policy.md`: complete, lines 1–70.
  - `src/apps/molresponse/calc/calc_manager.hpp`: relevant search matches and lines 567–605.
  - `src/apps/molresponse/tests/test_calc_manager.cpp`: relevant search matches and lines 702–780.
  - `src/madness/chem/StageCheckpoint.h`: relevant search matches and lines 1–18, 292–325, 361–403.
  - `src/madness/chem/test_stage_checkpoint.cc`: relevant search matches and lines 72–111, 197–216.
- No prior reviews, progress files, conversation history or other scratch artifacts inspected. No edits, tests or agent delegation performed.
