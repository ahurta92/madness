# Final design review — Stage 8: Intermediate Zoom Out

**Date:** 2026-11-12\
**Target:** [design_draft.md](design_draft.md), all 2445 lines.\
**Feature range:** 1–1894; **implementation range:** 1895–2445.\
**Status:** resolved — timing correction approved and applied; Stage 9 authorized.\
**Result:** 0 Blocker, 0 Major, 1 Nit.

Both approved Stage 7 wording edits are applied. The critique is closed after
two consecutive clean rounds. Stage 8 found no architectural scope drift or
unnecessary new framework in the intervening clarifications. Source inspection
supports reuse of existing drivers, metrics, response planning, archive bindings
and checkpoint publication.

## Proposed correction before Stage 9

**zoom8-1 — Nit — PROPOSE FIX.** At [line 2105](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/design_draft.md:2105),
the implementation plan says to reuse the existing wall-clock source. The POSIX
`madness::wall_time()` uses `gettimeofday`, which can jump with system-clock
adjustments. This does not meet the monotonic elapsed-time requirement at lines
1344–1351. The source evidence and complete reviewer rationale are below.

Replace:

```text
Reuse the existing wall-clock source; keep function-bearing work out of JSON records.
```

With:

```text
Use `std::chrono::steady_clock` for elapsed intervals; keep function-bearing work out of JSON records.
```

This changes one implementation instruction. Measurement boundaries, driver-rank
scope, the no-extra-fence rule and controlled-clock tests remain as specified.
It does not require changing the runtime's existing timer API.

[Exact proposed diff](/tmp/madqc-final-independent-review-yapnzvuv/stage8/proposed_changes.diff). The correction was approved with an explicit scope: use steady_clock inside the
shared collector and keep existing MADNESS timer APIs and callers unchanged.
That scope is now recorded in both the feature contract and implementation plan.
Stage 9 is authorized.

## Review evidence and limits

The independent reviewer read the full current draft, persona, standard, root
README, linked-file list, initial Stage 5 report and full diff from the Stage 5
snapshot. The reviewed target and 74 tracked context files (72 linked files,
root README and initial report) remained unchanged during review. The primary
agent also inspected the clock implementations before proposing the correction.

Target SHA-256: `69597c264e3a1a190f0fea244c798c0a8af073509190252d62638b591ccc58d1`.

This completes Stage 8's audit, not implementation readiness or numerical
validation. Stage 9 remains. No production code was edited, tests or builds run.

## Independent reviewer — verbatim report

**Status: completed.** Stage 8 Intermediate Zoom Out found **0 Blocker, 0 Major, 1 Nit**.

The design remains focused on consolidating competing lifecycle, convergence, restart and publication behavior while retaining engine-owned scientific state and algorithms. The additions since Stage 5 clarify existing boundaries without introducing another scheduler, recovery system or convergence framework. I assessed lines 1895–2445 at implementation altitude.

**zoom8-1 — Nit — Clarify the clock reuse instruction**

Affected current scope: [design_draft.md:2105](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/design_draft.md:2105), with the timing requirement at **lines 1344–1351**.

The implementation table directs the new collector to reuse the existing wall-clock source, while the contract requires monotonic elapsed timing. The existing `madness::wall_time()` uses `gettimeofday` on POSIX ([timers.cc:48–75](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/src/madness/world/timers.cc:48)). Response’s existing `StageTimer` also uses that source ([response_workflow.hpp:114–120](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/src/apps/molresponse/orchestrator/response_workflow.hpp:114)). System-clock adjustments can therefore distort intervals or make them negative.

This is a narrow suitability issue in the reuse proposal. Specify `std::chrono::steady_clock` for the new collector, preserving the existing measurement boundaries and no-extra-fence rule. The repository already uses that standard primitive in [timers.h:139–141](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/src/madness/world/timers.h:139). No broader timer framework or runtime-wide clock refactor is needed.

**Architectural assessment**

- Bundle validity remains an application prerequisite; the clarification adds no aggregate convergence metric or shared orthogonality policy.
- Unconverged observations reuse existing measurements and remain separate from completed property output. The clarification requires no additional scientific computation.
- Checkpoint trust, path and concurrency restrictions match the existing storage implementation. They do not expand recovery guarantees.
- Per-task progress consolidates existing response planning. The draft explicitly identifies the threshold-only matching that needs adaptation.
- The proposed records and thin bindings connect existing operations. Source inspection supports reusing workflow handoff, named-runner callbacks, numerical metric helpers, equation evaluators, archive serializers and checkpoint publication.
- The finite lifecycle test harness fits the repository’s existing distinction between numerical QC cases and scripted sequences.

**Read receipt**

Read completely:

- [Zoom Out persona](/home/ahurtado/.codex/plugins/cache/egm/egm/0.1.1/agents/zoom_out.md), **74 lines**.
- [Zoom Out Standard](/home/ahurtado/.codex/plugins/cache/egm/egm/0.1.1/resources/document_standards/zoom-out-standard.md), **60 lines**.
- [Current target](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/design_draft.md), **all 2445 lines**: feature specification **1–1894**, implementation plan **1895–2445**.
- [Root README](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/README.md), **86 lines**.
- [Supplied reading list](/tmp/madqc-final-independent-review-yapnzvuv/stage8/reading_list.txt), **72 paths**.
- [Authorized initial Stage 5 report](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/zoom_out_report.md), **86 lines**.
- Complete diff between the authorized [Stage 5 snapshot](/tmp/madqc-final-independent-review-yapnzvuv/design_draft.snapshot.md) and the current target.

Target SHA-256, confirmed unchanged after inspection:

`69597c264e3a1a190f0fea244c798c0a8af073509190252d62638b591ccc58d1`

Actual source inspection covered:

- Workflow composition and handoff: [Applications.hpp](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/src/madness/chem/Applications.hpp), [Drivers.hpp](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/src/madness/chem/Drivers.hpp), [WorkflowBuilders.hpp](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/src/madness/chem/WorkflowBuilders.hpp), and [RunNamedStages.hpp](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/src/apps/madqc/RunNamedStages.hpp).
- Metric arithmetic and policy/planning: residual and density helpers, scalar metric summaries, kinetic operations, convergence settings, explicit criteria and ladder planning.
- Response scheduling and evaluation: [calc_manager.hpp](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/src/apps/molresponse/calc/calc_manager.hpp), [fd_solver.hpp](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/src/apps/molresponse/solvers/fd_solver.hpp), root-identity and kinetic-evidence helpers, with focused searches of seed and response-space operations.
- Storage and bindings: [StageCheckpoint.h](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/src/madness/chem/StageCheckpoint.h), [CheckpointRetention.h](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/src/madness/chem/CheckpointRetention.h), response/VBC archive membership, Nemo save/load, result types and Dirac binding/build searches.
- Timing implementations and the [QC test boundary](/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/src/examples/qc/README.md:388). Relevant baseline policy and application documentation were also inspected.

This was a bounded static review. No files were edited, tests or builds run, additional agents spawned, or other prior review reports retrieved. The proposed behavior remains implementation work; this report does not certify implementation readiness or numerical correctness.

Resolution: zoom8-1 fixed with the author’s explicit existing-timer preservation constraint.
