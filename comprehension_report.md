## What this proposes

The proposal consolidates chemistry calculations around the retained workflow interface and a shared lifecycle for settings, convergence, restart, checkpoints and reporting, while each engine keeps its scientific state and numerical methods. Applications measure the prepared state through shared helpers; a common controller applies the selected criteria, permits intermediate-protocol advancement and requires final-protocol convergence. Restarts assess compatible saved measurements against the current request, preserve each task’s progress and fail strictly when an expected checkpoint is invalid. Response keeps its internal dependency manager, treats excited roots as an active coupled bundle and releases direct computations only when their inputs qualify. The implementation plan requires adoption by E01–E12, preserves E13–E15 as deferred standalone programs and explicitly identifies scientific definitions and numerical validation still needed.

## How the existing system works, per this document

Madqc currently has separate named-stage and unnamed-workflow execution paths that share numerical engines but make some lifecycle, restart and publication decisions separately. Engines differ in stopping rules, measurement conventions and measurement timing; some excited-state paths freeze roots or perform separate warmup work. Molresponse already has iteration and protocol drivers, an internal dependency manager, root-tracking machinery and numerical persistence, while SCF and Nemo retain their own loops. Existing checkpoint infrastructure already provides validated publication and retention operations, but the proposal still needs to connect these pieces to consistent current-state evidence, dependency identities and output records.

No comprehension findings. After reading the complete target, I could paraphrase every section and connect the proposed responsibilities, execution policy, restart evidence, dependency handling, record layout and implementation sequence. The unresolved scientific definitions and validation work are identified explicitly rather than presented as settled behavior; they do not prevent understanding the proposal.

Read receipt:

- **Coverage:** All lines 1–2427 of `/home/ahurtado/.codex/worktrees/madqc-simple-policy/madness-contract/design_draft.md`, including FeatureRange 1–1876 and DetailedRange 1877–2427. Truncated output was reread.
- **SHA-256:** `053ebf01ed36d6d98b5732600f761dd15b74a9c4c2edc7874e05a422cb8da4fa`
- **Persona:** Read the complete `/home/ahurtado/.codex/plugins/cache/egm/egm/0.1.1/agents/comprehension.md` before the document; applied Mode REVIEW.
- **Context isolation:** No linked files, source code, README files, prior reviews, private progress files, conversation history or memory were consulted. No external standards or reading list were used.
- **Limitations:** This was a document-only comprehension review. It does not verify source accuracy, scientific correctness, numerical choices or implementation readiness. No edits, report files, builds, numerical runs or subagents were used.
