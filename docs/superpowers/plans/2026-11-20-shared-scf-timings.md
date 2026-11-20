# Shared SCF timing publication implementation plan

**Status — November 20, 2026:** implemented and verified. Six native checks passed, followed by five final serial/MPI/public checks and the existing response wiring check. The one important review finding was reproduced and fixed. No solver equations, update order, references or archive locations changed.

> **For agentic workers:** Use superpowers:executing-plans inline; one independent review of this increment.

**Goal:** Save SCF and response timings together using the existing shared timing schema, including fresh, reused and failed calculations.
**Architecture:** Forward the existing SCF iteration measurements to an optional application-owned observer. SCFRun uses the existing steady-clock collector and checkpoint store to save calculation/protocol/task intervals in the same timings.json that response already updates. Each calculation keeps its own existing invocation identity; no new workflow identity or numerical work is introduced.
**Tech Stack:** C++, MADNESS World/native stores, JSON, CMake/Ninja and Python/CTest.
**Spec:** design_draft.md, Timing breakdown; docs/design/madqc-response-milestone.md, Shared timing conventions. The user's request to continue authorizes closing this observed reporting gap.

## Global constraints
- Reuse the current worktree/build and preserve all other uncommitted edits. No commit, push or reference regeneration.
- Reuse existing measured iteration intervals, including the combined prepare/evaluate/decide interval. No new fences, numerical kernels, metric evaluations, or changed iteration budgets.
- Keep earlier timing records on restart; a reused calculation records only current loading/preparation/publication work, never historical iterations.
- Keep collective publication/error handling and preserve the original failure if writing failure timings also fails.
- Retain current archive paths/IDs. Relocating shared files from task_0/moldft requires coordinated consumer/test migration and is the next separate increment.

## Review focus
- Both ordinary and frozen-virtual SCF must forward timing measurements; final checks have no advance step.
- Quiet runs must still persist timing data.
- SCF reuse and response publication must retain earlier entries without copying historical work.
- Preparation and protocol changes must have correct scope and membership; frozen occupied qualification remains separate.
- Failure publication must be collective, retain previously valid checkpoints/records and rethrow the original exception.

### Task 1: publish SCF measurements
**Files:** src/madness/chem/SCF.h, SCF.cc, SCFRun.cc, test_scf_run.cc; existing SCF iteration/frozen tests as needed.
**Interfaces:** optional IterationControl timing observer receives existing IterationInfo and IterationTiming after a completed evaluation/update interval. SCFRun owns the collector; CheckpointStore::publish_timings appends/upserts the existing record format.
- [x] Extend native SCFRun tests for fresh terminal timing, quiet runs, reuse without iterations, failure records and preserved earlier entries; observe the missing timing file failure.
- [x] Forward existing measurements on ordinary and frozen paths; collect setup, protocol/task, iteration, load and result-publication timings in SCFRun, including partial failure scopes.
- [x] Run native tests in serial and two ranks; test source failure with valid retained native checkpoint.

### Task 2: workflow evidence and documentation
**Files:** src/apps/madqc/test_shared_workflow_timings.py, CMakeLists.txt; docs/design/madqc-test-index.md, docs/design/madqc-stdout-preservation.md.
**Interfaces:** public SCF-to-response run uses one timings.json with ground_state and response records. Console iteration values equal saved values to display precision.
- [x] Add a public fresh/reuse/tighten fixture that requires both engines' records, preserved history, no reused iterations, actual step timings and no update after terminal evaluation. Observe the original binary fail.
- [x] Build current madqc and affected native targets; run focused existing SCF/restart/output/response checks plus the new serial/MPI fixture. Do not rerun the scientific matrix for reporting-only changes.
- [x] Run clang-tidy on changed translation units (exclude known crashing argument-comment check); one independent scoped review and focused fixes if needed.
- [x] Update design/test index/vault with actual evidence and clearly identify layout migration and cold-water root selection as separate remaining work.
