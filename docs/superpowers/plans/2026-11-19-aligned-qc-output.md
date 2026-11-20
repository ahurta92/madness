# Aligned chemistry output implementation plan

> **For agentic workers:** Use superpowers:executing-plans for this approved increment.

**Goal:** Use common aligned metric, timing and restart tables for adopted SCF and driven-response paths, preserving scientific detail.

**Architecture:** Applications supply labels and already-computed measurements. IterationOutput.h owns formatting only; no numerical evaluation, MPI calls or new synchronization. Model-specific reports and existing runtime timers remain available.

**Tech Stack:** C++20, ostream, steady_clock, Python public workflow tests, Ninja/CTest.

**Spec:** docs/design/madqc-stdout-preservation.md and the user's approved console example, including aligned columns and consistent application output where possible.

## Global constraints
- Do not change solver decisions, defaults, checkpoint formats, retention, numerical references or output destinations.
- Use measured diagnostics even when unselected, clearly marked reported only. Missing history is unavailable, never zero. Driven kinetic is amplitude-dependent, not ordinary ground-state kinetic energy.
- Show actual effective limits; distinguish intermediate protocol completion from final convergence.
- Only report successful checkpoint saving after publication returns. Reuse reports saved evidence and current load time, never historical timings or invented iterations.
- Keep rank-zero and existing print-level gates. New timing uses steady_clock with no added fences; existing timers remain unchanged.
- Existing dirty isolated worktree is the authorized workspace. No commits, resets, unrelated cleanups or reference regeneration.

## Review focus
1. Long state/task labels and extreme numbers must not shift columns; unit tests cover both.
2. Unavailable/not-applicable/invalid or unselected measurements must not look converged; unit tests cover these.
3. Reuse and tightened restarts must describe saved evidence honestly without new evaluation; public SCF/response tests verify this.
4. Quiet output, relative tolerances and multi-state grouping remain correct; targeted unit and public tests cover these.
5. Final orbital preparation must keep saved-state metrics and show existing work accurately; preserve numerical tests and inspect print-level 3 smoke output.

### Task 1: Shared tables
**Files:** IterationOutput.h and molresponse/tests/test_shared_iteration.cpp.
**Interfaces:** optional MetricLabels supplied by applications; common numeric tables; metric-table routine reused by iteration and restart display; no schema additions.
- [x] Add behavior tests for column positions, unselected diagnostics, missing/inapplicable data, relative limits, stream flags and timing scope; observe failure.
- [x] Implement the shared printers and readable decision/timing labels; run the unit test.

### Task 2: Application integration
**Files:** SCF.h/cc, SCFCheckpoint.h/cc, molresponse/solvers/shared_fd_checkpoint.hpp, calc/calc_executor.hpp, madqc/test_shared_scf.py and test_shared_response.py.
**Interfaces:** app labels include density and kinetic definitions; final-protocol flag; existing assessment of saved metrics reused only when valid; restart summary follows successful loading.
- [x] Update public assertions before integration. SCF energy table uses the shared value-table layout; drop the stale legacy conv banner and redundant metric line on shared paths.
- [x] Restore elapsed intervals around final overlap/eigensolve/transformation work without changing legacy clocks or adding fences; identify final reevaluation phase.
- [x] Wire common tables and post-publication notices in SCF and driven response; display saved metrics/current limits on reuse or continuation when comparable.
- [x] Build and run public serial/restart/two-rank and curated helium numerical tests. Do not regenerate references.

### Task 3: Verification and documentation
**Files:** stdout preservation audit, madqc README, active vault session, this increment's evidence/ledger.
- [x] Run changed-TU clang-tidy and inspect real console output; record practical coverage and remaining legacy engines.
- [x] Produce a scoped diff against the pre-increment snapshot; obtain one fresh independent final review, fix important findings once and verify affected checks.
- [x] Record results and show a concise real output excerpt.

Plan self-review: reporting-only boundaries are explicit; no engine-specific stopping logic moves into formatting. Task 2 consumes Task 1's optional-label interface; Task 3 consumes verified real output. No conflicting producers.


Verification: 14 distinct checks pass (formatter/controller/SCF evidence tests,
public SCF and response in serial/two ranks, static response, five numerical
reference cases). All 74 reference files retain their original hashes. One
independent review, one fix pass; no new changed-line lint diagnostics. Evidence:
`.superpowers/sdd/2026-11-19-aligned-qc-output/verification.json`. Existing legacy
engine migration and controlled performance validation remain outside this
reporting increment. Changes remain uncommitted in the ongoing shared-core branch.
