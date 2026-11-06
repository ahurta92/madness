# Fixed intermediate work and resumable final checkpoints Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to execute this plan inline.

**Goal:** Deliver the first bounded implementation slice of the accepted madqc iteration/restart policy: fixed intermediate work and a fresh attempt for an unconverged final checkpoint.

**Architecture:** Extend the shared rung sample with an explicit final-protocol flag. Pass that flag through the response drivers and from the named SCF ladder. Keep the numerical engines responsible for measurements and divergence. Remove the planner's special refusal to resume a stalled final state. Preserve existing records, protocol scheduling, fresh per-run allowances, and checkpoint publication.

**Tech Stack:** C++20, pure header tests, existing MADNESS CMake/Ninja build and named QC decks.

**Spec:** `docs/design/madqc-simplified-policy.md`, based on the accepted `/home/ahurtado/repos/concepts_brief.md` and the drop-6 alignment review.

## Global Constraints

- Baseline: madqc-drop-6, 9179f5d32076a2ba555906dd146dc330ffc040d7. Work only in the isolated codex/madqc-simple-policy worktree.
- No numerical kernel, metric definition, binary archive format, or release-checkout change. Additive progress metadata distinguishes fixed intermediate work from older records.
- Intermediate work may terminate on divergence; convergence and plateau observations cannot end it before its configured count. At the boundary, report the actual convergence verdict; an unmet count limit is not convergence.
- Explicit diagnostic iteration counts keep their existing precedence.
- Final-stage criteria selection, removal of excited-state freezing, changed-setup validity, expanded roots, refined-reference reuse and persistence/history auditing are subsequent slices, explicitly tracked in the amendment.
- Existing final-stage stopping behavior remains in this slice; a stalled final saved state is allowed another attempt on the next invocation.

## Review Focus

- The response scheduler invokes one threshold at a time: role must come from the node's entire requested protocol list, not the one-element executor vector.
- A fresh count is per invocation, independent of state/history iteration counters.
- Final protocol remains mandatory; divergent states never advance.
- Check both TDA and full excited-state adapters, FD and SCF. Preserve legacy four-argument solver calls where practical with a default final-protocol argument.

### Task 1: Shared policy and restart planner

**Files:** `src/madness/chem/ConvergenceSettings.h`, `LadderPlanner.h`, `test_convergence_settings.cc`, `test_plan_ladder.cc`.

- [x] Establish the baseline by compiling/running the two existing pure tests against this checkout.
- [x] Add cases for intermediate convergence/plateau before the cap, met/unmet verdicts at the cap, diagnostic override, cap one, and stalled final resume with recorded historical counts.
- [x] Add `RungSample::final_protocol` defaulting true and first observe behavioral test failures, then implement fixed intermediate exit ordering and remove stalled-final load-only handling.
- [x] Run both tests, inspect diff, commit the shared policy change.

### Task 2: Wire protocol roles through solver drivers

**Files:** `src/apps/molresponse/solvers/{iterate.hpp,iterate_protocol.hpp,fd_solver.hpp,es_solver.hpp}`, `src/apps/molresponse/calc/calc_executor.hpp`, `src/madness/chem/SCF.cc`, new `src/apps/molresponse/tests/test_iteration_protocol.cpp`, `src/apps/molresponse/CMakeLists.txt`. Progress evidence also requires `fd_explicit_accuracy.hpp`, both named checkpoint adapters, their tests and `named-stages.md`. The all-locked ES branch must retain eigenvalue-change evidence; the H2 acceptance test covers both TDA and RPA.

- [x] Add pure driver scenario tests: early convergence/plateau at intermediate stage, mandatory final stage, single-threshold scheduler call marked intermediate, fresh allowance from nonzero history, divergence, one-iteration cap, diagnostic override.
- [x] Add explicit final-protocol fields to driver policies and optional final-protocol arguments to solver stop_reason methods; run tests to observe wrong early exit before implementing forwarding.
- [x] Wire the complete requested node ladder into the response executor's role flag; named SCF uses its existing supplied-ladder lookup. Preserve reported convergence at the completed intermediate boundary.
- [x] Register the pure test in CMake, run focused tests and compile the changed production paths. Commit adapters/tests.

### Task 3: Validate and record the slice

**Files:** `docs/design/madqc-simplified-policy.md`, this plan, validation artifacts under the ignored plan workspace.

- [x] Build affected production targets in isolation (reuse installed dependencies). Run focused CTest and small named-stage scenarios if build/runtime resources permit; record any limits honestly.
- [x] Run clang-tidy on changed translation units; distinguish pre-existing findings.
- [x] Request one fresh whole-branch review as required by executing-plans. Fix substantive findings and rerun affected checks.
- [x] Update the amendment with actual evidence and remaining alignment work; leave the branch reviewable without changing the other agent's drop or vault history.
