# Simple restart rules implementation plan

> **For agentic workers:** Use superpowers:executing-plans; execute inline and obtain one independent final review.

**Goal:** Apply the agreed model/numerical restart distinction to shared SCF and driven response.
**Architecture:** Compare model projections of the existing scientific snapshots before accepting native functions. Keep full snapshots for convergence-evidence qualification. Response may continue against a refined, same-model SCF source with fresh evidence and existing per-state protocol planning.
**Tech Stack:** C++, native MADNESS archives, existing CTest/Python drivers.
**Spec:** docs/design/madqc-response-after-scf-refinement.md and the user's approved restart table.

## Milestone placement

The parent response plan governs validation: focused interface checks during integration, then metric and representative converged-value checks in task 8. Historical iteration-count agreement and reference recomputation are not first-milestone requirements; broader reference/performance calibration is later validation work. The verification list below is an inventory, not a requirement to run the full matrix in this increment.

The [full response plan](2026-11-19-full-response-unification.md) now sequences this work as task 1. Reuse this plan's existing regression edits and evidence; do not restart the increment from scratch.

This is the first increment of [complete molresponse unification](../../design/madqc-response-milestone.md).
After driven response, cover coupled excited roots and quadratic-response tasks
and properties. Shared timing conventions are part of the milestone alongside
metrics, convergence, restart, checkpoints, results and printing.

## Global constraints
- Existing worktree and build; preserve unrelated dirty changes. No commit/push or reference regeneration.
- No SCF orbital-preparation changes, automatic orbital matching, new input controls, checkpoint-format redesign, or optimizer restart.
- Model changes reject ordinary restart. Explicit imported seeds remain application-owned.
- Optimizer integration and new live-geometry seeding behavior are deferred. Do not make an optimizer change a prerequisite for this increment; record any existing geometry-workflow interactions separately.
- Fresh checks after operator/ground-state changes; missing history stays unavailable; preserve per-state protocol progress and fresh maxiter allowance.
- Native archives, SafeMPI/World collectives; serial CTest with two MAD workers. Scope numerical validation to affected workflows per developer instruction.

## Review focus
- Geometry/field/XC/localization and unknown-setting mismatches must not silently become guesses.
- Numerical smoothing, operator accuracy and XC screening changes must retain usable functions but invalidate old completion.
- Changed SCF sources must not reuse old response properties or history; invalid requests must not replace response provenance.
- Restart preserves each independent state's saved protocol and qualification, including mixed progress.
- Explicit seed loading must remain separate from ordinary restart. Record geometry-workflow compatibility findings without expanding into optimizer implementation.

### Task 1: SCF model boundary
**Files:** SCFScientificSettings.h/.cc, SCFCheckpoint.cc, SCFRun.cc, existing SCF settings/checkpoint/run tests and public shared SCF test.
**Interfaces:** Add scf_model_settings(snapshot) returning the model-only JSON projection; unknown fields stay model-significant. Full snapshots remain unchanged for evidence comparison.
- [ ] Add changed-geometry rejection assertions to existing native/public tests; run and observe rejection tests fail against existing behavior.
- [ ] Implement explicit numerical exclusions (nuclear/core smoothing, lo, derivative discretization, orbital shift, XC screening), preserving physical fields, occupations, localization, functional coefficients and unknown settings. Test projections using real configured SCF snapshots.
- [ ] Require model equality in ordinary SCF planning. Preserve explicit imported-seed loading; do not add optimizer or live-geometry APIs in this increment.
- [ ] Run settings/checkpoint/run/public tests, serial and relevant two-rank tests. Expected: all pass, ordinary mismatches retain selected payload.

### Task 2: Response continuation
**Files:** shared_fd_checkpoint.hpp, new test_shared_response_refinement.py, madqc CMake registration and related rejection tests.
**Interfaces:** Consume scf_model_settings; existing changed context selects guess-only planning and fresh comparisons.
- [ ] Add real two-state response workflow covering changed SCF source, fresh zero-update checks, repeated reuse, same-model numerical settings and model rejection; observe existing changed-source rejection.
- [ ] Validate saved/current ground models and response identity/layout before registering changed inputs. Allow valid changed sources; preserve exact recorded provenance and existing conversion/planner behavior.
- [ ] Run serial/two-rank refinement plus existing response, ladder and rejection tests. Expected: original-source metrics never qualify changed equations; no lost state progress.

### Task 3: Verification and review
- [ ] Run affected CTest suite and unchanged He HF, water HF and He alpha numerical references. Save results and reference hashes.
- [ ] Run clang-tidy on changed translation units; exclude the previously identified crashing argument-comment check if necessary.
- [ ] One fresh independent review of this increment against its before snapshots; fix important findings with regressions, then rerun affected checks.
- [ ] Update design/status/vault with actual evidence, scoped patch and remaining limitations. Keep existing checkpoint generations/retention implementation unchanged.
