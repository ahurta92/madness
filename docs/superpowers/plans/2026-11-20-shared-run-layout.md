# Shared run directory layout

**Status:** implemented and verified, November 20, 2026. Twenty focused CTest checks pass, including serial/MPI layout, restart, seeding, subgroup/VBC, native inspection, 2PA, coupled roots and timing checks. The independent review finding has an artifact-reader RED/GREEN check; all 74 reference files are unchanged.

Execute inline using superpowers:executing-plans; request one independent review of this increment after verification.

The agreed SCF/response layout puts results.json, timings.json and checkpoints/ directly under the run prefix. Diagnostics live in analysis/ground_state and analysis/response. Explicitly separate the shared store root from engine working directories; archive records, identities, retention and numerical work stay unchanged. Preserve existing worktree edits and numerical references. Optimizer and unadopted workflows retain their existing paths.

- [x] Add a public fresh/reuse/continue/old-layout recovery fixture; observe the current binary fail the new layout assertion.
- [x] Thread an absolute shared run root through the single-point workflow, SCF and response adapters. Update seed detection, inspection/continue, response history/shards and plot discovery together.
- [x] Reject old layouts collectively before creating output. Provide an explicit offline development-run relocation tool with preflight, no overwrites and recoverable rename steps; test conflicts/interruption and byte-preserving native relocation. No persistent solver compatibility fallback.
- [x] Migrate active consumers/tests with separate shared versus analysis paths. Run focused serial/MPI restart, seed, coupled/direct response, subgroup, native inspection and timing checks. Do not regenerate numerical references or retest iteration-count parity.
- [x] Lint changed C++ translation units; independently review the scoped diff, address findings, then update design, test index and vault with evidence.

Evidence: .superpowers/sdd/2026-11-20-shared-run-layout/. Reference hashes remain those of the response milestone. Existing default cold-water root-selection gap is unrelated and remains open.
