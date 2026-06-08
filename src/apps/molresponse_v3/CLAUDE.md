# CLAUDE.md — molresponse_v3 (active response module)

You are working in `molresponse_v3`, the ground-up rebuild of the MADNESS
response pipeline and the **current active front** of response development
(v2 = previous production, `molresponse_legacy` = frozen reference).

**Read first, in this order — do not re-explore the tree to rediscover this:**
1. `docs/00_status.md` — active workstreams, hot-file conflict map, parallel-agent
   protocol, standing contracts. The anti-context-loss page.
2. `/gpfs/scratch/ahurtado/madness_es_bench/README.md` — the `cm.sh` build/run/
   validate harness (the canonical way to test v3).
3. The design docs in `docs/` (numbered): `06_refactor_plan`, `07_increment_plan`,
   `15_calc_manager_design`, `13_unified_persistence_schema`, `TPA_DESIGN`.

## Build / run / validate loop
- On the allocation, via `cm.sh` (see its README):
  `cm_build` → `cm_unit` → the workstream's `cm_*` command (`cm_es`, `cm_beta`,
  `cm_check_raman`, `cm_idem`, …) → `cm_record <calc_dir> <mol>`.
- Run MRA solves on a compute node via the `run-on-allocation` skill, never the
  login node. The no-MPI `cm_unit` / `test_calc_manager` runs anywhere.
- Build environment: GCC 13.2.0 — run binaries with its libstdc++ on
  `LD_LIBRARY_PATH` (`/gpfs/software/gcc/13.2.0/lib64`).
- New `test_*.cpp` → register in `CMakeLists.txt` AND `cm_build`'s ninja line.

## Hard contracts (do not regress)
- **Propose a concrete diff + get approval** before non-trivial solver/runtime
  edits.
- **Collective ops** (`Function::size()/norm2()/trace()`, `inner`, `matrix_inner`)
  run on **all ranks** — gate only printing on `if (world.rank()==0)`; a collective
  inside a rank guard → `MPI_ERR_TRUNCATE`.
- **Never write `response_metadata.json` directly** — go through the metadata
  layer (`solvers/response_metadata.hpp`).
- **ES is closed-shell only.** Open-shell ES + Full-RPA open-shell are out of
  scope. ES / diffuse-state Gaussian references use **d-aug-cc-pVQZ**
  (`madness_studies/refs/dalton_tdhf.json`); single-augmentation manufactures a
  phantom ~3% error.
- After any `kernels/` (esp. `two_electron.hpp`) edit, run `cm_equiv`.
- Commit with `git -c core.hooksPath=/dev/null` (a repo hook corrupts
  `.git/index`). Branch: `molresponse-feature-next`.

## Parallelism note
`calc/calc_executor.hpp` is touched by every active workstream — highest merge-
conflict risk. When running parallel agents, give writers of shared files their
own git worktree (`isolation: "worktree"`). See `docs/00_status.md` for the map.
