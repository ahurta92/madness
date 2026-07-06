# 18 — Context brief: R2 export / visualization agent

Status: **handoff brief, 2026-06-09.** Onboarding context for an agent starting
R2 (the L4 export / viz / ML layer). Read this + `16_architecture.md` (L4),
`14_property_architecture.md` (Tier A/B), `13_unified_persistence_schema.md`
(how state lives on disk). The harness is `cm.sh` (its README). Start with
`00_status.md` and the molresponse_v3 `CLAUDE.md` for the standing contracts.

## What R2 is (and is NOT)

R2 = **L4**: export converged response/excited-state data for **visualization**
and **ML** (the near-term goal is multiresolution-ML models targeting excited
states). It is a **Tier-B sibling** (doc 14): post-processing that **consumes
already-converged states from persistence and NEVER triggers a solve**. It runs
off the critical path, after property assembly. It is NOT a solver, scheduler, or
property-tensor assembler.

## The contract it plugs into

The single entry point is `run_response(world, ResponseWorkflowInput)` →
`ResponseWorkflowOutput` (`orchestrator/response_workflow.hpp`). The Output has a
**reserved `exports` JSON slot** (currently `{}`). R2 adds an `ExportManager`
stage in `run_response`, after assembly:
```
out.exports = ExportManager(world, in.settings, in.plan, calc_dir).run(...);
```
`exports` is a **manifest** of artifacts written — one descriptor per artifact:
`{object, state_id, omega, protocol_key, target, path, grid_spec?}` — so every
file is reproducible and joinable to `response_metadata.json`. Export is
**additive**: a knob turns it on; off changes nothing in the solve.

## How to GET the data (load, don't solve)

Converged states live on disk as MADNESS archives, keyed by **state identity +
`protocol_key`** in `response_metadata.json` (doc 13). Load them with the
existing persistence calls — all **collective, call on every rank**:
- FD points: `try_load_fd_state<Type,Shell>(world, calc_dir, pert, freq)` →
  `ResponseStateXY` (`solvers/fd_save_load.hpp`).
- ES bundles: `try_load_es_bundle<Type,Shell>(world, calc_dir)` → roots
  (`solvers/es_save_load.hpp`).
- VBC sources: `load_vbc<Shell>(...)` (`solvers/vbc_save_load.hpp`).
- Ground orbitals: `GroundState` (`orbitals_alpha()`, `V_local()`, …).

Densities are built from loaded states (the kernels have `compute_density`):
- closed-shell response density `ρ⁽¹⁾ = Σ_i (φ_i·x_i + φ_i·y_i)` (×2 spin),
- ES transition density `ρ_{0n} = Σ_i φ_i · X_i^{(n)}`.

## What to export (physics objects)

- **Excited states (ML PRIORITY, now unblocked):** ES eigenvectors `X` (and `Y`),
  transition densities `ρ_{0n}`, optionally NTOs. The ES guess now converges to
  the *correct* states (doc 17 Path A, `ESGuessMode::VirtualAO`), so ES-density
  export is finally viable.
- **Linear response:** `ρ⁽¹⁾` + response orbitals `x`(`,y`) per `(pert, ω, protocol)`.
- **Ground:** occupied orbitals (ML reference).

## The four output targets (all on MADNESS primitives — nothing exotic needed)

| Target | MADNESS API | Consumer | Fidelity |
|---|---|---|---|
| **A. Native MRA coefficients** | `Function::save(world, path)` (parallel archive) + a JSON descriptor | multiresolution ES-ML (the differentiator) | full |
| **B. Sampled grids → numpy/json** | regular-grid sampler / `madness::plot_line` (`mra/funcplot.h`) → `.npy`/`.json` | gecko + Python ML | lossy |
| **C. Volumetric cube/DX/VTK** | `plotvtk_begin/plotvtk_data/plotvtk_end`, `plotdx` (`mra/funcplot.h`) | VMD / ParaView | viz |
| **D. Structured property/state JSON** | generalize `solvers/es_analysis.hpp` (ω, transition dipole, oscillator strength, transition quadrupole, dominant-occ weights, Mulliken pop) | analysis + ML features | scalar/tensor |

Suggested order: **A (ES transition densities)** + **D (generalize `es_analysis`)**
first (ML priority), then C, then B. funcplot.h is the one-stop plotting header.

## Contracts / gotchas (do NOT regress)

- **Collective ops on all ranks.** `Function::save`, `norm2`, `inner`,
  `matrix_inner`, density builds are collective — run on every rank; gate only
  *printing* on `if (world.rank()==0)`. A collective inside a rank guard →
  `MPI_ERR_TRUNCATE`.
- **Closed-shell ES only** (open-shell ES is out of scope).
- **Never write `response_metadata.json` directly** — go through the metadata
  layer (`solvers/response_metadata.hpp`). The `exports` manifest goes in the
  `run_response` Output (and optionally an `exports/` metadata subtree via that layer).
- **`es_analysis.hpp` (target D) is built but DISABLED** behind `if(false)` in the
  ES solve path: there is a **heap-OOB in the ES-bundle load path**
  (`load_es_roots` / np-mismatch `ParallelInputArchive`) that aborts at teardown
  when an ES bundle is loaded standalone (repro: `--es-analyze-only --es-load-only`).
  Loading ES bundles for export hits the same path — **expect to fix that load-path
  bug** to ship ES export. (Memory: `project_v3_es_analysis_parked`.)
- **`calc/calc_executor.hpp` is the all-workstream conflict hotspot.** Wiring the
  export stage near the ES/FD solve touches it — coordinate / use a git worktree
  to avoid colliding with other agents.
- **Don't solve.** If producing something needs a new electronic solve, it's a
  Tier-A property (planner + calc manager), not an export.

## Build / test / division of labor

- Build only the cm targets (`cm_build`); the export is additive, so validate via
  a run that produces states then checks the exported artifacts + the `exports`
  manifest. `cm.sh` README is the command catalog.
- **Division of labor:** propose a concrete diff + `cm_build` (agent); the **user
  runs the MRA solver** on the full 64/8-node allocation and verifies. Multi-node
  is NOT agent-drivable here (OMPI has no SLURM PMI) — see the
  `reference_run_on_allocation` memory. Agent may self-serve single-node smoke
  via `mpirun --mca plm isolated --oversubscribe -np N`.
- Commit with `git -c core.hooksPath=/dev/null` (a repo hook corrupts `.git/index`).
  Branch `molresponse-feature-next`.

## Key files

`orchestrator/response_workflow.hpp` (the seam + `Output.exports`),
`solvers/{fd,es,vbc}_save_load.hpp` (load states), `solvers/es_analysis.hpp`
(transition properties — target D, currently disabled), `kernels/*` (`compute_density`),
`solvers/response_state.hpp` + `ResponseFunctions.hpp` (state types),
`src/madness/mra/funcplot.h` (plotvtk/plotdx/plot_line), `Function::save/load`.
