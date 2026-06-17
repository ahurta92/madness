# `dump_mra_trees --amr` — Overlapping-AMR function-value export (PROTOTYPE)

**Status:** prototype on branch `feat/amr-export` (off `feat/htg-writer`). Not yet
built/validated — needs a `-DMADNESS_ENABLE_VTK=ON` build on the cluster and a
ParaView pass. This note pins the design so the build shake-out + iteration start
from a documented base.

## Why a third export format

Three complementary views of one MADNESS `Function`:

| format | what it stores | isosurfaces | scales? |
|---|---|---|---|
| `.cube` (have) | function on a **dense uniform grid** | smooth | no — uniform fine grid everywhere |
| `.htg` (have) | **one scalar per leaf box** (skeleton: level/norm/rms/owner) | blocky (cell-resolution) | yes (compact) |
| **`.vthb` (this)** | function sampled on an **m³ grid per leaf box, across levels** | **smooth + adaptive** | yes (fine blocks only where refined) |

`.vthb` (`vtkOverlappingAMR`) is the adaptive replacement for the dense cube: ParaView's
AMR Contour gives smooth isosurfaces at sub-box resolution, but only the refined regions
carry fine blocks — so it stays small for large molecules where a uniform cube blows up.
This is the natural home for the "k³ quadrature values per box" idea (HTG, being one-value-
per-cell, cannot hold sub-box grids).

## MADNESS → OverlappingAMR mapping (the load-bearing logic — verified)

- **Leaf walk:** rank 0 walks the reconstructed tree (`collect_live_tree_r`), keeps leaves
  (`!has_children`), broadcasts the flat `(level,lx,ly,lz)` list to all ranks so the
  sampling loop is collective + in lockstep.
- **Per-leaf sampling:** every rank calls `f.eval_cube(box_cell, {m,m,m})` for each leaf
  (collective; rank 0 receives the m³ grid). `box_cell` (bohr) = `[lo+l·h, lo+(l+1)·h]`,
  `h = cell_width / 2^level` per axis.
- **AMR assembly (rank 0):** refinement ratio **2** (octree). Global origin = `cell_lo`.
  Level-n cell spacing = `(cell_width/2^n)/m`. Each leaf box → one `vtkUniformGrid` block
  (origin = box lower corner bohr, m³ cells) at AMR level n, with index-space
  `vtkAMRBox` lo = `l·m`, hi = `l·m + m − 1`. Values attached as cell data `"psi"`.
- Output `<stem>.vthb` via `vtkXMLUniformGridAMRWriter`. Ground MOs only in this
  prototype (extending to FD response orbitals / ρ⁽¹⁾ is the same threading as
  `--cube`/`--htg`).

## Known prototype simplifications (intentional — fix during iteration)

1. **Per-box `eval_cube` is O(#leaves) collective calls** — fine for h2o/c2h4, slow at
   scale (each call fences). Production path: transform each leaf's scaling coefficients
   to point values **locally** (no per-box collective) — `node.coeff()` → values via the
   quadrature transform. Prototype favors the simple `eval_cube` reuse.
2. **Cell-vs-point half-cell offset:** the m³ `eval_cube` samples (grid points, inclusive
   endpoints) are attached as **cell** data on an m-cell block — a half-cell offset.
   Acceptable for a first look; tighten by sampling at cell centers (offset the cell by
   h/2m) or by using point data with matching `vtkAMRBox` dims.
3. **Blanking** (coarse cells covered by finer blocks) is left to ParaView's AMR overlap
   resolution from the box geometry; we do not pre-blank. Verify it renders correctly;
   if not, call `vtkOverlappingAMR::GenerateParentChildInformation()` / set blanking.
4. **Ground-state only** — FD response/ρ⁽¹⁾ not yet wired (trivial follow-up).

## VTK API to verify at the first build (version-sensitive)

The MADNESS-side mapping is solid; these VTK calls are a best-effort first cut:
- `vtkOverlappingAMR::Initialize(int nLevels, const int* blocksPerLevel)` — confirm signature.
- `SetGridDescription(VTK_XYZ_GRID)` — the constant may need `#include <vtkStructuredData.h>`.
- `SetSpacing(level, double[3])`, `SetRefinementRatio(level, ratio)` — confirm they exist on
  `vtkOverlappingAMR` in the cluster VTK (some versions route via the AMR metadata object).
- `SetAMRBox(level, id, vtkAMRBox)`, `SetDataSet(level, id, vtkUniformGrid*)` — confirm.
- `vtkXMLUniformGridAMRWriter` writes `.vthb` — confirm; alternatively the HDF AMR writer.
- **CMake components:** the existing `find_package(VTK COMPONENTS CommonCore CommonDataModel
  IOXML)` may need **`FiltersAMR`** (and possibly `IOAMR`) added for the AMR classes/writer.
  Update the `MADNESS_ENABLE_VTK` find_package in `CMakeLists.txt` if the link fails.

## Build + validate (on the alloc)

```bash
cmake -B build -DMADNESS_ENABLE_VTK=ON -DVTK_DIR=$VTK_DIR <usual flags>
cmake --build build --target dump_mra_trees
mpirun -np 4 dump_mra_trees --archive=.../h2o/mad.restartdata --out=amr --amr --amr-m=8
```
Then in ParaView: open `amr/mo_0.vthb` → **Contour** on `psi` → expect a **smooth**
isosurface (vs the blocky `.htg`), with **finer blocks near the nuclei**. Cross-check the
lobe shape against `mo_0.cube`'s isosurface (should match) and compare file size +
ParaView load time vs the `.cube` (the adaptivity payoff grows with molecule size).
```
