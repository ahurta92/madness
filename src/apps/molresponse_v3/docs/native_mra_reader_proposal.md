# Proposal: native MRA coefficient export → native ParaView reader

**Branch:** `feat/amr-export` (this branch owns the `dump_mra_trees --amr/--htg`
Overlapping-AMR export).
**Depends on:** **`io-hdf5`** — pull it before implementing (see "What to pull").
**Status:** proposal / design note. No code here yet.

## Why
The `--amr` export (`be7ebc9e4` prototype, `933261a6f` VTK 9.1 build fix) writes
a `vtkOverlappingAMR` (`.vthb` + per-leaf `.vti` blocks). It renders beautiful
adaptive isosurfaces in the ParaView GUI — but it is a **point-sampled** image of
each leaf box, not the function's native multiresolution representation. Every
export we have is lossy or metadata-only:

| Export | Branch | Representation | Limitation |
|---|---|---|---|
| `_recon.json` / `_comp.json` | base | per-node **metadata only** (`has_coeff` 0/1, `norm`, `snorm`, `dnorm`) | **no coefficient values** — cannot evaluate offline |
| `.cube` | molresponse-feature-next | dense uniform grid sample | adaptivity discarded; large; fixed resolution |
| `.vthb` (this branch) | feat/amr-export | each leaf **sampled** onto a uniform sub-block | smooth + adaptive, but still a *sample*, not the basis |

A **native MRA reader plugin** would instead evaluate the MADNESS expansion on
demand (any point, any zoom) from the stored coefficients — exact, compact, no
resampling, no MADNESS at view time.

## What the reader needs, per leaf box
- **`n`** — refinement level.
- **`l = (lx, ly, lz)`** — translation index (box identity / Morton position at
  level `n`); with `n` it fixes the box's physical extent in the cell.
- **`s[k][k][k]`** — the **scaling-coefficient tensor** at that leaf
  (reconstructed form). `k` = polynomial/wavelet order (k6 here).
- Optionally **`d`** (wavelet/difference coeffs) for the compressed tree, plus a
  header (`cell`, `k`, `thresh`).

`dump_mra_trees` already walks every node (that's where `has_coeff` is emitted)
and the leaf `FunctionNode` already holds its `coeff()` tensor — so this is
"serialize the values we already touch," not new tree traversal.

## What to pull from `io-hdf5`  ← the point of this note
**Do not roll a new serializer.** `io-hdf5` already built the HDF5 Function-I/O
layer this needs:
- `c8b8f994c` — *opt-in HDF5 Function I/O — structured + archive backend* [P2 Layer A]
- `cb8fc38cc`, `5ab0525df` — opt-in HDF5 restart wired into `ResponseStateX/XY<ClosedShell>`
- `98d7c3041` — `test_function_hdf5` round-trip (legacy vs archive-HDF5, exact)
- HEAD of `io-hdf5` at time of writing: **`5ab0525df`**

That backend already serializes MADNESS `Function` data through HDF5 with a
structured layout. The native-coefficient dump should write per-leaf `(n, l, s)`
**through that same writer**, so the `.h5` it produces is the durable native
artifact the reader consumes.

**Action:** merge `io-hdf5` into `feat/amr-export` (or cherry-pick the P2 Layer A
HDF5-IO commits) first:
```
git checkout feat/amr-export
git merge io-hdf5            # brings in the HDF5 Function-I/O backend
```
Then add a `--coeffs` (native) path to `dump_mra_trees` that emits per-leaf
`(n, l, s[k^3], k)` via the io-hdf5 HDF5 writer, alongside the existing `--amr`
sampled path.

## Suggested split
- **`io-hdf5`** = the bulk coefficient payload + HDF5 layout (one `.h5` per
  function; group per level; datasets keyed by `l` holding `k^3` tensors; header
  `cell`, `k`, `thresh`). It is the native dump the reader loads.
- **`feat/amr-export`** = either (a) carry coefficients as field data on the
  existing `.vthb` blocks so one file serves sampled + native, or (b) keep
  `.vthb` as the quick-look sampled export and ship a native VTK reader that
  targets the io-hdf5 `.h5` and reconstructs on load (can reuse this branch's
  AMR machinery to build a `vtkOverlappingAMR` from the coefficients).

## Open questions
1. Reconstructed `s` only, or also compressed `d` (multi-scale / progressive load)?
2. HDF5 layout: per-level group vs flat `(n,l)`-keyed table? chunking/compression?
3. Reader: reconstruct to `vtkOverlappingAMR` on load (reuse feat/amr), or a true
   on-demand evaluator?
4. Precision: store `double` `k^3` as-is, or truncate by `thresh`?
