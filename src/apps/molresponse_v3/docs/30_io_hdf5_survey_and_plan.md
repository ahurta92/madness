# 30 — MADNESS I/O survey + HDF5 prototype plan  (thread: `io-hdf5`)

**Status:** P0 (read-only audit + plan) — this document. No core IO edits yet.
**Mandate:** survey where MADNESS spends I/O (checkpoint/restart, function
load/store, archive), prototype an HDF5-backed path, quantify vs the current
serialization. Legacy archive stays the **default**; HDF5 is **opt-in** until
proven. Coordinate the HDF5 layer with `feat/amr-export` (VTKHDF) — **one shared
HDF5 stack**, not two. (See the `io-hdf5` brief in
`madness_studies/RELEASE_STATUS.md`.)

---

## 1. The legacy I/O stack (what exists today)

### 1.1 Archive core (`src/madness/world/`)
- `archive.h` — `BaseArchive` + `ArchiveStoreImpl`/`ArchiveLoadImpl`; the `&`
  operator dispatches `serialize(ar,t)` with a runtime type cookie (preamble/
  postamble). Every concrete archive is a thin sink/source under this.
- `binary_fstream_archive.h` — `BinaryFstream{Output,Input}Archive`: raw binary
  `store(const T*, long n)` / `load(T*, long n)` to one `std::fstream`. **The
  production on-disk primitive.**
- `text_fstream_archive.h`, `vector_archive.h`, `buffer_archive.h`,
  `mpi_archive.h`, `cereal_archive.h` — text / in-memory / MPI-message / optional
  Cereal variants. Not the persistence path.
- `parallel_archive.h` — `Parallel{Output,Input}Archive<Archive>` with an `nio`
  (number of I/O nodes) field. `io_node(rank)=rank%nio`; opens `nio` files named
  `<stem>.00000`, `<stem>.00001`, …; non-I/O ranks gather to their I/O node.
- `parallel_dc_archive.h` — `ContainerRecord{Output,Input}Archive` bridges a
  local `VectorOutputArchive` to the parallel archive (per-container-record).

### 1.2 Function persistence (`src/madness/mra/`)
- `mra.h:1533/1569` — `Function::load/store`: magic cookie `7776769`, type id,
  `NDIM`, `k`, cell, then delegates to `FunctionImpl`.
- `mra.h:2901/2907` — free `madness::save/load(Function)`:
  `ParallelOutputArchive<BinaryFstreamOutputArchive>` with **`nio=1`**.
- `funcimpl.h:1323/1305` — `FunctionImpl::store/load`: metadata first
  (`k`, thresh, initial/max level, truncate mode, autorefine, **`tree_state`**),
  then `ar & coeffs` — the distributed
  `WorldContainer<Key<NDIM>, FunctionNode<T,NDIM>>`.
- `worlddc.h:2069/2339` — the WorldContainer (de)serialization: each rank
  serializes its local nodes into a `BufferOutputArchive`, sizes go up by
  `MPI_Gather`, buffers by `MPI_Gatherv`, **rank 0 writes one file** with header
  `(magic=-5881828, 1, -magic, count, data…)`. Load mirrors it with
  `MPI_Scatterv`. **`nio=1` ⇒ a single rank-0-gathered file** — exactly the
  model an HDF5 prototype mirrors first.

### 1.3 molresponse_v3 checkpoint/restart call sites (`src/apps/molresponse_v3/`)
| Site | File:line | What it persists | Archive |
|---|---|---|---|
| Response vectors | `solvers/response_state.hpp:153/161` (X), `:210/220` (XY) | `x_alpha`(+`x_beta`/Y) `Function` vectors | `ParallelOutputArchive<BinaryFstream…>` `nio=1` |
| FD point | `solvers/fd_save_load.hpp:90/193` | one FD state (pert×freq×protocol) + JSON meta | collective `.save`/`.load` + rank-0 JSON |
| ES roots | `solvers/es_save_load.hpp:90/~180` | bundle: `root_<s>` archives + `roots.json` | per-root collective `.save`/`.load` |
| VBC source | `solvers/vbc_save_load.hpp:42/86` | quadratic source + JSON meta | collective `.save`/`.load` + rank-0 JSON |
| Metadata | `solvers/response_metadata.hpp:57/133` | all status/timing/identity | **JSON only** (atomic tmp+rename, rank 0) |

**Metadata is JSON and stays JSON** — do not route it through HDF5 (the
"never write `response_metadata.json` directly" contract owns that file).

### 1.4 SCF ground-state restart
`GroundState::from_archive` (the v3 fix in
[[project_v3_validation_baseline]] lineage) resolves the moldft SCF archive from
`scf_calc->work_dir` + prefix `.restartdata` and reloads the MO `Function`s via
the same `ParallelInputArchive` path. So the ground state is the **largest, most
frequently re-read** function payload — a prime HDF5 target.

---

## 2. Existing HDF5 prototype state — `src/examples/writecoeffs/`

What's there, and what it actually does:

| File | Reality |
|---|---|
| `FunctionIO.h` | per-leaf **grid-point values** (`coeffs2values`), **text** round-trip; rank-0 tree-walk via remote `coeffs.find(key).get()`. Validated by `norm2`/error. |
| `FunctionIOHDF5.h` | includes `<h5cpp/h5cpp.hpp>` + adds physical `coords`, **but body is still text/JSON** — no HDF5 dataset write of the tree. |
| `writecoeff{,2,3}.cc`, `*_json.cc`, `h2_write*.cc` | text / nlohmann-JSON drivers of the values format. |
| `writecoeff_hdf5.cc` (in build list) | `#include <FunctionIO.h>` — the **text** path; "hdf5" in name only. |
| `h5cpp_test.cc`, `core.hpp`, `app.cpp` | toy `h5cpp` demos (std::vector / ChunkedDataset). **Not in `EXAMPLE_SOURCES` ⇒ never built.** |

**Provenance:** these files come from `upstream/pr_writecoeff`
(m-a-d-n-e-s-s/madness; already merged into our trunk, ours slightly newer).
That branch adds **no HDF5 build wiring** — only `add_subdirectory(writecoeffs)`
+ an `EXAMPLE_SOURCES` list that *excludes* the h5cpp files. So it gives us a
validated *values* serialization to copy, but **zero HDF5 plumbing**.

**Furthest-developed path = values + JSON.** `writecoeff3.cc` round-trips
`FunctionIOData` (FunctionIO2.h) → JSON → file → `create_function`, validated by
`(f-f2).norm2()`. `FunctionIO2.h` also emits per-leaf **physical grid coords**
(write-only — reconstruction needs only `k`+`nl`+`cell`). ⚠️ **Bug:**
`FunctionIO2.h:242-249` (NDIM==3 coords) drops the outer x-loop and never sets
`c[0]` — emits k² pts, x missing. Fix before porting the interchange writer.

**Conclusion:** function-tree → HDF5 is **not wired**. The examples prove (a) a
values-based per-leaf serialization round-trips, and (b) someone intended the
`h5cpp` wrapper — which **isn't installed**. The HDF5 access layer is a fresh
decision.

Two representations are in play and they differ:
- **Native archive** stores the raw coeff tensors **+ `tree_state`** (works for
  compressed *and* reconstructed trees, GenTensor, any `T`). General.
- **writecoeffs values** stores reconstructed **leaf grid-point values** only
  (`reconstruct()` first; `MADNESS_EXCEPTION` under GenTensor). Narrower, but
  human-meaningful and the basis the viz/cube tools already use.

---

## 3. Build reality (probed 2026-06-19)

- **No `find_package(HDF5)`** anywhere in MADNESS CMake (`cmake/`, top
  `CMakeLists.txt`, `external/`). HDF5 is brand-new to the build.
- **System HDF5 that matches the v3 toolchain (GCC 13.2.0 / Milan):**
  `/gpfs/software/hdf5/gcc13/milan/1.14.3` — full `bin/include/lib`, static
  `libhdf5.a` + shared, `h5dump`/`h5ls`/`h5cc` present. **Serial**
  (`H5_HAVE_PARALLEL` undefined in `H5pubconf.h`).
- Other builds exist (`hdf5-1.{10.5,12.1}-parallel`, …) but they're **Intel**-
  compiler builds — wrong ABI for the GCC13 static MADNESS build.
- **`h5cpp` / `HighFive` are NOT installed** anywhere on the cluster.

**Implication:** a first prototype using the rank-0-gather model (mirroring
`nio=1`) needs only **serial HDF5** — the gcc13/1.14.3 build is a direct fit and
static-linkable into the `BUILD_SHARED_LIBS=OFF` MADNESS build. Parallel HDF5 /
MPI-IO is a later optimization that would require a **GCC13 parallel HDF5 build**
(not currently present).

---

## 4. Plan (phased; each phase gated on the prior)

**P0 — audit + plan. (this doc) ✅**

**P1 — build wiring (opt-in, zero-effect when off).**
- Add `MADNESS_ENABLE_HDF5` CMake option (default OFF). When ON,
  `find_package(HDF5)` (point `HDF5_ROOT` at gcc13/milan/1.14.3) and define
  `MADNESS_HAS_HDF5`. Guard **all** new code behind it so a no-HDF5 build is
  byte-for-byte unchanged (the legacy `.vtk`/archive paths are the model).
- **Validate:** configure with the flag, compile a trivial `H5Fcreate` →
  dataset write/read → `H5Fclose` unit, run it. No MADNESS function code yet.

**P2 — `Function` HDF5 round-trip prototype (the thread's Test recipe).**
Two layers (per the representation decision, §6.2) — start with **Layer A**:

- **Layer A — restart/checkpoint: raw coeff tensors (the archive alternative).**
  New isolated header (e.g. `solvers/function_hdf5_io.hpp`), **not** touching the
  archive core. Rank-0 gather (reuse the `worlddc.h` gather) → one `.h5`:
  - group `/meta`: attrs `k, thresh, ndim, cell, tree_state, T-id` (mirror
    `FunctionImpl::store`);
  - group `/tree`: a `keys` dataset (level + `NDIM` translations) + a **chunked**
    `coeffs` dataset (one leaf's `k^NDIM` raw coeff tensor per chunk).
  - Assert `(f - f2).norm2() < 10·thresh` **and** `norm2(f)==norm2(f2)`, vs the
    legacy `madness::save/load` path. **Time both.** This is the Test recipe.

- **Layer B — interchange/plotting: reconstructed values + coords (separate
  writer/version).** Port the validated `FunctionIOData` (FunctionIO2.h) JSON
  path to HDF5 datasets: `values` (chunked, `k^NDIM` per leaf) + write-only
  `coords` for plotting / MRChem hand-off. Reconstructed-only; **fix the
  NDIM==3 coords bug first**. Independent of Layer A — different payload, same
  libhdf5 + cell/chunk conventions.

**P3 — quantify.**
- File size + wall-time, HDF5 vs legacy, on: h2o ground MOs (5 orbs), one FD
  response state, at k=6 and k=8/k=10. Report a table back to
  `RELEASE_STATUS.md`. This is the go/no-go for P4.

**P4 — wire as an opt-in alternative at ONE call site (gated on P3).**
- Likely `response_state.hpp::save/load` behind a runtime/CMake gate, legacy
  default. If rank-0 gather is the bottleneck, evaluate Parallel HDF5 (needs the
  GCC13 parallel build first). No metadata-layer changes.

---

## 5. Coordination with `feat/amr-export` (the shared-stack contract)

Viz writes **VTKHDF** (`.vtkhdf`/HTG) via VTK's own writer; io-hdf5 writes
**function/checkpoint HDF5** via libhdf5 directly. To honor "one shared HDF5
stack, not two":
1. **Same libhdf5:** both link the system `gcc13/milan/1.14.3` libhdf5 (VTK's
   VTKHDF uses libhdf5 underneath — pin VTK's `HDF5_DIR` to the same build).
2. **Same conventions:** coordinates in **bohr** / `FunctionDefaults::get_cell()`
   (already pinned in `madness_studies/scripts/HTG_VTKHDF_EXPORT_PLAN.md` §3);
   chunk = one leaf node's tensor; `owner` from the live container
   (`coeffs.owner(key)`), not recon-JSON.
3. **Pin both in `docs/operator_contracts.md`** once P1/P2 land.

(Note: VTKHDF's HTG container is *grid/topology* export for ParaView; io-hdf5's
is *coefficient/checkpoint* I/O. Different payloads, same libhdf5 + same
cell/chunk conventions — that's the shared stack.)

---

## 6. Decisions
1. **Access layer** — *recommendation, pending confirm:* **native HDF5 C API**
   for the prototype. Zero vendoring, statically links the gcc13/1.14.3 libhdf5,
   matches MADNESS's existing low-level archive style, lowest build risk, cleanest
   upstream-merge story; payloads are flat double arrays + a few scalar attrs, so
   the C API isn't painful. **HighFive** (header-only, vendor into `src/external/`)
   is the fallback if we want nicer C++ once the prototype proves out — a P4
   reconsideration, not now. **h5cpp** (what the examples reach for) is rejected:
   not installed, needs building.
2. **Representation — DECIDED (user 2026-06-19): both, as two layers.**
   Layer A = **raw coeff tensors + `tree_state`** for *restart* (like the native
   archive). Layer B = **reconstructed values + coords** for *plotting / MRChem
   interchange*, as a **separate version** of the writer. Build A first.
3. **Concurrency — DECIDED:** **rank-0 gather, single file first** (serial HDF5
   is all the GCC13 build has; mirrors `nio=1`). Parallel HDF5/MPI-IO deferred to
   P4 (needs a GCC13 parallel build).

---

## 7. Don't-touch (this thread)
- The archive/serialization format stays the **default**; HDF5 is opt-in.
- `response_metadata.json` — go through `response_metadata.hpp`, never direct.
- The writecoeffs examples are reference; don't make MADNESS depend on `h5cpp`.
- Zero-effect-when-off: a no-HDF5 build must be unchanged.
