# Agent Instructions for MADNESS (`molresponse-feature-next`)

This file is the repo-local quick reference for agents working in this worktree.
Read [`CLAUDE.md`](./CLAUDE.md) first for architecture, branch status, and
response-specific design constraints.

For cross-project Gecko + MADNESS environment setup, molecule-library workflow,
helper commands from `~/.bashrc`, HBM launch recipes, and shared MCP server
details, use `$gecko-madness-workspace`.

## Read This First

- Active development happens in `src/apps/molresponse_v2/` and `src/madness/chem/`.
  Do not add new solver logic to the legacy `src/apps/molresponse/` tree.
- Workflow dispatch is centralized in `src/madness/chem/WorkflowBuilders.hpp`
  and `src/madness/chem/MolresponseLib.hpp`.
- Many response bugs only reproduce under MPI. Prefer parallel smoke tests over
  single-rank-only validation.
- Never move `madness-worktrees/builds/`. CMake caches in this workspace use
  absolute paths.

## Environment Setup

Intel oneAPI (MKL, TBB) must be loaded before any build or run.

```bash
source setenv.sh
```

`setenv.sh` lives at the repository root and sources `~/load_xeonmax.sh`.
No `cmake` or `ninja` invocation will work reliably without it.

## Build Layout

| Item | Path |
|------|------|
| Worktree source | `/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next` |
| `./build` symlink target | `/gpfs/projects/rjh/adrian/development/madness-worktrees/builds/molresponse-feature-next/96core` |
| Active debug build | `/gpfs/projects/rjh/adrian/development/madness-worktrees/builds/molresponse-feature-next/debug` |
| Available configs | `40core`, `96core`, `debug`, `warn-check` |

Notes:
- `./build` is the safe default for `ninja -C build ...` commands from the repo root.
- The debug build is the active development build used by the shared shell helpers and `madqc` path overrides.

## Building

Preferred application target:

```bash
source setenv.sh
ninja -C build madqc
```

Other useful targets:

```bash
source setenv.sh
ninja -C build molresponse2
ninja -C build moldft
ninja -C build testsuite
ninja -C build check-short-madness
```

To configure a fresh build:

```bash
source setenv.sh
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_NEVER_SPIN=ON \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_TESTING=ON
```

For debug builds, [`CLAUDE.md`](./CLAUDE.md) also recommends:

```bash
-DCMAKE_CXX_FLAGS_DEBUG="-O0 -g -Wall"
```

## Testing

Quick validation:

```bash
source setenv.sh
ninja -C build check-short-madness
```

Numerical library tests:

```bash
source setenv.sh
ninja -C build testsuite
./build/src/madness/src/testsuite
MAD_NUM_THREADS=2 mpiexec -np 2 ./build/src/madness/src/testsuite
```

DFT smoke test:

```bash
./build/src/apps/moldft/moldft --geometry=water --dft="xc=lda"
```

Python response smoke tests live in `src/apps/madqc_v2/`:

```bash
python3 src/apps/madqc_v2/test_molresponse_excited_metadata_smoke.py
python3 src/apps/madqc_v2/test_molresponse_excited_restart_reuse.py
python3 src/apps/madqc_v2/test_molresponse_excited_protocol_projection.py
```

Guidance from [`CLAUDE.md`](./CLAUDE.md):
- Build incrementally after edits. Prefer the narrow target first, then broader validation.
- Test in parallel with `mpiexec -np 2` when touching response scheduling, futures, metadata, or restart code.
- Keep `ENABLE_NEVER_SPIN=ON` in debug builds to avoid spin-wait hangs.

## Response-Code Guardrails

- `molresponse_v2` is the active response frontend.
- Keep workflow wiring in `WorkflowBuilders.hpp` and `MolresponseLib.hpp`.
- Do not reintroduce standalone solver logic into `molresponse2.cpp`.
- Metadata should flow through `ResponseRecord2` and the
  `StateSolvePersistence` / `JsonStateSolvePersistence` interfaces.
  Do not write metadata JSON ad hoc.
- Any type sent through MADNESS tasks or stored in `WorldContainer` must be
  serializable with the MADNESS archive system.

## Dependencies

Ubuntu packages used by the project:

```bash
sudo apt-get update && sudo apt-get install -y \
  build-essential gcc gfortran libopenblas-serial-dev \
  cmake git ninja-build openmpi-bin openmpi-common libopenmpi-dev
```
