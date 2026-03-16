# CLAUDE.md — AI Assistant Guide for MADNESS

MADNESS (Multiresolution Adaptive Numerical Environment for Scientific Simulation) is a high-performance C++ framework for solving integral and differential equations in many dimensions using adaptive multiresolution analysis. It includes a full quantum chemistry stack (DFT, HF, post-HF, response properties) and a distributed task-based runtime (MADWorld).

---

## Repository Layout

```
madness/
├── src/
│   ├── madness/               # Core libraries
│   │   ├── world/             # MADWorld distributed runtime
│   │   ├── tensor/            # N-dimensional tensor/BLAS/LAPACK wrappers
│   │   ├── mra/               # Multiresolution analysis (wavelets, operators)
│   │   ├── chem/              # Quantum chemistry (SCF, DFT, CC, response)
│   │   ├── misc/              # Shared utilities
│   │   └── external/          # Bundled: gtest, nlohmann_json, muParser
│   ├── apps/                  # Executables
│   │   ├── moldft/            # DFT/HF solver
│   │   ├── molresponse/       # Frequency-dependent response (legacy/stable)
│   │   ├── molresponse_v2/    # Next-gen response module (active development)
│   │   ├── madqc_v2/          # Workflow orchestration entry point
│   │   ├── cc2/, mp2/, nemo/  # Post-HF methods
│   │   ├── cis/, zcis/        # Configuration interaction
│   │   └── oep/, pno/         # OEP and PNO methods
│   └── examples/              # 60+ tutorial programs
├── doc/                       # Doxygen + Sphinx documentation
├── cmake/                     # CMake modules and toolchains
├── external/                  # Dependency version configs
├── admin/                     # Test data and FCI reference data
├── .github/workflows/         # CI/CD (cmake.yml, make_doxygen.yml)
├── AGENTS.md                  # Agent build/test quick reference
├── GEMINI.md                  # Getting started guide (Ubuntu)
└── INSTALL.md                 # Full installation reference
```

---

## Build System

**CMake ≥ 3.12**, Ninja preferred, C++17 minimum (C++20 also supported).

### Fast debug build (standard for development)

```bash
mkdir build && cd build
cmake .. -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS_DEBUG="-O0 -g -Wall" \
  -DENABLE_NEVER_SPIN=ON \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_TESTING=ON \
  -DLAPACK_LIBRARIES="-L/usr/lib/x86_64-linux-gnu/openblas-serial -lopenblas -llapack"
```

### Key CMake options

| Option | Default | Purpose |
|--------|---------|---------|
| `CMAKE_BUILD_TYPE` | Debug/Release | Build optimization level |
| `ENABLE_MPI` | ON | MPI distributed parallelism |
| `BUILD_TESTING` | ON | Build test targets |
| `BUILD_SHARED_LIBS` | OFF | Use static libs (required for ASLR) |
| `MADNESS_ENABLE_CEREAL` | OFF | Enable cereal serialization |
| `ENABLE_NEVER_SPIN` | OFF | Disable spin-wait (useful for debugging) |

### Required Ubuntu packages

```bash
sudo apt-get install -y build-essential gcc gfortran libopenblas-serial-dev \
  cmake git ninja-build openmpi-bin openmpi-common libopenmpi-dev
```

---

## Building Targets

```bash
# All targets
ninja

# Specific targets
ninja testsuite      # Numerical library tests
ninja moldft         # DFT solver
ninja molresponse2   # Response property solver (v2)
ninja check-short-madness  # Quick validation suite
```

---

## Testing

```bash
# Quick validation (build + run)
ninja check-short-madness

# Numerical library test suite (build)
ninja testsuite
# Run single-process
./src/madness/src/testsuite
# Run with MPI + threads
MAD_NUM_THREADS=2 mpiexec -np 2 ./src/madness/src/testsuite
# Success: prints "testsuite passed:  true" with exit code 0

# Smoke test DFT
./src/apps/moldft/moldft --geometry=water --dft="xc=lda"
```

Test infrastructure uses **Google Test** (bundled) with CMake helpers: `AddUnittests`, `AddMPITests`, `AddScriptedTests`.

---

## Code Conventions

### Language & style
- **C++17** minimum; some new code targets C++20.
- All MADNESS library code lives in the `madness` namespace.
- `madness::World` must be initialized before using any MADNESS functionality.
- Header guards use the pattern `MADNESS_<COMPONENT>_<FILE>_H__INCLUDED`.
- Use `#include <madness/component/file.h>` style includes, not relative paths.

### Parallelism model
- **MADWorld** provides the task-based distributed runtime. All parallel tasks go through `World::taskq`.
- `madness::Future<T>` is the primary async primitive.
- `WorldGOP` provides collective operations (broadcast, reduce, barrier).
- MPI process groups are first-class; `World` wraps an MPI communicator.
- Thread count is controlled at runtime via `MAD_NUM_THREADS` environment variable.

### Numerical conventions
- `madness::Function<T,NDIM>` is the core MRA type. `T` is the value type (usually `double`), `NDIM` is dimension.
- Functions are represented adaptively; do not assume a fixed grid.
- Operators (`DerivativeOperator`, `SeparatedConvolution`) act on `Function` objects.
- Use `compress()` / `reconstruct()` explicitly when needed before operations.
- Threshold/precision is set via `FunctionDefaults<NDIM>::set_thresh(eps)`.

### Serialization
- MADNESS archive system (`madness/world/archive.h`) handles serialization for MPI, checkpointing, and I/O.
- Objects passed to tasks or stored in `WorldContainer` must be serializable.

---

## Active Development: molresponse_v2

The `molresponse_v2` module (`src/apps/molresponse_v2/`) is under active development and is the **primary focus of recent commits**. Key points:

### Entry points
- **Preferred:** `madqc --wf=response [options] [input_file]`
- **Compatibility wrapper:** `molresponse2 [options] [input_file]` (delegates to same workflow)

### Workflow architecture
Response workflow wiring is centralized in:
- `src/madness/chem/WorkflowBuilders.hpp` — dispatch via `add_workflow_drivers()` / `add_response_workflow_drivers()`
- `src/madness/chem/MolresponseLib.hpp` — `response_calculation()` orchestration

**Do not add standalone solver logic back into `molresponse2.cpp`.**

### Three-stage pipeline
1. **Planning** — `StateGenerator.hpp`, `StateParallelPlanner.hpp`, `DerivedStatePlanner.hpp`
   - Converts response knobs → `PlannedStates`
   - Handles state ownership, frequency canonicalization, derived-state dependencies
2. **Solving** — `FrequencyLoop.cpp`, `ResponseSolver.cpp`
   - Serial or subgroup (state-parallel) execution paths
   - Restart-aware: loads from exact checkpoint, coarser protocol, frequency continuation, or fresh init
   - Timing recorded via `StateSolvePersistence::record_timing()`
3. **Properties** — `PropertyManager.hpp`, `MolecularProperty.hpp`
   - Assembles polarizability, hyperpolarizability, Raman from solved states

### Persistence and metadata
- `JsonStateSolvePersistence` wraps `ResponseRecord2` (status/timing) + `ResponseDebugLogger` (iteration-level debug)
- Subgroup mode: sharded `response_metadata.group<gid>.json` files, merged by rank 0
- Validate correctness by comparing `*.calc_info.json` from `madqc --wf=response` vs `molresponse2`

### Key files in molresponse_v2
| File | Role |
|------|------|
| `FrequencyLoop.cpp/hpp` | Per-state iterative solve loop |
| `ResponseManager.cpp/hpp` | Orchestration and state dispatch |
| `ResponseSolver.cpp/hpp` | Core iterative solver |
| `StateParallelPlanner.hpp` | Subgroup scheduling policy |
| `ExcitedStateBundleSolver.cpp/hpp` | Excited-state bundle solver (in progress) |
| `ResponseRecord.hpp` | Status/timing metadata |
| `ResponseDebugLogger.hpp` | Iteration-level debug logging |
| `MOLRESPONSE_TUTORIAL.md` | Developer pipeline walkthrough |
| `EXCITED_STATE_REINTEGRATION_PLAN.md` | Roadmap for excited-state features |

---

## CI/CD

Workflows in `.github/workflows/`:

- **`cmake.yml`**: Matrix build on Ubuntu 24.04 and macOS-latest
  - Compilers: GCC 13, Clang
  - Task backends: Threads, OneTBB, PaRSEC
  - Both Debug and Release
  - MPI tests with 3 threads
- **`make_doxygen.yml`**: API documentation generation

Static analysis: `.clang-tidy` with `bugprone-*`, `readability-*`, and other checkers.

---

## Documentation

- `INSTALL.md` — Comprehensive build reference
- `AGENTS.md` — Quick build/test commands for agents
- `GEMINI.md` — Ubuntu quick-start guide
- `doc/` — Sphinx (ReadTheDocs) + Doxygen sources
  - `doc/numerical_library.md`, `doc/quantum.md`, `doc/runtime.md`
  - Tutorials: `doc/tutorial/`, `doc/getting_started/`

---

## Tips for AI Assistants

1. **Read before modifying.** The codebase is large (~317K lines). Always read relevant files before suggesting changes.
2. **Parallelism is pervasive.** Changes to task submission, futures, or global ops have distributed consequences.
3. **molresponse_v2 is the active front.** New response features go in `src/apps/molresponse_v2/` and `src/madness/chem/`, not in the legacy `src/apps/molresponse/`.
4. **Workflow wiring is centralized.** `WorkflowBuilders.hpp` and `MolresponseLib.hpp` are the single sources of truth for workflow dispatch; do not duplicate logic elsewhere.
5. **Build incrementally.** After edits, run the specific target (`ninja moldft`, `ninja molresponse2`) before the full suite.
6. **Test in parallel.** Many bugs only appear under MPI. Use `mpiexec -np 2` for response tests.
7. **Serialization requirements.** Any type stored in `WorldContainer` or sent via MPI tasks must implement MADNESS archive serialization.
8. **Never spin.** Use `ENABLE_NEVER_SPIN=ON` in debug builds to avoid spin-wait hangs during debugging.
