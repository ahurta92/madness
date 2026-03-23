# CLAUDE.md — AI Assistant Guide for MADNESS

MADNESS (Multiresolution Adaptive Numerical Environment for Scientific Simulation) is a high-performance C++ framework for solving integral and differential equations in many dimensions using adaptive multiresolution analysis. It includes a full quantum chemistry stack (DFT, HF, post-HF, response properties) and a distributed task-based runtime (MADWorld).

For excited-state reintegration contracts, hot files, and the focused validation matrix, use `$madness-excited-state-dev`. For environment/setup and cross-project run workflow, use `$gecko-madness-workspace`.

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
│   │   │   ├── MolresponseLib.hpp      # Response orchestration (Stage 1-3)
│   │   │   ├── WorkflowBuilders.hpp    # Workflow dispatch (single source of truth)
│   │   │   └── ResponseParameters.hpp # All response input knobs
│   │   ├── misc/              # Shared utilities
│   │   └── external/          # Bundled: gtest, nlohmann_json, muParser
│   ├── apps/                  # Executables
│   │   ├── moldft/            # DFT/HF solver
│   │   ├── molresponse/       # Frequency-dependent response (legacy/stable)
│   │   ├── molresponse_v2/    # Next-gen response module (active development)
│   │   ├── madqc_v2/          # Workflow orchestration entry point + Python tests
│   │   ├── cc2/, mp2/, nemo/  # Post-HF methods
│   │   ├── cis/, zcis/        # Configuration interaction
│   │   └── oep/, pno/         # OEP and PNO methods
│   └── examples/              # 60+ tutorial programs
├── doc/                       # Doxygen + Sphinx documentation
├── docs/                      # Extended molresponse design docs and gap analyses
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

Python-based integration tests for molresponse (in `src/apps/madqc_v2/`):
```bash
python3 test_molresponse_excited_metadata_smoke.py
python3 test_molresponse_excited_restart_reuse.py
python3 test_molresponse_excited_protocol_projection.py
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
- **MacroTaskQ** provides subworld-based worker pools for state-parallel execution.

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

The `molresponse_v2` module (`src/apps/molresponse_v2/`) is the **primary focus of active development**, with major additions landing via `molresponse-feature-next`. This branch adds excited-state support, enhanced subgroup scheduling, and a new persistence/metadata layer.

### Entry points
- **Preferred:** `madqc --wf=response [options] [input_file]`
- **Compatibility wrapper:** `molresponse2 [options] [input_file]` (delegates to same workflow)

### Workflow architecture — single sources of truth
- `src/madness/chem/WorkflowBuilders.hpp` — full dispatch via `add_workflow_drivers()` / `add_response_workflow_drivers()`
- `src/madness/chem/MolresponseLib.hpp` — `response_calculation()` top-level orchestration

**Do not add standalone solver logic back into `molresponse2.cpp`.**

---

### Four-stage pipeline

#### Stage 1 — Planning
Files: `StateGenerator.hpp`, `StateParallelPlanner.hpp`, `DerivedStatePlanner.hpp`

- Converts response input knobs → `PlannedStates` struct
- Canonicalizes frequencies to avoid floating-point near-duplicates
- Builds `StateParallelPlan` with channel ownership assignments, effective mode (`serial`/`parallel`), and `point_parallel_start_protocol_index`
- Builds `DerivedStatePlan` of VBC-driven quadratic requests with dependency tracking
- Builds `ExcitedStateBundlePlan` if `response.excited.enable = true`

Key types:
- `PlannedStates` — aggregates linear, derived, excited plans + parallel scheduling
- `PerturbationChannelAssignment` — ownership map per channel
- `DerivedStateRequest` / `DerivedStateDependency` — quadratic task with linear prereqs
- `DerivedStateGateEntry` / `DerivedStateGateReport` — dependency gate output

#### Stage 2a/2b — Linear State Solves
Files: `FrequencyLoop.cpp/hpp`, `ResponseSolver.cpp/hpp`, `ResponseSolverUtils.hpp`

Two execution paths (selected by `StateParallelPlan::effective_mode`):
- `execute_serial_state_solve(...)` — single communicator
- `execute_subgroup_state_solve(...)` — MacroTaskQ-based worker pools with sharded metadata

Ownership model per protocol:
- **Protocol 0 (channel-series mode):** all groups solve their assigned channels
- **Protocol 1+ (channel-point mode):** fine-grained frequency-block distribution per channel

Inside `computeFrequencyLoop()`:
- Restart precedence: exact checkpoint → coarser-protocol snapshot → frequency continuation → fresh init
- Convergence diagnostics via `ResponseSolveDiagnostics` (residual norms, alpha step, stall detection, explosive-growth detection)
- Wall-time stall detection via `MADQC_POINT_STALL_TIMEOUT_S` env var
- Per-point timing recorded via `StateSolvePersistence::record_timing()`

#### Stage 2c — Excited-State Bundle (NEW)
Files: `ExcitedStateBundleSolver.cpp/hpp`

- Runs between linear state solves (2b) and derived states (2d)
- Executes per-protocol bundle solve for excitation energies / trial states
- Restart precedence: `current protocol snapshot → nearest lower protocol snapshot → generic guess archive → carryover → fresh guess` (`LocalizedGaussianGuess` functor for fresh)
- Records per-protocol results (convergence, per-root energies, slot permutations) through `ResponseRecord2`
- Abstract interface: `ExcitedStateBundleSolver`; placeholder: `ExcitedStateBundleNoopSolver`
- Factory: `make_excited_state_bundle_solver_adapter()`

Key types:
- `ExcitedRootDescriptor` — stable root identity with JSON serialization
- `ExcitedBundleProtocolInput` / `ExcitedBundleProtocolResult` — per-protocol I/O
- `RestartSnapshot` — typed protocol checkpoint (restricted static/dynamic; unrestricted: guess-only until phase 6)

Phase completion status:
- **Phases 1–3 complete:** root descriptors, archive/restart system, solver initialization
- **Phase 4 in progress:** iteration/convergence loop, per-root accelerators, step-restriction policies

#### Stage 2d — Derived-State Execution (ENHANCED)
- Dependency gate filters `DerivedStateRequest` list to only ready requests (all linear dependencies converged)
- Ready requests dispatched via round-robin owner assignment
- Results assembled by `PropertyManager`

#### Stage 3 — Property Assembly
Files: `PropertyManager.hpp`, `MolecularProperty.hpp`

- Assembles polarizability, hyperpolarizability (SHG, OR, all triplets), Raman from solved states
- `try_claim_property_component_task()` — file-based distributed task claiming
- Partial-result support for incomplete state sets
- `PropRow` / `PropKey` — serializable result rows

---

### Persistence and metadata

All persistent metadata flows through `ResponseRecord2` → `response_metadata.json`.

| Scope | File(s) |
|-------|---------|
| Linear state status/timing | `response_metadata.json` keyed by `freq/protocol` |
| Excited state results | `excited_states/protocols/<pkey>/` subtree |
| Debug iteration logs | `response_log.json` |
| Subgroup shards | `response_metadata.group<gid>.json`, merged by rank 0 |

`JsonStateSolvePersistence` (defined in `MolresponseLib.hpp`) implements the `StateSolvePersistence` interface and wraps `ResponseRecord2` + `ResponseDebugLogger`.

Validate correctness by comparing `*.calc_info.json` from `madqc --wf=response` vs `molresponse2`.

---

### Response parameters (new knobs)

| Parameter | Purpose |
|-----------|---------|
| `response.state_parallel` | `off`/`auto`/`on` — subgroup scheduling mode |
| `response.state_parallel_groups` | Number of processor subgroups |
| `response.state_parallel_min_states` | Auto-activation threshold |
| `response.state_parallel_point_start_protocol` | Protocol index for point-mode fanout |
| `response.excited.enable` | Enable excited-state bundle |
| `response.excited.num_states` | Number of excitation roots |
| `response.excited.tda` | Tamm-Dancoff approximation |
| `response.excited.owner_group` | Dedicated subgroup for excited bundle |
| `response.beta.shg`, `.or`, `.all_triplets` | Beta frequency triplet selection |
| `response.force_retry_removed_frequencies` | Retry policy for failed points |

---

### Key files in molresponse_v2

| File | Role |
|------|------|
| `FrequencyLoop.cpp/hpp` | Per-state/per-point iterative solve loop |
| `ResponseManager.cpp/hpp` | Orchestration and state dispatch |
| `ResponseSolver.cpp/hpp` | Core iterative solver |
| `ResponseSolverUtils.hpp` | Step restriction, iteration diagnostics |
| `StateParallelPlanner.hpp` | Subgroup scheduling policy |
| `DerivedStatePlanner.hpp` | VBC-driven quadratic state planning + dependency gates |
| `ExcitedStateBundleSolver.cpp/hpp` | Excited-state bundle solver |
| `ResponseRecord.hpp` | `ResponseRecord2` status/timing/excited metadata |
| `ResponseDebugLogger.hpp` | Structured iteration-level debug logging |
| `ResponseVector.hpp` | Typed variant: static/dynamic × restricted/unrestricted |
| `PropertyManager.hpp` | Property storage, distributed claim-lock, assembly |
| `StateGenerator.hpp` | Linear state list from input knobs |
| `StateParallelPlanner.hpp` | Scheduling policy and ownership plan |
| `MOLRESPONSE_TUTORIAL.md` | Developer pipeline walkthrough |
| `STATE_PARALLEL_DESIGN.md` | Terminology, execution model, data flow, failure handling |
| `EXCITED_STATE_REINTEGRATION_PLAN.md` | Roadmap for excited-state phases |

Design and reference documentation under `src/apps/molresponse_v2/docs/reintegration/`:
- `context.md` — architectural guardrails and standing contracts
- `status.md` — phase dashboard
- `active_contracts.md` — identity contracts, restart support matrix
- Phase contracts/notes/prompts/validation for phases 3 and 4

---

## CI/CD

Workflows in `.github/workflows/`:

- **`cmake.yml`**: Matrix build on Ubuntu 24.04 and macOS-latest
  - Compilers: GCC 13, Clang
  - Task backends: Threads, OneTBB, PaRSEC
  - Both Debug and Release; includes `GENTENSOR=1` compilation check
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
- `docs/` — Extended molresponse design docs
  - `excited_state_reintegration_gap_analysis.md` — feature roadmap and gap analysis
  - `current_excited_state_implementation_report.md` — phase completion status
  - `legacy_excited_state_report.md` — legacy code reference
  - `legacy_molresponse_revival_analysis.md` — historical context

---

## Tips for AI Assistants

1. **Read before modifying.** The codebase is large (~317K lines). Always read relevant files before suggesting changes.
2. **Parallelism is pervasive.** Changes to task submission, futures, or global ops have distributed consequences. MacroTaskQ subworlds are used for state-parallel scheduling — be aware of universe vs subgroup rank.
3. **molresponse_v2 is the active front.** New response features go in `src/apps/molresponse_v2/` and `src/madness/chem/`, not in the legacy `src/apps/molresponse/`.
4. **Workflow wiring is centralized.** `WorkflowBuilders.hpp` and `MolresponseLib.hpp` are the single sources of truth for workflow dispatch; do not duplicate logic elsewhere.
5. **Excited-state phases 1–3 are complete; phase 4 is in progress.** Do not regress the restart/archive contracts defined in `docs/reintegration/active_contracts.md`.
6. **All metadata flows through `ResponseRecord2`.** Never write `response_metadata.json` directly; use the `StateSolvePersistence` / `JsonStateSolvePersistence` interface.
7. **ResponseVector is typed.** Use `StaticRestrictedResponse`, `DynamicRestrictedResponse`, etc. — not raw vector containers. Unrestricted restart is guess-only until phase 6.
8. **Build incrementally.** After edits, run the specific target (`ninja moldft`, `ninja molresponse2`) before the full suite.
9. **Test in parallel.** Many bugs only appear under MPI. Use `mpiexec -np 2` for response tests. Use the Python smoke tests in `src/apps/madqc_v2/` to validate excited-state metadata and restart behaviour.
10. **Serialization requirements.** Any type stored in `WorldContainer` or sent via MPI tasks must implement MADNESS archive serialization.
11. **Never spin.** Use `ENABLE_NEVER_SPIN=ON` in debug builds to avoid spin-wait hangs during debugging.

---

## Scaling Goal and Known Memory Bottleneck

**Project goal:** extend the response solver to larger molecular systems (>20 occupied
orbitals). The three primary constraints are memory per MPI task, wall time, and
strong-scaling efficiency.

### Empirical data (mul_sparse benchmark, March 2026)

In `state_parallel "on"` (subworld) mode every MPI task holds a **complete replica of
all occupied ground-state orbitals** at the current protocol accuracy. Per-task memory:

```
memory ≈ n_occupied × n_leaves(k) × k³ × 8 bytes  +  response_overhead
```

At the final protocol (`k=10`, `thresh=1e-7`), `n_leaves` is ~4.6× the k=6 count.
Measured MaxRSS per MPI task (6 nodes × 8 tasks/node, Xeon Max ≈576 GB/node):

| Molecule    | Occupied orbs | MaxRSS/task | Per-node (×8) | Result  |
|-------------|--------------|-------------|----------------|---------|
| H2O         | 5            | ~14 GB      | ~112 GB        | ✓       |
| C2H4        | 8            | ~32 GB      | ~256 GB        | ✓       |
| CH3OH       | 9            | ~43 GB      | ~344 GB        | ✓       |
| C6H6        | 21           | ~67 GB      | ~536 GB        | ✗ OOM   |
| naphthalene | 34           | >67 GB      | >576 GB        | ✗ OOM   |

C6H6 and naphthalene fail with exit 137 (SIGKILL/OOM) immediately after `PROTOCOL_POLICY`
lines, before the first response iteration. The current workaround is to halve
`--ntasks-per-node` (and double `--cpus-per-task`) to reduce the per-node multiplier.

The ground-state orbital replication across subworlds is the primary blocker.

### Priority improvements

**Memory (high priority):**
- **Pre-flight memory estimate** — before Stage 2 subworld allocation in `MolresponseLib.hpp`,
  compute and print `n_occupied × n_leaves × k³ × 8 bytes` per task and compare against
  available memory. `worldmem.h` already has RSS query infrastructure; `/proc/meminfo`
  gives node-level available memory. A clean abort with a useful message is far better
  than a silent SIGKILL.
- **Per-task RSS at protocol boundaries** — use `worldmem.h` to emit a
  `MEMORY_HWM  rank=N  protocol=N  rss_GB=X` line inside `FrequencyLoop` at each
  protocol transition. This replaces the current post-mortem sacct workflow.
- **Lazy/on-demand orbital loading** — investigate whether `GroundStateData` can load
  orbitals on demand and evict when not needed, rather than keeping a full replica
  resident throughout the run. `ResponseIO.hpp` archive infrastructure could be reused.
- **Shared node-local orbital replica** — rather than one full copy per MPI task,
  a single copy per node (via MPI shared-memory windows) would reduce the memory
  multiplier from `ntasks_per_node` to 1 for the ground-state contribution.

**Logging (medium priority):**
- **Machine-readable protocol lines** in `FrequencyLoop`:
  `PROTOCOL_START index=N thresh=X k=N` and `PROTOCOL_DONE index=N iters=N converged=true`
- **Input path echo** in `madqc.cpp`: `INPUT_FILE  given=X  resolved=/absolute/path`
  (madqc currently lowercases its input argument; callers must pass basename not
  absolute path as a workaround — echoing the resolved path surfaces this immediately)
- **SIGTERM handler** — a signal handler that flushes a
  `MADQC_EXIT reason=signal last_action=X protocol=N` line before the process dies
  is far more useful than a silent kill entry in the SLURM log

**Scaling (longer term):**
- **Delay full subworld expansion** — `protocol0_owner_groups` already reduces active
  groups at protocol 0; investigate extending this so the expensive k=10 protocol
  allocates subworlds only as needed rather than all upfront
- **Per-task timing breakdowns** in `ResponseRecord` (orbital projection, BSH apply,
  orthogonalization) to enable molecule-by-molecule bottleneck identification

### Reference data and operator guide

```
Benchmark study:  /gpfs/scratch/ahurtado/mul_sparse_study/
Run log:          RUNLOG.md               (empirical memory table, job history)
Diagnosis guide:  agents/madqc_error_handling.md  (OOM + error triage steps for agents)
```

When suggesting memory or scaling changes, verify they do not regress H2O/CH3OH/C2H4
and aim to bring C6H6 and naphthalene within budget at 8 tasks/node
(target: <72 GB/task so that 8 × 72 = 576 GB fits within the Xeon Max node capacity).
