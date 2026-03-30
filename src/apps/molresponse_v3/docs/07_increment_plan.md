# molresponse_v3 — Incremental Build Plan

## Strategy

Create a new app `src/apps/molresponse_v3/` that builds the refactored
response pipeline from scratch, one piece at a time. It coexists with
`molresponse_v2` (current production) and `molresponse_legacy` (frozen
reference). Nothing in v2 or legacy changes.

Each increment is a working, testable program. We never have a phase
where "it'll work once all the pieces are in place." Every step compiles,
runs, and can be validated against v2 or legacy output.

## Directory Structure

```
src/apps/molresponse_v3/
├── CMakeLists.txt
├── README.md
├── main.cpp                    # Entry point — grows with each increment
│
├── types/                      # Response type definitions
│   ├── ResponseTypes.hpp       # Static, Full, TDA type traits
│   └── ResponseVector.hpp      # Storage types (x-only, x+y)
│
├── ground/                     # Ground-state context (thin wrapper)
│   └── GroundState.hpp         # Load from archive, provide orbitals/energies
│
├── operators/                  # Building blocks (single-vector operations)
│   ├── Density.hpp             # compute_density_i per type
│   ├── Coulomb.hpp             # apply_coulomb (type-independent)
│   ├── Exchange.hpp            # apply_exchange_i per type
│   ├── XCKernel.hpp            # apply_xc_i per type (no-op for HF)
│   ├── Potential.hpp           # assemble_potential_i per type
│   ├── BSH.hpp                 # apply_bsh_i per type
│   └── Projection.hpp          # project_occupied_i
│
├── solvers/                    # Solver skeletons
│   ├── FDSolver.hpp            # Frequency-dependent (one vector)
│   └── ESSolver.hpp            # Excited-state (coupled bundle)
│
├── naming/                     # State identity and persistence
│   ├── StateIdentity.hpp       # FD keys, ES root_id
│   └── Persistence.hpp         # save/load/query interface
│
├── properties/                 # Property computation
│   ├── ComponentStore.hpp      # Tensor component I/O + computation
│   ├── PropertyEvaluator.hpp   # Physics formulas (pure arithmetic)
│   └── PropertyFormatter.hpp   # Human + JSON output
│
├── orchestrator/               # Pipeline coordination
│   ├── SerialOrchestrator.hpp  # Stage A: simple serial pipeline
│   └── ParallelOrchestrator.hpp # Stage B: added later
│
└── tests/                      # Per-piece tests
    ├── test_types.cpp
    ├── test_operators.cpp
    ├── test_fd_solver.cpp
    ├── test_es_solver.cpp
    ├── test_naming.cpp
    ├── test_properties.cpp
    └── test_orchestrator.cpp
```

## Incremental Build Sequence

### Increment 0: Skeleton App

**What:** Empty app that compiles and links against MADchem. Loads a
ground-state checkpoint from a moldft archive and prints basic info
(molecule, orbital count, energies).

**Why:** Proves the build infrastructure works. Establishes the pattern
for how v3 consumes ground-state data.

**Files:**
- `CMakeLists.txt` — build target linking MADchem
- `main.cpp` — load archive, print info
- `ground/GroundState.hpp` — thin wrapper around archive loading

**Test:** Run on a pre-computed H2 or H2O checkpoint. Verify it prints
correct orbital count and energies.

**Reuses from v2:** `GroundStateData` loading logic (can adapt or
call directly for now — replace later with v3's own implementation
if needed).

---

### Increment 1: Type System + Single-Vector Building Blocks

**What:** Implement the three response types (Static, Full, TDA) and
the single-vector building blocks: density, Coulomb, exchange, BSH,
projection. No solver yet — just the building blocks, tested in isolation.

**Why:** This is the foundation from the type system design doc. Every
building block gets its own test before any solver uses it.

**Files:**
- `types/ResponseTypes.hpp` — type traits / tag dispatch for Static/Full/TDA
- `types/ResponseVector.hpp` — x-only and x+y storage
- `operators/Density.hpp` through `operators/Projection.hpp`
- `tests/test_types.cpp` — storage layout tests
- `tests/test_operators.cpp` — per-operator tests

**Tests (Level 1 from refactor plan):**
- Density: known orbitals + known response → verify density matches
  analytic expression for each type
- Coulomb: apply to known density → verify potential
- Exchange: apply to known input → compare to legacy output
- BSH: apply to known theta → compare to legacy output
- Projection: verify projected result is orthogonal to occupied space

**Validation:** For each operator, run with same inputs as legacy code
and compare outputs. This is where you catch any misunderstanding of
the type definitions.

---

### Increment 2: FD Solver (Static)

**What:** Implement the FD solver skeleton for static response (the
simplest case: x-only, y=x, omega=0). Given a perturbation and ground
state, iterate to convergence and produce a response vector.

**Why:** The static case has the fewest moving parts. If this works,
the solver skeleton is correct and the static building blocks are correct.

**Files:**
- `solvers/FDSolver.hpp` — iteration loop calling building blocks
- `tests/test_fd_solver.cpp`

**Tests (Level 2 from refactor plan):**
- Single-iteration test: verify one step produces expected output
- Convergence test: solve H2 static dipole, verify convergence
- Compare converged result to legacy and v2 output

**Property check:** Compute static polarizability from the converged
state (simple inner product). Compare to legacy / v2 / known value.

---

### Increment 3: FD Solver (Full Dynamic)

**What:** Extend FD solver to handle Full type (x+y) at nonzero
frequency. This adds the coupled BSH update with mu+ and mu-.

**Files:** Update `solvers/FDSolver.hpp` to dispatch by type.

**Tests:**
- Solve H2O dynamic dipole at omega=0.02
- Compare alpha(omega) to v2 output

---

### Increment 4: Naming + Persistence

**What:** Implement the naming convention for FD states and the
persistence interface (save/load/query). The solver now saves its
results and can restart.

**Files:**
- `naming/StateIdentity.hpp` — FD key construction
- `naming/Persistence.hpp` — save/load interface
- `tests/test_naming.cpp`

**Tests:**
- Save a converged state, reload it, verify identity
- Restart a solve partway, verify it produces same result as
  straight-through
- Protocol ramp: solve at low protocol, reload at higher protocol

---

### Increment 5: Serial Orchestrator (FD only)

**What:** The simple three-stage pipeline for FD response: plan states,
solve over protocols, compute alpha. This is Stage A from the refactor
plan.

**Files:**
- `orchestrator/SerialOrchestrator.hpp`
- `main.cpp` updated to run the orchestrator
- `tests/test_orchestrator.cpp`

**Tests (Level 3 from refactor plan):**
- End-to-end: run H2O polarizability, compare to v2
- Restart: run partway, kill, restart, verify same result
- Protocol ramp: verify low→high protocol produces expected accuracy

**This is the first milestone where v3 can replace v2 for simple
polarizability calculations.**

---

### Increment 6: Property Layer

**What:** Implement ComponentStore, PropertyEvaluator, PropertyFormatter.
Add beta/raman component computation and property assembly.

**Files:** `properties/` directory

**Tests:**
- Component computation against known values
- Property evaluation against known tensors
- JSON and text output verification

---

### Increment 7: ES Solver (TDA)

**What:** Implement the ES solver for TDA (simplest ES case: x-only,
single BSH, but coupled bundle with rotation).

**Files:**
- `solvers/ESSolver.hpp`
- `tests/test_es_solver.cpp`

**Tests:**
- Solve H2 TDA excited states
- Compare excitation energies to legacy code
- Verify root naming and identity tracking

---

### Increment 8: ES Solver (Full)

**What:** Extend ES solver to full response (x+y, coupled BSH).

**Tests:**
- Solve H2 full TDDFT excited states
- Compare to legacy

---

### Increment 9: ES Naming + Restart

**What:** Add ES-specific naming (root_id, bundle archives,
slot_permutation) and restart.

**Tests:**
- Protocol ramp with excited states
- Verify root identity stability across protocols

---

### Increment 10: madqc Integration

**What:** Wire v3 into the madqc workflow system so `madqc --wf=response`
can optionally use the v3 pipeline. Add structured timing output.

**Tests:**
- Run same calculation through v2 and v3, compare calc_info.json
- Verify timing report structure

---

### Increment 11: State-Parallel (Stage B)

**What:** Add the parallel orchestrator on top of the serial one.

**Tests:**
- Serial vs parallel on same calculation, verify identical results
- Scaling test: measure speedup

---

## Build Integration

### CMakeLists.txt (initial)

```cmake
cmake_minimum_required(VERSION 3.12)

if(NOT TARGET MADchem)
  message(STATUS "molresponse_v3: MADchem not available, skipping")
  return()
endif()

# Start minimal — just the skeleton
set(MOLRESPONSE_V3_SOURCES
    main.cpp
)

add_executable(molresponse_v3 ${MOLRESPONSE_V3_SOURCES})
target_link_libraries(molresponse_v3 MADchem)

target_include_directories(molresponse_v3 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/madness
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_BINARY_DIR}/src/madness
)

# Can also link MADresponse2 to reuse pieces during transition
# target_link_libraries(molresponse_v3 MADresponse2)

add_dependencies(applications-madness molresponse_v3)

install(TARGETS molresponse_v3
        DESTINATION "${MADNESS_INSTALL_BINDIR}")
```

### Relationship to v2

During development, v3 can optionally link against MADresponse2 to
reuse pieces (GroundStateData, response vector I/O, etc.) while
building its own implementations. As v3 matures, those dependencies
are replaced one by one until v3 is self-contained.

v2 continues to be the production code. v3 is validated against v2
at every increment. Once v3 passes all tests and produces identical
results, the switch happens.

---

## What to Build First (Today)

Increment 0: the skeleton app.

1. Create `src/apps/molresponse_v3/`
2. Write CMakeLists.txt
3. Write main.cpp that loads a ground-state archive and prints info
4. Build and run on your H2 checkpoint
5. Verify it works

That's a 30-minute task. Once it compiles and runs, you have the
workspace and you know the build infrastructure is right. Then
Increment 1 starts the real work.
