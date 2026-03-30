# MADNESS Development — Project Overview

**Timeline:** Next MADNESS release in ~6 months

---

## What is molresponse?

molresponse is a module within MADNESS for computing frequency-dependent
molecular response properties using TDDFT/TDHF. It computes independent
response vectors which are then combined to produce linear and nonlinear
response properties.

The code has four layers:

1. **Numerical solvers** — Iterative solvers for coupled response equations
   (static and dynamic). The same algorithmic steps apply across model types
   (closed-shell/open-shell, static/dynamic), but each variant defines its
   own building blocks within that shared algorithm. There are two solver
   skeletons: FD (frequency-dependent, one vector at a time) and ES
   (excited-state, coupled bundle with rotation/diagonalization). Both
   solvers call the same single-vector building blocks; ES loops them
   over the bundle.

2. **Naming conventions** — Response vectors are named to encode the
   state(s), frequency/frequencies, and protocol. For FD states, identity
   is (perturbation, frequency) — both are inputs and stable. For ES
   states, identity is a permanent root_id (stable_index) because the
   excitation energy is an output that changes during iteration. ES
   bundles are archived as a unit per protocol, with a slot_permutation
   tracking root reordering.

3. **Property computation** — Two layers: a ComponentStore that loads
   converged response vectors and computes tensor components via MRA
   inner products, and a PropertyEvaluator that takes those tensors and
   applies physics formulas to produce physical properties. This
   separation means the expensive MRA work is decoupled from the cheap
   arithmetic.

4. **State-parallel execution** — Distributes independent response state
   computations across parallel subworlds, managing the tradeoff between
   memory and parallelization.

### What is protocol?

Protocol controls the accuracy of a calculation. It is defined by the
truncation threshold and polynomial order used in the MRA representation.
Lower protocol means coarser trees (faster, less accurate). Higher
protocol means finer resolution (slower, more accurate). In practice,
SCF and response equations start at low protocol and ramp up — it is
the convergence metrics (dconv, econv, BSH residual) that measure the
achieved accuracy, but protocol is the knob we turn to get there.

### Response Type System

The fundamental distinction is how many independent function sets the
response state carries and what symmetry relationship holds between them:

- **Static (x-only, y=x):** Response density = 2·Σ[φᵢ·xᵢ]. Single BSH
  equation per orbital, no frequency shift.
- **Full (x+y, independent):** Response density = Σ[φᵢ·(xᵢ+yᵢ)]. Coupled
  BSH equations with frequency-shifted Green's functions.
- **TDA (x-only, y=0):** Response density = Σ[φᵢ·xᵢ]. Single BSH equation
  but different density (and therefore different operators) than static.

The density definition determines every operator (Coulomb, exchange, XC).
Two types with the same storage layout can have completely different
operators because their densities differ.

### Core design principle

The numerical solver and property computation algorithms should each be
written once. Each response type (static/full/TDA) supplies its own
definitions of the foundational building blocks (density, exchange, BSH,
etc.). Each property type supplies its own combination rules. The
algorithm skeleton contains no type-specific or property-specific logic.
Adding a new type or property means supplying new building block
definitions — the algorithm doesn't change.

### Current capabilities

- Dipole frequency response: polarizability and hyperpolarizability
- Linear response to nuclear displacement: Raman properties

---

## Projects

### T0 — molresponse Refactor

The code has grown organically and the four layers are not cleanly
separated. Model-specific logic leaks into what should be generic
algorithms, and the orchestrator (MolresponseLib.hpp) mixes scheduling,
restart, metadata, and subgroup logic in thousands of lines.

**Goal:** Rebuild from the ground up as molresponse_v3. Clean separation
of the four layers with well-defined interfaces. Start with a simple
serial orchestrator, add state-parallel on top.

**Strategy:** Create `src/apps/molresponse_v3/` alongside v2 (production)
and legacy (frozen reference). Build incrementally — each increment is a
working, testable program validated against v2. Use the legacy code as
the structural reference for what the clean algorithm should look like.

**First step:** Increment 0 — skeleton app that loads ground state.
Then Increment 1 — type system and single-vector building blocks.

**Unblocks:** T1, T2, T5

---

### T1 — Excited-State Response

A legacy version of MADNESS had working excited-state code, but it is
incompatible with the current molresponse framework — different naming
conventions, different structure.

**Goal:** Port the legacy excited-state code into the v3 framework so
we can compute two-photon absorption and resonant Raman response.

**Approach:** The ES solver in v3 uses the same single-vector building
blocks as the FD solver, looped over the bundle. The type system
(Full or TDA) supplies operator definitions. ES naming uses root_id
for permanent identity, bundles as archive units.

**Depends on:** T0 (v3 type system and solver skeletons)

---

### T2 — State-Parallel Resource Optimization

The state-parallel layer works but choosing the number of subworlds is
currently manual. The right choice depends on molecule size, property
type, total work volume, and available resources.

**Goal:** An automated, general strategy for choosing the number of
subworlds that balances parallelization against memory constraints.

**Approach:** Profile efficiency and memory usage as a function of system
size, then develop a heuristic or algorithm. In v3, this is Increment 11
— added after the serial orchestrator is proven correct.

**Depends on:** T0 (clean serial orchestrator as reference for validation)

---

### T3 — Protocol-Accuracy Workflow

Different properties may have different relationships between protocol
settings and achieved property accuracy.

**Goal:** A reusable methodology that, given any response property type,
produces its precision-accuracy curve and a recommended protocol ramp-up
strategy. The deliverable is the workflow itself — when a new property is
developed in the future, this workflow is ready to characterize it.

**Currently being developed alongside:** T4 (Raman benchmarking is the
first test case)

---

### T4 — Raman Benchmarking

**Goal:** Compute Raman response properties for simple molecules and
benchmark accuracy versus protocol, comparing MADNESS results to
Gaussian basis set results via Gecko.

**Depends on:** T3 (protocol methodology), T6 (Gecko for comparison)

---

### T5 — TDDFT / XC Kernel Reintegration

**Goal:** Ensure molresponse has correct XC kernels for both linear and
nonlinear response within the TDDFT framework. In the v3 type system,
the XC kernel is a building block with a no-op implementation for HF.

**Depends on:** T0 (clean building block interfaces)
**Enables:** TDDFT variants of all properties (T1, T4)

---

### T6 — Gecko

Gecko is the platform for comparing MADNESS calculations against
Gaussian basis set calculations.

**Goal:** A reliable, general comparison engine that supports validation
for any property type.

**Consumed by:** T3, T4, and all future benchmarking.

---

## How everything connects

```
T0 (refactor/v3) ──────unblocks──▶ T1, T2, T5
T1 (excited states) ──consumes──▶ T3 (protocol workflow), T2 (state-parallel)
T2 (state-parallel) ──enables───▶ T3 (at scale), T4 (throughput)
T3 (protocol workflow)──consumes──▶ T6 (Gecko for comparison)
T4 (Raman bench) ─────consumes──▶ T3 (protocol workflow), T6 (Gecko)
T5 (TDDFT/XC) ────────enables───▶ T1, T4 (TDDFT variants)
T6 (Gecko) ────────────produces──▶ validation for T3, T4, and all future work
```

Key insight: these threads produce reusable infrastructure for each other.
The protocol-accuracy workflow (T3) is developed while doing Raman
benchmarking (T4), but once it exists, it's ready for excited-state
properties (T1). Gecko (T6) validates everything. The state-parallel
optimization (T2) determines whether protocol studies can run at scale.
The v3 refactor (T0) is what makes independent progress on any of this
tractable.

---

## What's next

- [ ] Run Increment 0 agent task — create molresponse_v3 skeleton app
- [ ] Run Increment 1 — type system and single-vector building blocks
- [ ] Run agent inventory tasks on legacy and current code
- [ ] Fill in current status and recent progress for each thread
- [ ] Identify the most immediate blockers and decide what to work on first
