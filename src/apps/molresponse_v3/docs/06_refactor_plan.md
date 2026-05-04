# molresponse Refactor Plan

## Strategy

Use the legacy excited-state code as the structural reference for refactoring
molresponse. The legacy code is clean and correct but lacks naming conventions,
restart, and protocol handling. The current molresponseLib has those features
but is messy. The refactored version should combine the clarity of the legacy
code with the full capabilities of the current code.

The refactor is being built as `molresponse_v3` in `src/apps/molresponse_v3/`.
It coexists with v2 (production) and legacy (frozen reference). Each increment
is a working, testable program validated against v2. See `increment_plan.md`
for the full build sequence (Increments 0–11).

Design documents guiding the refactor:
- `type_system_design.md` — solver skeletons, response types, building blocks
- `naming_conventions.md` — FD and ES state identity and persistence
- `property_design.md` — component store + property evaluator separation
- `madqc_interface_timing.md` — clean interface and structured timing

---

## Phase 1: Freeze and Document Legacy Code

Minor printing fixes only — no new features.

Produce a document that covers:
- Every function: name, purpose, inputs, outputs
- The algorithmic flow from start to finish (the solver skeleton)
- What building blocks the algorithm calls at each step
- What is hardcoded vs what is parameterized
- What is missing (naming, restart, protocol handling, state-parallel)

**Agent task:** Walk through the legacy codebase and produce this inventory
automatically. Output should be a structured list that can be reviewed
and corrected.

---

## Phase 2: Document Current molresponseLib

Same inventory, but with additional questions:
- Every function: name, purpose, inputs, outputs
- The algorithmic flow (may be harder to trace if logic is tangled)
- Where does naming convention logic live?
- Where does restart logic live?
- Where does protocol handling live?
- Where does state-parallel logic live?
- Where are concerns mixed? (e.g., solver function that also handles naming,
  or property computation that contains model-specific logic)

**Agent task:** Walk through molresponseLib and produce the same structured
inventory, flagging every location where more than one concern is present
in the same function or code block.

---

## Phase 3: Map Legacy to Current

Compare the two inventories side by side:
- For each step in the legacy algorithm skeleton, find the corresponding
  code in molresponseLib
- Identify what molresponseLib adds at each step (naming, restart, protocol)
- Identify where those additions are cleanly separated vs entangled with
  the core algorithm
- Produce a mapping: legacy function → current function(s) → what needs
  to be preserved vs restructured

This mapping IS the refactor plan. It tells you:
- What the clean algorithm skeleton looks like (from legacy)
- What building blocks need to be factored out as interfaces (from current)
- Where the current code deviates from the clean structure and why

---

## Phase 4: Incremental Build with Tests

Build the refactored version incrementally, testing at three levels.

### Level 1 — Building Block Tests (unit level)

Test each foundational building block in isolation before composing them
into the full algorithm. These tests verify that individual pieces are
correct independent of the rest.

Solver building blocks:
- **Initial guess construction:** Given ground-state orbitals and a
  perturbation type, does it produce a reasonable starting point?
  Compare to legacy output for the same input.
- **Right-hand side assembly:** Given orbitals and a perturbation, does
  it produce the correct RHS vector? Compare to legacy and/or analytic
  expressions for simple cases.
- **Coupling operator application:** Apply to a known input vector,
  compare to reference.
- **XC kernel application:** Apply to a known input, compare to reference.
  (Critical for T5 — TDDFT reintegration.)
- **BSH operator application:** Apply to known input, compare to reference.
- **Residual computation:** Given current and previous iterates, does
  the residual match expected values?
- **Convergence check:** Given a residual and threshold, does it correctly
  report converged/not converged?
- **Protocol update:** Does it correctly tighten thresholds and update
  polynomial order according to the protocol schedule?
- **State I/O (naming):** Does save/load round-trip correctly? Does the
  name encode the right metadata? Does restart recover the correct state?

Property building blocks:
- **State requirement lookup:** Given a property type, does it correctly
  identify which response states (and at what frequencies) are needed?
- **Combination rule:** Given known response states, does the contraction
  produce the correct property value? Compare to legacy output and/or
  analytic results for simple systems (e.g., hydrogen atom, helium).
- **Symmetry:** Does it correctly identify redundant components and
  skip or deduplicate them?

### Level 2 — Algorithm Flow Tests (integration level)

Test that the solver skeleton correctly orchestrates the building blocks
through a full solve. These don't test whether the building blocks are
right — that's Level 1. These test that the algorithm calls them in the
right order, handles the iteration loop correctly, and respects the
protocol ramp-up.

- **Single-iteration test:** Run one iteration of the solver with known
  inputs. Verify it calls the right building blocks in the right order
  and produces the expected output.
- **Convergence test:** Run the solver on a trivial system (very small
  molecule, low protocol) and verify it converges in the expected number
  of iterations and reaches the expected accuracy.
- **Protocol ramp-up test:** Start at low protocol, verify the solver
  correctly detects when to ramp up, and that higher protocol produces
  tighter convergence.
- **Restart test:** Run partway, save, reload, continue. Verify the
  result matches a straight-through run.
- **Model-agnosticism test:** Run the same algorithm with two different
  model type definitions (e.g., closed-shell static and closed-shell
  dynamic). Verify that the algorithm code path is identical — only the
  building block calls differ.

### Level 3 — End-to-End Property Tests (validation level)

Test that the full pipeline (solver → converged states → property
computation) produces correct physical results. This is where you
compare against the legacy code and against Gecko.

- **Legacy parity:** For every property the legacy code computes correctly,
  the refactored code should produce the same result (within numerical
  tolerance) at the same protocol.
- **Gecko comparison:** For key test molecules, compare MADNESS results
  to Gaussian basis set results at high enough protocol that the
  comparison is meaningful.
- **Known values:** For systems with known analytic or high-accuracy
  reference values (e.g., hydrogen atom polarizability), verify
  agreement.

### Test Molecule Set

Choose a small set of molecules that covers the important variations:
- A very small system for fast iteration and debugging (e.g., H2, He)
- A small polyatomic for realistic but manageable tests (e.g., H2O, HF)
- A medium system to test scaling behavior (e.g., formaldehyde, benzene)

Use the same molecules consistently across all tests so results are
comparable and regressions are obvious.

---

---

## Phase 5: Rebuild the Orchestrator from Scratch

The current MolresponseLib.hpp is ~3000+ lines of scheduling, restart,
metadata persistence, subgroup management, and shard merging. Rather
than refactoring it incrementally, build a new orchestrator in two stages.

### Stage A: Serial Orchestrator (build first)

A clean, minimal orchestrator that runs on one world with no
state-parallel complexity. It should read like the three-stage pipeline:

```
run_response(world, params):
    # Stage 1: Plan
    ground   = load_ground_context(world, params)
    states   = plan_required_states(world, params, ground)

    # Stage 2: Solve
    for protocol in protocol_thresholds:
        set_protocol(world, ground, protocol)
        for state in states.linear:
            solve_fd_response(world, ground, state, protocol)
        if states.excited.enabled:
            solve_es_bundle(world, ground, states.excited, protocol)
        solve_derived_states(world, ground, states.derived, protocol)

    # Stage 3: Properties
    properties = compute_properties(world, ground, states, params)
    return properties
```

Each function called above is a clean, testable unit. The orchestrator
itself is just sequencing — no scheduling logic, no ownership mapping,
no shard merging.

**Naming and persistence** are handled by a thin persistence layer
that the solver calls through:
- `persistence.is_saved(point)` → bool
- `persistence.is_converged(point)` → bool
- `persistence.record_status(point, converged)` → void
- `persistence.save_state(point, response_vector)` → void
- `persistence.load_state(point)` → response_vector

The orchestrator doesn't touch JSON directly. It doesn't know about
file formats. It asks the persistence layer "is this done?" and
"save this."

**Restart** is handled by checking persistence before solving:
```
for state in states.linear:
    if not persistence.needs_solving(state, protocol):
        skip
    guess = persistence.best_available_guess(state, protocol)
    result = solve_fd_response(world, ground, state, protocol, guess)
    persistence.save(state, protocol, result)
```

This gives you a working, testable, readable orchestrator that
produces correct results. It is the reference implementation.

### Stage B: State-Parallel Orchestrator (add on top)

Once Stage A works and is validated, add state-parallel as a wrapper:

```
run_response_parallel(world, params):
    ground = load_ground_context(world, params)
    states = plan_required_states(world, params, ground)
    plan   = plan_parallel_execution(world, states, params)

    # Same three stages, but work is distributed
    for protocol in protocol_thresholds:
        if plan.use_subgroups(protocol):
            with_subworld(world, plan.num_groups) as subworld:
                my_work = plan.work_for(subworld.group_id, protocol)
                local_ground = load_ground_context(subworld, params)
                set_protocol(subworld, local_ground, protocol)
                for state in my_work:
                    solve_fd_response(subworld, local_ground, state, protocol)
                # merge metadata shards
            merge_metadata_shards(world, plan.num_groups)
        else:
            # fall back to serial
            set_protocol(world, ground, protocol)
            for state in states.linear:
                solve_fd_response(world, ground, state, protocol)

    properties = compute_properties(world, ground, states, params)
    return properties
```

The solver functions are identical — they don't know whether they're
running on the universe or a subworld. The parallel layer only handles:
- Deciding how many groups and which work goes where
- Creating/destroying subworlds
- Merging metadata after each protocol step

**Testing:** Run the same calculation with both orchestrators. Results
must match within numerical tolerance. Any discrepancy is a bug in
the parallel layer, not the solver.

### What This Eliminates from Current Code

The following structures in the current MolresponseLib become unnecessary
or are replaced by much simpler equivalents:

- ProtocolExecutionPolicy → replaced by simple serial loop (Stage A)
  or a thin plan object (Stage B)
- RuntimePointOwnershipPolicy → same
- StateSolveScheduleContext → same
- PendingPointWorkItem / PendingProtocolManifest → same
- build_protocol_execution_policy → not needed in Stage A
- compute_runtime_point_ownership_policy → not needed in Stage A
- build_pending_work_manifest → not needed in Stage A
- build_pending_manifest_from_metadata → not needed in Stage A
- All the claim-file / done-file coordination → simplified persistence
- FilteredLineStreambuf / ScopedRankLogRedirect → deferred to Stage B
- merge_state_metadata_json → deferred to Stage B
- The 500+ line execute_subgroup_state_solve → deferred to Stage B

Stage A should be achievable in hundreds of lines, not thousands.

---

## What the Agent Does vs What You Do

**Agent work (tedious, well-defined, high-volume):**
- Phase 1: Inventory every function in the legacy code
- Phase 2: Inventory every function in molresponseLib, flag mixed concerns
- Phase 3: Produce the side-by-side mapping
- Write boilerplate test scaffolding once test strategies are defined
- Run tests and report pass/fail

**Your work (judgment, design, domain knowledge):**
- Review and correct the agent's inventories
- Make design decisions about interfaces and building block boundaries
- Decide what the clean algorithm skeleton should look like
- Choose test molecules and reference values
- Evaluate whether test results are physically meaningful
