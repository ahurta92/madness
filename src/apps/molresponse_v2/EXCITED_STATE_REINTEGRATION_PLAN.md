# Excited-State Reintegration Plan (molresponse_v2)

## Objective

Reintegrate excited-state solving into the new `molresponse_v2` execution flow so
linear frequency-dependent response states and excited-state bundles can coexist
in one stage-2 solver pipeline, with metadata/restart support and future support
for mixed properties (first target: two-photon absorption).

## Constraints

1. Linear response states are independent per `(channel, frequency, protocol)` and
   remain state-parallel eligible.
2. Excited states are solved as a coupled bundle of `M` roots and **must run as one
   unit on a single subgroup** (no root-level splitting).
3. We still want protocol-aware restart behavior and metadata parity with the
   existing linear-state path.

## Legacy Code To Reuse / Port

Primary old implementation:
- `src/apps/molresponse/ExcitedResponse.hpp`
- `src/apps/molresponse/ExcitedResponse.cpp`

Legacy parameter model:
- `src/apps/molresponse/response_parameters.h`

Legacy dispatch context:
- `src/apps/molresponse/molresponse.cc`
- `src/apps/molresponse/ResponseBase.hpp`

## New Architecture Fit

Current stage orchestration lives in:
- `src/madness/chem/MolresponseLib.hpp` (`solve_all_states`)

Current linear planning/state model:
- `src/apps/molresponse_v2/StateGenerator.hpp`
- `src/apps/molresponse_v2/ResponseState.hpp`
- `src/apps/molresponse_v2/ResponseRecord.hpp`
- `src/apps/molresponse_v2/StateParallelPlanner.hpp`

### Proposed Model Extension

Add two work classes to stage 2:
- `LinearChannelPoint` (existing behavior)
- `ExcitedStateBundlePoint` (new coupled solve unit)

A bundle point is keyed by:
- protocol
- excited configuration (e.g. `M`, `tda/full`, guess mode)

A bundle point result includes:
- excited energies `omega[0..M-1]`
- coupled state vectors (X/Y-like orbital response pair space)
- convergence/timing metadata

## Phase 1 (Scaffolding, no numerical behavior change)

1. Add excited-state request knobs to `ResponseParameters` (new format):
   - `response.excited.enable` (bool, default false)
   - `response.excited.num_states` (int)
   - `response.excited.tda` (bool)
   - `response.excited.guess_max_iter` (int)
   - `response.excited.maxiter` (int)
   - `response.excited.maxsub` (int)
   - `response.excited.owner_group` (int, default 0)
2. Extend planning output to include an optional excited bundle request set.
3. Extend metadata schema (`ResponseRecord2`) with `excited_states` section:
   - protocol-indexed `saved`, `converged`, `timings`, `energies`.
4. Add empty/stub `ExcitedStateSolver` interface under `molresponse_v2`.

## Phase 2 (Single-group excited bundle solve integration)

1. Port core iterative excited-state kernel into
   `src/apps/molresponse_v2/ExcitedStateSolver.[hpp/cpp]`.
2. Integrate into `solve_all_states` protocol loop:
   - run only on designated owner subgroup,
   - do not split roots across groups,
   - record status in metadata shards and merged metadata.
3. Add restart lookup for bundle points at same/nearest protocol.

## Phase 3 (Mixed property support: Two-Photon Absorption)

1. Add `two_photon_absorption` to requested property set.
2. Implement property stage that mixes:
   - linear dipole frequency-response states,
   - excited-state energies/vectors.
3. Keep component-first design so TPA components can later reuse dynamic
   claiming behavior similar to current beta/raman component precompute.

## Phase 4 (Open-shell follow-on)

After closed-shell path stabilizes:
- generalize excited bundle solve to unrestricted/open-shell variants,
- validate metadata/persistence compatibility in mixed-spin workflows.

## Scheduling Policy (Explicit)

- Linear stage: keep current protocol-aware state/point scheduling.
- Excited stage: one bundle task per protocol, pinned to one subgroup
  (`excited.owner_group`), with serial fallback when subgroup mode is off.
- Property stage: keep current property-group assembly behavior.

## Initial Implementation Slice (recommended first PR on this branch)

1. Parameter + metadata schema additions.
2. Planning and orchestration hooks in `solve_all_states`.
3. Stub solver with clear TODOs + compile-time plumbing.

This gives a safe incremental merge point before porting the full excited
iterative kernel.

## Acceptance Criteria For Reintegration

1. Existing linear-response workloads remain numerically unchanged.
2. Excited-state-only run can execute through `madqc --wf=response` using new
   knobs and produces protocol-tagged excited metadata.
3. Restart of excited bundle works for partially completed protocol ladders.
4. Subgroup execution honors single-owner constraint for excited bundle tasks.
