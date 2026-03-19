# Excited-State Reintegration Checklist (Concrete, File-Mapped)

## Goal

Restore legacy excited-state capability inside `molresponse_v2` without
reintroducing `X_space` all-state container design.

## Legacy Algorithm Sketch (What must be preserved)

Primary legacy path:
- `src/apps/molresponse/ExcitedResponse.cpp`
- `src/apps/molresponse/ResponseBase.cpp`

Per protocol, per iteration (legacy intent):
1. Prepare/orthonormalize trial space and project out occupied space.
2. Build response potentials:
   - `T0X`, `V0X`, `E0X`, `gamma`
   - (`compute_response_potentials`, `compute_V0X`, gamma builders)
3. Rotate/deflate excited space (`rotate_excited_space`, generalized eigen
   update) to obtain updated root ordering and `omega`.
4. Build `theta` residual-driving terms and apply BSH update operators
   (`bsh_update_excited`).
5. Form residual (`update_residual`), perform KAIN update
   (`kain_x_space_update`), optional step restriction, normalize.
6. Check convergence using density + residual criteria, then save/restart.

## New-Architecture Mapping

- Legacy monolithic driver:
  `ExcitedResponse::iterate`
  -> `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`

- Legacy protocol/restart orchestration:
  `calc_runner` + `ResponseBase`
  -> `src/madness/chem/MolresponseLib.hpp` stage-2c (`execute_excited_state_bundle_stage`)

- Legacy checkpoint content:
  `save/load` in `ExcitedResponse`
  -> protocol restart snapshots in
  `ExcitedStateBundleSolver.cpp` + metadata in `response_metadata.json`

- Legacy root identity:
  implicit by vector index/order
  -> explicit protocol root manifest (`state_names`, `roots`) in metadata.

## Execution Checklist

### Phase A: Metadata and Identity (protocol boundary only)

- [x] Add protocol root naming payload to solver result:
  - file: `src/apps/molresponse_v2/ExcitedStateBundleSolver.hpp`
- [x] Persist/load root names in protocol restart snapshots:
  - file: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- [x] Assign root names at protocol boundaries (not every iteration), with
      grouping tolerance `10 * threshold` and suffix scheme `a/b/c...`:
  - file: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- [x] Store protocol root manifest (`roots`: `root_index`, `name`, `energy`)
      in stage metadata and protocol events:
  - file: `src/madness/chem/MolresponseLib.hpp`
- [x] Merge root manifest fields across subgroup shards:
  - file: `src/madness/chem/MolresponseLib.hpp`

### Phase B: Solver Correctness (legacy numerical core)

- [ ] Replace scaffold per-state update with legacy-equivalent coupled terms:
  - file: `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
  - port from:
    - `src/apps/molresponse/ExcitedResponse.cpp` (`update_response`, `bsh_update_excited`)
    - `src/apps/molresponse/ResponseBase.cpp` (`compute_response_potentials`, residual/KAIN)
- [ ] Introduce explicit convergence gates mirroring legacy criteria
      (density and residual targets) with protocol-aware thresholds.
- [ ] Add per-iteration diagnostics for energy/BSH/density residuals in
      machine-readable metadata fields.

### Phase C: State Model and Downstream Compatibility

- [ ] Introduce explicit excited-state descriptor objects (state identity decoupled
      from array index) in `molresponse_v2` state model.
  - files:
    - `src/apps/molresponse_v2/ResponseState.hpp`
    - `src/apps/molresponse_v2/StateGenerator.hpp`
- [ ] Ensure property workflows can consume finalized excited-state roots from
      metadata/manifests (Raman/derived-state integration path).
  - files:
    - `src/apps/molresponse_v2/PropertyManager.hpp`
    - `src/madness/chem/MolresponseLib.hpp`

### Phase D: Validation

- [ ] Fixture comparison against Dalton for `H2 d-aug-cc-pVQZ`, `32` roots.
- [ ] Verify protocol ladder behavior (`1e-4 -> 1e-6`) with restart.
- [ ] Confirm root naming stability across restart and protocol transitions.

## Current Status

The architecture/plumbing side (protocol root identity + manifest metadata) is
implemented. The remaining critical gap is numerical parity with the legacy
excited-state update kernel.
