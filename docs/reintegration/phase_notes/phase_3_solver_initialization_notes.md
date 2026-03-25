# Phase 3 Notes - Solver Initialization Path

Date:

- 2026-03-13

Related roadmap phase:

- Phase 3 — Solver Initialization Path

Related contract:

- [`../contracts/phase_3_solver_initialization_contract.md`](../contracts/phase_3_solver_initialization_contract.md)

## What Was Implemented

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_trial(...)` now always refreshes
  `trial.num_states` and `trial.num_orbitals` after refinement.
- Added explicit solver-local helpers:
  `sort_trial_space(...)`,
  `select_lowest_trial_roots(...)`,
  and `seed_active_bundle_from_trial_space(...)`
  in `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`.
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_fresh_guess(...)` now:
  - computes a requested root count explicitly
  - refines the generated trial space with `iterate_trial(...)`
  - sorts and truncates to the lowest requested roots
  - retries with a larger random trial if the first generator under-fills the request
  - seeds `active_response_bundle_` from the selected trial-space roots
  - emits `EXCITED_FRESH_GUESS_SELECT ... snapshot_kind=guess_bundle|guess_trial_seed`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: initialize_protocol_guess(...)` now reports
  fresh-start `snapshot_kind`, `bundle_state_present`, and `restart_capable`
  consistently with the seeded active bundle.
- `src/apps/madqc_v2/test_molresponse_excited_guess_fallback.py` now validates:
  - the explicit fresh-start selection marker
  - the persisted guess-archive snapshot contract
  - stable root ids and state names across forced guess-archive reuse
  - parity between workflow metadata and `response_metadata.json` after fallback reuse

## Files Changed

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`
- `src/apps/madqc_v2/test_molresponse_excited_guess_fallback.py`

## Key Symbols

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_fresh_guess`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_trial`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: sort_trial_space`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: select_lowest_trial_roots`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: seed_active_bundle_from_trial_space`
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: initialize_protocol_guess`
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::initialize`

## Planned Focus

- Compare the current fresh-start initialization path with legacy `ExcitedResponse::initialize(...)`.
- Make trial-bundle construction and truncation behavior explicit and reproducible.
- Preserve Phase 1 stable-root identity and Phase 2 restart/archive semantics.

## Actual Design Decisions

- Reused the current generator family
  (`create_trial_functions2`, `make_nwchem_trial`, `create_trial_functions`, `make_random_trial`)
  instead of forcing a literal legacy `2 * num_states` trial-bundle size in every case.
- Made fresh-start selection explicit by layering `iterate_trial(...)`, energy sorting, and lowest-root
  truncation on top of the existing generators rather than replacing them.
- Seeded `active_response_bundle_` directly from the selected trial roots through the existing
  `ResponseVector` variant constructors so the fresh-start path produces a real `guess_bundle`
  restart seed for restart-capable variants.
- Kept the Phase 2 restart support boundary unchanged:
  restricted variants remain `full_bundle_resume`,
  unrestricted variants remain explicitly non-full-restart.

## Deviations From The Original Plan

- The implementation did not force exact legacy oversubscription semantics for every case.
  Instead, it keeps the existing generator heuristics and only expands with `make_random_trial(...)`
  when the first guess does not provide enough candidate roots.
- This phase stopped at initialization-path parity and did not change:
  - the main iteration kernel
  - convergence gates
  - `maxsub` / KAIN behavior

## New Constraints Discovered

- Restart artifact filenames are driven by `config_.output_prefix`, so the authoritative archive names are
  `response.excited_bundle...`, not `molresponse.excited_bundle...`.
- In the current interactive 40-core allocation, `load_xeonmax.sh` is suitable for serial smoke tests
  but does not provide a usable MPI runtime for the full scripted excited-state test slice.
  The MPI-backed slice still needs `load_40core.sh`.

## Remaining Follow-Up

- Phase 4 still needs restricted-shell iteration, convergence, and `maxsub` / KAIN parity work.
- Cross-protocol naming policy and degeneracy handling remain Phase 5 work.
- If future legacy parity testing shows the exact oversized trial-bundle heuristic matters,
  revisit the generator sizing strategy then rather than widening this phase retroactively.

## Validation Summary

- Passed the focused excited-state scripted slice in the interactive 40-core job:
  `test_molresponse_excited_metadata_smoke.py`,
  `test_molresponse_excited_restart_reuse.py`,
  `test_molresponse_excited_protocol_projection.py`,
  `test_molresponse_excited_guess_fallback.py`,
  `test_molresponse_excited_restart_mode_compat.py`.
- Passed a serial excited-state metadata smoke inside the same interactive job using `load_xeonmax.sh`.

## Reconciliation Checklist

- [x] roadmap updated
- [x] status dashboard updated
- [x] active contracts updated
- [x] validation record written
