# Molresponse Refactor Progress

Date: 2026-03-23

## Scope of this checkpoint

This checkpoint captures the recent `molresponse_v2` cleanup work around the restricted response and excited-state code paths. The main goals were:

- split the new restricted excited-state bundle workflow into model-specific ops headers
- improve runtime diagnostics for the restricted excited-state path
- move solver code away from raw `.flat` manipulation toward an `all()/x()/y()` API
- remove duplicated x-only storage for `StaticRestrictedResponse` and `TDARestrictedResponse`

## Completed changes

### 1. Restricted excited-state bundle split

The new `ExcitedResponse` path now uses a restricted bundle-ops layout under `src/apps/molresponse_v2/excited_ops/`:

- `RestrictedBundleCommon.hpp`
- `RestrictedTDABundleOps.hpp`
- `RestrictedFullBundleOps.hpp`

This split pulls the restricted TDA and restricted full bundle algorithms out of the old monolithic flow so the model-specific excited-state math lives close to the model type. Shared orchestration stays in `ExcitedResponse` and shared bundle helpers stay in `RestrictedBundleCommon.hpp`.

Related integration changes landed in:

- `src/apps/molresponse_v2/ExcitedResponse.cpp`
- `src/apps/molresponse_v2/ExcitedResponse.hpp`
- `src/madness/chem/MolresponseLib.hpp`
- `src/apps/molresponse_v2/CMakeLists.txt`

### 2. Restricted excited-state diagnostics

Additional debug logging was added so the restricted excited-state solver reports where it is inside the bundle iteration. The logging changes include:

- named phase/timing markers in the restricted bundle workflow
- state-indexed phase markers for `compute_potentials`, `build_theta`, `bsh_update`, residual construction, KAIN, and projection/normalization
- extra potential diagnostics around `<x|gamma>`, `<x|v0>`, and `<x|lambda>`
- test-side wiring so `test_tda_h2` writes a debug log JSON file

Related files:

- `src/apps/molresponse_v2/ResponseDebugLogger.hpp`
- `src/apps/molresponse_v2/test_tda_h2.cpp`
- `src/apps/molresponse_v2/excited_ops/RestrictedBundleCommon.hpp`

This instrumentation was used to narrow the 32-rank H2 TDA failure boundary to the first diagnostic inner product in `print_potential_diagnostics(...)`, after `compute_potentials` and before the space rotation.

### 3. ResponseVector accessor migration

The restricted solver code now has a clearer API for channel semantics versus solver/storage representation. `ResponseVector.hpp` now exposes:

- `all()` for the concatenated solver/storage view
- `x()` for the x-channel view
- `y()` for y-capable restricted types
- generic helpers `response_all(...)`, `response_x(...)`, `response_y(...)`
- `assign_all_and_sync(...)` as the preferred replacement for direct `flat` assignment

Shared kernel helpers now include:

- `ensure_initialized_all(...)`
- `clone_with_all(...)`

The following restricted code paths were migrated onto the new accessor layer:

- `src/apps/molresponse_v2/ResponseKernels.hpp`
- `src/apps/molresponse_v2/ResponseBundle.hpp`
- `src/apps/molresponse_v2/FrequencyLoop.hpp`
- `src/apps/molresponse_v2/ops/DynamicRestrictedOps.hpp`
- `src/apps/molresponse_v2/excited_ops/RestrictedBundleCommon.hpp`
- `src/apps/molresponse_v2/excited_ops/RestrictedTDABundleOps.hpp`
- `src/apps/molresponse_v2/excited_ops/RestrictedFullBundleOps.hpp`

The archive layout and external `flat`-named surface were preserved for compatibility. `assign_flat_and_sync(...)` remains as a compatibility wrapper during the transition.

### 4. X-only restricted storage dedup

`StaticRestrictedResponse` and `TDARestrictedResponse` no longer store duplicated `x_alpha` and `flat` vectors. They now inherit from a shared alias-backed storage helper so both names refer to the same underlying vector.

Consequences:

- x-only restricted states still present `x_alpha` and `flat` to existing code
- `sync()` and `flatten()` are now intentional no-ops for those two types
- the all-channel representation for x-only restricted types is still exactly `[x_alpha]`

One legacy helper in `ExcitedStateBundleSolver.cpp` had to be updated because it used member pointers like `&ResponseType::x_alpha`, which do not work with reference-backed alias storage. That helper now takes accessors instead.

## Verification performed

The following verification was run from the debug build tree:

```bash
source setenv.sh
ninja -C /gpfs/projects/rjh/adrian/development/madness-worktrees/builds/molresponse-feature-next/debug madqc test_tda_h2
```

Focused runtime smoke test:

```bash
source /gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/setenv.sh
/gpfs/projects/rjh/adrian/development/madness-worktrees/builds/molresponse-feature-next/debug/src/apps/molresponse_v2/test_tda_h2 \
  --num-states 2 --max-iter 1 --print-level 3 .
```

Working directory used:

- `molresponse_testdata/results/h2_tda_test`

Observed result:

- build succeeded
- serial TDA H2 smoke run exited successfully
- phase-timing output and `tda_excited_debug_log.json` were produced
- the x-only alias-backed storage change did not change the observed one-iteration serial behavior

## Remaining follow-up

The most obvious next cleanup steps are:

1. finish migrating the remaining restricted x-only helpers and ops to `response_x(...)` / `response_all(...)` rather than direct field access
2. add focused regression coverage for x-only alias semantics in addition to the end-to-end H2 smoke test
3. decide how far to push the accessor cleanup into the legacy `ExcitedStateBundleSolver.cpp` path versus keeping that file in compatibility mode
