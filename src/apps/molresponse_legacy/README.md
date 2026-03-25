# molresponse_legacy

## Status

`src/apps/molresponse_legacy/` is a frozen compatibility baseline imported from
`bsundahl/madness` (`src/apps/molresponse`). Its purpose is to preserve a
working reference for closed-shell excited-state response while equivalent
behavior is re-established in `src/apps/molresponse_v2/`.

This directory is not the place for new feature work. New development belongs
in `molresponse_v2`.

## Freeze Contract

- Do not add new algorithms, new features, or new workflow plumbing here.
- Keep edits to the minimum needed to make the imported code compile and run in
  the current MADNESS tree.
- Every local deviation from the imported source must be marked `LEGACY_PATCH`
  and should be explainable as one of:
  - build/include compatibility
  - current MADNESS API compatibility
  - current archive-format compatibility
  - correctness fix needed to make the imported solver run at all
- If a numerical behavior change is desired, implement it in `molresponse_v2`
  and validate against this baseline instead of rewriting the legacy copy.

## Supported Scope

The checked-in target is a partial but runnable legacy import aimed at the
closed-shell excited-state path.

Included in the build:
- `TDDFT.cc`
- `iterate_excited.cc`
- `iterate_gamma.cc`
- `iter_freq2.cc`
- the supporting density, property, operator, and driver files needed by that
  path

Explicitly excluded from the build because required upstream pieces are missing
from the imported source set:
- `property_functions.cc`
- `iterate_freq.cc`
- `iterate_xy.cc`

Because of those exclusions, `molresponse_legacy` should be treated as a
validated reference subset, not a complete historical re-host of all legacy
molresponse functionality.

## Local Compatibility Patches

The current baseline already contains a small set of documented local patches:
- `legacy_preamble.h`: forced namespace preamble for the imported code style
- `response_parameters.h`: current `QCCalculationParametersBase` shim and
  public input wrapper
- `ground_parameters.h`: current moldft archive read order (`version 4`)
- collective-call fixes around `inner()` / `norm2()` debug printing
- a load-balance out-of-bounds fix in `TDDFT.cc`
- the collective input/archive fix in `molresponse.cc`

This means the baseline is a frozen, reviewed compatibility port, not a byte-
for-byte fossil.

## Validation

Build in a configured tree after sourcing the environment:

```bash
source setenv.sh
ninja -C /gpfs/projects/rjh/adrian/development/madness-worktrees/builds/molresponse-feature-next/debug molresponse_legacy
```

Current binary path after a successful debug build:

```bash
/gpfs/projects/rjh/adrian/development/madness-worktrees/builds/molresponse-feature-next/debug/src/apps/molresponse_legacy/molresponse_legacy
```

Recommended next lock-down step:
- calibrate at least one small closed-shell excited-state case against this
  binary and wire the result into an automated regression check
- `src/apps/madqc_v2/test_tda_h2_smoke.py` already has a placeholder hook via
  `LEGACY_MOLRESPONSE`, but the legacy reference energies are still unset

## Relationship To molresponse_v2

Use this directory to answer questions like:
- what did the legacy restricted excited-state update sequence do?
- what diagnostics and convergence gates existed?
- what output should a small closed-shell excited-state case approximately
  reproduce?

Do not reintroduce the legacy container/layout model wholesale into
`molresponse_v2`. Port the numerical behavior, not the public representation.
