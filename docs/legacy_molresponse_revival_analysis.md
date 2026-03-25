# Legacy `molresponse` Revival Analysis

Line references in this note refer to the current checkout in `/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next` as inspected on 2026-03-16.

## Scope

This note answers a narrower question than the legacy excited-state report:

- What would be required to get the standalone legacy `src/apps/molresponse` code compiling and running again in the current tree?
- What must be fixed before it can serve as a usable excited-state comparison executable?

It does not propose a redesign. It records the current blockers and the minimum practical recovery path.

## Current status

The legacy source tree is still present:

- [src/apps/molresponse](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse)

The standalone target is not currently part of the top-level apps build:

- [src/apps/CMakeLists.txt:8](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/CMakeLists.txt#L8)

```cmake
#add_subdirectory(molresponse)
```

The legacy target definition itself still exists:

- [src/apps/molresponse/CMakeLists.txt](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/CMakeLists.txt)

It defines:

- `MADresponse`
- `molresponse`
- `MADresponse_base`
- `MADall_response`

So the first blocker is not source deletion. It is that the build has been intentionally disabled.

## Compile blockers confirmed in the current tree

I checked compile compatibility with syntax-only compiles of the legacy translation units using the current build include paths. The first hard failures are:

### 1. `ResponseParameters` no longer satisfies the base-class interface

`ResponseParameters` derives from `QCCalculationParametersBase`:

- [src/apps/molresponse/response_parameters.h:24](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/response_parameters.h#L24)

But it does not override the now-required pure virtual:

- [src/madness/mra/QCCalculationParametersBase.h:300](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/madness/mra/QCCalculationParametersBase.h#L300)

This makes every use of `ResponseParameters` ill-formed, including:

- [src/apps/molresponse/global_functions.h:17](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/global_functions.h#L17)
- [src/apps/molresponse/ResponseBase.hpp:656](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.hpp#L656)
- [src/apps/molresponse/ResponseBase.hpp:673](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.hpp#L673)

Minimum fix:

- add `std::string get_tag() const override { return tag; }` to `ResponseParameters`

### 2. The legacy exchange-operator API usage is stale

Legacy code uses `Exchange<double, 3>::Algorithm::...`:

- [src/apps/molresponse/ResponseBase.hpp:339](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.hpp#L339)
- [src/apps/molresponse/ResponseBase.cpp:219](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L219)
- [src/apps/molresponse/FrequencyResponse.cpp:551](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/FrequencyResponse.cpp#L551)
- [src/apps/molresponse/FrequencyResponse.cpp:1827](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/FrequencyResponse.cpp#L1827)

The current API exposes `ExchangeAlgorithm` values directly on `Exchange`:

- [src/madness/chem/SCFOperators.h:117](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/madness/chem/SCFOperators.h#L117)
- [src/madness/chem/SCFOperators.h:181](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/madness/chem/SCFOperators.h#L181)

Minimum fix:

- replace `Exchange<double, 3>::Algorithm::multiworld_efficient` with `Exchange<double, 3>::multiworld_efficient`
- same for `multiworld_efficient_row`, `large_memory`, `small_memory`

### 3. The standalone executable is disabled at the top-level build

Even if the compile issues above are fixed, CMake still will not build the legacy app until:

- [src/apps/CMakeLists.txt:8](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/CMakeLists.txt#L8)

is restored or otherwise made conditional.

Minimum fix:

- re-enable `add_subdirectory(molresponse)` or gate it behind a CMake option

## Runtime blockers for excited-state use

Compile fixes alone are not enough for the legacy excited-state path to run correctly.

### 4. `ExcitedResponse` never initializes `response_context`

The common operator infrastructure lives behind `ResponseBase::response_context`:

- [src/apps/molresponse/ResponseBase.hpp:429](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.hpp#L429)
- [src/apps/molresponse/ResponseBase.hpp:634](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.hpp#L634)

Those strategies are explicitly installed in `FrequencyResponse`:

- [src/apps/molresponse/FrequencyResponse.hpp:564](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/FrequencyResponse.hpp#L564)
- [src/apps/molresponse/FrequencyResponse.hpp:576](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/FrequencyResponse.hpp#L576)

But `ExcitedResponse` has only a forwarding constructor:

- [src/apps/molresponse/ExcitedResponse.hpp:19](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ExcitedResponse.hpp#L19)

The excited-state code still calls helpers that require the context to be set:

- `response_context.compute_j1(...)`
  - [src/apps/molresponse/ResponseBase.cpp:478](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L478)
  - [src/apps/molresponse/ResponseBase.cpp:720](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L720)
- `response_context.compute_k1(...)`
  - [src/apps/molresponse/ResponseBase.cpp:488](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L488)
  - [src/apps/molresponse/ResponseBase.cpp:744](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L744)
- `response_context.compute_VXC1(...)`
  - [src/apps/molresponse/ResponseBase.cpp:504](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L504)
  - [src/apps/molresponse/ResponseBase.cpp:733](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L733)
- `response_context.inner(...)`
  - [src/apps/molresponse/ResponseBase.cpp:2562](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L2562)

Without strategy installation, these paths throw runtime exceptions from `Context`.

Minimum fix:

- initialize `response_context` in the excited-state path, by analogy with `FrequencyResponse`
- for TDA excited-state runs, use the static strategy set
- for full excited-state runs, use the full strategy set
- use `LoadExcitedXSpace`, not `LoadFrequencyXSpace`

### 5. The full excited-state branch uses the wrong `Y`-channel switch

Several shared helpers decide whether to operate on both `x` and `y` channels by checking:

- [src/apps/molresponse/ResponseBase.cpp:1466](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L1466)
- [src/apps/molresponse/ResponseBase.cpp:1515](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L1515)
- [src/apps/molresponse/ResponseBase.cpp:1582](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L1582)
- [src/apps/molresponse/ResponseBase.cpp:1619](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L1619)

```cpp
bool compute_y = r_params.omega() != 0.0;
```

That works for frequency response, but it is not a correct excited-state discriminator. A full TDDFT/TDHF excited-state run is selected by `tda == false`, not by nonzero external driving frequency.

Minimum fix:

- replace the excited-path `compute_y` logic with something derived from `!r_params.tda()` or `calc_type`
- do not use `omega != 0.0` as the excited-state branch selector

Without this, the legacy full excited-state solver is internally inconsistent even if it compiles.

### 6. Restart/load metadata is read into temporary getter values

Legacy restart I/O writes and reads parameter metadata through getter return values:

- [src/apps/molresponse/ExcitedResponse.cpp:2767](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ExcitedResponse.cpp#L2767)
- [src/apps/molresponse/ExcitedResponse.cpp:2805](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ExcitedResponse.cpp#L2805)
- [src/apps/molresponse/ResponseBase.hpp:61](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.hpp#L61)
- [src/apps/molresponse/FrequencyResponse.cpp:390](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/FrequencyResponse.cpp#L390)

For output archives this is harmless. For input archives, the archive layer uses `const_cast` internally:

- [src/madness/world/archive.h:748](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/madness/world/archive.h#L748)

So these loads succeed syntactically but write into temporaries, not back into `r_params`.

Practical effect:

- archive metadata is not actually restoring `archive`, `tda`, `num_orbitals`, or `num_states`
- the restart logic is therefore relying on the already-loaded input parameters rather than the archived values

Minimum fix:

- read archive header fields into local variables and validate them explicitly
- or serialize/deserialize `ResponseParameters` itself instead of serializing getter return values

This is not the first thing blocking a fresh run, but it matters immediately if restart-based testing is part of the plan.

## Operational limitations that are not hard blockers, but matter for testing

### 7. The legacy CLI is standalone and file-name driven

The executable hardcodes `response.in`:

- [src/apps/molresponse/molresponse.cc:71](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/molresponse.cc#L71)
- [src/apps/molresponse/global_functions.cc:14](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/global_functions.cc#L14)

It expects:

- a local `response.in`
- a compatible `moldft.restartdata` archive path referenced from that file
- local file outputs like `response_base.json`
  - [src/apps/molresponse/ResponseBase.cpp:2145](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ResponseBase.cpp#L2145)
- a hard-coded initial guess archive `guess_restart`
  - [src/apps/molresponse/ExcitedResponse.cpp:111](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ExcitedResponse.cpp#L111)

For one-off testing this is acceptable. For systematic comparison against `molresponse_v2`, it is awkward.

Nice-to-have fix:

- accept `--input` or a positional input file
- make JSON and restart output prefix-aware instead of hard-coded

### 8. The legacy driver depends on the old standalone moldft-style archive contract

Ground-state loading is done through `GroundStateCalculation`:

- [src/apps/molresponse/global_functions.cc:20](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/global_functions.cc#L20)
- [src/apps/molresponse/ground_parameters.h:73](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/molresponse/ground_parameters.h#L73)

So a runnable comparison setup needs a compatible archive in the old format. This is likely still possible because the loader is still present, but it should be verified with an actual `moldft.restartdata` produced by the current tree before spending time on solver debugging.

## Minimum recovery path for excited-state comparison use

If the goal is specifically:

- “build legacy `molresponse` again”
- “run the old excited-state solver as a comparison executable”

then the minimum practical sequence is:

1. Re-enable the legacy app build in [src/apps/CMakeLists.txt](/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/apps/CMakeLists.txt).
2. Fix compile compatibility:
   - add `ResponseParameters::get_tag()`
   - update all stale `Exchange<double,3>::Algorithm::...` usages in `ResponseBase` and `FrequencyResponse`
3. Initialize `response_context` in `ExcitedResponse`
   - TDA path: static strategies
   - full path: full strategies
4. Fix excited-state `compute_y` selection so full excited-state runs do not depend on `omega != 0.0`
5. Repair restart header loading so archive metadata is read into real storage, not temporaries
6. Build and smoke-test the standalone executable with:
   - a minimal TDA excited-state input
   - a known-good ground-state archive
7. Only after TDA is running cleanly, attempt full excited-state runs

That ordering matters. Steps 1-2 are compile recovery. Steps 3-5 are runtime correctness recovery. Trying to compare legacy and `molresponse_v2` numerics before doing step 3 and step 4 would not be meaningful.

## Recommended validation ladder

To avoid mixing compile recovery with solver debugging, validate in this order:

1. `molresponse --help` / `--print_parameters`
2. static frequency-response smoke test
3. TDA excited-state smoke test with `states=1`
4. multi-root TDA excited-state test
5. full excited-state test (`tda=false`)
6. restart test for excited states

Reason:

- the frequency-response path already contains the needed `response_context` initialization and is therefore the easier runtime smoke test
- the full excited-state path is the least trustworthy part of the legacy code because of the current `compute_y` coupling to `omega`

## Bottom line

The legacy code is not far from being revivable, but it is not just “comment in the CMake line and build.”

There are three classes of work:

- Build enablement: re-enable the subdirectory.
- Mechanical compile repair: `get_tag()` plus exchange enum updates.
- Excited-state runtime repair: initialize `response_context`, fix `compute_y` semantics, and clean up restart header loading.

If those are done in that order, the old standalone `molresponse` should be recoverable as a comparison executable for excited-state testing without having to redesign it.
