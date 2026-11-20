# mp2 component

## Overview

Standalone MP2 correlation-energy driver.

## Interface

- MP2::get_hf().value -> MP2::value; help, print_parameters, TT=TT_2D or TT_TENSORTRAIN. [mp2.cc:main](./mp2.cc)
- Active installed executable mp2 linked to MADchem. [CMakeLists.txt:add_mad_executable](./CMakeLists.txt)

## Invariants

- USE_GENTENSOR required; sets 6D tensor type and apply randomization; disabled builds only print guidance. [mp2.cc:main](./mp2.cc)

## Side Effects

- Prints HF, correlation and sum; caught exceptions still return zero. [mp2.cc:main](./mp2.cc)

## Source files

| File | Role |
| --- | --- |
| [mp2.cc](./mp2.cc) | Executable entrypoint; initialization, CLI dispatch, solver sequencing, reporting. |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

```sh
cmake --build <build-dir> --target mp2
```

No enabled local CTest registration; use qctest cases where applicable.
