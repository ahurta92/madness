# zcis component

## Overview

Complex-orbital CIS executable.

## Interface

- Znemo::value -> Zcis::value; help and print_parameters. [zcis.cc:main](./zcis.cc)
- Active installed executable zcis linked to MADchem. [CMakeLists.txt:add_mad_executable](./CMakeLists.txt)

## Invariants

- Complex ground-state reference evaluated before excited-state construction. [zcis.cc:main](./zcis.cc)

## Side Effects

- Prints calculation/exception information; delegates calculation to MADchem. [zcis.cc:main](./zcis.cc)

## Source files

| File | Role |
| --- | --- |
| [zcis.cc](./zcis.cc) | Executable entrypoint; initialization, CLI dispatch, solver sequencing, reporting. |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

```sh
cmake --build <build-dir> --target zcis
```

No enabled local CTest registration; use qctest cases where applicable.
