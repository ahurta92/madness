# pno component

## Overview

Pair-natural-orbital MP2/F12 executable.

## Interface

- Nemo reference, PNOParameters(pno), F12Parameters(f12), PNO::solve(all_pairs). [pno.cpp:main](./pno.cpp)
- Active installed executable pno linked to MADchem. [CMakeLists.txt:add_mad_executable](./CMakeLists.txt)

## Invariants

- Reference is solved before pair correlation; total energy printed only when MP2 correlation energy is nonzero. [pno.cpp:main](./pno.cpp)

## Side Effects

- Reports SCF and pair energies, elapsed times, mean/maximum PNO rank, runtime statistics. [pno.cpp:main](./pno.cpp)

## Source files

| File | Role |
| --- | --- |
| [pno.cpp](./pno.cpp) | Executable entrypoint; initialization, CLI dispatch, solver sequencing, reporting. |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

```sh
cmake --build <build-dir> --target pno
```

No enabled local CTest registration; use qctest cases where applicable.
