# znemo component

## Overview

Complex regularized Hartree–Fock executable.

## Interface

- Znemo::value or analyze using read_reference; help, print_parameters, analyze. [znemo.cc:main](./znemo.cc)
- Active installed executable znemo linked to MADchem. [CMakeLists.txt:add_mad_executable](./CMakeLists.txt)

## Invariants

- Requires USE_GENTENSOR disabled; forces no_orient before molecule construction; optimization replacement not supported. [znemo.cc:main](./znemo.cc)

## Side Effects

- analyze reads saved reference and copies alpha/beta orbitals; prints final energy/time. [znemo.cc:main](./znemo.cc)

## Source files

| File | Role |
| --- | --- |
| [znemo.cc](./znemo.cc) | Executable entrypoint; initialization, CLI dispatch, solver sequencing, reporting. |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

Use a configured build with ENABLE_GENTENSOR=OFF and BUILD_TESTING=ON. [CMakeLists.txt:NOT ENABLE_GENTENSOR](./CMakeLists.txt)

```sh
cmake --build <build-dir> --target znemo
```

```sh
ctest --test-dir <build-dir> -R "^madness/test/(znemo/build|scripted_tests/znemo/test_energy.py/run)$"
```
