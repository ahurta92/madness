# interior_bc component

## Overview

Dormant embedded Dirichlet-boundary convergence example for four sphere problems.

## Interface

- embedded_dirichlet positional arguments k thresh problem eps penalty-expression mask. [embedded_dirichlet.cc:main](./embedded_dirichlet.cc)
- EmbeddedDirichlet selects SURFACE, DIRICHLET_RHS, EXACT or DOMAIN_MASK; ConstantSphere, CosineSphere, Y20Sphere, InhomoConstantSphere implement analytic problems; DirichletCondIntOp supplies GMRES operator. [test_problems.h:EmbeddedDirichlet](./test_problems.h) [test_problems.h:DirichletCondIntOp](./test_problems.h)

## Invariants

- Requires k>=4, thresh<=1e-4, problem 1..4, eps>0, positive parsed penalty, mask 1 or 2. [embedded_dirichlet.cc:main](./embedded_dirichlet.cc)
- GMRES uses 30 max iterations and 1e-5 residual/update thresholds; computes interior error and exact-point ratios. [embedded_dirichlet.cc:main](./embedded_dirichlet.cc)

## Side Effects

- Mutates FunctionDefaults, redistributes adaptive trees, prints surface/volume/error diagnostics and z-axis solution samples. [embedded_dirichlet.cc:main](./embedded_dirichlet.cc)

## Source files

| File | Role |
| --- | --- |
| [embedded_dirichlet.cc](./embedded_dirichlet.cc) | Command-line solve/report pipeline. |
| [test_problems.h](./test_problems.h) | Required production problem/operator declarations despite test_* filename. |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

No active root CMake target: src/apps/CMakeLists.txt comments out add_subdirectory(interior_bc). Local declaration exists but is not configured.
