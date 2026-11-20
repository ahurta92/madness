# cosine_potential component

## Overview

Dormant MATLAB direct-sum Coulomb experiment.

## Interface

- fd_coulomb_conv.m projects a Gaussian density on a 20^3 grid and explicitly sums density/distance pairs. [fd_coulomb_conv.m:nsize](./fd_coulomb_conv.m)

## Invariants

- No singular self-distance guard is implemented; this is an experimental script, not a validated Coulomb solver. [fd_coulomb_conv.m:dist](./fd_coulomb_conv.m)

## Side Effects

- Clears MATLAB workspace and figures; prints progress and populates rho/rv arrays. [fd_coulomb_conv.m:clear](./fd_coulomb_conv.m)

## Source files

| File | Role |
| --- | --- |
| [fd_coulomb_conv.m](./fd_coulomb_conv.m) | Direct Coulomb sum experiment. |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

No CMake or CTest registration.
