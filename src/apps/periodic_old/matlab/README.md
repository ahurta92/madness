# Periodic MATLAB experiments

## Overview

MATLAB scratch calculations accompanying the dormant periodic solver.
This directory groups the child experiment; it has no executable entry point.

## Interface

- [Cosine-potential experiment](./cosine_potential/README.md): direct-sum
  Coulomb calculation in [fd_coulomb_conv.m](./cosine_potential/fd_coulomb_conv.m).

## Invariants

The script is a historical experiment without registered numerical validation;
its child contract records the unguarded self-distance singularity.

## Side Effects

Running the child script clears the MATLAB workspace and figures, prints progress,
and populates density and potential arrays. This directory adds no wrapper around it.

## Verification

No CMake target or CTest registration exists for this subtree.
