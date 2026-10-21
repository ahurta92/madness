# MADQC convergence evidence report

## Scope

This report summarizes the first controlled evidence slice. The source worktrees and numerical references were not modified. Forced trajectories remain diagnostic; they are not convergence calibration references.

## Completed evidence

- H₂: existing forced and HDF5-enabled response evidence, including coarse/fine ground-state orbital exports and a static plus excited-state response run.
- LiH: ordinary static response rerun with the HDF5 backend enabled. The static dipole response converged in 9 iterations with `max_bsh_res = 4.693172e-4`, `dρ = 3.999584e-4`, and `alpha_zz = 22.15485`.
- H₂O: an ordinary HDF5-enabled run reached the fine `1e-6`, `k=8` response phase and produced a native ground-state checkpoint, but was stopped after a time-bounded exploratory window before all response channels completed. Its log and partial artifacts are preserved and marked `interrupted-time-budget`.

## Field comparison

The H₂ coarse/fine structured orbital exports were read with `pymra` on a common `64^3` midpoint grid over `[-20,20]^3` bohr. The closed-shell density proxy `rho = 2*psi_alpha^2` had a sampled L2 difference of `4.7901983e-5`, pointwise RMS `1.8934921e-7`, and maximum absolute difference `2.8245135e-5`. The inner `r<3` contribution was `4.7865053e-5`; the outer contribution was `1.8806088e-6`.

This is a ground-state orbital-derived comparison. It does not establish a response-state density difference. Native application HDF5 files contain `madness_archive` byte streams and require an explicit structured export before `pymra` can read them.

## Interpretation

The evidence supports three conclusions so far:

1. The HDF5 response backend is operational for both H₂ and LiH. The LiH log explicitly shows response state reloads with `backend=hdf5`.
2. Ground-state endpoint differences are measurable and spatially localized; the H₂ coarse/fine difference is concentrated near the molecular core.
3. The current slice is insufficient to classify response convergence as reference-limited, response-solver-limited, or gate-limited. The response runs are not yet matched across multiple exhausted ground-state endpoints for the same molecule and response observable.

The water partial run shows that the realistic matrix is expensive. It should be rerun as a deliberately scoped experiment, likely one perturbation and one protocol comparison at a time, rather than treated as an all-at-once regression job.

## Next decision

Before changing solver gates, complete at least one matched closed-shell comparison (H₂ or water) and one open-shell comparison (LiH) where the same response observable is run from two exhausted ground-state endpoints. Compare response iteration count, residual trajectory, and final property. Keep total-density and spin-density reporting separate, retain per-function residual and kinetic contributions, and keep protocol progress distinct from explicit user criteria.
