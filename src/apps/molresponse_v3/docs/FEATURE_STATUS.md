# molresponse / madqc — feature status

The living status of the modern MADNESS response stack (`molresponse` solver +
`madqc` driver). This is the **user-facing** source of truth for what is ready,
what is opt-in/experimental, and what is on the roadmap. Update it as features
graduate — bump the badge, flip the default, add the doc link.

**Legend**
- ✅ **Stable** — validated, on by default (or a stable opt-in flag); safe for production.
- 🟡 **Experimental** — works, opt-in and OFF by default; interface/behavior may change.
- 🔬 **In development** — not in this release; tracked for a follow-up.

## Features

| Feature | Status | Enable | Why you'd use it | Guide |
|---|---|---|---|---|
| Response solver (α, β) | ✅ Stable | default (`engine v3`) | Modern unified solver; α/β validated (H₂O α_zz≈8.53, β_zzz≈7.76), restart-safe, structured diagnostics | Getting started with molresponse |
| madqc orchestration | ✅ Stable | `engine v3` in the response block | One SCF→response driver; input decks unchanged; output is a v2 superset | `madqc/README.md`, `MIGRATION_FROM_V2.md` |
| State-parallel FD | ✅ Stable | `--fd-subworlds=G` (0 = single-world ref) | 2–2.56× faster + ~8× less per-rank memory → runs c6h6/naphthalene at k=10 that v2 OOMs | Scaling FD across nodes |
| FD exchange tensor | ✅ Stable (opt-in) | `--fd-tensor` | Shared convolution tensors + per-protocol g0 cache → fewer Poisson convolutions on the FD path | Exchange tensor path, `operator_contracts.md` |
| HDF5 restart / I/O | ✅ Stable (opt-in, closed-shell) | build `-DMADNESS_ENABLE_HDF5`; run `MADRESPONSE_IO_HDF5=1` | 15–40% faster restart reads; backward-compatible; optional gzip | Accelerating restart with HDF5 |
| MRA visualization export | ✅ Stable (NP=1) | build `-DMADNESS_ENABLE_VTK`; `dump_mra_trees --htg/--amr/--coeffs` | Inspect trees + response orbitals/ρ⁽¹⁾ in ParaView; native `.mad.h5` coefficient archive | Visualizing MRA trees & response orbitals |
| Performance cost model | 🟡 Experimental | build `-DENABLE_WORLD_PROFILE`; `perf_model_fit.py` | Per-phase profile → predict wall-time / find bottlenecks for large runs | Profiling & the cost model |
| State-parallel auto-select | 🟡 Experimental | — (manual `--fd-subworlds` is stable) | Pick the subworld count automatically from problem shape + nodes | (roadmap) |

## Response properties

| Property | Status | Notes |
|---|---|---|
| Polarizability α (static + dynamic) | ✅ Stable | validated vs v2 / reference |
| First hyperpolarizability β | ✅ Stable | 2n+1 / VBC path; all VBC contributions must climb protocols |
| Single-component Raman | 🟡 Experimental | β(dipole; dipole, nuclear) for one (atom, axis); runs end-to-end, reference validation pending |
| Excited states (TDA / Full RPA) | 🟡 Experimental | eigenpairs compute; convergence/climb stabilizing (VirtualAO guess) |
| Full-tensor per-atom Raman | 🔬 In development | v2-only today; gated on further v3 work |
| Two-photon absorption (2PA) | 🔬 In development | kernel designed; gated on ES convergence |
| Resonance Raman | 🔬 In development | gated on ES convergence |

## Roadmap (deferred to follow-up PRs)

- Exchange Inc-3 batching/tiling/parallelization (doc 33).
- State-parallel: auto-selector, multi-node-per-state, sub-node packing, ES state-parallel.
- HDF5: interop value+coords (Layer B), parallel HDF5 / MPI-IO, open-shell states.
- Visualization: ParaView reader plugin for `.mad.h5` (native format is frozen).
- Properties: full-tensor Raman, 2PA, resonance Raman.
