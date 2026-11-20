# kernels component

## Overview

Header-only linear-response operators, response-subspace algebra, quadratic sources, and property contractions shared by FD and ES solvers.

## Interface

- Static/TDA select X-only storage and Full selects XY storage, independently tagged ClosedShell/OpenShell [tags.hpp:StorageOf](./tags.hpp).
- Kernels specializations provide density, gamma, potential/Fock, BSH and residual operations; TDA in both shells and Full closed shell additionally satisfy ESKernel [kernel_interface.hpp:FDKernel](./kernel_interface.hpp) [kernel_interface.hpp:ESKernel](./kernel_interface.hpp) [tda.hpp:MV3_ASSERT_ES_KERNEL](./tda.hpp) [full.hpp:MV3_ASSERT_ES_KERNEL](./full.hpp).
- ResponseGroundState supplies occupied orbitals/energies, local potentials, Fock blocks, Coulomb/projectors, optional closed-shell XC response and exchange caches [tda.hpp:ResponseGroundState](./tda.hpp).
- rs provides Euclidean inner products, RPA metric, kinetic Gram matrices, slot-preserving diagonalization, transformations and orthonormalization [response_space_ops.hpp:diagonalize](./response_space_ops.hpp) [response_space_ops.hpp:metric](./response_space_ops.hpp).
- source_spec assembles Apply and OccupiedMatrix entries; vbc computes quadratic sources; beta and tpa contract them; raman projects tensors onto modes and calculates rotational invariants [source_spec.hpp:SourceEntry](./source_spec.hpp) [source_spec.hpp:assemble_source](./source_spec.hpp) [vbc.hpp:compute_vbc](./vbc.hpp) [beta.hpp:beta_abc](./beta.hpp) [tpa.hpp:observables](./tpa.hpp) [raman_spectrum.hpp:project_on_modes](./raman_spectrum.hpp).

## Invariants

- Static means Y=X while TDA omits Y: closed-shell density factors are 4 for Static, 2 for TDA and 2 multiplying X+Y for Full [static.hpp:compute_density](./static.hpp) [tda.hpp:compute_density](./tda.hpp) [full.hpp:compute_density](./full.hpp).
- The Full response metric is <X|X>-<Y|Y>, with the same signs on beta blocks; kinetic Gram uses a positive gradient form over stored blocks [response_space_ops.hpp:metric_inner](./response_space_ops.hpp) [response_space_ops.hpp:kinetic_gram](./response_space_ops.hpp).
- Subworld callers must pass their World to rs::diagonalize; its compatibility default uses World::get_default. Output slot order preserves diagonal dominance, phase and near-degenerate continuity; ascending energy is separate [response_space_ops.hpp:diagonalize](./response_space_ops.hpp).
- Exchange leg {bra,ket} means |ket><bra|. Source entries carrying exchange legs must also supply Coulomb density with the closed-shell factor already included; assemble_source leaves final truncation to the caller [source_spec.hpp:LegPair](./source_spec.hpp) [source_spec.hpp:detail::eval_fock_action](./source_spec.hpp) [source_spec.hpp:assemble_source](./source_spec.hpp).
- VBC and beta throw for open-shell instantiations; their second XC kernel is not implemented, and DFT quadratic calculations are refused by the workflow [vbc.hpp:compute_vbc_spec](./vbc.hpp) [beta.hpp:beta_abc](./beta.hpp) [../orchestrator/response_workflow.hpp:run_response_with_ground](../orchestrator/response_workflow.hpp).
- Cached ground exchange and g0 tensors depend on the current protocol and must be rebuilt after ground-state re-preparation [common_ops.hpp:make_ground_exchange](./common_ops.hpp) [tda.hpp:ResponseGroundState](./tda.hpp).

## Side Effects

- Numerical operations create, compress, project, truncate and transform distributed Function trees and invoke collective fences/reductions; rs::inner compresses its shared underlying handles even through const views [response_space_ops.hpp:inner](./response_space_ops.hpp) [response_space_ops.hpp:transform](./response_space_ops.hpp) [assembly.hpp:assemble_theta](./assembly.hpp).
- Exchange configuration reads MADRESPONSE_HFEXALG (actual default multiworld_row), MADRESPONSE_HFEX_GRANULARITY, MADRESPONSE_HFEX_ACCUMULATION and MADRESPONSE_HFEX_COST_AWARE [common_ops.hpp:configure_exchange](./common_ops.hpp).
- Diagnostic residue contractions and rank-deficient diagonalization print timing/warnings on rank 0; this directory does not write checkpoint files [tpa.hpp:tpa_moment_residue](./tpa.hpp) [tpa.hpp:tpa_e3_residue](./tpa.hpp) [response_space_ops.hpp:diagonalize](./response_space_ops.hpp).

## Source files

| File | Role |
| --- | --- |
| [assembly.hpp](./assembly.hpp) | State-generic theta/lambda linear combinations |
| [beta.hpp](./beta.hpp) | Closed-shell 2n+1 beta contraction |
| [common_ops.hpp](./common_ops.hpp) | Shared kinetic, BSH and configurable exchange operators |
| [exchange_ctx.hpp](./exchange_ctx.hpp) | Optional cached pair-tensor exchange path |
| [full.hpp](./full.hpp) | XY full-response kernels; closed-shell ES and both-shell FD |
| [kernel_interface.hpp](./kernel_interface.hpp) | Compile-time kernel shape validation |
| [raman_spectrum.hpp](./raman_spectrum.hpp) | Pure normal-mode projection and Raman rotational invariants |
| [response_space_ops.hpp](./response_space_ops.hpp) | Response subspace metrics, rotations, kinetic Gram and diagonalization |
| [source_spec.hpp](./source_spec.hpp) | Declarative perturbed-Fock source evaluation |
| [static.hpp](./static.hpp) | X-only static FD kernels in both shells |
| [tags.hpp](./tags.hpp) | Response/shell tags and storage mapping |
| [tda.hpp](./tda.hpp) | Prepared kernel ground state and TDA kernels in both shells |
| [tda_batch.hpp](./tda_batch.hpp) | Closed-shell TDA bundle batching and optional tensor gamma |
| [tpa.hpp](./tpa.hpp) | TPA contractions and orientational observables |
| [tpa_source_spec.hpp](./tpa_source_spec.hpp) | State-independent photon-pair residue specifications |
| [two_electron.hpp](./two_electron.hpp) | Raw or projected Coulomb-minus-exchange application |
| [vbc.hpp](./vbc.hpp) | Closed-shell quadratic source used by beta and Raman |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

```sh
cmake --build <build-dir> --target test_bsh_shift test_diagonalize_slots test_raman_spectrum test_exchange_ctx
```

Registration evidence: [../CMakeLists.txt:BUILD_TESTING](../CMakeLists.txt)

```sh
ctest --test-dir <build-dir> -R "^madness/test/molresponse/(test_bsh_shift|test_diagonalize_slots|test_raman_spectrum|test_exchange_ctx)/run$"
```

Fixture-dependent: test_kernel_equivalence, test_vbc_spec_equivalence, test_tpa_pq_vs_vbc, test_vbc_pq_terms, test_tpa_regroup_identities (MADNESS_RESPONSE_FIXTURE_DIR) [../CMakeLists.txt:MADNESS_RESPONSE_FIXTURE_DIR](../CMakeLists.txt)
