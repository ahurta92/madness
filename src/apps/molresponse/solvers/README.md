# solvers component

## Overview

Iterative FD/ES solvers, convergence evidence, state archives, metadata, imported seeds and post-solve analysis.

## Interface

- FDSolver iterates driven states; ESSolver iterates coupled root bundles; iterate and iterate_protocol supply the rung loop and threshold ladder [fd_solver.hpp:FDSolver](./fd_solver.hpp) [es_solver.hpp:ESSolver](./es_solver.hpp) [iterate.hpp:iterate](./iterate.hpp) [iterate_protocol.hpp:iterate_protocol](./iterate_protocol.hpp).
- ResponseStateX/XY hold alpha and optional beta function blocks; handles are shallow while copy operations make independent trees [response_state.hpp:ResponseStateX](./response_state.hpp) [response_state.hpp:ResponseStateXY](./response_state.hpp).
- States are saved as nodes (states/response/<name>/) and starting guesses from elsewhere as seeds (seeds/response/<name>/, same archive form, kind "seed"); ResponseMetadata writes protocols, fd_states, excited_states, vbc_states, properties and provenance as output; nothing reads it to decide [shared_response_checkpoint.hpp:DrivenState](./shared_response_checkpoint.hpp) [shared_es_checkpoint.hpp:ExcitedSet](./shared_es_checkpoint.hpp) [response_seed.hpp:save_seed](./response_seed.hpp) [response_metadata.hpp:ResponseMetadata](./response_metadata.hpp).
- Analysis computes closed-shell transition properties and optional AO descriptions; full Raman assembles a nuclear Hessian and projects complete derivative rows onto normal modes [es_analysis.hpp:analyze_es_bundle](./es_analysis.hpp) [raman_analysis.hpp:assemble_raman_spectrum](./raman_analysis.hpp).
- run_dalton_import validates geometry and HF provenance and seeds requested FD/ES states; it never invokes Dalton [dalton_import.hpp:run_dalton_import](./dalton_import.hpp).

## Invariants

- Final criteria remain split: FD supplies requested accuracy plus independent saved-output reevaluation to shared rung_exit; ES bundle_rung_sample always selects engine checks and leaves named stage criteria to its caller [fd_solver.hpp:stop_reason](./fd_solver.hpp) [fd_solver.hpp:evaluate_saved_state](./fd_solver.hpp) [convergence_policy.hpp:bundle_rung_sample](./convergence_policy.hpp).
- FD engine checks require at least two iterations and all-channel BSH/density limits. Named BSH accuracy uses maximum raw per-function residual, not aggregate quadrature; genuine input and saved-output phases are separate [fd_solver.hpp:engine_check_passed](./fd_solver.hpp) [fd_explicit_accuracy.hpp:valid_fd_explicit_accuracy](./fd_explicit_accuracy.hpp) [function_residuals.hpp:measure_function_residuals](./function_residuals.hpp).
- ES engine convergence requires density and absolute omega-change limits after at least two iterations; BSH is explosion sanity only. Optional lock_converged retains freezing: consecutive-pass roots leave active subspace, BSH and KAIN work until protocol reset; locked and converged metadata are distinct [es_solver.hpp:es_root_converged](./es_solver.hpp) [es_solver.hpp:step](./es_solver.hpp).
- Root indices preserve bookkeeping identity through discrete sorting. Continuous mixing is not validated physical correspondence, so temporal kinetic deltas can be unavailable [es_root_identity.hpp:ExcitedRootDescriptor](./es_root_identity.hpp) [es_kinetic_evidence.hpp:es_kinetic_evidence](./es_kinetic_evidence.hpp) [es_solver.hpp:step](./es_solver.hpp).
- Rung tolerances scale from max(threshold,user dconv); intermediate requested-criterion rungs are nonfinal [convergence_policy.hpp:effective_for_thresh](./convergence_policy.hpp) [iterate_protocol.hpp:iterate_protocol](./iterate_protocol.hpp).
- A state starts from its own node, else a saved node at the nearest frequency, else its seed when the seed's settings (ground-state geometry and model, the state's physical settings) match, else the default guess. A seed with other settings is ignored with one line [response_seed.hpp:find_seed](./response_seed.hpp).
- Structured Function HDF5 requires HDF5, double functions and one rank; load restores cell and validates coefficient/tree shape. Opaque archive HDF5 supports serialized multirank state payloads [function_hdf5_io.hpp:load_function_hdf5](./function_hdf5_io.hpp).
- Dalton seeds are written only for states with no save of their own; a seed with matching settings is left, one with other settings replaced. Missing requested records warn and skip; nearest-frequency guesses require tolerance and are computed when the leg starts, not written. ES defaults zero y unless requested [dalton_import.hpp:seed_fd_from_dalton](./dalton_import.hpp) [dalton_import.hpp:dalton_fd_guess](./dalton_import.hpp) [dalton_import.hpp:seed_es_from_dalton](./dalton_import.hpp).
- Seed overlap guards are informational, not convergence or lowest-N guarantees; the seed is reloaded from its folder in the active representation [es_seed_guard.hpp:evaluate_es_seed_guard](./es_seed_guard.hpp).
- Full Raman needs all 9*3N derivative rows and nuclear legs; missing legs return no analysis and open-shell skips. Row convergence is reported separately from completeness [raman_analysis.hpp:assemble_raman_spectrum](./raman_analysis.hpp).
- Closed-shell DFT attaches f_xc; open-shell DFT response with nonunit exact exchange throws. ES Full builder is compile-time closed-shell-only; guess adapters reject near-zero trials [build_response_ground_state.hpp:build_response_ground_state_closed_shell](./build_response_ground_state.hpp) [build_response_ground_state.hpp:build_response_ground_state_open_shell](./build_response_ground_state.hpp) [build_response_ground_state.hpp:build_es_problem_full](./build_response_ground_state.hpp).

## Side Effects

- Iterations mutate trees, KAIN histories and protocol defaults, perform collectives, print diagnostics and checkpoint [iterate.hpp:iterate](./iterate.hpp) [fd_solver.hpp:step](./fd_solver.hpp) [es_solver.hpp:step](./es_solver.hpp).
- Node and seed archives are replaced with their description in one Node::save; response_metadata.json lives under calc_dir and is written through temporary files [response_seed.hpp:save_seed](./response_seed.hpp) [response_metadata.hpp:save](./response_metadata.hpp).
- Shard merge unions group FD/VBC entries, quarantines malformed JSON as .corrupt and removes merged shards. Replacing FD data invalidates dependent VBC entries [response_metadata.hpp:ResponseMetadata::merge_state_shards](./response_metadata.hpp) [response_metadata.hpp:set_fd_state](./response_metadata.hpp).
- Metrics collect global coefficient storage, worst-rank RSS, caller wall/iteration counts, and unavailable reasons [state_metrics.hpp:StateMetrics](./state_metrics.hpp) [function_metric_diagnostics.hpp:valid_metric_diagnostics](./function_metric_diagnostics.hpp) [es_kinetic_evidence.hpp:es_kinetic_evidence](./es_kinetic_evidence.hpp).
- Node subworlds group by hostname using SafeMPI split, clamped to minimum ranks/host; returned Worlds must die before universe/finalization [node_subworlds.hpp:make_subworld_pool](./node_subworlds.hpp).
- Rank-zero discovery may shell out to tar extracting RSPVEC and molden.inp into calc_dir/dalton_import. Prepared path registry prevents subworld extraction; per-leg fallback is read-only [dalton_import.hpp:locate_dalton_dir](./dalton_import.hpp) [dalton_import.hpp:register_prepared_dalton_seed](./dalton_import.hpp) [dalton_import.hpp:dalton_fd_guess](./dalton_import.hpp).
- Ground-state seed writing projects/orthonormalizes occupied Molden orbitals, changes defaults and writes restartdata plus optional HDF5; active geometry avoids Molden rounding [dalton_gs_seed.hpp:write_gs_seed_from_molden](./dalton_gs_seed.hpp).

## Source files

| File | Role |
| --- | --- |
| [convergence_policy.hpp](./convergence_policy.hpp) | Threshold scaling and distinct FD/ES stopping policies |
| [iterate.hpp](./iterate.hpp) | One-rung solver loop and checkpoint cadence |
| [iterate_protocol.hpp](./iterate_protocol.hpp) | Threshold ladder and rung finality |
| [fd_problem.hpp](./fd_problem.hpp) | Driven source/frequency descriptor |
| [fd_solver.hpp](./fd_solver.hpp) | FD fixed-point, KAIN, diagnostics and output reevaluation |
| [es_solver.hpp](./es_solver.hpp) | Bundle ES iteration, root tracking, locking and kinetic evidence |
| [fd_explicit_accuracy.hpp](./fd_explicit_accuracy.hpp) | Named-criterion input/saved-state evidence |
| [response_subspace_kain.hpp](./response_subspace_kain.hpp) | Response-state KAIN vector arithmetic |
| [response_state.hpp](./response_state.hpp) | Shell/storage state layouts and backend selection |
| [state_metrics.hpp](./state_metrics.hpp) | Coefficient/RSS/wall/iteration metrics |
| [function_residuals.hpp](./function_residuals.hpp) | Per-function residual reductions |
| [function_kinetics.hpp](./function_kinetics.hpp) | Gradient-form kinetic metrics |
| [function_metric_diagnostics.hpp](./function_metric_diagnostics.hpp) | Function summaries and provenance diagnostics |
| [es_kinetic_evidence.hpp](./es_kinetic_evidence.hpp) | Validation and unavailable reasons for root kinetic evidence |
| [es_root_identity.hpp](./es_root_identity.hpp) | Stable bookkeeping root indices and correspondence |
| [function_hdf5_io.hpp](./function_hdf5_io.hpp) | Structured Function and opaque archive HDF5 persistence |
| [response_metadata.hpp](./response_metadata.hpp) | JSON registry, property upserts and shard merging |
| [response_seed.hpp](./response_seed.hpp) | Seeds: starting guesses from elsewhere, in the node layout under seeds/ |
| [saved_states.hpp](./saved_states.hpp) | Lists and loads a run folder's saved response nodes through the node module, read-only |
| [node_subworlds.hpp](./node_subworlds.hpp) | Hostname-aware task partition and World lifetimes |
| [es_analysis.hpp](./es_analysis.hpp) | Transition properties and AO analysis |
| [raman_analysis.hpp](./raman_analysis.hpp) | Hessian, modes and Raman spectrum |
| [es_seed_guard.hpp](./es_seed_guard.hpp) | Informational imported-root overlap evidence |
| [dalton_gs_seed.hpp](./dalton_gs_seed.hpp) | Molden occupied orbitals to restartdata |
| [dalton_import.hpp](./dalton_import.hpp) | Dalton discovery, validation and seeding |
| [build_response_ground_state.hpp](./build_response_ground_state.hpp) | Ground kernel builders, guess adapters, warmup/downselection and TDA-to-Full promotion |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

Targets and CTest registrations are in the parent manifest; BUILD_TESTING is needed for registered tests [../CMakeLists.txt:BUILD_TESTING](../CMakeLists.txt).

```sh
cmake --build <build-dir> --target molresponse test_fd_gate test_es_plateau test_response_metadata test_fd_explicit_accuracy test_iteration_protocol test_es_kinetic_evidence
```

```sh
ctest --test-dir <build-dir> -R "^madness/test/molresponse/(test_fd_gate|test_es_plateau|test_response_metadata|test_fd_explicit_accuracy|test_iteration_protocol|test_es_kinetic_evidence)/run$"
```

Function residual, kinetic and ES archive tests have MPI variants when MPI/launcher enabled; HDF5 tests require MADNESS_ENABLE_HDF5 [../CMakeLists.txt:ENABLE_MPI](../CMakeLists.txt) [../CMakeLists.txt:MADNESS_ENABLE_HDF5](../CMakeLists.txt).
