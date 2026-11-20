# calc component

## Overview

Pure response DAG reconciliation and World-bound execution, restart, parallel waves and property assembly.

## Interface

- calc_manager.hpp supplies CalcNode/WorkItem/StatePlan and pure build_dag, reconcile_protocol, planned_action, schedule and accuracy/drop audits; ICalcExecutor is the numerical execution seam [calc_manager.hpp:CalcNode](./calc_manager.hpp) [calc_manager.hpp:ICalcExecutor](./calc_manager.hpp).
- ExecutorSettings carries caller policy/seeding/parallel settings; ExecutorContext adds borrowed World/GroundState; FdResponseExecutor dispatches Static/Full FD in both shells and closed-shell TDA/Full ES [calc_executor.hpp:ExecutorSettings](./calc_executor.hpp) [calc_executor.hpp:ExecutorContext](./calc_executor.hpp) [calc_executor.hpp:FdResponseExecutor](./calc_executor.hpp).
- CalcManager builds and repeatedly runs the next wave; assemble_alpha/beta/tpa contract saved states and publish property rows [calc_executor.hpp:CalcManager](./calc_executor.hpp) [calc_executor.hpp:assemble_alpha](./calc_executor.hpp) [calc_executor.hpp:assemble_beta](./calc_executor.hpp) [calc_executor.hpp:assemble_tpa](./calc_executor.hpp).

## Invariants

- Only schedule(...).front() has final reconcile actions: later waves use the same metadata snapshot and must be rescheduled after completed work. Ordering is deterministic because every rank must issue matching collectives [calc_manager.hpp:schedule](./calc_manager.hpp) [calc_executor.hpp:CalcManager::run](./calc_executor.hpp).
- FD identity includes perturbation plus frequency key f%.5f; stored frequency is checked at tighter tolerance to reject aliased exact states. ES identity checks method/root count where recorded [calc_manager.hpp:fd_entry](./calc_manager.hpp) [calc_manager.hpp:es_entry_matches](./calc_manager.hpp).
- Hard prerequisites require same-rung convergence or an explicitly met planner load; soft seed hints do not gate correctness. All-atom nuclear sentinels expand by atom count; derived FD wildcard expands only after matching ES bundle converges at its finest rung [calc_manager.hpp:build_dag](./calc_manager.hpp) [calc_manager.hpp:prerequisites_converged](./calc_manager.hpp) [calc_executor.hpp:expand_converged_es](./calc_executor.hpp).
- Budget-exhausted or stalled FD/ES rungs can skip forward while remaining unconverged; a larger budget enables Resume. Repeated wave plus unchanged state metadata quarantines those nodes and continues independent work; wall/RSS measurements are excluded from progress fingerprints [calc_manager.hpp:reconcile_protocol](./calc_manager.hpp) [calc_executor.hpp:strip_volatile_metrics](./calc_executor.hpp) [calc_executor.hpp:CalcManager::run](./calc_executor.hpp).
- FD fresh guesses are deep-copied zeros; restart, planner frequency seed, optional ES-root and Dalton guesses are explicit separate routes. Written requested criteria never receive best-effort accept-at-maxiter [calc_executor.hpp:solve_fd_protocol](./calc_executor.hpp).
- Current ES policy remains engine-decided and production executor defaults root locking on. Warmup is KAIN-free, oversampled and downselected; Full promotes X to independent Y=0. Own warmup clears default KAIN hold-off unless explicitly written [calc_executor.hpp:ExecutorSettings](./calc_executor.hpp) [calc_executor.hpp:es_iteration_policies](./calc_executor.hpp) [calc_executor.hpp:solve_es_full_closed_shell](./calc_executor.hpp).
- ES warmup performs additional steps outside the main stage allowance. Downselection and TDA-to-Full promotion reset the main iteration count; checkpoint `iterations_this_rung` counts main-solve iterations. Bundle `wall_s` starts before ground preparation and warmup and ends just before saving, so it includes preparation, warmup and main work but excludes that save. Dividing it by the main iteration count does not measure main-iteration cost alone [../solvers/build_response_ground_state.hpp:slice_state_lowest](../solvers/build_response_ground_state.hpp) [calc_executor.hpp:solve_es_tda_closed_shell](./calc_executor.hpp) [calc_executor.hpp:solve_es_full_closed_shell](./calc_executor.hpp).
- Open-shell dipole FD is solved/saved but alpha/beta/TPA assembly and ES execution are skipped; nuclear beta-spin source throws. Magnetic source is unimplemented [calc_executor.hpp:pert_source_beta](./calc_executor.hpp) [calc_executor.hpp:pert_source_alpha](./calc_executor.hpp) [calc_executor.hpp:FdResponseExecutor](./calc_executor.hpp) [calc_executor.hpp:assemble_alpha](./calc_executor.hpp).
- FD/NuclearFD/VBC may run in subworlds; ES remains on the universe. Original universe pmap object must be restored before destroying subworlds. Per-group metadata shards merge after the universe fence [calc_executor.hpp:CalcManager::run](./calc_executor.hpp).
- Property accuracy records the actual loaded source rung and distinguishes strict convergence, best-effort acceptance and planner-met loads. Alpha/beta accept planner-met rows; TPA currently checks strict FD convergence without passing state_plans [calc_manager.hpp:fd_leg_accuracy](./calc_manager.hpp) [calc_executor.hpp:assemble_alpha](./calc_executor.hpp) [calc_executor.hpp:assemble_beta](./calc_executor.hpp) [calc_executor.hpp:assemble_tpa](./calc_executor.hpp).

## Side Effects

- Executors change global FunctionDefaults, prepare/reproject ground and response states, run collectives, log protocol start/end and persist states/metrics [calc_executor.hpp:solve_fd_protocol](./calc_executor.hpp) [calc_executor.hpp:solve_es_tda_closed_shell](./calc_executor.hpp).
- Named checkpoint observers run after saved protocol boundaries; explicit FD accuracy requires saved-state validation. Loaded converged ES bundles missing transition analysis can generate that analysis without iterations [calc_executor.hpp:solve_fd_protocol](./calc_executor.hpp) [calc_executor.hpp:analyze_loaded_es_bundle](./calc_executor.hpp).
- Run summary persists dropped VBC/two-photon work and exposes scheduler trace, pass count and stalled/dropped stop reasons; empty schedule alone is not a successful numerical verdict [calc_manager.hpp:audit_dropped_work](./calc_manager.hpp) [calc_executor.hpp:CalcManager::run](./calc_executor.hpp).
- Beta/Raman and TPA property rows are rank-zero metadata writes. Production TPA uses residue P,Q composition; diagnostic-only skips contraction and property writes, root filters skip property writes, while legacy TPA source archives are still produced by that legacy branch [calc_executor.hpp:assemble_beta](./calc_executor.hpp) [calc_executor.hpp:assemble_tpa](./calc_executor.hpp).

## Source files

| File | Role |
| --- | --- |
| [calc_manager.hpp](./calc_manager.hpp) | Pure DAG/state-plan/reconcile scheduler and audits |
| [calc_executor.hpp](./calc_executor.hpp) | Numerical execution, parallel scheduling and properties |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

Pure scheduler test is unconditional as a target; CTest registration requires BUILD_TESTING [../CMakeLists.txt:test_calc_manager](../CMakeLists.txt).

```sh
cmake --build <build-dir> --target molresponse test_calc_manager
```

```sh
ctest --test-dir <build-dir> -R "^madness/test/molresponse/test_calc_manager/run$"
```

Integration registration requires MADNESS_RESPONSE_FIXTURE_DIR with h2o/mad.restartdata [../CMakeLists.txt:MADNESS_RESPONSE_FIXTURE_DIR](../CMakeLists.txt).

```sh
cmake --build <build-dir> --target test_calc_manager_run test_run_response
```

```sh
ctest --test-dir <build-dir> -R "^madness/test/molresponse/(test_calc_manager_run|test_run_response)/run$"
```
