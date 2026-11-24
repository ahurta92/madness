# calc component

## Overview

Pure response DAG scheduling from the session's plans, and World-bound execution through the response session, parallel waves and property assembly.

## Interface

- calc_manager.hpp supplies CalcNode/WorkItem/StatePlan and pure build_dag, planned_action, schedule and accuracy/drop audits; ICalcExecutor is the numerical execution seam [calc_manager.hpp:CalcNode](./calc_manager.hpp) [calc_manager.hpp:ICalcExecutor](./calc_manager.hpp).
- ExecutorSettings carries the response session and caller policy/guess/parallel settings; ExecutorContext adds borrowed World/GroundState; FdResponseExecutor dispatches Static/Full FD in both shells, closed-shell TDA/Full ES and VBC through the session [calc_executor.hpp:ExecutorSettings](./calc_executor.hpp) [calc_executor.hpp:ExecutorContext](./calc_executor.hpp) [calc_executor.hpp:FdResponseExecutor](./calc_executor.hpp).
- CalcManager builds and repeatedly runs the next wave; assemble_alpha/beta/tpa contract saved states and publish property rows [calc_executor.hpp:CalcManager](./calc_executor.hpp) [calc_executor.hpp:assemble_alpha](./calc_executor.hpp) [calc_executor.hpp:assemble_beta](./calc_executor.hpp) [calc_executor.hpp:assemble_tpa](./calc_executor.hpp).

## Invariants

- Every schedulable node needs a plan from the session; schedule() refuses one without. Only schedule(...).front() has final actions: later waves use the same snapshot and must be rescheduled after completed work. Ordering is deterministic because every rank must issue matching collectives [calc_manager.hpp:schedule](./calc_manager.hpp) [calc_executor.hpp:CalcManager::run](./calc_executor.hpp).
- FD identity includes perturbation plus frequency key f%.5f; stored frequency is checked at tighter tolerance to reject aliased exact states. ES identity checks method/root count where recorded [calc_manager.hpp:fd_entry](./calc_manager.hpp) [calc_manager.hpp:es_entry_matches](./calc_manager.hpp).
- Hard prerequisites require same-rung convergence or an explicitly met planner load; soft seed hints do not gate correctness. All-atom nuclear sentinels expand by atom count; derived FD wildcard expands only after matching ES bundle converges at its finest rung [calc_manager.hpp:build_dag](./calc_manager.hpp) [calc_manager.hpp:prerequisites_converged](./calc_manager.hpp) [calc_executor.hpp:expand_converged_es](./calc_executor.hpp).
- A rung the run has scheduled once is skipped afterwards; the session's node decides when a state stops and saves. Repeated wave plus an unchanged snapshot quarantines those nodes and continues independent work; wall/RSS measurements are excluded from progress fingerprints [calc_manager.hpp:planned_action](./calc_manager.hpp) [calc_executor.hpp:strip_volatile_metrics](./calc_executor.hpp) [calc_executor.hpp:CalcManager::run](./calc_executor.hpp).
- An FD state starts from its own node, from the nearest saved frequency of the same perturbation, from an imported Dalton guess, or from deep-copied zeros [calc_executor.hpp:solve_fd_protocol](./calc_executor.hpp).
- An ES bundle starts from its node, an imported Dalton bundle, or the lowest Ritz roots of an oversampled TDA guess with no warmup updates; Full promotes X to independent Y=0 [calc_executor.hpp:solve_shared_excited](./calc_executor.hpp) [calc_executor.hpp:es_iteration_policies](./calc_executor.hpp).
- Open-shell dipole FD is solved/saved but alpha/beta/TPA assembly and ES execution are skipped; nuclear beta-spin source throws. Magnetic source is unimplemented [calc_executor.hpp:pert_source_beta](./calc_executor.hpp) [calc_executor.hpp:pert_source_alpha](./calc_executor.hpp) [calc_executor.hpp:FdResponseExecutor](./calc_executor.hpp) [calc_executor.hpp:assemble_alpha](./calc_executor.hpp).
- FD/NuclearFD/VBC may run in subworlds; ES remains on the universe. Original universe pmap object must be restored before destroying subworlds. Per-group metadata shards merge after the universe fence, for the output file only [calc_executor.hpp:CalcManager::run](./calc_executor.hpp).
- Each pass, rank 0 lists the completely saved nodes under states/response/, judges each against the current request (the restart's reuse judgement, from stored metrics), builds the snapshot the scheduler functions read and broadcasts it; response_metadata.json is written for people and scripts and never read to decide [node_snapshot.hpp:node_snapshot](./node_snapshot.hpp) [calc_executor.hpp:session_snapshot](./calc_executor.hpp) [shared_response_session.hpp:Session::node_meets_request](../solvers/shared_response_session.hpp).
- Property assembly loads only states whose nodes meet the final request; other components are published as unavailable [calc_executor.hpp:assemble_alpha](./calc_executor.hpp) [calc_executor.hpp:assemble_beta](./calc_executor.hpp) [calc_executor.hpp:assemble_tpa](./calc_executor.hpp).

## Side Effects

- Executors change global FunctionDefaults, prepare/reproject ground and response states, run collectives, log protocol start/end and save states through their session nodes [calc_executor.hpp:solve_fd_protocol](./calc_executor.hpp) [calc_executor.hpp:solve_shared_excited](./calc_executor.hpp).
- The excited-state analysis is computed whenever a complete bundle is assembled, including a bundle reused from its node with no updates [calc_executor.hpp:assemble_shared_excited](./calc_executor.hpp).
- Run summary persists dropped VBC/two-photon work and exposes scheduler trace, pass count and stalled/dropped stop reasons; empty schedule alone is not a successful numerical verdict [calc_manager.hpp:audit_dropped_work](./calc_manager.hpp) [calc_executor.hpp:CalcManager::run](./calc_executor.hpp).
- Beta/Raman and TPA property rows are rank-zero metadata writes. TPA uses the residue (P,Q) composition only [calc_executor.hpp:assemble_beta](./calc_executor.hpp) [calc_executor.hpp:assemble_tpa](./calc_executor.hpp).

## Source files

| File | Role |
| --- | --- |
| [calc_manager.hpp](./calc_manager.hpp) | Pure DAG/state-plan scheduler and audits |
| [calc_executor.hpp](./calc_executor.hpp) | Numerical execution, parallel scheduling and properties |
| [node_snapshot.hpp](./node_snapshot.hpp) | The scheduler's snapshot of the saved states, from node descriptions |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

Pure scheduler test is unconditional as a target; CTest registration requires BUILD_TESTING [../CMakeLists.txt:test_calc_manager](../CMakeLists.txt).

```sh
cmake --build <build-dir> --target molresponse test_calc_manager
```

```sh
ctest --test-dir <build-dir> -R "^madness/test/molresponse/test_calc_manager/run$"
```

End-to-end coverage is the madqc response qctests (`ctest -L qctest -R "madness/test/qc/response_"`).
