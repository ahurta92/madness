# orchestrator component

## Overview

Workflow seam that, given a prepared ground state and a response session, builds and drives the calculation DAG, assembles response properties, and returns JSON on rank 0.

## Interface

- ResponseWorkflowInput carries archive_file, a nonempty protocol ladder, a prebuilt ResponsePlan and ExecutorSettings (which hold the session); ResponseWorkflowOutput contains property, metadata, vibrational, timing, diagnostics and reserved export/debug slots [response_workflow.hpp:ResponseWorkflowInput](./response_workflow.hpp) [response_workflow.hpp:ResponseWorkflowOutput](./response_workflow.hpp).
- run_response_with_ground drives CalcManager and property assembly; it refuses settings without a session. The madqc adapter is its caller [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).

## Invariants

- Protocol ladders must be nonempty; run_response_with_ground requires the supplied ground state prepared at the first protocol [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).
- DFT ground states with VBC or derived-FD requests are rejected because the second XC kernel required by beta/Raman/2PA is absent [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).
- Saved response states are matched to the ground state by the session's nodes; seeds (seeds/response/) by their settings, the ground state's geometry and model [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).
- Requested FD subworld fan-out needs an archive path: without one the workflow prints a warning and keeps the single-World path [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).
- Collected output property/metadata/timing/diagnostic objects are authoritative on rank 0; other ranks retain their default output objects [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).

## Side Effects

- Merges stale per-state metadata shards before restart reconciliation; shard-sweep failures become collective failures [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).
- Changes response protocol defaults, runs solver/checkpoint writes, optionally loads separate ground states in FD subworlds, assembles alpha/beta/Raman/2PA as requested, and records effective native/HDF5 backend [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).

## Source files

| File | Role |
| --- | --- |
| [response_workflow.hpp](./response_workflow.hpp) | Workflow input/output seam and collective ground-state restart gate |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

```sh
cmake --build <build-dir> --target molresponse
```

End-to-end coverage is the madqc response qctests (`ctest -L qctest -R "madness/test/qc/response_"`).
