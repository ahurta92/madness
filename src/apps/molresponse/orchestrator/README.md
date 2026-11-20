# orchestrator component

## Overview

Public workflow seam loading a ground-state archive, building and driving the calculation DAG, assembling response properties, and returning JSON on rank 0.

## Interface

- ResponseWorkflowInput carries archive_file, a nonempty protocol ladder, a prebuilt ResponsePlan, ExecutorSettings and optional Molecule; ResponseWorkflowOutput contains property, metadata, vibrational, timing, diagnostics and reserved export/debug slots [response_workflow.hpp:ResponseWorkflowInput](./response_workflow.hpp) [response_workflow.hpp:ResponseWorkflowOutput](./response_workflow.hpp).
- run_response loads and prepares the archive ground state, then run_response_with_ground drives CalcManager and property assembly [response_workflow.hpp:run_response](./response_workflow.hpp) [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).
- molecule_from_archive_dir reads moldft.calc_info.json or mad.calc_info.json; fock_json_from_archive_dir selects the corresponding Fock JSON if present [response_workflow.hpp:molecule_from_archive_dir](./response_workflow.hpp) [response_workflow.hpp:fock_json_from_archive_dir](./response_workflow.hpp).

## Invariants

- Protocol ladders must be nonempty; run_response_with_ground requires the supplied ground state prepared at the first protocol [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).
- DFT ground states with VBC or derived-FD requests are rejected because the second XC kernel required by beta/Raman/2PA is absent [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).
- When archive_file is nonempty, the ground-state archive fingerprint gate runs before restart reuse; a mismatched stamp fails collectively unless MADRESPONSE_ALLOW_GS_MISMATCH starts with 1. A direct run_response_with_ground caller with an empty archive path bypasses this gate and must establish compatibility with saved response states before reusing calc_dir. Named madqc supplies an archive path [response_workflow.hpp:gs_fingerprint_gate](./response_workflow.hpp) [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).
- Requested FD subworld fan-out needs an archive path: without one the workflow prints a warning and keeps the single-World path [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).
- Collected output property/metadata/timing/diagnostic objects are authoritative on rank 0; other ranks retain their default output objects [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).

## Side Effects

- Creates calc_dir; on the archive-backed path, writes the ground-state fingerprint to response_metadata.json on rank 0. Missing older stamps warn and are restamped [response_workflow.hpp:gs_fingerprint_gate](./response_workflow.hpp) [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).
- Merges stale per-state metadata shards before restart reconciliation; fingerprint and shard-sweep failures become collective failures [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).
- Changes response protocol defaults, loads archives, runs solver/checkpoint writes, optionally loads separate ground states in FD subworlds, assembles alpha/beta/Raman/2PA as requested, and records effective native/HDF5 backend [response_workflow.hpp:run_response](./response_workflow.hpp) [response_workflow.hpp:run_response_with_ground](./response_workflow.hpp).

## Source files

| File | Role |
| --- | --- |
| [response_workflow.hpp](./response_workflow.hpp) | Workflow input/output seam and collective ground-state restart gate |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

```sh
cmake --build <build-dir> --target molresponse
```

Integration test registration requires MADNESS_RESPONSE_FIXTURE_DIR with h2o/mad.restartdata [../CMakeLists.txt:MADNESS_RESPONSE_FIXTURE_DIR](../CMakeLists.txt).

```sh
cmake --build <build-dir> --target test_run_response
```

```sh
ctest --test-dir <build-dir> -R "^madness/test/molresponse/test_run_response/run$"
```
