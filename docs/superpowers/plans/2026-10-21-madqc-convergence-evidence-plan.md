# MADQC Convergence Evidence Study Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Produce a reproducible H₂/LiH/water evidence set that identifies whether response convergence is limited by the ground-state reference, response solver, or stopping gates before shared convergence behavior is changed.

**Architecture:** Keep source worktrees and numerical references unchanged. Use isolated analysis runs, portable JSON manifests, native MADNESS HDF5 checkpoints, structured `dump_mra_trees --coeffs` exports, and a matplotlib notebook as the analysis interface. Promote conclusions into a reviewed metric and gate contract before any solver implementation work.

**Tech Stack:** MADNESS `madqc`, `dump_mra_trees`, JSON, native HDF5, `pymra`, `h5py`, matplotlib, Jupyter/nbclient, Python 3.

**Spec:** `docs/superpowers/specs/2026-10-21-madqc-convergence-evidence-design.md`

## Global Constraints

- Do not modify the Fable or Astra source worktrees, application branches, numerical references, or source commits during the evidence phase.
- Forced runs are diagnostic trajectories and are not calibrated convergence references.
- Missing or incomparable measurements must be recorded as unavailable, never as zero.
- Native application HDF5 archives and structured `pymra` exports must be identified separately.
- All response results must retain molecule, ground-state endpoint, protocol, frequency/state, phase, and run classification.
- Do not change solver gates or defaults until the evidence phase and metric/gate contract are reviewed.

## Review Focus

- Ground-state endpoints must be matched to response runs; test the endpoint identifier and protocol metadata in every response record.
- Total-density and spin-density metrics must not be conflated; test an open-shell record with both values present.
- Kinetic trace must be `sum_i T_ii`, while per-function contributions remain available; test scalar and contribution consistency.
- Native HDF5 archives are not automatically `pymra` structured fields; test the unavailable/exported paths separately.
- Forced progress must not be labeled as calibrated convergence; test ordinary, forced, and configured-progress classifications in the manifest.

---

### Task 1: Create the experiment manifest and portable evidence schema

**Files:**
- Create: `analysis/madqc_convergence_manifest.json`
- Create: `analysis/madqc_evidence_schema.py`
- Test: `analysis/test_madqc_evidence_schema.py`

**Interfaces:**
- Produces `RunRecord` dictionaries with `run_id`, `molecule`, `workflow`, `source_branch`, `ground_state_endpoint`, `protocol`, `phase`, `sampling_phase`, `classification`, and `artifacts`.
- Produces `MetricRecord` dictionaries with `value`, `available`, `aggregation`, `units`, and optional `per_function` values.

- [ ] **Step 1: Write failing schema tests**

  Test that a response record without a ground-state endpoint, a record with an unknown classification, or a metric represented as numeric zero when marked unavailable is rejected. Test that an open-shell record can carry both total-density and spin-density metrics.

- [ ] **Step 2: Run the schema tests and verify failure**

  Run: `/home/ahurtado/miniforge3/envs/qca070/bin/python -m pytest analysis/test_madqc_evidence_schema.py -q`

  Expected: FAIL because the schema validator is not implemented.

- [ ] **Step 3: Implement the validator**

  Implement small validation functions in `analysis/madqc_evidence_schema.py`; keep the schema JSON-oriented and dependency-light. Encode unavailable measurements explicitly as `available: false` and `value: null`.

- [ ] **Step 4: Run the schema tests and verify success**

  Run the same pytest command. Expected: PASS.

- [ ] **Step 5: Populate the initial manifest**

  Add the existing H₂ HDF5 run, coarse/fine orbital exports, and current portable Astra trajectories to `analysis/madqc_convergence_manifest.json`. Include exact artifact paths and source/build provenance.

### Task 2: Run the controlled ground-state and response matrix

**Files:**
- Create: `analysis/run_madqc_evidence.py`
- Create: session-folder run manifest and logs under `parallel-implementations/2026-10-20/evidence-runs/`
- Modify: `analysis/madqc_convergence_manifest.json`

**Interfaces:**
- Consumes: existing QC decks under `src/examples/qc/`, isolated HDF5-enabled MADNESS build, Task 1 record schema.
- Produces: ordinary and forced run records with ground-state endpoint IDs, protocol metadata, iteration histories, response histories, and native HDF5 artifact paths.

- [ ] **Step 1: Define the matrix in the runner**

  Define H₂, LiH, and water cases; ordinary protocol execution; fixed-iteration ground-state runs at each selected rung; and matched static plus frequency-dependent or excited-state response runs where the deck supports them. Make the iteration budget and protocol list explicit command-line/configuration values.

- [ ] **Step 2: Add dry-run validation**

  Add a `--dry-run` mode that prints every planned run, source deck, protocol, iteration budget, output directory, and classification. Verify that every response run names a ground-state endpoint.

- [ ] **Step 3: Execute the matrix in isolated output directories**

  Use `MAD_NUM_THREADS` and BLAS thread limits appropriate to the host. Keep outputs outside the source worktrees. Preserve stdout/stderr, input decks, calc-info JSON, native HDF5 archives, and exit status.

- [ ] **Step 4: Validate the generated records**

  Run Task 1 validation against every generated record. Reject incomplete response-to-ground-state linkage or missing phase/classification metadata.

- [ ] **Step 5: Record unavailable cases explicitly**

  If water or a response type cannot run within the available build/dependency configuration, add an unavailable record with the reason and preserve the attempted command/log.

### Task 3: Export and compare HDF5 fields

**Files:**
- Create: `analysis/compare_madqc_hdf5.py`
- Create: session-folder `hdf5-analysis/<molecule>/<endpoint>/`
- Modify: `analysis/madqc_convergence_manifest.json`

**Interfaces:**
- Consumes: native restart archives, `dump_mra_trees --coeffs`, `pymra.read_function`, endpoint metadata.
- Produces: machine-readable field comparison JSON containing grid, cell, representation, endpoint IDs, L2 differences, RMS, max absolute difference, and limitations.

- [ ] **Step 1: Add export discovery tests**

  Test that a native `madness_archive` HDF5 file is classified as an application archive and that a structured `.mad.h5` file is classified as a `pymra` field export.

- [ ] **Step 2: Implement matched-endpoint checks**

  Require compatible cell, state/spin identity, representation, and endpoint metadata before comparison. Return an unavailable result with a reason when compatibility is not established.

- [ ] **Step 3: Implement common-grid comparison**

  Use a documented midpoint grid and calculate sampled L2, RMS, maximum absolute difference, and optional spatial-region contributions. Label orbital-derived closed-shell density proxies distinctly from response densities.

- [ ] **Step 4: Run comparisons for exhausted endpoints**

  Export matched H₂, LiH, and water endpoints where possible. Save native and structured artifacts under the session folder and update the manifest.

- [ ] **Step 5: Verify reproducibility**

  Re-run the comparison script from the manifest and confirm identical numeric results within floating-point serialization tolerance.

### Task 4: Expand the analysis notebook and interpretation note

**Files:**
- Modify: `analysis/madqc_trajectory_review.ipynb`
- Create: `analysis/madqc_convergence_report.md`
- Modify: session `trajectory-analysis-2026-10-21.md`

**Interfaces:**
- Consumes: Task 2 portable trajectories, Task 3 field comparison JSON, manifest metadata.
- Produces: plots and an evidence-backed classification for representative closed-shell and open-shell cases.

- [ ] **Step 1: Add endpoint-indexed response plots**

  Plot response iteration count, final response residuals, and final properties against the ground-state endpoint ID and protocol rung.

- [ ] **Step 2: Add full ordinary/forced/configured-progress histories**

  Preserve protocol shading, switch markers, phase labels, sampling labels, and unavailable measurements. Do not collapse forced runs into endpoint-only bars.

- [ ] **Step 3: Add HDF5 field plots**

  Plot available density differences and report the exact grid and field construction. Render an unavailable panel/message when `pymra` or archives are absent.

- [ ] **Step 4: Execute the notebook in `qca070`**

  Run with `nbclient` and verify all code cells complete. Save the executed notebook only if outputs remain reasonably sized and reproducible.

- [ ] **Step 5: Write the interpretation report**

  For each representative case, classify evidence as reference-limited, response-solver-limited, gate-limited, or inconclusive. Include protocol values, endpoint metadata, and limitations.

### Task 5: Review and freeze the metric/gate contract

**Files:**
- Create: `docs/superpowers/specs/2026-<date>-madqc-metric-gate-contract.md`
- Modify: session decision notes

**Interfaces:**
- Consumes: Task 4 report and all machine-readable evidence.
- Produces: reviewed contract for the subsequent implementation phase.

- [ ] **Step 1: Summarize observed behavior**

  Record whether total-density, spin-density, BSH, energy, kinetic trace, response residuals, and properties agree or diverge across endpoint quality.

- [ ] **Step 2: Decide aggregation and progress rules**

  Resolve maximum aggregation, quadrature-gate treatment, explicit user criteria, protocol continuation, stall reporting, and unavailable measurements using observed trajectories.

- [ ] **Step 3: Define implementation acceptance tests**

  Specify exact JSON keys, scalar/per-function relationships, protocol transition behavior, and tests required before solver behavior changes.

- [ ] **Step 4: Obtain review before implementation**

  Do not modify moldft/molresponse convergence behavior until the contract is reviewed and approved.

## Verification Commands

- Schema tests: `/home/ahurtado/miniforge3/envs/qca070/bin/python -m pytest analysis/test_madqc_evidence_schema.py -q`
- Notebook execution: `/home/ahurtado/miniforge3/envs/qca070/bin/python -m nbclient analysis/madqc_trajectory_review.ipynb --execute --output analysis/madqc_trajectory_review.executed.ipynb`
- JSON validation: `/home/ahurtado/miniforge3/envs/qca070/bin/python -m json.tool analysis/madqc_convergence_manifest.json`
- Source isolation: `git -C /home/ahurtado/repos/madness-step4-astra status --short` and `git -C /home/ahurtado/repos/madness-step4-fable status --short`; confirm no analysis command edits source files.

## Commit Strategy

Commit each completed task separately with plain messages. Keep analysis artifacts and source changes in separate commits. Do not merge or change application source until Task 5 is approved.
