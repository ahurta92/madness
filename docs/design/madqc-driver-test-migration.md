# Driver and test migration — 2026-11-19

`madqc` is the public executable for supported engines. This index records the
purpose of each retired test independently of record format. Numerical reference
files and comparison tolerances are unchanged. A listed unit test covers its
own component; it does **not** certify an application adapter that is not wired.

| Retired test | Purpose / disposition | Retained checks | Remaining acceptance gap |
|---|---|---|---|
| `test_named_stages.cc` | Shared policy, identity, serialization and malformed-input checks; named grammar retired. | test_driver_contract; test_iteration_metrics; test_iteration_restart; test_checkpoint_store; test_scf_settings | Response adapter integration remains pending. |
| `test_named_replay.cc` | Offline convergence replay against v1/v2 named records retired; backward format compatibility is not a release requirement. | test_iteration_restart; test_checkpoint_store | — |
| `test_response_properties.cc` | Named property selector/formatter retired. Engine property planning and convergence gating remain. | test_response_property_planner; test_calc_manager; test_raman_spectrum; test_workflow_response | Grouped shared response output needs its own adapter tests. |
| `test_named_checkpoint_e2e.py` | Publication, restart, corruption, incomplete prerequisites and geometry changes. | test_shared_scf; test_scf_recovery; test_workflow_completion; test_checkpoint_store | — |
| `test_named_checkpoint_cadence.py` | Endpoint/periodic saves and bounded retention. keep_all and arbitrary older-selection modes retired. | test_shared_scf; test_scf_recovery; test_checkpoint_retention; test_iteration_metrics | — |
| `test_named_checkpoint_hdf5.py` | Named HDF5 generation mirror retired; shared restart uses native archives. | test_scf_checkpoint; optional molresponse test_state_archive_hdf5 | Shared HDF5 is deferred. |
| `test_named_convergence.py` | Measured residual/density/energy evidence and final representation. | test_scf_final_evidence; test_iteration_metrics; test_scf_evaluation; test_function_residuals | Response shared publication remains pending. |
| `test_named_density_changes.py` | Charge/spin diagnostics and density differences. | test_density_changes; test_metric_helpers; test_scf_evaluation; test_function_metric_summary | — |
| `test_named_econv.py` | Old econv engine-versus-request modes retired; explicit energy tolerance is the contract. | test_scf_inputs; test_iteration_metrics; test_shared_scf; test_scf_final_evidence | — |
| `test_named_es_econv.py` | Legacy report-only econv for excited roots retired. | test_iteration_metrics; test_es_kinetic_evidence | Shared excited-state energy-change criteria must be wired and tested during response adoption. |
| `test_named_es_restart.py` | Aliased protocol keys and missing transition analysis on loaded roots. | test_calc_manager; test_response_metadata; test_protocol_key; response_h2_es_tda/rpa | End-to-end shared bundle expansion/analysis recovery is still a response release gate. |
| `test_named_from_store.py` | External immutable named fixture store retired, including ancestor-only compatibility checks. | test_scf_recovery; test_shared_scf; test_workflow_response | No shared response external-fixture replay claim. |
| `test_named_moldft_seed.py` | Cold/seeded energy agreement, iteration benefit, restart precedence and geometry rejection. | test_workflow_seed; test_workflow_seed_accuracy; test_dalton_import | Copying every seed source into each checkpoint is not retained; seed provenance adoption is separate. |
| `test_named_dalton_leg.py` | DALTON FD seed import and seed-versus-solved distinction. | test_dalton_import; response_h2o_tpa_seeded; test_calc_manager | Shared per-state publication/reuse of seeded FD functions remains pending. |
| `test_named_reference_from_generation.py` | Exact qualified ground-state handoff and changed-reference rejection; arbitrary named standalone dependency decks retired. | test_shared_reference; test_workflow_completion; test_workflow_requalification | — |
| `test_named_response_checkpoint.py` | Scientific response solve/restart remains in the engine; old named publication wrapper retired. | test_calc_manager; test_response_metadata; test_fd_kfallback_load; test_workflow_response; response_he_alpha; response_h2_es_tda/rpa | Shared FD retention/cadence, metric-based tightening, per-state additions and coupled-bundle restart still require response adapter acceptance tests. |
| `test_named_response_eligibility.py` | Missing density comparisons must not be false zero. | test_iteration_metrics; test_fd_explicit_accuracy; test_function_metric_summary | Shared response controller integration remains pending. |
| `test_named_stages_e2e.py` | Fresh numerical chemistry retained; arbitrary stage naming/isolation language retired. | qctest SCF, response alpha/beta, Raman, TDA/RPA, 2PA, Nemo/CIS/OEP; test_driver_contract | OEP model coverage remains a release gap, not implied by its smoke deck. |
| `test_named_two_photon.py` | 2PA tensor identities, half-frequency legs and saved-function reuse. | test_workflow_response; test_tpa_regroup_identities; test_tpa_pq_vs_vbc; response_h2o_tpa(_seeded) | Shared mixed-request per-state reuse remains pending. Named-only output field names and TDA parser restriction are retired; engine capabilities govern workflows. |
| `test_convergence_modes.py` | Engine/request dual-verdict modes retired in favor of one shared controller. | test_iteration_metrics; test_scf_inputs; test_shared_scf | — |
| `test_response_convergence_modes.py` | Named engine/request dual verdict retired. | test_fd_explicit_accuracy; test_iteration_metrics | Response shared controller integration remains pending. |
| `test_configured_progress.py` | Experimental stall policy retired; intermediate thresholds follow the reviewed heuristic. | test_iteration_metrics; test_scf_run | — |
| `test_response_configured_progress.py` | Named experimental stall policy retired. | test_iteration_protocol; test_calc_manager | Shared response protocol decisions remain pending. |
| `test_explicit_accuracy.py` | SCF criterion selection, saved metrics, tighter requests and corruption. | test_shared_scf; test_scf_recovery; test_scf_inputs; test_iteration_restart | — |
| `test_response_explicit_accuracy.py` | Scientific FD criterion behavior stays in engine tests; old named wrapper retired. | test_fd_explicit_accuracy; test_function_residuals; test_iteration_metrics | Shared FD output eligibility remains pending. |
| `test_fixed_intermediate_work.py` | Fixed-only intermediate policy superseded by threshold-dependent tolerances; final representation still mandatory. | test_iteration_metrics; test_shared_scf; test_scf_run | Response shared controller integration remains pending. |
| `test_forced_rungs.py` | Diagnostic forced-iteration knob retired; maxiter is a fresh update allowance. | test_iteration_metrics; test_shared_scf; test_scf_recovery | — |
| `test_moldft_restart_e2e.py` | Legacy calc_info/raw-archive fallback replaced by strict selected native checkpoint recovery. | test_shared_scf; test_scf_recovery; test_workflow_completion | — |

## Input mapping

Only SCF blocks for `scf` and `response` workflows are migrated here. Each now
selects `bsh_residual`, `density_change` and `energy_change`. The former `dconv`
value becomes an explicit absolute density-change tolerance. The BSH tolerance
keeps the old effective `5 * max(final_threshold, dconv)` gate; `econv`
becomes the energy-change tolerance (1e-5 when omitted). This is deliberately
stricter through the added density-change check: the old loop
applied an additional update after checking that residual. New checks describe
the functions actually saved. Numerical references and comparison tolerances
remain fixed; discrepancies remain visible, never hidden by updating an oracle.

Nemo/CIS/OEP/CC2 still use their existing engine inputs until their separate
shared-layer adoption. Response-specific `dconv`, guess and iteration controls
also remain engine inputs in this increment. Removing their executables does
not claim that all engines already implement the shared record contract.

## Executables

Removed duplicate targets/mains: `moldft`, `nemo`, `molresponse`, `cis`, `oep`,
`cc2`. Their algorithms remain in MADchem and the response engine. Tests run
`madqc --wf=...`. Separate MP2 uses the `MP2` implementation, whereas the madqc
MP2 route uses CC2; retain `mp2`. Also retain `pno`, `znemo`, `zcis`, `DFdriver`,
periodic DFT, TDSE, import/analysis utilities and private test binaries.

Legacy CIS/OEP `--test` helpers remain callable in the library (`TDHF::test`,
`OEP::selftest`); their removed standalone CLI branches were not registered tests.
Workflow numerical tests remain the release checks. OEP coverage beyond OAEP
is still required before release.

## Known science gates

The existing OAEP smoke deck caps its Nemo reference at three iterations and
previously consumed the unfinished result. The workflow completion gate now
stops before OEP. Keep this numerical test visible as a release gap; do not
weaken the gate or regenerate its partial-reference numbers to obtain a pass.

See [madqc-test-index.md](madqc-test-index.md) for the broader scientific
coverage matrix and [madqc-acceptance-checklist.md](madqc-acceptance-checklist.md)
for the release contract.
