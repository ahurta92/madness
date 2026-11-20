# Preserve useful chemistry output during madqc unification

Source and captured-output audit, November 19, 2026. The author requires the
useful MOLDFT/SCF progress, energy diagnostics and timing information to remain
available on screen as applications adopt the shared layer. Shared records and
summary lines do not replace that information. This audit records what is
present, what was lost, and the focused reporting work still needed. It changes
no numerical code, print levels, output destinations or timer implementations.

## Evidence and scope

Compared `SCF.cc` at drop-6 commit
`9179f5d32076a2ba555906dd146dc330ffc040d7` with the current working tree based on
`76255552e16e542f6b5ae2d0d7e2b6b034c7ed48`. Inspected the current driver, shared
printers, response banners, timing helpers and existing output tests.

Captured examples include an older H2 fixed-work run (its banner identifies a
November 6 build based on `dc2062490`, dirty), the current helium SCF regression
log, an SCF reuse log, and fresh/reused driven-response logs. These show output
content, not performance equivalence: the molecules, requests and execution
sizes differ. The old source comparison is the authority for which SCF print
sites existed, rather than assuming every optional branch ran in these samples.

The complete SCF print/timer call-site inventory, distinct timer labels, five
captured logs, their original paths and verification evidence are retained in
`.superpowers/sdd/design_draft/stdout-preservation-audit/` in this worktree.
No new chemistry calculation was run for this audit. The user's example
directory was only read.

## Where the output actually goes

Engine diagnostics go to **stdout**, normally the terminal or the Slurm job
log. madqc also prints its final summary to stdout and appends that summary to
`<prefix>.out`. It does not automatically copy all stdout into that file.
Consequently, an `.out` file can contain only a short summary even though the
engine printed detailed timings during the calculation.

If a launcher already directs stdout to `<prefix>.out`, both paths write there.
The current append behavior prevents the driver from truncating the job log;
it can also duplicate the summary. Keep that no-truncation guarantee. For
current runs, capturing stdout/stderr in a separate job log is unambiguous.
Choose any future full-log/summary naming change explicitly, without silently
changing existing command-line behavior.

Source: [`WorkflowDriver::run`, summary publication](../../src/apps/madqc/WorkflowDriver.hpp#L184).

## SCF output inventory

The levels below describe an ordinary fresh SCF process. Optional model and
diagnostic paths may add output; `print_level 0` is not a promise of a completely
silent program.

| Information | Current status and visibility | Preserve or clarify |
| --- | --- | --- |
| Build/runtime identity, geometry, resolved parameters and working directory | Present in startup/workflow output | Keep enough context to interpret and reproduce the calculation. Preserve readable parameter tables. |
| Protocol threshold and polynomial order | Present; SCF banner at level 2+, shared metric rows at level 1+ | Keep. Replace the banner's ambiguous legacy `conv` value with the actual resolved shared criteria. |
| Initial guess, imported orbitals and archive loading | Existing engine messages; shared load/check elapsed time at level 1+ | Keep source and starting protocol clear. Explain reuse versus continued iteration, including which requested checks remain unmet. |
| Iteration number and elapsed time | Legacy iteration header at level 2+; shared iteration/update index at level 1+ | Keep. Explain that a terminal evaluation can occur without another update. |
| Energy breakdown | Present at level 2+: kinetic, nonlocal pseudopotential, nuclear attraction, Coulomb, PCM, exchange-correlation, nuclear repulsion, optional D3, total | Preserve all applicable terms and the total. These quantities are useful even when energy/kinetic change is not selected as a stopping criterion. |
| Per-orbital BSH norms, RMS and maximum | Present at level 2+ | Preserve detailed norms and label the maximum as the shared convergence reduction. RMS is diagnostic. |
| Density change | Shared maximum spin-density change at level 2+ and when selected at level 1+; legacy separate alpha/beta line remains on the fallback path | Preserve the aggregate and make separate spin detail available at a suitable detailed level. Do not substitute total charge-density change without changing its label and convention. |
| Required metric values, limits and pass/fail | Shared rows at level 1+ | Keep. Missing history says unavailable. Measured but unselected changes are not currently all shown; distinguish them as reported diagnostics rather than extra stopping requirements. |
| KAIN, step restriction and orthonormalization | Step-change RMS/max at level 3+, coefficients/subspace details at level 4+, existing warnings and relevant timers | Preserve numerical diagnostics; distinguish update size from the raw BSH residual. |
| Orbital spectrum and analysis | Alpha/beta spectra and MO analysis at level 2+, occupations, centers/radii and coefficient counts | Preserve. Keep localized Fock expectations distinct from the canonical spectrum used for DOS. |
| Frozen-occupied virtual solve | Virtual-block residuals, eigenvalues, occupied energy labelled as coming from the archive | Preserve that distinction; do not present the occupied energy as a newly converged virtual energy. |
| Requested physical properties | Dipole and gradient/property printers remain in application paths | Preserve labels, units and requested detail. Optional models need their own coverage beyond the helium example. |
| Completion, cap and failure | Shared decision reasons during iteration; completion and representation in the final summary | Keep an unmet final request visibly incomplete. A finished intermediate protocol is not final convergence. |
| Calculation and process totals | Application wall interval/execution size in the summary; MADNESS process footer and runtime statistics | Preserve both scopes and distinguish them. A reused result must report this run's actual work, not replay earlier iteration timings. |

Sources: [`SCF::set_protocol`](../../src/madness/chem/SCF.h#L327),
[energy/progress block](../../src/madness/chem/SCF.cc#L3397),
[BSH details](../../src/madness/chem/SCF.cc#L2230),
[step restriction](../../src/madness/chem/SCF.cc#L2880),
[final analysis](../../src/madness/chem/SCF.cc#L3515),
[shared printers](../../src/madness/chem/IterationOutput.h#L31),
[summary timing](../../src/madness/chem/ResultsSummary.hpp#L427).

## What the output looks like

This excerpt from the older H2 log shows the useful style to preserve. Other
operation timers and energy terms are omitted here only to keep the example
short:

```text
Iteration 0 at time 1.4s
timer:              Coulomb     0.16s     0.16s
timer:          HF exchange     0.12s     0.12s
              kinetic       1.20630737
   nuclear attraction      -3.71599219
              coulomb       1.34190684
 exchange-correlation      -0.67095445
    nuclear-repulsion       0.71428571
                total      -1.12444672
BSH residual: rms 3.757231e-02    max 3.757231e-02
```

The current helium log still contains detailed timers and the energy table.
The shared lines add an explicit verdict tied to the evaluated orbitals:

```text
              kinetic       2.86166808
   nuclear attraction      -6.74911394
              coulomb       2.05153320
 exchange-correlation      -1.02576661
                total      -2.86167928
SCF scf  threshold=1.0e-06 k=8  iteration=2  final_criteria_met
  bsh_residual      2.6040e-06  <= 1.0000e-04  met
  density_change    4.6215e-06  <= 1.0000e-04  met
  energy_change     3.4781e-09  <= 1.0000e-05  met
```

The same endpoint reports `prepare/evaluate/decide=12.827 s`,
`checkpoint/report=0.038 s`, and `total=12.864 s (driver rank)`, with no
`advance` entry. That preparation interval includes final preparation and
reevaluation; it is not just one ordinary BSH update. An existing SCF reuse
log instead reports `load/check=0.006 s` and `saved_criteria_met`, with no
new iteration timing. These numbers are recorded examples, not benchmarks.

## Timing inventory and identified gaps

The legacy SCF timer uses a stack, so nested intervals remain nested. Its lines
contain two columns: `cpu_time()` then `wall_time()`, printed without column
headings and rounded to two decimals. The label is truncated to 20 characters.
On this ARM host, `cpu_time()` uses `steady_clock`; on x86 it uses a cycle
counter. **It is not process CPU consumption summed over worker threads.** The
separate `process_cpu_time()` helper has that meaning. Preserve the existing
timer methods as previously agreed, but describe their scope honestly.

`START_TIMER` already fences; the shared collector adds no new fences and uses
`steady_clock` at the existing evaluation boundaries. Neither the detailed
rank-zero intervals nor shared driver-rank intervals are a sum over MPI ranks.
Do not add inclusive nested timers together or treat their sum as the
application's elapsed time. The process footer includes work outside the SCF
application interval.

Sources: [SCF timer implementation](../../src/madness/chem/SCF.cc#L118),
[`cpu_time` versus `process_cpu_time`](../../src/madness/world/timers.h#L123),
[application elapsed interval](../../src/madness/chem/Applications.hpp#L298),
[runtime footer](../../src/madness/world/world.cc#L424).

The source comparison found **51 distinct SCF timer labels before the split
and 47 now**. The retained families cover initial guesses/projection,
localization, load balancing, density construction, nuclear/Coulomb/PCM work,
HF exchange and XC, potential application, kinetic/Fock construction, BSH,
KAIN/subspace transforms, orthonormalization, properties and orbital analysis.

| Gap | Evidence | Recommended focused work |
| --- | --- | --- |
| Four final-orbital timing labels disappeared | `Overlap`, `focka eigen sol`, `fockb eigen sol`, `Transform MOs` are absent from current SCF timer sites. The corresponding work remains in `prepare_final_orbitals` and final-spectrum evaluation. | Restore equivalent phase timing around the existing work, retaining alpha/beta attribution. Do not repeat numerical operations merely to time them. |
| Old convergence banner remains | `SCF::set_protocol` prints `conv=max(thresh,param.dconv())` even for shared solves | Print the current protocol and actual shared metric limits; avoid displaying the dormant legacy default as the stopping rule. |
| Only selected shared metrics appear in the shared table | `write_iteration_metrics` iterates `decision.checks` | Report available requested diagnostic information separately from the selected stopping requirements. In particular, retain energy and kinetic quantities and make their measured changes visible without requiring them to stop the solve. |
| Broad timings hide final reevaluation and archive cost | `prepare/evaluate/decide` and `checkpoint/report` combine several operations | Keep the useful broad totals; add attribution for final preparation/reevaluation and archive I/O where already available or cheap to measure. Do not call report overhead pure archive time. |
| Timer columns and nesting are unclear | Unlabelled legacy columns, truncated operation names; inclusive nested intervals | Explain the clock/scope, use clear labels, and avoid implying process CPU accounting or additive totals. Very short intervals rounded to zero are not evidence that no work happened. |
| Restart explanation is incomplete | Shared SCF prints load/check time and final `saved_criteria_met`; it does not print the saved metrics against current limits on the reuse path | Show loaded protocol, reuse/continue decision and relevant saved/requested metrics without a new scientific evaluation. |
| Log destination is easy to misunderstand | Summary is separately appended to `<prefix>.out` | Document the difference from captured stdout; preserve existing logs and address duplicate summaries when destinations coincide. |
| Detailed output coverage is incomplete | Current tests assert broad timing/decision lines and reuse behavior, not every legacy energy term or lost final-phase timer | Add semantic coverage for preserved information; do not compare timing numbers or require exact whitespace. |

The timer visibility flag is currently process-global: an SCF constructor turns
it on at level 3+, but does not turn it off for a later quieter SCF instance.
Record this existing limitation when making verbosity reliable across repeated
calculations; do not claim print levels are fully isolated today.

Four optional `print_meminfo` call sites also disappeared from density,
load-balancing, Coulomb and PCM work. These write profiling files when
`WORLD_MEM_PROFILE_ENABLE` is enabled, not ordinary stdout. Keep them in the
profiling-preservation inventory without confusing them with the missing timing
lines. Source: [memory-profile helper](../../src/madness/world/worldmem.h#L166).

## Response output to retain during adoption

Legacy driven output identifies channels/frequencies, BSH and density progress,
distance from each tolerance, and (at higher verbosity) the source overlap
`<v|x>`. Excited-state output includes root frequencies in hartree and eV,
per-root residual/density details, final status, and optional subspace/debug
diagnostics. Preserve the information and make its metric conventions explicit;
do not equate every source overlap with a final assembled property.

Shared driven output already adds per-state evaluation/decision/checkpoint/update
times, native load/check intervals, and application planning/solve/property
assembly totals. These supplement detailed engine diagnostics. The optional
`ES_TIMING` instrumentation is in one non-locking iteration path and is off by
default; its existence does not establish universal per-root or per-step timing
coverage. Internal direct tasks such as VBC still need explicit timing coverage
as their shared adoption proceeds.

Response currently maps deck levels 1–3 to Normal, 4–9 to Verbose, and 10+ to
Debug. SCF's thresholds differ. Document what users get before attempting to
standardize that behavior. This is a first response-source inventory, not an
executed output audit of every response property/model.

Sources: [FD header and progress](../../src/apps/molresponse/solvers/fd_solver.hpp#L347),
[ES progress](../../src/apps/molresponse/solvers/es_solver.hpp#L257),
[optional ES phase timing](../../src/apps/molresponse/solvers/es_solver.hpp#L804),
[response level mapping](../../src/apps/molresponse/madqc_adapter.hpp#L90),
[shared loop](../../src/apps/molresponse/solvers/iterate.hpp#L36).

## Preservation checks for the next implementation

Keep the existing verbosity control. Level 1 provides compact shared progress
and broad times; level 2 preserves SCF energy/orbital detail; level 3 preserves
detailed SCF operation timers. Detailed/debug diagnostics remain available at
their established levels until a deliberate, documented change is made.

Check information rather than exact screen layout: energy components sum to the
reported total within printing precision; quantities and units are labelled;
selected metrics show their actual limits; absent history is not a false zero;
the final values describe the checkpointed functions; caps remain incomplete;
reuse prints actual current load time without fictitious iterations; stdout
capture is not truncated; and detailed final-orbital timing is present again.
Include an unrestricted and localized case for spin/spectrum distinctions, and
an optional-model case when its dependencies are available. Reuse existing
curated cases and do not regenerate numerical references to accommodate output
formatting. Timings must be finite/nonnegative, not pinned to reference seconds.

The earlier stdout-cleanup checks already cover broad timing, selected metrics,
missing-history handling, endpoint/no-update behavior and reused-checkpoint
timing. Their passing results do not close the newly identified preservation
gaps. No new numerical accuracy or performance claim is made by this audit.


## Approved aligned output implementation — November 19

The follow-up implementation uses `IterationOutput.h` for the adopted SCF and
driven-response paths. Applications supply model-specific labels; the shared
printer aligns numerical columns, keeps state identifiers on separate headings,
reports measured diagnostics even when unselected, and distinguishes unavailable,
inapplicable and invalid values. Limits are the controller's actual effective
limits, including relative scaling; the printer does not evaluate orbitals.

SCF energy terms use the same value-table layout as wall-time intervals. The
stale legacy `conv` protocol banner and duplicate shared density/residual line
are removed. Existing detailed engine timers remain unchanged. At print level 3,
chronological phase headings distinguish initial preparation and final
reevaluation, and new steady-clock intervals restore visibility of final
orbital overlap, eigensolve and transformation work without new fences. These
intervals are included in the broad iteration time, and asynchronous operations
can carry work into the next observing operation; they are not standalone
performance benchmarks. Localized and canonical preparation keep their own
labels and computations.

A successful publication prints `Checkpoint saved` only after it returns.
Restart output shows saved metrics/current limits when scientifically comparable,
with an explicit reuse or fresh-check statement. Reuse does not replay saved
iteration timings or claim another publication. Intermediate completion is
reported separately from final convergence.

This is the output convention for further application adoption. Legacy response
families and Nemo are not newly migrated by this reporting increment; their
existing detailed output remains available. The earlier audit excerpts are
historical evidence, not the current formatting. Verification and the scoped
diff are recorded in `.superpowers/sdd/2026-11-19-aligned-qc-output/`.


The independent review added three reporting clarifications: unavailable selected
metrics retain a known absolute-only limit and are marked `required`; completion
claims apply to applicable selected criteria; and extra-virtual preparation uses
its actual intermediate limits and workflow-selected bundle disposition on load.
The SCF density label explicitly names the maximum over spin channels. These
are presentation changes, not new qualification or preparation policies.


The regression coverage is organized by meaning: `test_shared_iteration.cpp`
checks alignment, diagnostic versus required values, missing comparisons,
inapplicable selections and relative limits without chemistry; `test_shared_scf.py`
checks current versus saved measurements, final-orbital timing and intermediate
preparation/restart reporting through the public driver; `test_shared_response.py`
checks fresh, continued and reused response states, including the property-only
load path. Their scientific references remain separate and unchanged. This
coverage should survive future wording changes: preserve the information and
numerical relationships rather than freezing an entire stdout transcript.

## Persist SCF timings alongside response — November 20, 2026

The output walkthrough exposed a remaining gap: SCF printed measured iteration
intervals but did not publish them to the shared `timings.json`. The focused
follow-up uses those same intervals for both destinations. Ordinary SCF and
frozen-virtual SCF forward completed iteration measurements even at print level
zero. The saved breakdown keeps the honest combined
`prepare/evaluate/decide` interval, checkpoint/report work and any actual update;
it does not repeat computation to obtain a finer breakdown.

SCFRun owns measured preparation, protocol/task, loading and result-publication
scopes. Successful runs and failures use the existing collective timing writer.
Each calculation retains its existing invocation ID, so one timing file can
contain both ground-state and response entries without introducing a second
workflow identity. Earlier entries remain unchanged. A reused SCF calculation
records current preparation/loading/publication work and no solver iterations.
The outer driver summary still includes work outside the engine collector;
those totals have distinct boundaries and are not expected to be identical.

Fresh/reused/tighter public workflows pass on one and two ranks; quiet SCF,
frozen virtuals, preparation recovery and collective failure preservation pass
the focused checks. The exact evidence is tracked in
[the focused plan](../superpowers/plans/2026-11-20-shared-scf-timings.md).
The subsequent [run-layout migration](../superpowers/plans/2026-11-20-shared-run-layout.md)
places shared results, timings and checkpoints at the run root and moves engine
diagnostics into analysis/ground_state and analysis/response. That relocation
preserves archive contents, identities and the existing retention rules.
