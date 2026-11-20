# madqc acceptance checklist

**Driver consolidation — 2026-11-19.** The sole supported application entry is
`madqc`; named stages and duplicate engine executables are retired. The
[driver/test migration index](madqc-driver-test-migration.md) records current
test replacements and remaining gaps. Older named-test references below describe
the pre-migration inventory, not active targets.

Reviewed for acceptance coverage, 2026-11-13. This checklist indexes the agreed
[design](../../design_draft.md); it does not report passing implementation tests.
Nemo definitions await Florian. CC2/6D work stays outside this walkthrough;
release scope is unchanged.

Each stable ID records a purpose and expected outcome. Map IDs to existing or
adapted tests during implementation. Preserve numerical references and their
scientific meaning independently of output filenames and JSON paths.

Release entry-point amendment, approved 2026-11-19: verify supported engines
through madqc. The release must contain no named-stage interface and no separate
application executables for those engines. Migrate their scientific checks before
retiring the old entry points. Preserve separate builds for unsupported engines;
private numerical test binaries are not additional public applications. This
does not add optimizer replacement or CC2/6D work to the current increment.

All 16 cases (A01–A06, B01–B03, C01–C04 and D01–D03) were reviewed and
accepted for coverage. On 2026-11-13 the author deferred optimizer replacement
until after workflow unification: D01/D02 remain reviewed follow-on requirements,
not gates for the current work. A01–C04 and standalone D03 remain in scope.
Scientific confirmations, implementation tests and independent review remain
separate requirements.

OEP release coverage applies the relevant measurement, restart, checkpoint and
output cases to all supported OAEP/OCEP/DCEP/MRKS calculations, including fresh,
resumed and refined runs. Preserve the OAEP case and supply the missing model
migration evidence. Record pre-existing limitations separately with affected
conditions and supporting evidence. Missing tests remain release gaps; they
cannot be closed by relabeling an existing model unsupported.

## 1. SCF: inputs, evaluation and restart

Contract: design sections 1, 2, 4 and 7.

| ID / purpose | Scenario and expected outcome |
| --- | --- |
| A01 — consistent input | Equivalent supported file/CLI settings enter the same lifecycle. Explicit `excited.maxiter` gets the generic unknown-input error, including alongside valid response `maxiter`; no custom migration message. |
| A02 — trustworthy measurements | Evaluate before updating. Passing criteria cause no KAIN update or step restriction; continued iteration reuses raw BSH work. Saved metrics describe the saved functions. |
| A03 — final representation | Intermediate targets or caps permit advancement; a plateau alone does not. Always reach the requested final threshold and k, with unchanged final limits. |
| A04 — evidence-based reuse | Qualifying final saved evidence permits reuse without numerical evaluation, update or archive rewrite. Unmet criteria resume work; missing supported evidence requires fresh measurements. Never manufacture zero change from a loaded state. |
| A05 — iteration allowance | Restart gives a fresh per-protocol main-update allowance. Terminal evaluation permits no extra update. Warmup is off; KAIN history restarts empty. |
| A06 — metric meaning | Verify declared reductions, units, weights and sampling. One large residual cannot disappear in an average. Selected applicable metrics alone decide convergence using agreed absolute/relative scales. Missing/nonfinite required evidence cannot pass; unsupported selections fail validation. |

**A02/A05 clarification approved 2026-11-18.** Final SCF rotations, truncation
and normalization precede the definitive evaluation and stopping decision.
Cover a provisional pass that fails after rotation: resume from the rotated
orbitals with their fresh residual if updates remain, resetting KAIN when the
orbital basis changed; otherwise save a valid unconverged final endpoint.
Retain the intermediate-cap rule. A post-rotation pass performs no further
update or rotation. Change metrics use the preceding comparable iteration,
not the pre-rotation state, and missing comparisons remain unavailable.
Verify that loaded checkpoint functions match the published metric evidence;
the extra evaluation changes timing, not the main-update count. This is required
SCF migration coverage, not a claim about existing engine tests.

**A03/A05/A06 protocol clarification approved 2026-11-18.** Carry orbitals
forward and test both automatic `1e-4, k=6` to `1e-6, k=8` conversion and a
threshold change with explicit fixed k. Verify the new representation and
operator settings, fresh BSH work, empty KAIN and cleared comparison history.
The first evaluation has zero main updates and unavailable change metrics;
the next comparable iteration uses the new baseline. Preserve old analysis
and checkpoint measurements without treating them as fine-protocol evidence.
Preparation consumes no update; the new protocol gets its own allowance.
Smoothing does not change implicitly. Numerical adapter coverage remains pending.

A01 includes repeated CLI groups/options in both orders: reject every retired,
unknown or malformed occurrence before applying overrides, including an invalid
first group followed by a valid one. Preserve file-then-CLI, last-group replacement,
path case and winning-input provenance. Valid earlier payloads are checked but
not merged or applied; cross-setting checks use the resolved request. Validation
failures perform no solver work or checkpoint writes. Both occurrence validation
and final CLI reading distribute the same error type and message to every MPI
rank before throwing. A bounded two-rank test covers retired keys (including
masked earlier occurrences) and malformed accepted values, requiring nonzero
exit, the intended diagnostic and no solver/checkpoint activity. Timeout or an
unrelated launcher failure cannot satisfy the case. Keep serial/non-MPI coverage.

A04 example: a valid final residual of `5e-7` passes an absolute limit of
`1e-6`, but fails `1e-7`. Every other selected criterion and compatibility
check must also pass. These are illustrative values, not defaults.
A supported conversion from an unmatched saved protocol supplies a guess at
the first requested protocol, with fresh comparisons and no convergence claim.

**A04 clarification approved 2026-11-18.** Separate archive usability, matching
calculation settings and current tolerance checks. With unchanged settings,
passing saved evidence permits reuse and insufficient evidence resumes work at
the preserved protocol position. Changing only stopping criteria/tolerances or
maxiter must not change the scientific context. Changed geometry, model/reference,
actual frequency or operator/smoothing settings must prevent reuse of old
convergence claims even when their numbers pass; supported functions may still
start the continued solve with fresh measurements. Check representation and
metric definitions explicitly. Unsupported continuation preserves the checkpoint
and reports failure. No automatic nearby-state search or discretionary policy is
introduced. Application settings-comparison tests remain migration work.

SCF virtual-path coverage extends A02–A06 without adding a new workflow. Check
both spin blocks, exact frozen-occupied references, occupied/virtual membership,
alpha-only and beta-only padding, and an interrupted extra-virtual schedule.
A complete final-membership checkpoint can skip preparation; an interrupted one
resumes its saved phase. Added virtuals retain matching SCF protocol progress but
need new convergence evidence. `nv_its` is a separate preparation-update cap.
Fixed occupied-density/energy/kinetic changes cannot certify virtual convergence;
virtual residuals and eigenvalues must describe the saved functions after any
rotation. Include bound-virtual success and honest iteration-limit outcomes for
the PR's documented continuum case. Preserve explicit occupations and KAIN
numerics; inspection never pads missing functions. The
[source record](madqc-metric-validation.md#scf-virtual-orbitals--pr-798-and-source-review-2026-11-13)
separates PR/source findings from implementation validation still required.

A01/A03/A04 also cover seed entry: a new basis-set seed with a single fine
protocol starts there with fresh MRA checks; an explicit multi-stage protocol
runs from its first entry without hidden ladder rewriting. Matching checkpoint
restart retains its progress. Reject `seed.start_rung` and the test driver's
`--seed-start-rung` for every explicit value, including masked earlier arguments.

## 2. Response: independent tasks and coupled roots

Contract: design sections 3 and 4.

| ID / purpose | Scenario and expected outcome |
| --- | --- |
| B01 — independent progress | Reuse satisfactory final states; resume unmet final or interrupted intermediate work; advance completed intermediate work. Matched 2PA responses with refined frequencies retain their functions and requested saved protocol, with fresh convergence checks. Genuinely new/unmatched states start first. Unrelated seeds do not transfer progress. |
| B02 — dependencies | Response requires an eligible final ground state. VBC requires eligible B/C; the final quadratic property also requires A. Direct completion has no iterative metrics. Changed source generations invalidate dependent evidence. |
| B03 — coupled roots | Reuse requires a valid bundle and all required roots passing, with no extra aggregate gate. Otherwise all roots stay active. Change metrics share valid root correspondence; ambiguous comparisons cannot pass selected change criteria. Added roots require renewed bundle convergence. |

**B02 exact-input clarification approved 2026-11-18.** VBC/property builders
load the exact published generations qualified by the current final-protocol
metric checks. Coarse/partial guesses and legacy cap acceptance cannot release
direct work. VBC records B/C/reference inputs and completion; properties bind
all sources their formulas use. Changing A alone preserves valid VBC; changing
B or C rebuilds affected VBC/properties. Unchanged sources satisfying tighter
limits reuse without rebuilding. Missing/corrupt bound sources fail strictly;
newer `latest` pointers never rebind an existing result. Molresponse retains
dependency ownership. Curate these pending cases under
[RESPONSE-DIRECT-INPUTS](madqc-test-index.md#response-direct-inputs--qualified-sources-for-direct-computations).

**B01–B03 extension approved 2026-11-17: refined 2PA frequencies.** Match the
saved response by the same tracked root and perturbation, with unchanged ground
reference and other relevant scientific settings. A changed photon frequency
retains the requested saved protocol but forces new equation work, metrics and
comparison baselines, empty KAIN and a fresh allowance. Even a previously
completed intermediate protocol is repeated at the new frequency before normal
advancement. Exact unchanged-frequency/context matches retain ordinary
saved-metric reuse. Root reordering preserves the association; an ambiguous or
unmatched root cannot receive another root's progress. A removed protocol uses
the existing supported-conversion fallback, including exact k matching.

Check that old checkpoint identities and published result sources remain
unchanged; updated 2PA sources/properties use the new root and response
generations. Missing/corrupt selected sources fail strictly, and interruption
during replacement preserves the required old payloads. Curated coverage is
[RESPONSE-2PA-CONTINUATION](madqc-test-index.md#response-2pa-continuation--refined-photon-frequencies).
This is approved coverage to implement, not a claim of passing engine tests.

**Replacement/retirement extension approved 2026-11-17.** Preserve existing
response task IDs and record exact `replaces` references. Temporary continuation
inputs protect old functions; replacement history alone must not pin archives
forever. Retirement of an obsolete task waits for current results and every
other active consumer, including another property sharing the old frequency.
Persist retirement before deletion and leave other tasks in the calculation
active. Same-ID frequency refinements require a new generation and fresh checks
at the actual frequency, never retirement of the continuing task. Verify these
cases with interrupted replacement and immutable historical result sources.

Nemo CPHF extends A02/A04–A06 and B01/B02: atom/axis tasks retain independent
progress under the exact final ground reference. Check reuse-or-iterate for both
CPHF restart inputs, fresh/inspection conflicts and a changed reference. Residuals
and saved iterated functions share a sample; corrected full responses carry their
sources. No extra update occurs at convergence/cap, and no automatic LDA warmup
runs. A missing or unconverged perturbation blocks the requested Hessian and
vibrational properties. Validate owner-confirmed density/kinetic/correction
conventions and the Hessian-row check against curated numerical references.

Full molresponse acceptance additionally requires the
[test-index coverage matrix](madqc-test-index.md#full-molresponse-acceptance-goal):
all supported response properties, one-/multiple-electron systems, closed/open
shells, static/dynamic and equation variants, exchange/correlation selections
and terms, functional components and complete workflows. Require both accuracy
and controlled performance regression evidence. The selected helium HF and water
LDA alpha pair is the first validation slice, not full response certification.
Existing unsupported combinations retain explicit tested outcomes.

## 3. Checkpoints and output

Contract: design sections 5 and 6.

| ID / purpose | Scenario and expected outcome |
| --- | --- |
| C01 — recoverable publication | Interrupt before/after publication: current pointers select complete matching native payloads/records. Expected missing, corrupt or unsupported checkpoints fail strictly. No first publication differs from a missing expected checkpoint. Adopted paths reject HDF5 requests without fallback or checkpoint writes; shared HDF5 support is deferred. |
| C02 — bounded storage | Save endpoints/controlled stops by default; optional intervals work. Protect latest non-retired tasks, current-result sources and active dependencies transitively. Reclaim obsolete payloads only after successful publication, complete reference checks and explicit task retirement where needed. Historical maps alone do not pin archives; failed publication or unreadable live roots prevent deletion. |
| C03 — understandable results | Group requested tensors; record state metrics once and availability explicitly. Separate properties, runtime evidence and timings. Incomplete work remains visibly incomplete; reuse preserves measurement identity. |
| C04 — analysis and regression | Timings expose application steps and each state's iteration steps. Inspection does no solver work or new convergence assessment. New output compares the same physical quantities against preserved references. |

The November 19 [stdout preservation audit](madqc-stdout-preservation.md)
adds source evidence and pending coverage under C03/C04: preserve the engine's
energy/orbital diagnostics and detailed timers alongside shared progress,
restore missing final-orbital timing attribution, and keep log destinations and
timing scopes clear. These gaps are not closed by the earlier broad-timing tests
or by marking the original acceptance checklist reviewed.

CORE-CHECKPOINTS in the test index curates the shared native implementation of
C01/C02, including numerical restart when scientific output is missing and
retention of manifest identities after payload reclamation. Engine adoption and
inspection integration remain later acceptance work.

For C03/C04, [RESULT-SOURCES](madqc-test-index.md#result-sources--immutable-result-attempt-binding)
verifies the root run attempt and each calculation's explicit or inherited
source-map attempt. Two electronic attempts under one evaluation must leave an
older snapshot bound to its original map, even while the newer attempt is
selected but unfinished or after it publishes. Cover mixed-attempt collections,
older reused numerical samples and strict missing/mismatched-map errors. Result
inspection uses the snapshot binding, never current selectors. Retention follows
that binding for the current result; historical metadata keeps its original
sources after unneeded payloads are reclaimed. Loading reclaimed functions
reports unavailability without substituting another generation. RESULT-SOURCES
also covers C02 publication ordering, transitive protection and safe reclamation.

For C04, [INSPECT-NATIVE](madqc-test-index.md#inspect-native--saved-function-inspection)
checks the existing visualization writer and native saved-function reader against
a published unfinished sample. Shared records must preserve its exact source,
geometry, representation, available metrics and original incomplete status.
Readable functions remain discoverable without an available scientific property.
Inspection performs no solver evaluation/update, new convergence assessment,
protocol advance, checkpoint publication or completion change. Requested plot
files are allowed; solver records and native payloads remain unchanged. Activate
this check with each supported SCF/response reader migration, preserving existing
limitations and leaving shared HDF5/new export features deferred.

## 4. Optimization follow-on and finite differences

Contract: design section 6 for standalone D03; the
[optimizer follow-on requirements](madqc-optimizer-follow-on.md) for D01/D02.
Preserve the existing optimizer during unification; its SCF/Nemo calls must still
meet the electronic acceptance cases, including exact-geometry and current-request
checks before reuse. Verify independent defaults `gtol=1e-4` and
`gradient_precision=1e-5`, explicit optimization-group file/CLI overrides and
unchanged defaults when electronic criteria change. Keep the existing precision
clamp and other optimizer behavior. This does not promise saved optimizer history.

[OPT-ELECTRONIC](madqc-test-index.md#opt-electronic--retained-optimizer-electronic-records)
checks invocation-scoped evaluation IDs and exact-geometry/model reuse. Cover two
geometries, a repeated geometry, insufficient measurements and missing gradients;
keep published request/source records immutable and apply RESULT-SOURCES to
old snapshots after electronic selection changes. Test selection publication
failure and interruption before/after selection, including a selected null task.
Recovery uses only the selected electronic calculation and strict checkpoint
rules. It does not restore an accepted geometry, step, Hessian or line search.
Verify replacement selection and reference checks precede archive retirement.
These electronic checks are current scope; D01/D02 remain deferred.

| ID / purpose | Scenario and expected outcome |
| --- | --- |
| D01 — resume once (deferred) | Saved evaluation numbers match geometry-specific SCF/response work. Restart neither repeats committed moves nor applies Hessian updates twice. Energy/gradient share geometry, target and compatible sources. |
| D02 — valid history (deferred) | Preserve valid optimizer history; invalidating settings reset it at the saved current geometry. Restart gives a fresh geometry-step allowance. Trials, electronic iterations and resets do not replenish it. |
| D03 — valid differences | Changed point geometries get new identities; unchanged points retain theirs. Tightening reuses qualifying points and continues others. Assemble only a complete eligible set for the saved formula/displacement/target; protect required generations. |

For D03, [FD-PLAN](madqc-test-index.md#fd-plan--finite-difference-plan-publication)
checks the single versioned `checkpoints/finite_difference_plan.json` and atomic
replacement. The file owns revision/counter/point/source assignments; the run
inventory records only that it is required. Test failed replacement and first
registration interruption, changed-point allocation, child publication before
source binding, strict missing/corrupt/unsupported-plan errors and retirement
ordering. No alternate plan selection or optimizer history is introduced.

D03 can first use deterministic providers; use the same approach for D01/D02
when replacement work begins. New scientific gradient paths still require
validated providers and consistent target selection. Electronic tolerances alone
do not bound finite-difference error.
