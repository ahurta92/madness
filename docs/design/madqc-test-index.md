# madqc test index

**First response integration milestone — author clarification.** Connect all
supported response paths and check representative reasonable convergence, values,
metric consistency and records. Historical iteration-count/cost guards below
remain recorded for later calibration; count differences alone are not this
milestone's acceptance target. Recomputing reference runs and completing the
broader scientific/performance matrix are later validation work. Existing
numerical oracles remain unchanged and clear numerical failures still require
investigation. The release-wide inventory below remains in force for release
validation, distinct from this first integration milestone.

**Driver consolidation — 2026-11-19.** The sole supported application entry is
`madqc`; named stages and duplicate engine executables are retired. The
[driver/test migration index](madqc-driver-test-migration.md) records current
test replacements and remaining gaps. Older named-test references below describe
the pre-migration inventory, not active targets.

**Initial alpha curation — 2026-11-16.** The author selected helium HF plus one
DALTON-reference water case. Use water LDA for this first pair. The entries below
record existing tests and planned migration; no new numerical run or shared-layer
acceptance is claimed. This is the beginning of the release-wide index, not its
completion. Stable IDs describe scientific intent and survive input/output changes.

## Full molresponse acceptance goal

### RESPONSE-SHARED-LADDER — independent progress across protocols

- Purpose: adding or refining independent response states must preserve completed
  work and qualify the actual functions used for each property.
- Case: `response_shared_ladder/response_shared_ladder.in`, helium HF alpha at
  0 and 0.05 au, protocols `[1e-4,1e-5,1e-6]`. The last transition keeps k=8.
- Harness: `test_shared_response_ladder.py`, serial and two ranks; planner/floor
  checks in `test_shared_fd_planning.cpp`. Partial qualification, exact reuse
  without evaluation or payload rewriting, deleted analysis, tighter requests,
  new frequencies, an interrupted intermediate checkpoint, and all tensor axes.
- Review regressions: `test_shared_response_ladder_review.py` checks frequency
  collisions across runs, the static archive boundary at 1e-13, removal of stale
  alpha rows, preservation of unrelated analysis, and exact source bindings after
  all-state reuse. A stale results snapshot fails the reuse guard.
- Scientific reference: unchanged `response_he_alpha` alpha values at each
  matching frequency, absolute comparison tolerance 1e-3. Tensor diagonals also
  obey helium isotropy to 1e-3. Format changes must preserve these claims.
- Scope: this does not establish open-shell, XC-family, near-resonance or broad
  speed/accuracy coverage; those remain release gaps.

### RESPONSE-SHARED-ALPHA — driven state restart

**Purpose:** keep the numerical meaning of alpha separate from the mechanics of
restarting it. The public `test_shared_response.py` runs one helium HF z state
through the shared controller and native checkpoint store, at zero and nonzero
frequency. Helium has two electrons and one occupied spatial orbital.

**Contract:** evaluation precedes KAIN, a capped endpoint has current-state
metrics, missing history stays unavailable, satisfying tighter requests perform
zero evaluations, unmet requests get a fresh update allowance, periodic saves
are reclaimed, and damaged selected archives fail. A zero-update reload checks
the saved residual and kinetic quantity directly. The final converged alpha is
compared to the unchanged `response_he_alpha` reference with absolute tolerance
`1e-3`; large and zero gates in the lifecycle steps make no scientific accuracy
claim. The same test runs on two MPI ranks. `test_shared_iteration` separately
checks update/evaluation order and retained proposal use.

`test_shared_response_review` rejects malformed property/direction selections and
unsupported subworlds before response preparation. A changed SCF source must be
refused without changing the registered response inputs or selected generation.
`test_response_property_metadata` forces malformed reads and failed alpha writes
on two MPI ranks; both ranks must receive the same error promptly and a failed
write must preserve the previous record.

**Evidence:** `.superpowers/sdd/design_draft/response-driven-alpha/`. This is
initial shared-path coverage, not full response or performance qualification.
Water DALTON-reference migration, open shells, other functionals/properties,
multi-state/protocol workflows and detailed shared timings remain release gaps.

### Remaining release coverage

The initial alpha pair is a starting point. Before accepting the full molresponse
migration, inventory and test all supported functionality across the dimensions
below. Every supported path needs linked evidence; missing coverage remains a
release gap. Record existing limitations separately and test their expected
rejection/unavailability behavior. Source code determines support; a guide or a
missing test alone cannot establish it. Earlier engine-scope decisions remain.

| Dimension | Required coverage |
| --- | --- |
| Requested properties | All exposed response properties: alpha, beta, excitation energies and supported transition observables, Raman/polarizability derivatives, 2PA tensors and derived observables, and other supported selections. Cover selected components and full tensors, identity, units, availability and exact contributing states. |
| Electron count and spin | One-electron and multiple-electron systems; closed and open shells; all supported spin representations and channels. Helium is a two-electron case despite having one occupied spatial orbital; it does not cover the one-electron requirement. |
| Frequency and equation variants | Static and frequency-dependent calculations, supported quadratic frequency combinations, and excited-state TDA/full-response branches. Include documented frequency/pole limitations and supported refusal paths. |
| Exchange and correlation | Every supported exchange/correlation selection and distinct implemented contribution: Coulomb, exact exchange, XC functionals/kernels and their derivative/correlation terms. Inventory supported HF, local, gradient-dependent, hybrid and other implemented variants with their build requirements; do not infer response support from an SCF-only functional. |
| Solvers and scientific components | Fixed-point/BSH proposal, metric arithmetic and conventions, KAIN and plain updates, step restrictions, spin/X/Y membership, normalization, projection, root matching/mixing, source construction and property contractions. Test individual terms as well as assembled results. |
| Dependencies and execution | Independent states and coupled roots; ground-reference eligibility; B/C to VBC to final-property ordering; single-process, MPI and supported subworld/state-parallel paths. |
| Inputs, seeds and persistence | File/CLI agreement, protocol entry, fresh/seeded/restarted/refined work, partial-state reuse, added states/roots, tightening, exact sample/checkpoint agreement, strict recovery, retention and read-only inspection. Seeds obey the explicit protocol. |
| Outputs and timings | Grouped properties and state metrics, unavailable values, source identities, analysis records, and stage/task/iteration/step timings. Exercise the shared functional interfaces and application-specific callbacks. |

**Accuracy and performance are both acceptance requirements.** Preserve scientific
reference values, provenance and justified numerical comparison tolerances. Check
that requested convergence measurements and published results agree. The evidence
must distinguish a MADNESS reproducibility reference from an independent external
calculation, and identify a paper/table only when that association is established.

Compare speed at comparable achieved numerical accuracy using controlled builds,
hardware, MPI ranks, thread budgets, protocol settings, seed/restart state and
checkpoint cadence. Record iterations/evaluations, total/stage/state/step time and
checkpoint I/O. Use repeated measurements to establish noise before selecting
performance bounds; do not invent a universal percentage or treat machine-specific
historical timings as portable limits. Investigate regressions in either accuracy
or cost before accepting the affected path. Unit tests for timing records use a
controlled clock; performance evidence comes from measured runs. Keep expensive
chemistry in its measured long/nightly tier.

## Initial numerical cases

| Stable ID | Existing case | Purpose | Current tier and recorded cost |
| --- | --- | --- | --- |
| ALPHA-HE-HF | `response_he_alpha` | Small SCF-to-response regression; static/dynamic alpha_zz against an archived MADNESS calculation. | medium; about 10 s, historical Seawulf run |
| ALPHA-H2O-LDA | `response_h2o_lda_alpha` | Molecular, multi-direction LDA response against independent DALTON alpha diagonals at two frequencies. | verylong; about 512 s on dd, historical run |

The existing [QC README](../../src/examples/qc/README.md) and
[CTest registrations](../../src/examples/qc/CMakeLists.txt) record the timing
provenance. These costs were not remeasured here. Keep helium as the quick check
and water as the slower validation check; selecting water does not move it into
short/medium CI. Other existing tests remain required where they cover other
science or behavior.

### ALPHA-HE-HF — static and dynamic helium response

**Intent:** preserve the SCF reference and two independently solved driven states,
then associate each alpha_zz value with its requested frequency. Detect numerical
and result-association regressions through a small complete workflow.

**Setup:** closed-shell He, HF, z perturbation, omega = 0 and 0.05 in atomic angular
frequency, protocol `[1e-4]`, k = 6. Current input uses DFT maxiter 20 and response
maxiter 30; legacy dconv/econv are not proposed new convergence defaults.

| Observable | Existing reference | Comparison |
| --- | --- | --- |
| SCF total energy (hartree) | -2.861628622851094 | absolute tolerance 1e-5 |
| alpha_zz(0), bohr^3 | 1.3226741885537008 | absolute tolerance 1e-3 |
| alpha_zz(0.05), bohr^3 | 1.3261590297119599 | absolute tolerance 1e-3 |

**Evidence and limits:** archived MADNESS run, commit `80e89f829`, 2026-09-21,
RelWithDebInfo/TBB, 1 rank × 7 threads. The alpha tolerance is a reproducibility
allowance documented in the check file, not a bound on the physical error.
This case does not test a DALTON seed, open shells or a full tensor.

**Current implementation:** [input](../../src/examples/qc/response_he_alpha/response_he_alpha.in),
[runner](../../src/examples/qc/response_he_alpha/run.sh),
[checks](../../src/examples/qc/response_he_alpha/check.json),
[reference](../../src/examples/qc/response_he_alpha/reference/response_he_alpha.calc_info.json).
CTest: `madness/test/qc/response_he_alpha/run`; use the existing qctest runner.

**Migration record — all 14 current checks:**

| Current checks (1-based in check.json) | Meaning to preserve | Shared-record destination/action |
| --- | --- | --- |
| 1 | SCF energy | Requested ground-state total energy with the same reference and tolerance. |
| 2–3 | Old achieved-threshold/dconv metadata, both stored as 1e-4 | Replace with current-request convergence and final-representation evidence. These are numeric metadata, not Boolean convergence tests; do not carry the retired dconv field forward. |
| 4 | Response application identity | Resolve the response calculation explicitly, not by task array position. |
| 5 | Final k = 6 | Runtime/checkpoint representation evidence. |
| 6, 8, 11 | Complete workflow and two converged response legs | Current request/dependency assessments and available published properties; saved measurements remain separate from the assessment. |
| 7, 10 | Frequencies 0 and 0.05 | Match frequency and component identity. Preserve the explicit zero-frequency assertion. |
| 9, 12 | The two alpha_zz values | Grouped polarizability values/components, preserving both numeric oracles. |
| 13–14 | Each state's legacy iteration count <= 9 | Preserve the cost guard's intent. Define its mapping to main updates versus evaluations when adopting the new accounting; review changed semantics before updating this bound. Reference iteration 6 is not an exact target. |

Links to acceptance cases: A01–A03, A06, B02 and C03/C04. Restart/reuse and
interruption variants must be added separately; this fresh-run case is not proof
of A04 or B01.

### ALPHA-H2O-LDA — DALTON molecular comparison

**Intent:** preserve static/dynamic molecular alpha diagonals with independent
numerical evidence, covering three perturbation directions and LDA response.

**Setup:** closed-shell water at the exact stored geometry/orientation, LDA
(Slater + VWN5), xyz perturbations, omega = 0 and 0.05; protocol `[1e-4, 1e-6]`,
final k = 8. Keep geometry, orientation, physical model and units unchanged.
The existing response-block protocol duplicates the DFT ladder; the adopted
interface keeps the ladder in its owning DFT group.

| Frequency (atomic angular frequency) | alpha_xx | alpha_yy | alpha_zz |
| --- | --- | --- | --- |
| 0 | 10.49982 | 10.27436 | 10.34337 |
| 0.05 | 10.67091 | 10.35192 | 10.4557 |

Values are in bohr^3. Preserve the current relative comparison
`abs(actual - reference) <= 1e-3 * abs(reference)` for these six values.

**Evidence and limits:** `reference_source` identifies DALTON master `984913d`,
SVWN5/t-aug-cc-pVQZ. The case records basis-set comparison history and explains
the tolerance's allowance for finite-basis and MRA errors. The checked-out
`dalton_reference` contains the `.dal` and `.mol` inputs; a raw DALTON output is
not present there in this checkout. A specific reference-paper citation/table
has not been established by this curation. Off-diagonal values in the reference
are placeholders and are not numerical oracles. The case also does not establish
excitation-energy accuracy or seeded-import correctness.

**Current implementation:** [input](../../src/examples/qc/response_h2o_lda_alpha/response_h2o_lda_alpha.in),
[runner](../../src/examples/qc/response_h2o_lda_alpha/run.sh),
[checks](../../src/examples/qc/response_h2o_lda_alpha/check.json),
[reference](../../src/examples/qc/response_h2o_lda_alpha/reference/response_h2o_lda_alpha.calc_info.json),
[DALTON settings](../../src/examples/qc/response_h2o_lda_alpha/dalton_reference/response_h2o_lda_alpha.dal),
[DALTON geometry](../../src/examples/qc/response_h2o_lda_alpha/dalton_reference/response_h2o_lda_alpha.mol).
CTest: `madness/test/qc/response_h2o_lda_alpha/run`.

**Migration record — all 11 current checks:**

| Current checks | Meaning to preserve | Shared-record destination/action |
| --- | --- | --- |
| 1 | Response application identity | Explicit calculation identity. |
| 2 | Final k = 8 | Runtime/checkpoint representation evidence. |
| 3 | Dynamic frequency 0.05 (tolerance 1e-12) | Identify the dynamic property by frequency; also bind the static oracle to omega = 0 rather than assuming row order. |
| 4, 8 | Both frequency results converged | Current assessments and exact contributing-state evidence; verify all required xyz states. |
| 5–7, 9–11 | Three diagonal values at each frequency | Match named xx/yy/zz components at each frequency with the unchanged DALTON numbers and relative limits. Do not compare placeholder off-diagonals. |

Links to acceptance cases: A01/A03/A06, B01/B02 and C03/C04. The fresh-run case
covers multiple independent states; partial-restart B01 behavior still needs its
own test. Grouped tensor schema checks cover all requested components separately
from the six independently referenced numerical values.

## Common migration and reference rules

Use [the shared tensor contract](../../tensor_result_examples.md): select the
response calculation, then `properties.polarizability`, exact frequency and named
component. Read the immutable reference using its original field paths. Do not
rewrite it into current output merely to satisfy a new reader. The planned shared
reader maps identities; it does not decide convergence or recompute science.

Replace old dconv/econv input controls with explicitly reviewed metric criteria;
their coupled legacy meanings make a numeric alias inappropriate. Preserve each
numerical comparison tolerance above. Keep final representation checks separate
from convergence tolerances and iteration/performance guards.

Never overwrite independent DALTON values with a MADNESS run. Preserve the
existing `reference_source`/`--update` guard in
[run_qctest.py](../../bin/run_qctest.py) and its tests. No references, decks or
assertions were changed during this curation.

Use staged validation: focused helper/controller tests, the small helium run,
then the water comparison and controlled restart/seed variants as the new layer
becomes available. The full matrix above remains required for complete response
acceptance. Numerical baseline refreshes need scientific justification; timing
baselines need measured noise and environment information.

## Existing scripted consumers retained for migration

These are inventory entries, not completed assertion curation. Preserve their
current numerical references and all scientific checks; the initial alpha pair
does not replace their beta, Raman or symmetry coverage.

| Existing consumer | Owning migration and remaining work |
| --- | --- |
| [Water HF alpha/beta z script](../../src/apps/madqc/test_molresponse_h2o_alpha_beta_z.py) and [deck](../../src/apps/madqc/test_molresponse_h2o_alpha_beta_z.in) | Step 7, with SCF-reference effects in step 6. Curate energy, component/frequency key sets and every alpha/beta assertion; preserve the absolute/relative comparison rule and replace old record readers/inputs. |
| [LiH alpha/Raman/beta script](../../src/apps/madqc/test_molresponse_lih_alpha_raman_beta_xyz.py) and [deck](../../src/apps/madqc/test_molresponse_lih_alpha_raman_beta_xyz.in) | Step 7, with SCF-reference effects in step 6. Preserve all three property groups and their existing references; map each assertion before changing readers. |
| [CIS water symmetry script](../../src/apps/madqc/test_cis_symmetry_h2o.py) | Step 8, with SCF-reference effects in step 6. Preserve both irreps, excitation energies and oscillator strengths. Replace its reference-only `restart read_only` use with solver reuse; the second scientific solve is not inspection. |

## INSPECT-NATIVE — saved-function inspection

**Status:** planned acceptance case; no new test or numerical run yet.

**Intent:** keep existing native inspection and plotting usable when result and
checkpoint paths change, including saved work that has not converged. This covers
C04 and the output migration in steps 4, 6 and 7.

**Fixture and assertions:** publish an unfinished native sample with shared
calculation/task/source records. Exercise the existing visualization manifest
writer and native loader for supported ground-state, driven and excited-root
samples as their engines adopt the layout. Resolve the original geometry,
representation and reference generation from shared evidence. Preserve available
metrics and the original completion status; readable functions do not require an
available final property. Check the matching results-snapshot link when present.
Verify no solver evaluation/update, new convergence assessment, protocol advance,
checkpoint publication or changed completion claim. Requested plot files may be
written; compare solver records, pointers and payloads before and after.

**Planned implementation:** none at present; `test_native_inspection.py` was
deleted when the inspection tools moved to node folders. Production consumers are
[VizManifest.hpp](../../src/madness/chem/VizManifest.hpp) and the native-loading
parts of [dump_mra_trees.cpp](../../src/apps/molresponse/tools/dump_mra_trees.cpp).
Use the shared records and native loader, with no legacy metadata authority.
Preserve existing supported capabilities and limitations. No shared HDF5 backend
or new exporter is required. Measure the test's tier when implemented; this
functional check does not replace numerical accuracy or performance evidence.

## OPT-ELECTRONIC — retained optimizer electronic records

**Status:** planned adapter/record acceptance case; no implementation or run yet.

**Intent:** keep each SCF/Nemo result attached to its exact geometry and physical
model while preserving the existing optimizer. This covers the retained electronic
handoff in steps 5–6, not deferred optimizer restart cases D01/D02.

**Cases:** use deterministic providers first, then the actual SCF/Nemo boundary.
Two geometries, including distinct coordinates with equal squared norms, get
separate evaluation IDs. Repeating an exact geometry/model reuses qualifying work
under the same non-retired ID. Missing metrics require evaluation or continuation;
an energy-only result does not satisfy a gradient request. Changed requests or
further work after a concluded attempt use a new electronic attempt and cannot
rewrite an earlier request, source map or result association. Changed geometry
can use previous orbitals as a supported guess, never as convergence evidence.

Apply [RESULT-SOURCES](#result-sources--immutable-result-attempt-binding) when
one evaluation has multiple electronic attempts. The source-map owner in each
result entry is independent of the optimizer's current selected attempt.

Fail selector publication before numerical work; interrupt before and after a
successful selection and after checkpoint publication. Resolve only the durable
electronic selector, with a selected null task distinct from a missing/corrupt
checkpoint. A subsequent invocation uses its supplied geometry and normal
optimizer entry, reusing/continuing compatible selected electronic work or using
it only as a supported guess for different geometry. No directory-order recovery,
older-checkpoint fallback, accepted-step restoration or Hessian/line-search
recovery is implied. Preserve selected, actively used and referenced archives;
retirement follows replacement selection and dependency checks.

**Implementation locations:** extend
[test_workflow_builders.cpp](../../src/apps/madqc/test_workflow_builders.cpp) with
stub identity/selection checks and the planned `test_input_execution.py` with
adopted SCF/Nemo coverage. The production boundary is
[Drivers.hpp](../../src/madness/chem/Drivers.hpp) and
[SCFTargetAdapter.hpp](../../src/madness/chem/SCFTargetAdapter.hpp), using shared
registration/publication. Preserve existing optimization numerical references;
measure the implemented test's cost tier.

## FD-PLAN — finite-difference plan publication

**Status:** planned D03 record/workflow tests; no numerical-gradient implementation
or numerical run is claimed.

**Intent:** resume the right displaced calculations, pair only eligible results
from the selected plan and preserve their archives during publication/retirement.

**Cases:** the single versioned `checkpoints/finite_difference_plan.json` owns the
plan revision, next evaluation number, base/formula/displacements/request and
point/source assignments. A completed-minus/interrupted-plus fixture restores
those identities and resumes only needed work. Changed displacement allocates
new IDs for changed geometries; explicitly unchanged points retain their IDs.
Tighter criteria reuse qualifying points and continue others. A derivative is
assembled only after all required child results qualify for the selected plan.

Inject failed writes/closes/replacement and interruption before/after atomic
replacement, including first-plan publication before inventory registration and
child publication before plan source binding. Check that stale revisions fail,
allocation never reuses committed numbers, and registration/child execution and
cleanup occur in the specified order. The run inventory's marker makes a missing
required plan an error; malformed/unsupported plans also fail without child-
directory reconstruction or alternate-file fallback. Plan-dependent cleanup
stops on these errors. Old exact sources remain protected until replacement
bindings and dependency checks permit retirement; published derivative snapshots
retain their original interpretation after plan replacement.

**Planned implementation:** `src/madness/chem/FiniteDifferencePlan.h` and
`test_finite_difference_plan.cc`; extend
[test_checkpoint_retention.cc](../../src/madness/chem/test_checkpoint_retention.cc)
and [test_workflow_builders.cpp](../../src/apps/madqc/test_workflow_builders.cpp).
Use the existing checked publication and unit-test registration mechanisms;
workflow stubs supply the point values. Measure the test tier when implemented.
Optimizer D01/D02 and new scientific gradient providers remain separate work.

## RESULT-SOURCES — immutable result-attempt binding

**Status:** planned synthetic serializer/reader and native-inspection cases; no
implementation or numerical run is claimed.

**Intent:** published results keep their exact original sources when run attempts,
electronic selections or numerical checkpoints change. Covers C03/C04, step-4
records/readers and OPT-ELECTRONIC/INSPECT-NATIVE integration.

**Cases:** root `results.json.attempt_id` identifies the run writing the snapshot.
A calculation's `attempt_id` identifies its immutable source-map owner; omission
inherits the root, and a differing owner must be explicit. Test inherited and
explicit bindings, mixed owners in one results collection, and an older entry
carried into a later run without changing its source-map owner. A source map may
refer to numerical samples from still older producing attempts; preserve those
attempt/sample/generation IDs rather than inheriting the map owner's attempt.

Publish snapshot A for one evaluation, then select electronic attempt B without
a result. Resolve/inspect A and verify its original map and state. Publish B's
checkpoint while A remains the current result: both payloads remain loadable.
Atomically publish B's result: historical A still resolves the same source map,
sample and generation, but its payload may be reclaimed if no active dependency
needs it. Historical maps or copied old results alone do not retain archives.
An inspection of reclaimed functions reports unavailable payloads; it never
substitutes B or launches a solve.

For C02, test the recursive closure of current-result sources, latest non-retired
task checkpoints and active calculation/plan dependencies. An active dependency
on A protects it after B replaces the result; releasing the last dependency
permits reclamation. Omitting a task from a request does not retire it. Explicit
retirement is persisted before deletion. Failed publication, unreadable live
records and interrupted replacement preserve existing archives; cleanup is
repeatable after restart. Storage follows live work rather than historical
attempt count. These tests add no optimizer history or historical-pin registry.

Reject invalid attempt values and missing, unsupported or mismatched source
maps. Readers never fall back to the root when an explicit override is invalid,
or use optimizer selectors, latest pointers or directory searches to guess a
published result's sources. Merely resolving/inspecting a snapshot causes no
solver work or checkpoint publication.

**Planned implementation:** `src/apps/madqc/test_shared_records.py` and
`shared_records.py`; extend
[test_results_summary.cc](../../src/madness/chem/test_results_summary.cc) and
[test_response_properties.cc](../../src/apps/madqc/test_response_properties.cc)
for serialization/association. Integrate with the existing
[test_workflow_builders.cpp](../../src/apps/madqc/test_workflow_builders.cpp).
Use synthetic records first, then the owning engine fixtures. Preserve all
scientific comparison tolerances and measure integration-test cost tiers.

## CORE-METRICS — measurement-based decisions

**Status:** first shared-core implementation, 2026-11-16; scalar contract tests.
No engine adoption or scientific tolerance calibration is implied.

**Intent:** every application can submit measurements for its own prepared state
and receive one explained decision. The controller must never infer a zero from
missing evidence, skip a required root, relax final limits or request an extra
update after convergence.

**Setup and assertions:** declared model support/defaults and plain sample records
exercise absolute/relative criteria, relative-only overrides, intermediate
heuristic floors, negative energy scales, zero or missing scales, missing or
self comparisons, mismatched contexts/representations and coupled membership.
Application-owned loop fixtures count proposals, updates and terminal evaluations;
endpoint and optional interval publication are requests, not archive writes.

**Evidence and limits:** hand-checked scalar fixtures validate control flow and
arithmetic. They do not establish model-specific metric definitions, function
norms, MPI behavior, archive integrity or calibrated scientific accuracy.

**Acceptance:** inclusive limits, every selected applicable metric on every
required subject, exact sample binding, final-stage requirement, fresh update
budgets and no extra update after passing. Invalid or unavailable required
evidence cannot produce convergence. Reusing qualifying saved evidence returns
completion with its original identity and no publication request.

**Current implementation:**
[test_iteration_metrics.cc](../../src/madness/chem/test_iteration_metrics.cc),
CTest `madness/test/chem/test_iteration_metrics/run`. Uses
[IterationMetrics.h](../../src/madness/chem/IterationMetrics.h),
[ExplicitCriteria.h](../../src/madness/chem/ExplicitCriteria.h) and
[ConvergenceSettings.h](../../src/madness/chem/ConvergenceSettings.h).

**Migration record:** additive contract coverage. Existing legacy policy tests
remain unchanged; applications migrate with their owning engine phase.

**SCF endpoint extension — approved 2026-11-18; implementation pending.**
Extend A02/A05 coverage through the real SCF adapter: a pre-rotation pass cannot
certify a post-rotation failure. Test continued work with a fresh residual and
KAIN reset after a basis change, final-cap failure without another update, and
post-rotation success without another rotation or update. Verify the original
preceding-iteration baseline, unavailable comparisons, exact archive/metric
agreement and unchanged update count across the terminal evaluation. Exercise
supported SCF model/spin and occupied/virtual paths, including intermediate-cap
advancement. Existing scalar controller tests do not validate orbital rotations;
retain their scope and add the engine cases during SCF adoption.

**Protocol-transition extension — approved 2026-11-18; implementation pending.**
Exercise the real SCF projection/operator path for automatic `1e-4, k=6` to
`1e-6, k=8` and explicit fixed-k threshold changes. Verify carried orbitals,
rebuilt operators and fresh equation work, empty KAIN, cleared comparisons and
unavailable first change metrics. Subsequent comparisons must use the new
protocol's baseline. Preserve historical evidence, keep smoothing independent
and count no solver update for preparation. Bind these cases to A03/A05/A06;
scalar planner tests alone do not establish numerical conversion correctness.

## CORE-RESTART — independent progress and coupled restart

**Status:** first shared-core implementation, 2026-11-16; pure planner tests.

**Intent:** preserve each compatible task's progress without inheriting another
state's protocol position, old budget consumption or convergence label.

**Setup and assertions:** saved samples at completed/interrupted intermediate and
final protocols; tighter satisfied/unmet requests; changed numerical settings;
appended stages; removed protocols; same threshold with a different k; coupled
root expansion; new states; invalid selected records and unsupported conversion.

**Acceptance:** matching stages resume or advance appropriately, newly requested
states and converted removed protocols start first, every resumed stage has a
fresh allowance, and a saved final state is reused only through the shared
measurement assessment. Invalid selected records fail without fallback. Added
roots renew the coupled solve at the saved stage.

**Evidence and limits:** synthetic, already-loaded records exercise planning.
Native reading, conversion, saved/current calculation-settings comparison, root
construction and KAIN remain application/adapter responsibilities. The current
planner accepts `SavedEvidenceUse`; it does not perform the settings comparison.

**Approved adapter coverage — 2026-11-18; not yet implemented:** derive reuse or
fresh checks from explicit settings comparisons. Unchanged settings with passing
final metrics reuse; unmet/missing metrics continue with preserved progress.
Criteria/tolerance/maxiter changes alone must not invalidate scientific context.
Changed geometry, model/reference, actual frequency or operator/smoothing
settings require fresh measurements from supported saved functions. Cover
representation conversion and metric-definition changes separately. Unsupported
continuation and invalid selected archives fail without substitution. Add no
inferred equivalence or nearby-frequency search. These cases extend A04 without
claiming that the existing synthetic planner tests validate scientific adapters.

**Current implementation:**
[test_iteration_restart.cc](../../src/madness/chem/test_iteration_restart.cc),
CTest `madness/test/chem/test_iteration_restart/run`, using the additive shared
API in [LadderPlanner.h](../../src/madness/chem/LadderPlanner.h).

**Migration record:** existing `test_plan_ladder` assertions remain. The common
protocol value type is reused without changing the legacy matching policy.


## CORE-CHECKPOINTS — native publication and live-source retention

**Pending extension, 2026-11-17:** the current store retires whole calculations.
Individual task retirement and response `replaces` records are approved but not
implemented. Extend `test_checkpoint_store.cc` and `test_checkpoint_retention.cc`
with the replacement/retirement cases under RESPONSE-2PA-CONTINUATION; the
existing passing tests do not establish this new behavior.

**Status:** shared storage implementation, 2026-11-16; no engine activation.

**Intent:** a saved numerical sample, its measurements and its exact dependencies
stay together. Publishing new progress must preserve sources still used by
current results or active tasks. Heavy storage follows live work rather than
iteration history, while small records retain the original meaning of results.

**Setup:** temporary run directories with real native parallel archives of small
integer payloads. Use iterative samples and a direct VBC-shaped task, immutable
source maps, mixed attempt ownership and a workflow-owned plan fixture. Run in
serial and with two MPI ranks. Inject write errors, pointer/result replacement
failures, missing records, corruption and interrupted deletion.

**Assertions:** registration precedes writing; null means never published, while
missing or invalid selected data fails. A newer checkpoint does not rebind an old
result. Keep A while result A or a live dependency uses it; release A only after
result replacement and the last dependency ends. Preserve omitted tasks and
require explicit valid retirement. Historical maps alone do not pin payloads.
Keep manifests and reserve generation identities after reclamation. Repeated
cleanup finishes abandoned writes/deletions, without deleting on an incomplete
live-reference check. Missing results disables cleanup while numerical restart
and further saves remain possible. Reject unpublished manifests as exact sources.
Full checksum validation prevents a corrupt replacement from releasing the last
valid result source; keep historical publication receipts with the metadata.
A required workflow plan must have a validated root provider. Archive callbacks and failures follow one branch on every rank.

**Evidence and limits:** native serialization, filesystem failure injection and
pure graph cases test storage/control behavior. They do not validate chemistry,
FD formulas, scientific compatibility, process-loss recovery or performance at
production archive sizes. The workflow plan is a test fixture; its scientific
schema and writer belong to the later workflow adapter. Existing native tests
retain shard/cookie/checksum validation and reader process-count coverage.

**Current implementation:**
[test_checkpoint_store.cc](../../src/madness/chem/test_checkpoint_store.cc),
CTest `madness/test/chem/test_checkpoint_store/run` and
`madness/test/chem/test_checkpoint_store_mpi2/run` (measured below one second
each in the 2026-11-16 Release build; `short` tier);
[test_checkpoint_retention.cc](../../src/madness/chem/test_checkpoint_retention.cc)
adds transitive closure, release, cycle and missing-source cases.
[test_stage_checkpoint.cc](../../src/madness/chem/test_stage_checkpoint.cc)
continues exercising the underlying native primitives.

**Migration record:** the version-2 coordinator is additive. Existing engines
use their existing publication path until their owning migrations. Full property
formatting, application inspection and scientific reference comparisons remain
separate acceptance work; no reference values change for this storage increment.


### Iteration identity clarification — 2026-11-17

CORE-METRICS, CORE-RESTART and CORE-CHECKPOINTS now use `IterationInfo`,
`MetricValue`, `StateMetrics` and `IterationMetrics`. The identity is task,
attempt, protocol and completed-main-update number, scoped to the calculation;
zero is a valid initial/restarted state. Native records store that number once
as `updates_this_attempt`. Comparisons and source maps carry protocol and
iteration explicitly. Reuse preserves the producing identity.

Coverage distinguishes equal counters in different protocols/attempts, rejects
self/future/cross-task comparisons and stale live metrics, and validates stored
comparisons even when a metric is unavailable. It round-trips restart iteration
zero with its predecessor and rejects source maps or decisions bound to the
wrong iteration/protocol. Direct VBC records retain completion and exact inputs,
with no fabricated iteration or metrics. This clarifies identity without changing
scientific values, requested tolerances or checkpoint retention roots.


## RESPONSE-DIRECT-INPUTS — qualified sources for direct computations

**Status:** approved coverage, 2026-11-18; engine migration pending. Maps to B02
and complements RESULT-SOURCES and CORE-CHECKPOINTS.

**Intent:** VBC and property computations use the exact states whose measured
evidence satisfies the current final request. Molresponse keeps its dependency
manager; the shared layer supplies assessment, source binding and storage.

**Setup and assertions:** include an eligible final input alongside an available
coarse/partial guess, an unmet capped input and an old completion flag whose
metrics fail the current request. Only qualified sources release direct work.
Load their exact generations; an unusable bound archive fails without guess or
newer-generation substitution. VBC records B/C/reference inputs and completion,
without fabricated iteration metrics; properties bind every required source.
Changing only A preserves valid VBC; changing B or C rebuilds affected VBC and
properties. Tighter limits already satisfied by unchanged sources cause no
rebuild. Preserve original result sources after `latest` advances, and protect
bound inputs through publication using the existing retention rules.

**Implementation targets:** adapt `calc_manager.hpp`, `calc_executor.hpp`,
`vbc_save_load.hpp`, `fd_save_load.hpp` and `response_metadata.hpp`; extend
`test_calc_manager.cpp`, `test_response_metadata.cpp`,
`test_response_property_planner.cpp`, `test_fd_kfallback_load.cpp` and the migrated
native response-checkpoint tests. Keep solver guess-loading coverage separate
from direct-consumer eligibility. Existing scalar/store tests do not prove this
scientific boundary; preserve curated numerical references during adoption.

## RESPONSE-2PA-CONTINUATION — refined photon frequencies

**Status:** author-approved design extension, 2026-11-17; implementation and
engine acceptance evidence pending. Extends B01–B03 and C01–C03.

**Intent:** preserve useful response functions and their protocol after refinement
changes a tracked excited root's photon frequency, while ensuring convergence
and published properties describe the updated equation. This is continuation of
a matched task, not a general nearest-frequency seed search.

**Setup and assertions:**

- Same tracked root/perturbation, unchanged ground reference and other relevant
  settings, shifted photon frequency: start from the saved functions at their
  requested saved `(threshold, k)`, with fresh equation/metric/comparison work,
  empty KAIN and a fresh allowance. Old passing metrics cannot authorize reuse.
- Repeat a formerly completed intermediate protocol at the new frequency; do
  not skip ahead on its old completion flag. Other states retain their progress.
- Exact unchanged frequency/context: qualifying saved measurements still permit
  reuse without numerical evaluation, updates or archive writes.
- Reordered roots keep their tracked associations. Ambiguous/unmatched roots and
  genuinely added frequencies cannot inherit another task's protocol progress.
- Removed saved protocols use supported conversion at the first requested entry;
  changed k is not an exact match. Invalid selected checkpoints still fail
  strictly. No new frequency-distance input or frequency-rounding match is used.
- Old task/frequency/attempt/metric records and published result sources remain
  immutable. Updated source computations and tensors bind the actual new root
  and response generations. Test interruption before/after replacement and
  retirement; preserve payloads still required by current results or active
  continuation and reclaim superseded ones only through the existing rules.
- Preserve current frequency-based task IDs and record exact `replaces` sources.
  While the successor has no usable checkpoint, its pending continuation input
  protects the predecessor. After that need ends, replacement history alone
  must not keep the old payload alive.
- A still-requested old-frequency state, an active direct consumer, or a current
  result blocks retirement. Once those uses end, persist retirement of only the
  obsolete task and allow cleanup; other tasks in the calculation stay active.
- Tiny frequency changes can share the current rounded task key. Check the actual
  scientific frequency, publish a new generation with fresh evidence, and never
  retire the continuing task. Keep every older result bound to its original
  generation. Repeating registration/retirement after interruption is safe.

**Evidence and limits:** migrate the existing water 2PA cases and
`test_named_two_photon.py` to shared-record assertions; add focused planner and
DAG checks in `test_iteration_restart.cc`, `test_calc_manager.cpp` and
`test_calc_manager_run.cpp`. Compare continued and fresh solutions at identical
final settings using justified numerical tolerances; preserve existing references.
Record actual iterations, timings and achieved metrics before claiming a speed
benefit. Existing scalar planner tests do not yet establish this engine behavior.

## WORKFLOW-COMPLETION — unfinished producers do not release consumers

**Status:** verified boundary and transitional SCF/Nemo adapters, including
serial and two-rank numerical restart tests. This does not certify adoption of
the shared numerical metric or results.json contracts. Evidence is recorded in
`.superpowers/sdd/design_draft/workflow-completion/verification.json`.

**Intent:** a normal return from an electronic application must not be mistaken
for convergence. Preserve available results and valid restart work on a controlled
stop, while distinguishing numerical or publication failure. This covers the
workflow prerequisite in B02 and controlled-stop reporting in C02/C03.

**Cases:** the real workflow is exercised with small test applications for final
completion, iteration-limit stop and numerical failure. An incomplete producer
must publish no reference into StepContext and its consumer must not execute.
The incomplete summary retains available properties. A report replacement
failure must reach every MPI rank and preserve the previous report. Conversion
from the shared Decision accepts complete/stop and rejects update/finish_protocol.

The helium HF response deck supplies the numerical case: one SCF iteration
stops the workflow before response, leaves native orbitals and an incomplete
summary, and a later invocation with its ordinary allowance completes SCF and
executes response. A forced SCF report-write failure preserves the prior report
and prevents downstream execution. Run both serially and with two ranks.

**Checks:** [boundary test](../../src/apps/madqc/test_workflow_completion.cc),
[numerical test](../../src/apps/madqc/test_workflow_completion.py), and the unchanged
[helium alpha reference](../../src/examples/qc/response_he_alpha/check.json).
The boundary unit cases are short. Same-instance Nemo and the serial/two-rank
workflow cases are medium (11.26, 20.24 and 24.45 seconds on dd with seven
MADNESS threads and sequential BLAS); requalification cases are long (51.74
seconds serial and 61.89 seconds with two ranks).
These measurements select test tiers, not a solver performance baseline.

**Migration boundary:** SCF/Nemo currently bridge their existing engine outcomes;
moldft and Nemo record their actual solve or native-reuse outcome. Cached reports
are requalified through the native planner. Incomplete frozen-virtual work retains
the occupied source threshold while invalid dconv prevents combined completion;
a retry must finish successfully. No native archive fields or stopping criteria
change. Applications without an outcome keep
their existing path until adoption; absence is not shared convergence evidence.
Future shared metric and native publication adapters must retain these behavioral
assertions while replacing the legacy report/archive paths. No numerical
reference is regenerated for this behavior change.

The companion [requalification test](../../src/apps/madqc/test_workflow_requalification.py)
covers failed forced Nemo iteration followed by auto restart, changed-k SCF
reuse, changed-Hamiltonian read-only inspection followed by an ordinary run,
report-name collision and added virtuals saved incomplete then resumed with
frozen occupied orbitals. The numerical `--nemo-repeat` mode of the C++ boundary
test checks retrying the same application instance at the same geometry. Invalid
precision stamps and rank-zero JSON serialization failures are unit cases.

## SCF-EVALUATION — metrics belong to the evaluated and saved orbitals

**Intent:** separate orbital evaluation from advancement without changing the
HF/DFT equations. A stopping decision must describe the current functions; KAIN
must consume the residual already computed for those functions. Final orbital
preparation must be followed by a new evaluation before deciding or saving.

**Numerical coverage:** `test_scf_evaluation` checks canonical helium HF,
localized helium LDA, unrestricted lithium HF, canonical and localized water HF,
and immediate residual-only acceptance. Each observer sample is checked against
a fresh evaluation; closed-shell HF also uses the preexisting independent
evaluator. Reloading the native archive must reproduce final energy and residual.
Changing the supplied cached residual verifies that advancement uses that work.
The localized water DOS file is checked against the generalized Fock spectrum,
which differs from the localized orbital expectation values.

**Decision and publication coverage:** `test_scf_evaluation_policy` forces a
preliminary pass followed by failure after final preparation, both with an update
remaining and at the cap. `test_scf_final_evidence` exercises an engine pass with
an unmet explicit energy request, requiring final publication to retain the
exact deciding measurements and engine verdict. `test_named_stages` and
`test_named_density_changes` retain historical phase labels when loading existing
checkpoint histories. The first evaluation has a residual but no fabricated
zero density, energy or kinetic change.

**Accounting:** protocol iteration numbers count completed main updates; zero
is the initial evaluation. Endpoint evaluation consumes no update. The retained
legacy `scf_e_data.iterations` counts energy samples, including initial and terminal
samples; it is not interchangeable with the new update count. Its output-field
migration remains separate.

Measured on dd with the existing Ninja build: the numerical case takes 112 s
with one rank and seven threads (`long`), and 202 s with two ranks and two
threads each (`verylong`). The scalar endpoint policy stays `short`.

These tests accompany the existing workflow completion/requalification cases and
unchanged numerical references. They do not establish release-wide XC coverage,
Nemo adoption, or a performance bound. Test costs and final evidence are recorded
in `.superpowers/sdd/design_draft/scf-evaluation/verification.json`: 35 targeted
tests pass across the final runs, including the existing numerical references.
A prior changed-k test assumed that one allowed update could never converge;
the new terminal assessment proves it can. Its stale-stamp case now uses zero
updates. Both MPI configurations pass; the original failures remain recorded.

## SCF-SHARED-METRICS — one record for stopping and native restart evidence

**Intent:** validate ordinary SCF's shared measurement/controller adapter before
choosing numerical defaults or activating the new public inputs. A fresh run,
a tightened request and a protocol change must use the same metric definitions.

`test_metric_helpers` checks raw maximum residuals, component identities,
physical scalar changes, maximum field changes, relative scales, absent history
and invalid numerical evidence. It uses scalar fixtures without a chemistry
solve. These definitions survive changes in JSON layout.

`test_scf_iteration` drives the real shared-controller overload for helium HF,
unrestricted lithium HF and helium LDA. Independent reevaluation checks the
physical energy and kinetic quantities, both spin densities and residual blocks.
It checks intermediate threshold floors, strict final limits, fresh restart
allowances and history, fixed-k protocol changes, zero-update residual-only
completion, endpoint publication, native orbital reload and reassessment of
saved metrics against new tolerances. No reference value is regenerated.

The existing SCF-EVALUATION cases continue to cover the underlying equation
split, cached BSH advancement, endpoint reassessment and saved orbital rotations.
This adapter test is not evidence of public CLI migration, calibrated absolute
defaults, full workflow restart adoption, frozen-only/Nemo adoption or complete
XC/performance coverage. Executed checks and timings belong in the increment's
`scf-shared-adapter/verification.json` under the implementation ledger.

Final adapter validation: 30 targeted CTests pass, plus the standalone helper
compile/run and focused protocol-rejection regression. The adapter case takes
71.36 s on one rank/seven threads (`long`) and 136.77 s on two ranks/two
threads each (`verylong`). The single review finding, rejecting absent or empty
requested protocol sequences before evaluation, has failing/passing evidence.
Focused clang-tidy reports no diagnostics with the documented crashing checker
excluded. Native formats and numerical references are unchanged.

## SCF-DEFAULTS — scientific defaults and coarse-protocol limits

**Purpose:** keep chosen tolerances tied to actual prepared-state metrics and
unchanged scientific references, independently of output-format changes.
`src/madness/chem/test_scf_defaults.cc` exercises He/HF, He/PBE0, water/HF and
open-shell Li/TPSS using the existing QC decks. The scientific oracle is total
energy in hartree with the original `1e-5` absolute reference tolerance; the
existing qctests retain their other assertions. The default criteria are BSH,
per-spin density and total-energy change, with kinetic tested as an explicit
addition. The coarse Li case must stop with an unmet residual; refined Li must
complete, without comparing its different representation to a coarse oracle.

CTest entries `test_scf_defaults{,_hybrid,_water,_coarse,_refined,_kinetic,_invalid_stage}` keep
these contracts. The calibration also compares ten-times-tighter requests and
records metrics/update counts/times in the
[metric evidence](madqc-metric-validation.md#ordinary-scf--initial-default-validation-2026-11-18).
The hybrid/TPSS variants register only when Libxc is available. An injected
invalid-state controller decision verifies that a failed protocol cannot be
hidden by a later successful solve; intermediate allowance exhaustion still
advances normally. This is ordinary SCF evidence; frozen-only virtuals, other functionals and optional
models still need their migration coverage. No performance-regression bound is
inferred from these timings.

## INPUT-OCCURRENCES — validate before applying overrides

`src/madness/mra/test_cli_occurrences.cc` uses real repeated argv inputs. Its
serial and two-rank tests require invalid earlier values to fail even when
replaced, whole-group replacement without merging, preserved path case, and
identical error category/message on every rank. A private registration copy
prevents discarded values or partial failed reads from changing live settings.
Per-value hook cases cover a retired-key fixture and invalid bounds; original
arguments survive programmatic transformations without manufactured defaults.
A rank-zero-only parser tests collective participation, including explicit
missing-file and valid-file selection through the combined reader. A valid file
followed by a rejected override must roll back the whole call. The file case requires
the same validation and unchanged settings on failure. Owning application
retirement and standalone input migration remain separate coverage obligations.

`test_optimizer_settings` verifies independent `gtol=1e-4` and
`gradient_precision=1e-5`, retained protocol-derived energy/step defaults and
explicit-override precedence. Optimizer replacement and history restart remain
deferred; this test makes no new optimization-accuracy claim.

## SCF-FROZEN-METRICS — virtuals in a fixed occupied field

`src/madness/chem/test_scf_frozen_iteration.cc` checks the programmatic shared
virtual solver on restricted He/HF and He/LDA and unrestricted Li/HF with
unequal alpha/beta virtual counts. It first obtains an occupied solution under
explicit criteria, publishes it, and uses the existing native loader to add
virtuals. The virtual checkpoint names that exact occupied source.

The numerical assertions compare every recorded residual and virtual eigenvalue
with a fresh evaluation of the actual orbitals, and check that occupied functions,
occupied eigenvalues and source energy stay unchanged. Density, total-energy and
physical kinetic changes must be inapplicable, never measured zeros. Tests cover
inapplicable-only requests, invalid controls, incomplete spin membership, fresh
restart allowance, zero/one-update caps, intermediate advancement, endpoint
completion without KAIN and native metrics round-trip. Serial and two-rank runs
exercise the same contract.

Four focused review regressions supplement those ordinary numerical cases:
`test_scf_frozen_frozen` preserves a deliberately off-unit occupied function
through native padding and a KAIN update; `test_scf_frozen_empty` and
`test_scf_frozen_occupations` require unusable bundles to fail before changing
the previous outcome or convergence metadata; `test_scf_frozen_lifetime` uses
a weak reference at the existing beta-residual print boundary to verify that
the legacy solver has released the completed alpha KAIN history. The lifetime
check runs on one rank and adds no numerical hook to the solver. The off-unit
fixture tests preservation, not scientific eligibility of its source.

Measured on dd with sequential BLAS: the main numerical test takes about
58 seconds on one rank and 103 seconds on two ranks (`long`). The four focused
cases take 1.6–5.6 seconds each (`short`). These are harness cost tiers, not
solver performance-regression bounds.

These tests establish measurement and iteration behavior. Deliberately strict
caps and a measured loose completion request are not a physical virtual-spectrum
reference or a calibration of default virtual tolerances. Bound-virtual spectra,
broader functionals/corrections, exact workflow reference qualification and the
virtual-count preparation schedule retain their migration obligations.

## SCF-NATIVE-RESTART — measured evidence and the saved orbital bundle

`src/madness/chem/test_scf_checkpoint.cc` joins the real helium HF solver,
`SCFCheckpoint`, the shared restart planner and `CheckpointStore`. It verifies
that saved metrics can satisfy a tighter request without another evaluation,
while unmet requests receive fresh comparison history and a fresh update budget.
Before reevaluation it compares the reloaded orbitals, eigenvalues, occupations
and energy with the saved state. It also exercises final-protocol requirements,
changed scientific settings as guesses, unsupported cell conversion, a frozen
virtual task's exact occupied source and strict corruption rejection. The same
case runs on one and two MPI ranks.

Review regressions require current source qualification to override a historical
incomplete outcome, reject plans from another context or occupied source before
loading, and reject core-potential trimming before the legacy reader can remove
active occupied functions. The last case publishes a real helium active bundle
and uses a Be MCP model only to test conversion preflight; it makes no Be energy
or core-potential accuracy claim. Reconstructing the binding in one attempt must
also allocate a new immutable generation. The MPI case uses multiple I/O writers.
On dd the extended case takes about 23 seconds on one rank (`medium`) and
45 seconds on two ranks (`long`); the MPI harness reserves the full processor
allocation to avoid overlap with other numerical tests.
These numbers classify test cost, not solver performance.

This is a restart integration test, not a new numerical reference or a virtual
spectrum calibration. The binding captures scientific settings directly from SCF.
Additional cases reject changed geometry even when coordinate norms agree,
preserve reuse when only tolerances/budgets/reporting change, and reject a
rank-local settings change collectively before native loading. Public input migration,
virtual preparation scheduling and the legacy geometry-cache replacement remain
separate acceptance obligations.

Scientific-settings review regressions also cover native hole-to-default and
shortened-prefix occupation requests, checking the occupations during fresh
evaluation. Frozen occupation changes require a new qualified source. Two-rank
constructor cases reject differing orbital counts or frozen mode before branching.
An active PCM request without a matching prepared solver is rejected at binding.

## SCF-SCIENTIFIC-SETTINGS — what saved measurements describe

`src/madness/chem/test_scf_settings.cc` checks the scientific meaning of a
snapshot independently of its serialization layout. Exact atom positions,
effective electric fields, smoothing overrides, core/GTH coefficients, occupation
requests, derivative conventions, actual XC screening and active PCM settings
must distinguish calculations. Different input provenance with the same effective
geometry does not. Tolerances, update budgets, output and seed choices do not
change scientific identity; the shared planner checks protocol representation
separately. Unclassified new parameters participate until explicitly reviewed.

The fixtures exercise fields missing from the legacy molecule JSON and real
loaded MCP/GTH coefficient data. Unexposed external PCM reader/cavity data and
unloaded GTH data cannot qualify saved evidence. These are metadata/restart
checks, not numerical reference calculations for PCM, core potentials or
dispersion. Public workflow coverage must still establish correct operator
preparation and complete capture of supported external data before release.
When PCMSolver is enabled, the fixture also constructs a real host-configured
solver, changes its dielectric, and requires reconstruction before reuse. That
conditional path, and active DFTD3 construction, are not exercised by the current
`simple-policy` build, which disables both optional libraries.
The settings-only fixture takes about one second on dd (`short`).

## Public SCF adoption — 2026-11-19

Release amendment: madqc is the sole application executable for supported
engines. Standalone moldft checks below are transitional implementation evidence;
migrate their scientific assertions to madqc as that duplicate target is removed.
Named-stage assertions likewise migrate before retiring the old runner. Separate
application executables remain only for engines outside madqc support.

| Stable case | Purpose and oracle | Executable evidence |
| --- | --- | --- |
| SCF-PUBLIC-INPUT | Correct owner, explicit names/defaults, relative-only override, original CLI occurrences and invalid bounds. No old alias can silently affect a shared solve. | `test_scf_inputs`; public CLI/JSON failure cases in `test_shared_scf` |
| SCF-PUBLIC-RESTART | Madqc and the temporary moldft entry: capped endpoint, saved-metric reuse, tighter requests, same-norm changed coordinates, final k conversion and strict corruption. | `test_shared_scf` serial/two-rank; `test_scf_run` |
| SCF-INSPECT | Load saved native functions for density and Coulomb plotting; checkpoint hashes, run inventory and result records remain unchanged. New input settings do not certify the saved state. | `test_shared_scf`; `test_scf_recovery` |
| SCF-PREPARATION | Qualify occupied sources before frozen virtuals; cap preparation separately; resume after main-publication failure, including when only preparation exists; changed schedules/protocols require fresh checks; retire/reclaim preparation archives. | `test_shared_scf`; `test_scf_recovery`, serial/two-rank |
| SCF-SEED-OCCUPATIONS | Raw native functions are guesses. A removed source occupation hole must not survive the target's default occupations. | `test_scf_run` |
| SCF-RESPONSE-SOURCE | Qualified native archive carries verified charge/spin; virtual bundle size cannot become an electron count. Include charged open shell, a closed-shell bundle with virtuals, and rejection of unsupported beta core trimming. | `test_shared_reference` serial/two-rank; helium alpha workflow |
| SCF-EXACT-SOURCE | Publishing a newer task generation cannot redirect a consumer away from its already qualified, retained source. | `test_scf_run` |
| SCF-FROZEN-RECOVERY | A capped virtual bundle stays incomplete; a later solve must perform updates, publish a new generation and satisfy its explicit residual request. The coarse helium fixture uses 1e-2 for this lifecycle check; it does not calibrate physical virtual-spectrum accuracy or defaults. | `test_workflow_requalification`, serial/two-rank |

Measured on dd during integration: public restart/inspection coverage takes
about 116 s in serial and 202 s with two ranks. The focused recovery cases take
about 227 s and 150 s respectively. Recovery and the two-rank public test are
therefore `verylong`; they do not enter the short/medium CI budget. These
integration timings include build contention and are not performance benchmarks.

The He HF, He PBE0, water HF, Li TPSS and He alpha decks use the shared SCF
input names. Their numerical reference values and numerical tolerances are
unchanged. Legacy threshold-stamp comparisons are replaced by shared
controller/outcome/source tests, not relabeled as new metric measurements.
Li TPSS explicitly requests a 5e-4 BSH residual at its coarse final protocol;
the separate default-calibration case still requires the stricter shared
default and expects its documented cap. The response solver's own `dconv`
remains until that engine is adopted. Other historical decks still need explicit
input/schema migration; this focused set is not full release coverage.

## Full response milestone: alpha comparison sources

The initial integration checks use deliberately loose gates to isolate native
loading, task progress and output wiring. They do not establish scientific
accuracy. At the final numerical check, compare grouped
`response.properties.polarizability.values` by frequency and named Cartesian
component with these unchanged sources:

| Case | Preserved reference | Comparison |
| --- | --- | --- |
| `response_he_alpha` | Existing HF helium reference | Static/dynamic alpha_zz, existing absolute limit in `check.json`; record convergence and timings without requiring the old iteration count. |
| `response_he_lda_alpha` | DALTON t-aug-cc-pVQZ values in its reference directory | Static/dynamic alpha_zz; existing relative limit `1e-3`. |
| `response_he_pbe0_alpha` | DALTON t-aug-cc-pVQZ values in its reference directory | Static/dynamic alpha_zz; existing relative limit `1e-3`. |
| `response_h2o_lda_alpha` | Curated DALTON water references described above | Named xx/yy/zz at each stored frequency; preserve existing limits and geometry/orientation. |

Closed-shell linear response uses the existing LDA/GGA XC kernel through
`XCOperator::apply_xc_kernel`; hybrid exchange uses its existing exact-exchange
coefficient. Open-shell non-HF response is explicitly rejected by
`build_response_ground_state_open_shell`. Meta-GGA SCF support does not establish
meta-GGA response support: `XCfunctional::fxc_apply` explicitly rejects
functionals needing tau because the response kernel lacks perturbed kinetic
energy density. This remains unsupported, rather than an untested supported
variant. Quadratic DFT and open-shell
property contractions retain their documented limitations. These are inventory
statements, not newly implemented model support. External DALTON references must
never be replaced with values produced by this migration.


## RESPONSE-SHARED-INTEGRATION — public records and retained scientific oracles

`test_shared_response_numerics.py` runs a curated response deck and maps the
physical checks in its existing `check.json` to shared results. The original
reference file, absolute/relative tolerance and units retain their meaning.
The script records update counts and timing locations without making historical
iteration counts an integration gate. Original qctest format/cost checks remain
visible for later migration/calibration; this script does not regenerate data.

Mappings: SCF energy -> `total_energy.value`; alpha tensor entries -> frequency
and named Cartesian components; excited omega -> stable root energy; beta ->
frequency pair and Cartesian component; Raman derivative -> coordinate/frequency
and component; mode/spectrum quantities -> numeric arrays and spectrum rows;
2PA observables -> root values. DALTON alpha data remain external oracles.
The first representative batch is helium HF/LDA/PBE0, water LDA, H2 TDA/RPA
(including no-KAIN), and LiH beta. Full model/functional coverage, all reference
recomputations, defaults calibration and detailed cost studies remain later work.

Focused lifecycle tests have separate purposes:

| Test | Scientific or functional intent | Evidence / limits |
| --- | --- | --- |
| `test_shared_response_defaults.py` / `--families` | Ordinary alpha input converges with the default metric list; a mixed alpha/beta/Raman/excited/2PA request also uses shared records without selecting metrics explicitly. | The mixed case uses loose explicit absolute tolerances to isolate routing. Default accuracy remains provisional; no universal calibration claim. |
| `test_shared_response_seed.py` | Existing DALTON functions enter as guesses, receive new metrics, then reuse native records. | Loose BSH gate and zero updates intentionally isolate import/restart; not a DALTON accuracy check. |
| `test_shared_derived_response.py --full-raman` | Full coordinate selection reaches Hessian/modes/spectrum output and exact native reuse. | Numeric array contract; loose integration gates, no independent full-Raman oracle. Partial-coordinate mode is a separate test. |
| `test_shared_quadratic_response.py --ranks 2 --subworlds 2` | Independent Worlds publish/reuse states and exact VBC dependencies without collective or writer conflicts. | Includes participating-rank timing scope; numerical comparisons are separate. |
| `test_shared_es_iteration` | Actual native root functions agree with saved omega, density norms and positive kinetic sums; residual reduction/scales, update ownership and missing history are consistent. | TDA/full, batched/streamed, serial/two-rank fixtures; root-matching boundary cases in `test_es_comparison`. |

All these tests preserve existing solver kernels. Truncated ES analysis explicitly
marks uncomputed properties unavailable; it must not report default-initialized
zeros as measurements. Native inspection and seed conversion retain their own
public purposes: import tools produce guesses, while shared checkpoints own
ordinary solver restart.

### Response review boundary cases

| Test | Why it exists | Contract |
| --- | --- | --- |
| `test_shared_response_review` (one and two ranks) | Two physical frequencies can round to the same task name; an archive inspection can fail before the next root is saved. | Reject different frequencies, share identical requests, finish collective bundle writes before inspection, and export only current qualified properties at the requested protocol. |
| `test_shared_response_records` (one and two ranks) | A job may stop after publishing a replacement but before retiring its predecessor; Raman metadata writes may fail. | Prefer the exact current task, retire only its recorded predecessor, preserve exact live dependencies, report the failure collectively, retain existing data, and record failed-invocation timing. |
| Public quadratic and derived-response lifecycle tests | Existing analysis consumers still read `calc_info.json`. | Current beta, Raman and 2PA records remain available through the compatibility output, without exporting unrelated persisted results. |

These are interface and recovery checks, independent of reference formats and historical iteration counts. They do not certify new scientific methods or recalibrate numerical references.

### Cold excited-state root selection with warmup disabled

The water `response_h2o_tpa` solid-harmonic default reaches a converged pair at
0.3833018 and 0.4017967 Ha instead of the retained lowest-pair reference at
0.3223899 and 0.3832910 Ha. Restoring the legacy six-candidate trial space does
not resolve this ordering problem before the two-state selection. The existing
guess implementation documents that unrelaxed solid-harmonic trials are not
energy ordered; the former warmup helped select the desired roots.

This case remains a scientific reference gap. A convergence certificate says
that the saved functions satisfy the selected metrics; it does not prove that
they are the lowest states of the model. Keep the reference and the assertion.
The explicit `excited.guess virtual_ao` comparison and the DALTON-seeded case
are separate evidence, with their actual input choices recorded. They must not
silently replace or qualify the default solid-harmonic case. Warmup remains disabled, and the application default remains
`solid_harmonics`. Revising the initial guess policy is separate from shared
interface integration.

### Response verification — November 20, 2026

**Timing follow-up:** `test_shared_workflow_timings.py` checks that a public
SCF-to-response workflow retains both engines' entries in one timing file,
preserves earlier records and native sources on reuse, reports a real update
when tighter criteria require it, and omits an update at the terminal check.
Saved SCF iteration intervals must match printed values to display precision;
no absolute wall-time or speed threshold is used. `test_scf_run` additionally
covers quiet SCF publication and a collective model-rejection failure with the
previous checkpoint/history preserved. `test_scf_frozen_iteration` checks
forwarded timing for the frozen-virtual loop. Validation passed: six native checks, five final checks (including serial/MPI
workflow and error-preservation cases), and the existing response wiring check.
These sets overlap on native SCFRun. Evidence is under
`.superpowers/sdd/2026-11-20-shared-scf-timings/`.

The shared response integration, functional tests and final independent review
are complete. Full milestone acceptance retains the cold-water gap above.
The final scientific batch has ten passing cases, one passing explicitly
configured variant and one failing default case; no run remains pending.

| Evidence | Outcome and scope |
| --- | --- |
| Helium HF, LDA and PBE0 alpha; water LDA alpha | Four cases pass unchanged HF or DALTON numerical checks. |
| H2 TDA, RPA and RPA without KAIN | Three cases pass unchanged excitation-energy checks. |
| LiH beta and full Raman | Two cases pass retained checks. Full Raman uses the existing MADNESS reference; an independent full-Raman oracle remains a release gap. |
| DALTON-seeded water 2PA | Passes SCF energy, both excitation energies and both `D_linear` checks. Restart from six MPI ranks to one preserves all selected native generations and performs no solver iterations. |
| Water 2PA with explicit `excited.guess virtual_ao` | Passes the original water energy and 2PA checks after a six-rank restart. This is a separate variant, not a pass of the default input. |
| Default cold solid-harmonic water 2PA | Fails the intended root-energy check. The run was stopped after that failure, before unnecessary property assembly; warmup remains off. |

The final review's five important findings have focused fixes and regression
coverage: physical-frequency collisions, interrupted task replacement,
collective archive inspection, current compatibility output, and collective
property-write failures with timing records. Six focused serial/MPI checks and
five public beta/inspection/2PA/partial-Raman/full-Raman checks pass. The mixed
property-family test also passes without an explicit metric selector. All four
retained quadratic algebra checks pass. The current `madqc` and inspection tool
are built; changed translation units passed clang-tidy with the known crashing
`bugprone-argument-comment` check excluded and the intentional recursive timing
validator warning recorded.

All 74 baseline reference files remain byte-for-byte unchanged. Historical
iteration counts are recorded, not used as acceptance gates. Timings and counts
retain their original checkpoint source when a result is reused. This batch
does not establish performance equivalence or complete the release-wide model
and functional matrix.

Evidence is indexed by `scientific-verification.json`,
`references-verified.json` and `progress.md` under
`.superpowers/sdd/2026-11-19-full-response-unification/` in the active worktree.
The numerical comparison records are also retained there under
`scientific-records/`; full solver output and archives remain in the scratch
directories identified by the index. Earlier deliberately interrupted runs
retain their nonzero exit records and are not counted as passing runs.


## Shared run layout and offline development relocation — November 20, 2026

Stable intent: directory organization must not change which native functions,
measured metrics or exact source references a restart uses. Shared properties
and timings belong to the run; engine diagnostics belong to its calculation.

| Case | Evidence required |
| --- | --- |
| `test_shared_run_layout.py` (serial/MPI2) | Nested relative prefix; root results/timings/checkpoints; separate ground/response diagnostics; plot discovery with calculation identity; unchanged native sources on continue; collective old-layout rejection before mutation; offline relocation and native reuse; read-only preserves records. |
| `test_migrate_run_layout.py` | Rename preserves bytes; destination conflicts, locks, symlinks, unexpected tasks and duplicate diagnostic names are rejected; interrupted rename resumes; invalid journal cannot escape the run directory. |
| Existing SCF/requalification/recovery and workflow completion tests | Root discovery, frozen occupied/virtual store sharing, preparation, missing native state, incomplete upstream blocking and report-write failures still work with separate analysis paths. |
| Existing seed/response ladder/coupled/VBC/subworld/inspection/timing tests | Seeds never overwrite an existing native restart; state progress, exact direct inputs, subgroup records, native loading and both engines' timing history use the same root. |

The relocation tool is an explicit offline development aid. It does not validate
scientific compatibility or convert checkpoint formats; ordinary native readers
still perform those checks. Historical analysis can contain old path text and
is never restart authority. No numerical reference values or tolerances change.
Verification evidence: `.superpowers/sdd/2026-11-20-shared-run-layout/`.

Final verification for this increment: 20 focused CTest checks pass. The independent review also identified artifact globs in the no-KAIN H2 and cold/seeded water 2PA qctests. Their paths now use analysis/response; the actual qctest reader has RED/GREEN evidence against generated ES output. Numerical checks, tolerances and all 74 reference files are unchanged. This does not rerun or recalibrate those three scientific cases. See verification.json and review-artifact-paths.log in the evidence directory above.
