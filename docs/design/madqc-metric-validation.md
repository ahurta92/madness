# madqc metric definitions and validation evidence

This document records model definitions and their supporting evidence for the
shared convergence layer. Source inspection, scientific approval and numerical
validation are separate statuses. An empty or pending entry is not an approved
contract. Other required engines still need their entries.

The [legacy convergence audit](madqc-legacy-convergence.md) separately records
SCF and molresponse engine stopping rules and their shared-input mappings. It
is source evidence, not new calibration of tolerances or a performance result.

## Nemo — source review, 2026-11-12

**Status:** source behavior traced; the author selected the definitions below
as proposals in the design on 2026-11-12, pending Florian's confirmation.
Numerical validation has not been run. Stage 9 finding `readiness-1-1` remains open.

[Concise confirmation checklist for Florian](nemo-metric-confirmation.md).

### Proposed definitions recorded in the design

- Maximum regularized residual: `max_i norm2(F_i - G_F(F)_i)`, with no `R`
  multiplier, occupation weighting or function-count normalization.
- Physical per-spin density change; for unit-occupation closed-shell orbitals,
  the field is `R^2 sum_i abs(F_i)^2`, without an extra factor of two.
- Absolute change of the existing physical kinetic quantity `ke0`
  (`energies[1]`), retaining its closed-shell factors and nuclear-factor correction.

These are the author's proposed choices, not Florian's scientific approval.
Occupation support and weighting, derivative/smoothing conventions and exact
units still need confirmation. Total-energy change and the common prepared-state
measurement point remain as already designed. Numerical tolerances will be
validated during implementation.


**Source:** working tree at commit
`76255552e16e542f6b5ae2d0d7e2b6b034c7ed48`. The inspected implementation uses
real regularized spatial orbitals and closed-shell factors. Its occupation
support boundary needs confirmation; the existence of an occupation array is
not evidence that every energy and convergence path handles arbitrary values.

### Representation already defined

Nemo stores regularized functions `F_i = psi_i / R` in `calc->amo`; the physical
orbitals are `psi_i = R F_i`. The constructors mark the SCF-owned archive as
Nemo representation and record the nuclear correlation factor. Normalization
uses `||R F_i||_2`, and orthonormalization uses overlaps of those physical
orbitals. This is a defined implementation, not a new representation to design.

Sources: [constructors](../../src/madness/chem/nemo.cc#L155),
[normalization and orthonormalization](../../src/madness/chem/nemo.h#L90).

At protocol changes, `set_protocol` rebuilds required nuclear factors and
operators before reprojecting and normalizing orbitals with the new `R`. It
also updates the AO basis when needed. The metric context must identify these
settings; measurements from a changed representation cannot silently retain
comparison eligibility. See [protocol setup](../../src/madness/chem/nemo.h#L677).

### Metrics already computed

**BSH residual.** `solve` constructs the raw regularized proposal `G_F(F)`,
including the weak-XC flux contribution when active, then forms
`r_i = F_i - G_F(F)_i`. The scalar currently used for convergence is
`sqrt(sum_i ||r_i||_2^2 / N)`. It is an RMS of per-function L2 norms, with no
`R` multiplication or occupation weights. It is evaluated before the nonlinear
solver update. The shared maximum reduction is already agreed; the remaining
scientific question is which representation its norms describe.

Sources: [proposal and residual](../../src/madness/chem/nemo.cc#L701),
[vector norm implementation](../../src/madness/mra/vmra.h#L961).

**Density change.** The iteration field is
`d = R^2 sum_i |F_i|^2`. The loop reports `||d_new - d_previous||_2`.
`compute_density` sums all supplied functions with unit weights. For fully
occupied closed-shell spatial orbitals, this is one spin density, and the total
charge density is twice it. It is therefore not missing a factor of two if the
intended convergence field is the per-spin density.

The property paths instead construct
`rho_total = 2 R^2 sum_i occupation_i |F_i|^2` through `make_density`.
These paths differ in occupation weighting as well as the closed-shell factor.
With unit occupations the iteration field matches the one-spin convention used
by restricted SCF. For other occupations, the present two definitions must not
be described as interchangeable.

Sources: [iteration density](../../src/madness/chem/nemo.cc#L730),
[unit-weight helper](../../src/madness/chem/nemo.h#L172),
[property density](../../src/madness/chem/nemo.cc#L333),
[occupation-weighted density](../../src/madness/chem/nemo.cc#L999),
[SCF occupation weights](../../src/madness/chem/SCF.cc#L1476),
[restricted spin fields](../../src/madness/chem/SCF.cc#L2973).

**Energy change.** `energies[0]` is the assembled total energy. The current
check is `abs(energies[0] - oldenergies[0]) < econv` when selected. Assembly
includes the existing HF/DFT, nuclear, PCM and dispersion contributions when
applicable. No new total-energy definition is needed for the shared category.

The separate legacy `each_energy` check takes the largest change across the
returned component vector and compares it with `3*econv`. That vector includes
the physical kinetic diagnostic but is not equivalent to either total-energy
change or a dedicated kinetic-change criterion. The accepted input migration
already forbids silently translating `each_energy` into `energy_change`.

Sources: [energy assembly and returned vector](../../src/madness/chem/nemo.cc#L779),
[convergence checks](../../src/madness/chem/nemo.h#L327),
[input mapping](../../input_option_mapping.md).

**Kinetic quantity.** A physical kinetic diagnostic already exists in the live
energy path. With the current real, unit-occupation closed-shell expressions,
the code forms

- `Treg_ij = 0.5 sum_a <D_a(R^2 F_i), D_a F_j>`;
- `Kreg = 2 trace(Treg)`;
- `Kphysical = Kreg + <R^2 sum_i F_i^2, U1dotU1>`.

`U1 = -grad(R)/R`; `U1dotU1` is projected using the dedicated nuclear-factor
functor. `Kphysical` is `ke0`, printed as `kinetic only` and returned as
`energies[1]`. The implementation identifies it with the physical orbital
kinetic energy, including closed-shell occupancy. `Kreg` alone is not that
quantity. Both scalars are built from the same pre-update energy evaluation.

The header also contains product-rule and direct-orbital derivative versions
of physical kinetic energy. The current `solve` energy path uses the corrected
trace expression above, not those alternative helpers. The direct derivative
helper is labelled potentially imprecise in the source. Finite MRA truncation
and smoothing make numerical agreement a validation question, not something
proved by the algebraic identity alone. Do not replace the dedicated
`U1dotU1` construction with squared sampled components without scientific
review: the source records a difference at the nuclei.

Sources: [kinetic matrix construction](../../src/madness/chem/nemo.cc#L628),
[derivative matrix implementation](../../src/madness/chem/SCFOperators.cc#L81),
[physical correction](../../src/madness/chem/nemo.cc#L788),
[alternative helpers](../../src/madness/chem/nemo.h#L226),
[correction-factor projection](../../src/madness/chem/nemo.h#L201).

A separately selected `kinetic_change` measurement and comparison history are
not present in this loop. Reusing `ke0` would avoid inventing a new observable
or requiring a second full derivative evaluation solely for that metric.

### Implementation work already determined by the shared contract

Current ordering is: prepare/localize; build potentials and energy; optionally
canonicalize, truncate and normalize; build BSH residual; apply the nonlinear
update and step restriction; orthonormalize; measure the new density; check
convergence; save the updated orbitals. Energy therefore precedes the canonical
preparation, residual precedes the update, and density follows the update.
The saved orbitals are the updated ones while `current_energy` is assigned the
pre-update energy. They do not establish one evaluated checkpoint sample.

Implement the already-approved prepared-state evaluation boundary: all metrics
and saved evidence refer to the same prepared iterate before any next update.
Account for canonicalization/truncation when reusing earlier work; do not assume
finite-precision changes preserve every measurement exactly. Keep the raw BSH
proposal for the update if continuation is required. Stopping must not apply
an unnecessary nonlinear update. First-sample change metrics need a valid
comparison baseline: the current zero initial density and energy placeholders
are not evidence of a previous evaluated state.

The maximum reduction, measurement point, missing-comparison behavior and
checkpoint agreement do not need to be reopened as scientific choices. Nemo
owns functions and previous samples, uses shared helpers, and supplies scalar
measurements. Its new metric inputs still require implementation; `value`
currently rejects the old explicit BSH/kinetic accuracy controls and disabling
default criteria.

Sources: [iteration](../../src/madness/chem/nemo.cc#L572),
[current input restrictions](../../src/madness/chem/nemo.cc#L201).

### Narrow questions for Florian

1. **Residual representation.** Confirm the recorded proposal
   `max_i ||F_i - G_F(F)_i||_2`, preserving the regularized equation residual
   and replacing its RMS reduction with a maximum. Multiplying by `R` would
   instead measure the mapped physical orbital difference at fixed `R` and
   change the norm. The proposal therefore explicitly omits that multiplier;
   its scientific confirmation and tolerance validation are still pending.
2. **Spin and occupation scope.** Confirm that the convergence field should be
   the physical one-spin density, matching restricted SCF, and establish the
   supported occupation cases. Unit occupations reproduce the existing loop;
   occupation-weighted cases require consistent density and kinetic/energy
   handling, not merely adding weights to one diagnostic. Do not infer or
   change supported scope solely from this review.
3. **Physical kinetic estimator.** Confirm `ke0`/`energies[1]` as the intended
   quantity for `abs(K_n-K_previous)`, with its current closed-shell factors
   and dedicated correction functor. Confirm the applicable derivative and
   smoothing conventions and how to validate it against the product-rule
   evaluation. Energy and kinetic quantities are to be reported in the
   model's atomic energy units; record exact field/norm units in the approved
   contract. No default tolerance is being selected here.

These questions request confirmation of scientific meaning and supported
cases. Many operations already exist. They do not ask Florian to redesign the
shared controller or decide routine loop ordering. No message has been sent.

### Validation to perform during implementation

The nine current `src/examples/qc/nemo_*` case specifications were inspected.
They check total energies, selected orbital eigenvalues, model labels and some
convergence/optimization outputs. They do not directly assert the physical
kinetic diagnostic or the proposed maximum residual and per-spin density
measurements. These statements describe the inspected case specifications, not
a claim that no other repository tests cover related operators.

Preserve these references and add focused evidence, recorded in the curated
test index:

- Compare the selected residual norm on a fixed prepared state with explicit
  per-function measurements, including a nonuniform-residual case. Validate
  its target values after changing from RMS to maximum.
- Check the occupation-weighted one-spin density integral and closed-shell
  relation to the property density for the confirmed supported cases.
- Compare the corrected physical kinetic estimator with the approved
  product-rule evaluation on the same frozen orbitals at successively tighter
  representations. Record errors and cost; do not assert exact equality at
  finite truncation or prescribe an untested tolerance.
- Check matching energy, density, residual, kinetic and checkpoint samples;
  first-sample unavailability; terminal evaluation; and protocol/restart
  comparison resets through the shared contract tests.

Scientific approval, tested tolerances and measured validation results remain
pending. Completing this source review does not complete the full E01–E12
scientific-definition prerequisite.

## Nemo nuclear response (CPHF) — source review, 2026-11-13

**Status:** author-approved task/iteration/restart structure; scientific response
definitions and the Hessian-row acceptance role await Florian. No numerical tests
were run. Source checkout: `76255552e16e542f6b5ae2d0d7e2b6b034c7ed48`.
See the [confirmation questions](nemo-metric-confirmation.md#cphf-nuclear-response-metrics)
and design section 3 for the adopted contract.

| Existing code | Required shared-layer adoption |
| --- | --- |
| `compute_all_cphf` allocates up to `3N` response vectors and loops sequentially over atom/axis. Once the ground reference/fixed terms are known, a solve uses only its own response. | Nemo-owned independent perturbation tasks, initially sequential, with exact reference context and per-task progress. |
| `solve_cphf` forms raw `xi - G(xi)`, computes RMS and maximum, then applies KAIN/plain fixed-point update before its Hessian-row check. It returns only functions after at most ten evaluations. | Evaluate one sample before advancing, use the shared maximum/scales and resolved allowance, and return explicit convergence/limit/failure evidence. Preserve numerical update choices. |
| The residual gate is `RMS/norm(xi) < proto.dconv`; the row-change gate is `absmax(delta_h) < 1e-2`. Response-density terms include occupations, closed-shell factors and derivatives of the nuclear factor. | No automatic numerical equivalence to the new absolute/relative tolerances. Confirm response representation, density/kinetic definitions, units and the row-check role before adapter acceptance. |
| `read_cphf` loads `xi_` functions and returns; `restart_cphf` loads `xi_guess`. `save_function` stores a native vector length and functions without shared task metrics/reference identity. | Shared native checkpoint assessment and strict recovery; either restart request resolves through reuse-or-iterate, never an assumed-converged private reader. |
| A preparatory pass uses the `LDA` functional string on a shortened ladder; `drop_last_rung` keeps a one-rung ladder nonempty. The requested-model solve currently starts at the final rung. | Disable the separate iterative guess pass; keep direct guesses and common requested-model protocol stages. Validate the numerical/cost effect rather than assuming the old pass was redundant. |
| Finalization resets protocol from `econv`, reloads `xi_guess`, subtracts `parallel`, truncates and saves `xi_`. | Record the actual iterated representation and measured sample; preserve requested final representation. Derive the corrected full response with source identities instead of attaching raw residuals to a different unknown. |
| `hessian` then assembles electronic/nuclear/dispersion terms and vibrational properties. Both standalone Nemo and the madqc binding invoke it, without per-CPHF convergence outcomes. | Require eligible responses before direct property completion; publish supported requested quantities with units, source references, metrics once and timings. Incomplete dependencies leave properties unavailable. |

Sources: [Nemo equations/assembly](../../src/madness/chem/nemo.cc),
[controls and archive helpers](../../src/madness/chem/nemo.h),
[protocol helper](../../src/madness/chem/SCFProtocol.h),
[standalone dispatch](../../src/apps/nemo/nemo.cc),
[madqc binding](../../src/madness/chem/Applications.hpp).

Curate small supported HF/DFT Hessian and vibrational references with Florian,
including tighter calculations, units, mode conventions and costs. Existing
Nemo ground-state decks are not evidence for the CPHF contract. Add partial
perturbation restart, exact-reference invalidation, near-zero scaling, terminal
sampling, correction/representation, no-warmup, strict-checkpoint and property
eligibility tests under acceptance A02/A04–A06, B01/B02 and C01–C04. Keep scientific
confirmation separate from the fact that a native archive can be loaded.

## SCF virtual orbitals — PR #798 and source review, 2026-11-13

**Status:** the author approved adoption of both SCF virtual paths. PR metadata,
final patches, discussion and local source were inspected; no numerical tests
were run. The shared contract is in design section 3. This is source evidence,
not validation of the new metric defaults or checkpoint adapter.

[PR #798](https://github.com/m-a-d-n-e-s-s/madness/pull/798) merged on 2026-09-23
as `c51a79217c6469f848d4af7284d614abff285d9f`; its final head was
`bb161293c36a2f338504012dee68efeefef4289f`. The merge is an ancestor of the
inspected checkout `76255552e16e542f6b5ae2d0d7e2b6b034c7ed48`. Later changes
exist in SCF, restart planning, parameters and exchange; the PR is context for
adopting the current code, not a patch to reapply.

| Source fact | Consequence for adoption |
| --- | --- |
| `SCF::compute_residual` forms `psi - G(psi)` and returns the maximum function norm. `solve_virtuals` uses it per spin before KAIN, with a legacy `5*dconv` gate. | Reuse the raw residual construction; use the shared maximum over the actual iterated functions and controller decisions. Frozen occupied physical changes are inapplicable virtual criteria. |
| `solve_virtuals` builds a fixed occupied mean field, applies explicit occupations, solves alpha/beta virtual blocks and reports the archived occupied energy. | Keep exact occupied-reference and occupation context; do not present that inherited energy or unchanged density as convergence of the virtual task. Both spin blocks must be accounted for. |
| `canonicalize_virtuals` rotates functions and associated operator work. The coupled path also canonicalizes at termination. | Move/reconcile final transformations with the common sample boundary so residuals, eigenvalues and saved functions agree. Never carry pre-rotation measurements as evidence for a changed state. |
| `MolecularEnergy` applies `nv_extra`/`nv_step` on initial or padded work and overrides the intermediate cap with `nv_its`. | Save phase and current/target membership; preserve the schedule while adopting shared update counting. Enough archived orbitals alone does not prove preparation finished. |
| Padding checks loaded occupied counts and orthogonalizes new atomic-guess virtuals against loaded functions. Restart planning records alpha count; the loader handles both spin arrays. | Inventory both counts explicitly in shared evidence and test alpha-only/beta-only shortages. Added virtuals need fresh evidence; inspection cannot generate them. |
| Occupied-only exchange/localization and retained virtual–virtual coupling fix numerical behavior. KAIN solve/combine/trim were factored out. | Preserve these operations; changes to KAIN's near-dependence handling are separate numerical work. |

Sources: [SCF implementation](../../src/madness/chem/SCF.cc),
[SCF scheduler/reference checks](../../src/madness/chem/SCF.h),
[restart planner](../../src/madness/chem/RestartPlan.h),
[orbital loading](../../src/madness/chem/MolecularOrbitals.h),
[parameters](../../src/madness/chem/CalculationParameters.h),
[exchange](../../src/madness/chem/exchangeoperator.cc).

The PR reviews identified short-archive/read-only, occupation application,
reference validity and final-state consistency issues. The author reported fixes
before merge; current source includes guarded loading, explicit occupation
application, frozen-reference checks and no update after the terminal virtual
evaluation. This does not prove the new shared-layer guarantees already exist.
The [author's reply](https://github.com/m-a-d-n-e-s-s/madness/pull/798#issuecomment-5786077328)
and [approval](https://github.com/m-a-d-n-e-s-s/madness/pull/798#pullrequestreview-5294684286)
explicitly separate replacing KAIN with `XNonlinearSolver` from this work.

### Evidence to curate and validation still required

- Preserve [test_restart.cc](../../src/madness/chem/test_restart.cc)'s short/equal/long
  archive and continuation cases. Extend coverage to both spin counts and the
  shared prepared-state/progress contract; inspection remains free of solver work.
- Curate [h2o_hf_nvstep.in](../../src/apps/moldft/tests/h2o_hf_nvstep.in) as a schedule
  and capped-result case. Its own comments document unbound virtuals and a residual
  plateau; it is not a passing residual-convergence reference. Do not turn that
  plateau into stall acceptance or loosen final requirements.
- Preserve [ne_hf_corehole.in](../../src/apps/moldft/tests/ne_hf_corehole.in), explicit
  occupation validation and NWChem virtual seeding. Validate occupation/reference
  changes rather than inheriting energy or convergence from another model.
- Add curated bound-virtual comparisons for coupled versus frozen solves, HF/DFT,
  restricted/unrestricted spin and occupied localization choices. PR-reported
  manual comparisons motivate cases but are not newly executed or approved
  numerical reference files. Record costs and reviewed tolerances normally.
- Inject interruption during preparation and after an evaluated terminal sample;
  test phase/count/protocol recovery, added virtuals, changed reference/settings,
  no extra KAIN update, and archive/metric/eigenvalue identity. Preserve cadence
  and latest-required retention without an archive for every orbital count.

Residual units, absolute defaults and any additional virtual-specific
energy/density/kinetic definitions retain the normal scientific validation gate.
This revision does not invent a virtual eigenvalue-change convergence metric.

## Ordinary SCF — initial default validation, 2026-11-18

The shared C++ adapter now supplies absolute defaults of `1e-4` for the maximum
raw BSH residual and maximum per-spin density change, and `1e-5` hartree for
physical total-energy and kinetic-energy change. The initial requested criteria
are residual, density and energy; kinetic is diagnostic unless selected.
Relative defaults remain zero. Intermediate stages retain the reviewed threshold
floors; final requests have no floor. Public workflow input/restart activation
remains pending for both ordinary and frozen-occupied paths.

`test_scf_defaults.cc` drives real SCF evaluations and the shared controller from
unchanged curated QC decks. This calibration compares total energies with the
existing `tasks[0].properties.energy` references at their unchanged `1e-5`
hartree allowance. It does not replace the decks' eigenvalue/property checks.
Each run allows 25 main updates per protocol. The comparator request is ten
times tighter in all three selected metrics.

| Case | Protocol | Default / tighter updates | Default / tighter energy minus reference (hartree) | Outcome |
| --- | --- | --- | --- | --- |
| `scf_he_hf` | `[1e-4,1e-6]` | 6 / 6 | -2.967e-08 / -2.967e-08 | Both complete; reference agreement |
| `scf_he_pbe0` | `[1e-4,1e-6]` | 6 / 7 | -4.341e-08 / -4.385e-08 | Both complete; reference agreement |
| `scf_li_tpss` | `[1e-4]` | 25 / 25 | -2.383e-06 / -2.385e-06 | Both stop at the cap; residual remains unmet |
| `scf_h2o_hf` | `[1e-4,1e-6]` | 9 / 20 | -6.283e-07 / -5.713e-07 | Both complete; reference agreement |

The coarse open-shell Li/TPSS case reaches an energy within `2.4e-6` hartree of
its reference but retains a BSH residual near `2.95e-4` after 25 updates. It is
correctly **incomplete**. Changing its final protocol to `1e-6` meets the default
request in four updates. Its energy moves by about `1.46e-4` hartree from the
coarse reference: the finer representation is a convergence experiment, not a
claim of equality with a coarse oracle. At the same finer protocol, the tighter
request completes in six updates and changes the candidate energy by only
`1.79e-7` hartree. The regression keeps both representation cases, checking
an unmet residual at the coarse cap and completion after refinement.

Selecting kinetic change in He/HF adds four updates (ten total) and finishes
with a kinetic change of `7.64e-6` hartree. This verifies optional selection and
physical weighting in that case; it is not broad kinetic-default validation for
all systems or models. Measured solve times, all iteration metrics and commands
are retained under `.superpowers/sdd/design_draft/scf-input-defaults/`. Runs used
seven MADNESS threads, sequential BLAS and no rank binding on dd. Timings include
protocol preparation but exclude initial-guess construction, and some runs
shared the host with builds: they are study records, not a performance bound.
No reference values or deck settings were regenerated.

Final verification passed 25 targeted CTests, including serial/two-rank input
failures, shared SCF metrics/native evidence and unchanged SCF, Nemo, response
and optimizer references. Independent review fixes cover collective combined
input selection, whole-read rollback, optional-Libxc test registration and
propagation of an invalid-state protocol stop. Initial calibration histories
were checked: all their intermediate endpoints were valid, so those numerical
comparisons remain applicable. Final changed-line lint passed; the existing
clang-tidy 18 crashing-check workaround remains documented in the evidence.

## Frozen-occupied SCF — measurement and iteration validation, 2026-11-18

The typed virtual adapter reports one task maximum across the active alpha/beta
raw BSH residuals, with individual orbital norms and residuals retained as
evidence. Fixed occupied density, total-energy and kinetic changes are
inapplicable. Independent numerical reevaluation checks those residuals and the
prepared virtual eigenvalues on He/HF, He/LDA and unequal-spin Li/HF. The occupied
functions are compared with their exact native source before virtual iteration.

Focused regressions preserve slightly off-unit occupied functions during padding
and updates, reject incomplete arrays/empty active spans without erasing the
previous outcome, and verify that the legacy alpha KAIN history expires before
beta iteration. These address implementation correctness and resource lifetime.
The deliberately loose completion request is not a physical spectrum reference;
virtual default-tolerance, bound-spectrum, broader-model and performance
calibration remain open. Exact source eligibility is still the workflow caller's
responsibility until public input/native restart adoption.
