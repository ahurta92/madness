# Legacy SCF and molresponse convergence rules

Recorded November 19, 2026, before further restart simplification. This is a
source audit and a mapping to the shared metric inputs, not a change to solver
behavior or a new set of defaults. Preserve these rules as a reference for
reproducing comparable convergence without stall windows or automatic restarts.

“Legacy” here means the ordinary SCF engine check and the molresponse v3 engine
checks retained from drop 6. It does not mean the removed response v2 engine,
Nemo, or the separate TDHF implementation. The working tree is based on
`76255552e16e542f6b5ae2d0d7e2b6b034c7ed48`, with substantial uncommitted shared-layer
changes. The pre-split SCF loop was also checked at drop-6 commit
`9179f5d32076a2ba555906dd146dc330ffc040d7` to distinguish its old measurement order
from the current implementation. Function names below identify the source even
when line numbers move.

## The rules to preserve

Let `tau` be the current protocol's truncation threshold, `d` the resolved legacy
`dconv`, and `D = max(tau, d)`. These are empirical tolerances, not mathematical
bounds on the attainable residual or the error in a property. Polynomial order
affects the calculation, but does not enter these tolerance formulas directly.

| Solver | Legacy convergence test | Not required by this test |
| --- | --- | --- |
| Ordinary SCF | Each spin-density change `< max(5, natom) * D`, and maximum orbital BSH residual `< 5 * D` | Total-energy and kinetic-energy changes |
| Frozen-occupied SCF virtual solve | Maximum BSH residual of each iterated virtual block `< 5 * D` | Density, total-energy and kinetic-energy changes |
| Driven response, including static response | For every channel: combined BSH residual `<= 5*D + 5*D*response_norm`, and density change `<= 5*D + 5*D*density_norm` | Energy change, kinetic change, property change |
| Excited-state response, TDA or full/RPA | For every root: density change `< 5 * D` and excitation-energy change `abs(delta omega) < 5 * D` | BSH residual and kinetic change |

All required checks are combined with AND. A diagnostic may still be measured
and printed even when it does not determine convergence. In particular, the
excited-state BSH residual is reported and used for a divergence guard, but it
is not a convergence requirement in this legacy rule.

Sources: [`SCF::engine_check_limits`](../../src/madness/chem/SCF.h#L402),
[`SCF::solve_impl`](../../src/madness/chem/SCF.cc#L3363),
[frozen virtual check](../../src/madness/chem/SCF.cc#L2717),
[`ConvergencePolicy::effective_for_thresh`, `fd_leg_converged`, `es_root_within`](../../src/apps/molresponse/solvers/convergence_policy.hpp#L321).

## Defaults and measurement conventions matter

**SCF.** The registered legacy `dconv` default is `1e-4`; `econv` is `1e-5`,
but ordinary SCF's engine test does not use `econv`. Nor does the registered
`convergence_criteria` list select that engine test: other engines, notably
Nemo, consume it. The deprecated `conv_only_dens` switch removes the BSH test.

SCF uses `max_i ||psi_i - G(psi)_i||_2`, taking the maximum across iterated
alpha and beta functions. Density changes are the L2 differences of the
occupation-weighted spin densities. Restricted SCF compares the same one-spin
density in both slots. It does not compare the total charge-density change to
the same tolerance without a factor adjustment.

Sources: [legacy parameter registrations](../../src/madness/chem/CalculationParameters.h#L125),
[`SCF::compute_residual`](../../src/madness/chem/SCF.cc#L2230),
[`SCF::measure_iteration`](../../src/madness/chem/SCF.cc#L3045).

**Response.** The normal madqc parameter path derives an omitted response
`dconv` as **`100 * min(dft.protocol)`**. It uses `1e-6` as the finest-threshold
fallback if the list is empty. An explicit response `dconv` wins. Thus the
registered `ResponseParameters` default `1e-6`, and the bare `ConvergencePolicy`
default `1e-4`, are not reliable descriptions of a normal resolved madqc run.
The adapter copies the resolved value into `policy.dconv_user`.

For example, a finest threshold of `1e-6` gives `d=1e-4`, hence final absolute
factors `5e-4` and driven relative coefficients `5e-4`. A finest threshold of
`1e-8` gives `d=1e-6`, hence `5e-6`. A single `1e-4` protocol with no explicit
response `dconv` gives `d=1e-2`, hence `5e-2`. These describe legacy defaults;
they are not newly recommended shared defaults. The curated helium alpha and
H2 TDA decks explicitly set response `dconv=1e-4`.

The response factors for BSH, density, omega and the driven relative term are
all `5.0` in `ConvergencePolicy`. The relative coefficient also uses the
protocol-floored `D`; it is not just `5*d` when `tau > d`.

Sources: [`ParameterManager::set_derived_values`](../../src/madness/chem/ParameterManager.hpp#L201),
[adapter assignment](../../src/apps/molresponse/madqc_adapter.hpp#L553),
[policy factors](../../src/apps/molresponse/solvers/convergence_policy.hpp#L72),
[helium deck](../../src/examples/qc/response_he_alpha/response_he_alpha.in),
[H2 TDA deck](../../src/examples/qc/response_h2_es_tda/response_h2_es_tda.in).

Legacy response BSH uses `sqrt(sum_i ||r_i||_2^2)` across all stored functions,
including the applicable X/Y and spin blocks. There is no division by function
count. Driven `response_norm` is the corresponding combined norm **after** KAIN;
`density_norm` is the norm of the density calculated from the step's input.
The raw residual and density change are measured before KAIN.

The density is the existing kernel-induced response density, with its physical
spin and X/Y factors. For real closed-shell functions it is `4 sum(phi_i X_i)`
in static response, `2 sum(phi_i X_i)` in TDA, and
`2 sum(phi_i (X_i + Y_i))` in full response. Open-shell kernels sum the spin
contributions with their own factors. The density change is the L2 difference
between successive such fields. Excited-state comparisons use root slots;
mixing, phase and locking affect the trajectory and cannot be reproduced just
by copying tolerance numbers.

Sources: [static norms and densities](../../src/apps/molresponse/kernels/static.hpp#L50),
[full norms and densities](../../src/apps/molresponse/kernels/full.hpp#L60),
[TDA conventions](../../src/apps/molresponse/kernels/tda.hpp#L139),
[vector norm](../../src/madness/mra/vmra.h#L956),
[`FDSolver::step`](../../src/apps/molresponse/solvers/fd_solver.hpp#L772),
[ES density comparisons](../../src/apps/molresponse/solvers/es_solver.hpp#L657).

## When the checks run

It is not true of every current path that convergence is checked only at the
end. Separate *measuring/checking*, *leaving an intermediate protocol*, and
*declaring final convergence*:

| Path | Check and stopping behavior |
| --- | --- |
| Pre-split SCF, drop 6 | Checks spin-density change and the previous update's BSH residual on iterations after the first. The residual therefore need not describe the orbitals later saved. |
| Current SCF engine-check fallback | Uses the formulas above with freshly evaluated raw BSH work and a valid density predecessor; final orbital preparation is reevaluated. |
| Current legacy response loops | `step()` computes measurements and applies the update/KAIN; `iterate()` then asks for a stopping decision. FD and ES require at least two steps for their change checks. |
| Current unadopted response intermediate protocols | The November 6 fixed-work wrapper checks the engine verdict but suppresses early exits until the cap. At the final protocol, the engine check can stop the loop early. This is transitional behavior. |
| Old named excited-state explicit criteria | The loop uses the engine density/omega test; the separate requested verdict was assessed outside that loop. This machinery is not the released input contract. |
| Adopted shared SCF and driven-alpha paths | Evaluate the current functions and check selected metrics before every possible KAIN update. Intermediate criteria can finish a protocol early; final criteria can complete the solve early. Metrics are reevaluated after final SCF orbital preparation before saving. |

Legacy ES locking is a separate algorithmic choice. The bare policy defaults
to no locking, while the production executor enables it; locked roots retain
earlier density/residual measurements. This audit records that behavior without
restoring freezing to the shared design. Likewise, divergence and iteration-cap
exits are not convergence, and optional dormant stall machinery is not part of
the tolerance mapping.

Sources: drop-6 `SCF.cc::solve` at the commit above;
[current SCF assessment](../../src/madness/chem/SCF.cc#L3304),
[`iterate` and `iterate_shared`](../../src/apps/molresponse/solvers/iterate.hpp#L36),
[`rung_exit`](../../src/madness/chem/ConvergenceSettings.h#L230),
[`bundle_rung_sample`](../../src/apps/molresponse/solvers/convergence_policy.hpp#L376),
[ES eligibility](../../src/apps/molresponse/solvers/es_solver.hpp#L1209),
[executor locking default](../../src/apps/molresponse/calc/calc_executor.hpp#L229).

## Expressing comparable rules with the shared inputs

The four metric categories are available for reporting; they need not all be
required for stopping. Select the legacy subset explicitly when comparing with
a legacy calculation. No `dconv`/`econv` aliases or second “legacy mode” are needed.

| Legacy intent | Shared selection and absolute limits at a specified protocol |
| --- | --- |
| Ordinary SCF | `bsh_residual = 5*D`; `density_change = max(5,natom)*D`; relative tolerances zero |
| Density-only SCF | Select `density_change` alone, with the same density limit |
| Frozen virtuals | Select `bsh_residual` alone, with limit `5*D` |
| Driven response | Select `bsh_residual` and `density_change`; absolute tolerances `5*D` and relative tolerances `5*D` are comparison starting points, subject to the norm/scale differences below |
| Excited-state response | Select `density_change` and `energy_change`, both `5*D`, with `energy_change` meaning `abs(delta omega)`; shared ES adoption is still pending |

No legacy stopping heuristic for `kinetic_change` was found in these engine
checks. Selecting it adds a requirement. Selecting total-energy change for SCF,
or BSH convergence for excited states, likewise adds a requirement rather than
reproducing the legacy selection. Continue reporting these quantities when
available, independently of which ones stop the calculation.

For a molecule with at most five atoms, legacy SCF `dconv=1e-4` at final
`tau=1e-6` has the following equivalent tolerance magnitudes in a `dft` block:

```text
convergence.criteria [bsh_residual, density_change]
convergence.bsh_residual 5e-4
convergence.density_change 5e-4
```

This is an input fragment, not a complete deck. For twelve atoms the density
limit would instead be `1.2e-3`. Energy and kinetic changes are not selected in
this comparison. Current shared SCF defaults select BSH, density **and energy**,
with absolute defaults `1e-4`, `1e-4`, `1e-5`; copying the default input therefore
does not reproduce the old SCF gate.

With the same `D=1e-4`, the analogous driven-response fragment is:

```text
convergence.criteria [bsh_residual, density_change]
convergence.bsh_residual 5e-4
convergence.bsh_residual.relative 5e-4
convergence.density_change 5e-4
convergence.density_change.relative 5e-4
```

This reproduces the coefficient values, **not the old response norm or exact
iteration trajectory**. The adopted driven path currently requires explicit
absolute tolerances for selected applicable metrics. A fixed driven frequency
does not supply an energy-change criterion.

The shared comparison is `metric <= absolute + relative * scale`. Important
differences to preserve visibly:

- The shared BSH metric is the maximum per-function norm. Legacy response uses
  the combined norm. For `N` functions, `maximum <= combined <= sqrt(N)*maximum`;
  no single constant gives exact equivalence for arbitrary residuals. Do not
  silently add a function-count conversion or change the agreed metric.
- The shared BSH relative scale is the largest current per-function norm;
  legacy FD uses the combined post-update norm. Shared density scaling uses
  the maximum of the current and previous density norms; legacy FD uses the
  current step-input density norm. Equal coefficients need numerical comparison.
- Shared measurements describe the evaluated state before update. Retain that
  correction rather than reproducing legacy measurement lag. Missing comparison
  history is unavailable, not zero change.
- SCF and ES legacy tests use strict `<`; the shared evaluator uses `<=`.

Sources: [SCF defaults](../../src/madness/chem/SCF.cc#L2985),
[SCF selection](../../src/madness/chem/SCFInput.h#L26),
[driven input selection](../../src/apps/molresponse/solvers/shared_fd_input.hpp#L5),
[shared metric reductions/scales](../../src/madness/chem/MetricHelpers.h#L9),
[driven measurements](../../src/apps/molresponse/solvers/fd_solver.hpp#L665),
[shared tolerance comparison](../../src/madness/chem/ExplicitCriteria.h#L220).

## Protocol heuristic versus the final request

The accepted shared policy already keeps the simple heuristic at intermediate
protocols: `effective_absolute = max(requested_absolute, model_floor(tau))`.
Current SCF floors are `5*tau` for BSH and `max(5,natom)*tau` for density.
Current driven floors are `5*tau` for both. Energy and kinetic change have no
such floor in these adapters. Selected relative tolerances remain unchanged.

For a tightening ladder, choosing the SCF final absolute limits as
`5*max(d,tau_final)` and `max(5,natom)*max(d,tau_final)` reproduces the old SCF
limit magnitudes at every protocol through these intermediate floors. It does
not promise the same iteration count after correcting measurement timing.

Driven response has an additional difference: legacy relative coefficients grow
with `max(tau,d)`, while the shared relative inputs stay fixed. The shared
intermediate rule therefore does not exactly reproduce legacy relative gates
when the coarse threshold exceeds `max(d,tau_final)`. Record and validate this
difference; do not silently change the accepted relative-tolerance convention.

At the final protocol, use the written shared limits without automatic
relaxation. If the goal is to match an old run, write its actual resolved limits,
including the old factors, as the new request. If tighter limits are requested,
keep iterating up to `maxiter`, then report unmet criteria. This needs neither a
stall detector nor an automatic restart scheme.

Sources: [`resolve` and `decide`](../../src/madness/chem/ConvergenceSettings.h#L894),
[SCF floors](../../src/madness/chem/SCF.cc#L2985),
[driven floors](../../src/apps/molresponse/solvers/shared_fd_checkpoint.hpp#L60).
Nemo's distinct `0.1*sqrt(tau)` intermediate relaxation in `SCFProtocol` is
outside this audit; it must not be attributed to ordinary SCF or molresponse.

## Verification and the comparison still needed

Existing scalar tests cover the driven relative gates (`test_fd_gate`), the ES
density/omega test and separate stopping reasons (`test_es_plateau`), the legacy
engine/requested distinction (`test_convergence_settings`), and shared metric
assessment (`test_iteration_metrics`). These verify logic, not chemical accuracy
or equal runtime. This audit introduces no solver or numerical-reference edits.
All four tests were rerun successfully on November 19 in `build/simple-policy`;
the audit evidence is under
`.superpowers/sdd/design_draft/legacy-convergence-audit/` in this worktree.

For implementation validation, compare legacy and mapped shared requests using
the same molecule, model, seed, protocol/k, KAIN settings, and hardware/thread
configuration. Start with the agreed helium HF alpha plus DALTON-reference
water cases, include a multi-function case to expose the residual-norm change,
and include TDA/full roots when their shared adapter is adopted. Record each
protocol's resolved limits, measurements, stop reason, updates, time, and final
energy/property differences. Preserve existing numerical references. Similar
accuracy and useful iteration cost must be demonstrated, not inferred from
matching tolerance numbers. Broader release coverage remains in the test index.
