# Molresponse unification milestone

Status: integration implemented and functional verification complete,
November 20, 2026. Ten representative scientific cases and an explicit
virtual-AO water 2PA variant pass the retained numerical checks. The default
cold solid-harmonic water case selects higher roots and fails its retained
energy reference, so full milestone acceptance remains open. Guess defaults,
the no-warmup policy and reference data remain unchanged. See the
[verification record](madqc-test-index.md#response-verification--november-20-2026).

Implementation sequence and single-pass estimate: [full response plan](../superpowers/plans/2026-11-19-full-response-unification.md).

## Scope and order

1. Finish independent driven response, including the agreed SCF-to-response
   restart rules and per-state protocol progress.
2. Adopt coupled excited-state roots: per-root metrics, meaningful comparisons
   through mixing/reordering, collective stopping and bundle checkpointing.
3. Adopt quadratic-response tasks and properties, including VBC construction and
   any existing iterative quadratic-state paths. Cover the supported end-to-end
   combinations used by beta, two-photon absorption and Raman calculations.

Molresponse keeps ownership of its internal dependencies and numerical state.
The workflow owns the qualified ground-state handoff. This milestone covers the
existing supported response functionality; it does not add new scientific
methods to fill documented solver limitations.

## Preserve the numerical solver

This milestone connects existing components. Preserve the numerical kernels,
update formulas, operation ordering and optimized execution paths, except for
separately identified and verified bug fixes. Reuse available numerical work for
metrics and avoid duplicate solves, operator applications, unnecessary copies or
synchronization. The target is unchanged solver performance at equivalent
settings with minimal reporting/checkpoint bookkeeping overhead.

During integration, check interface wiring, persistence, records and collective
behavior. Perform numerical metric and representative converged-value checks after the response
paths are connected. Record iteration counts for analysis without making their
agreement with old runs an acceptance requirement. Reference recomputation and
detailed calibration belong to the later validation stage.

## Shared contract

- Use the shared metric helpers, records and convergence decisions for every
  iterative path. Report BSH residual, density, energy/excitation-energy and
  kinetic changes with their model-specific definitions and applicability.
  Retain the agreed intermediate heuristics and final requested criteria.
- Use the same simple restart rules and native checkpoint interface. Preserve
  independent-state progress; restart coupled roots as a bundle. Keep saved
  functions when compatible, use saved measurements only when valid, and reset
  unavailable comparison history when fresh checks are required.
- Direct tasks use the same checkpoint/completion interface and record their
  actual sources. Build VBC from qualified B/C states; evaluate a property only
  when every required state and direct input is ready. Recompute affected
  derived quantities when their inputs change.
- Use consistent input conventions, grouped property results, state metrics,
  partial-result availability, progress messages and aligned output. Numerical
  definitions remain application-owned and must be explained to the user.

## Shared timing conventions

Timing is a required part of this milestone across all response paths:

- The shared collector measures elapsed wall time with `std::chrono::steady_clock`
  and reports seconds. Retain existing timers as supplemental diagnostics.
- Keep the established driver-rank measurement scope explicit. Do not label
  elapsed time as summed CPU time or silently mix rank-local measurements with
  MPI maxima; do not add barriers merely to produce a timing row.
- Report application and internal-task totals, then iteration and named-step
  timings where iteration applies. Use common phase names/boundaries for shared
  operations and meaningful application-specific steps for distinct computations.
- Include preparation, evaluation, updates, checkpoint loading/publication,
  VBC construction and property evaluation when executed. Direct tasks have task
  and step timings without artificial solver iterations.
- Console and saved timing records use the same measurements and units. Clearly
  distinguish totals from their breakdowns; nested intervals are not extra work
  to sum into a total.
- Report the work performed by the current invocation. Reuse reports its actual
  load/check/property time; historical iteration timings remain separate.

## Completion evidence

Maintain a feature-to-test inventory covering all supported response properties,
static and frequency-dependent cases, one- and many-electron systems, open shells,
applicable XC functionals and supported correlation/kernel variants. Include
functional tests of metrics, stopping, restart, native archives, dependency
invalidation, grouped/partial results, inspection, printing and timings.

For this first milestone, verify that supported response paths use the common
interfaces and representative calculations converge reasonably to sensible
values with consistent metrics, saved functions and outputs. Existing references,
including the helium HF and DALTON-reference water alpha pair, remain useful
checks. Different iteration counts alone do not block acceptance; record them
without tuning the solver to reproduce them. Investigate clear numerical errors
and false completion.

The full scientific/model matrix, detailed performance/default calibration and
recomputation of MADNESS reference runs belong to a later validation stage.
Do not regenerate references for this milestone or replace independent external
oracles with new MADNESS output. Missing integration or representative evidence
remains a milestone gap; broader validation gaps remain release work. Record
existing scientific limitations separately.

Optimizer redesign, optimizer restart/history and new geometry-workflow behavior
remain deferred until workflow unification is complete. Broader engine adoption
remains release work; the outstanding response acceptance case is recorded above.
