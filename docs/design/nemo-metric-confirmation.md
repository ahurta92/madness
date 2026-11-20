# Nemo metric confirmation — Florian

**Pending scientific confirmation.** Please confirm or correct the ground-state
and CPHF definitions below, including units and supported occupation/spin cases.
Tolerance values will be validated during implementation.

Source checked on 2026-11-13 at `76255552e16e`. Nemo stores regularized functions
`F_i`; physical orbitals are `psi_i = R F_i`. See the [source evidence](madqc-metric-validation.md)
for equations and the longer audit.

## Ground-state metrics

- [ ] **Residual representation**

  **Code:** Forms `r_i = F_i - G_F(F)_i` before the nonlinear update and reduces
  the per-function L2 norms to an RMS. It applies no `R` multiplier or occupation
  weights. [Residual calculation](../../src/madness/chem/nemo.cc#L717).

  **Confirm:** The proposed shared metric is `max_i norm2(r_i)` in the
  regularized-function representation, without `R` or occupation weights.
  Is that the intended residual norm? The maximum reduction is already agreed;
  the scientific question is the representation and its interpretation.

- [ ] **Density and occupation weighting**

  **Code:** The iteration uses `R^2 sum_i abs(F_i)^2` with unit weights and
  reports its L2 change. For unit-occupation closed-shell orbitals this is one
  spin density. The property path instead uses
  `2 R^2 sum_i occupation_i abs(F_i)^2`.
  [Iteration](../../src/madness/chem/nemo.cc#L730),
  [property density](../../src/madness/chem/nemo.cc#L333).

  **Confirm:** Use physical per-spin density change. Which occupation/spin cases
  are supported, and what weights should each use? Flag any matching changes
  needed in energy/kinetic evaluation; adding weights only to a diagnostic would
  not establish consistent support.

- [ ] **Physical kinetic quantity**

  **Code:** `ke0`, returned as `energies[1]`, already includes closed-shell factors
  and the nuclear-factor correction. The correction uses a dedicated `U1dotU1`
  functor. A separate `kinetic_change` measurement/history is not yet implemented.
  [Kinetic quantity](../../src/madness/chem/nemo.cc#L779),
  [correction functor](../../src/madness/chem/nemo.h#L201).

  **Confirm:** Reuse `abs(ke0_current - ke0_previous)`, preserving those factors
  and correction. Which derivative/smoothing conventions must be recorded, and
  is the existing product-rule evaluator an appropriate validation reference?

## CPHF nuclear-response metrics

- [ ] **Residual and saved representation**

  **Code:** Forms `xi - G(xi)`, tests `RMS/norm(xi)`, and later subtracts a
  parallel-response term to form the final saved response.
  [Solve and correction](../../src/madness/chem/nemo.cc#L1678).

  **Confirm:** Is the shared maximum norm of the raw iterated regularized
  response residual appropriate, with what units and relative scale? Confirm
  the correction/sign convention and what must be stored to reconstruct that
  iterated unknown from a full response. The maximum reduction itself is agreed.

- [ ] **Perturbed density and kinetic quantity**

  **Code:** Builds perturbed density with occupation/closed-shell factors and
  nuclear-factor derivatives; it does not supply shared density/kinetic change
  histories. [Density construction](../../src/madness/chem/nemo.cc#L1631).

  **Confirm:** Which response density and kinetic quantity should be compared
  between iterations, with which weights, normalization and units? Confirm
  supported spin/occupation cases and missing evaluators. Ground-state total
  energy is not a CPHF response-change metric.

- [ ] **Hessian-row check and validation cases**

  **Code:** Requires a Hessian-row diagnostic change below `1e-2`, measured
  after updating, alongside the earlier residual. It also runs a preparatory
  pass labelled LDA. [Checks](../../src/madness/chem/nemo.cc#L1692),
  [preiterations](../../src/madness/chem/nemo.cc#L1791).

  **Confirm:** Should the row change remain a scientific acceptance check or
  become diagnostic? If required, define its sampling, units and tolerance.
  Which small reference cases validate corrected-response metrics and removal
  of the separate iterative guess pass? Warmup-off is already the design choice.

**Already agreed implementation work:** measure all metrics and checkpoint the
same prepared state before the next update; applications own comparison history.
Both the ground-state and CPHF loops need that refactor. CPHF tasks use shared
restart/outcomes and gate direct property assembly; their separate iterative
guess pass is disabled. These structural choices are already agreed.

**Confirmation record:** reviewer/date; confirmed items; corrections, supported
cases and exact units; any missing evaluator or validation requirement.
