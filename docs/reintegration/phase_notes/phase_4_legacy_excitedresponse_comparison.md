# Phase 4 Legacy `ExcitedResponse` Comparison

Date:

- 2026-03-16

Purpose:

- Provide a code-grounded comparison between the legacy excited-state solver in
  `src/apps/molresponse/ExcitedResponse.cpp` and the current Stage 2c bundle
  solver in `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp`.
- Identify the concrete numerical-parity gaps that still block Phase 4.

Scope:

- Restricted-shell excited-state flow first.
- Focus on `initialize`, `iterate_trial`, bundle rotation/diagonalization,
  `iterate`, `update_response`, and the BSH update path.
- This note does not propose a redesign. It records what the code actually does.

## Function Mapping

| Legacy path | Current path | Comparison |
| --- | --- | --- |
| `ExcitedResponse::initialize` | `RestartAwareExcitedScaffoldSolver::build_fresh_guess`, `initialize_protocol_guess` | Partial parity. Current code preserves the high-level structure, but not the same trial refinement algorithm. |
| `ExcitedResponse::iterate_trial` | `RestartAwareExcitedScaffoldSolver::iterate_trial` | Major mismatch. Legacy does a real TDA/BSH refinement loop; current code only reprojects, reorthogonalizes, and re-estimates heuristic energies. |
| `ExcitedResponse::rotate_excited_space` | `build_rotation_matrices`, `diagonalize_excited_bundle`, `rotate_excited_bundle_states` | Partial parity. Current code mirrors most of `excited_eig`, but not legacy `reduce_subspace`. |
| `ExcitedResponse::iterate` | `RestartAwareExcitedScaffoldSolver::iterate` | Partial parity. Current code now has the legacy-style threshold-derived convergence contract, but control-flow and fallback behavior still differ. |
| `ExcitedResponse::update_response` | `iterate_typed_bundle_legacy_sequence`, `iterate_state_from_potentials` | Partial parity. Current code preserves bundle rotation followed by per-root updates, but the update kernel is not mathematically identical. |
| `ExcitedResponse::bsh_update_excited` | `bsh_update_from_theta` | Close structurally on the restricted path. Shift construction, `-2` scaling, BSH application, projection, and normalization are present. |

## 1. Initialization And Trial-Space Construction

Legacy facts:

- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::initialize`
  allocates `trial` with `2 * r_params.num_states()` and selects the guess
  source from random, NWChem, derivative-based, or projected AO guess paths.
- It explicitly load-balances the trial functions, projects out the ground
  space with `QProjector`, runs two rounds of Gram-Schmidt plus normalization,
  calls `iterate_trial`, sorts by `omega`, then keeps the lowest
  `num_states` in `Chi.x` and initializes `Chi.y` to zero.
- It saves a hard-coded `guess_restart` immediately after guess construction.

Current facts:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_fresh_guess`
  also starts from an oversized trial bundle and then calls
  `iterate_trial`, `select_lowest_trial_roots`, and
  `seed_active_bundle_from_trial_space`.
- The current fresh-start logic resets root descriptors, names, slot
  permutation, and the active typed response bundle after selecting the trial
  roots.
- `initialize_protocol_guess` wraps this with restart precedence:
  current protocol snapshot, lower protocol snapshot, guess archive, carryover,
  then fresh guess.

Comparison:

- The control-flow role is similar: construct too many candidates, refine them,
  sort them, then seed the active excited-state container.
- The container model is intentionally different. Legacy writes directly into
  `Chi` as an `X_space`; current code writes into `trial_space_`,
  `omega_`, and then materializes typed `ResponseVector` objects.
- This stage is not the main Phase 4 blocker. The major mismatch starts inside
  `iterate_trial`.

## 2. Trial Refinement (`iterate_trial`)

Legacy facts:

- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate_trial`
  is not just an energy-sort loop.
- Each guess iteration does:
  - transition density construction via `transition_densityTDA`
  - projection and truncation
  - normalization for TDA
  - `compute_response_potentials(..., "tda")`
  - `rotate_excited_space(...)`
  - `theta_X = rotated_v_x - E0X + rotated_gamma_x`
  - `create_shift`, `apply_shift`, `create_bsh_operators`, `apply`
  - projection, masking, then two more Gram-Schmidt/normalize passes
- That means the trial space is already being pushed through a simplified
  excited-state update kernel before the main solve starts.

Current facts:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_trial`
  only:
  - projects and orthonormalizes the trial states
  - estimates energies with `estimate_state_energies`
  - damps the energy estimates
  - sorts by energy
- No response potentials, rotation, `theta`, shift, or BSH refinement appears in
  the current trial refinement path.

Comparison:

- This is a first-order parity gap.
- Legacy `iterate_trial` produces a numerically informed initial bundle.
  Current Stage 2c still starts from a much weaker trial space even after
  Phase 3.
- This is a likely contributor to the poor early-iteration behavior, but it is
  still distinct from the more serious mismatch in the main update kernel.

## 3. Potential Builder: Legacy `Lambda_X` vs Current `lambda`

Legacy facts:

- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::compute_response_potentials`
  explicitly builds:
  - `T0X = T(world, chi_copy)`
  - `E0X = chi_copy * hamiltonian`
  - `V0X = compute_V0X(...)`
  - `gamma = compute_gamma_*`
  - `Lambda_X = (T0X + V0X - E0X) + gamma`
- This is the quantity passed into `rotate_excited_space`.

Current facts:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: compute_state_response_potentials`
  builds:
  - `k0` and `gx` from exchange tasks
  - `v_local = ground_data_->V_local * current.flat`
  - `v0_flat = v_local - c_xc * k0`
  - `e0_flat = apply_hamiltonian_no_diag(world, current)`
  - `lambda_flat = v0_flat - e0_flat + gx`
- `src/apps/molresponse_v2/GroundStateData.cpp :: GroundStateData::computePreliminaries`
  shows `V_local = V_nuc + V_coul (+ V_xc)`; it does not contain the kinetic
  term.

Comparison:

- This is the strongest currently visible mathematical mismatch.
- Legacy `Lambda_X` includes an explicit `T0X` term before the subtraction of
  `E0X`.
- Current `lambda_flat` omits that explicit kinetic term and therefore is not
  the same quantity that legacy diagonalizes.
- Because bundle rotation is driven by `A = <chi|lambda>`, this mismatch affects
  every rotated root, not just the post-rotation BSH update.

Implementation target implied by the code:

- Reconstruct the legacy `Lambda_X = (T0X + V0X - E0X) + gamma` quantity in the
  current restricted bundle path before touching later convergence heuristics.

## 4. Bundle Rotation And Generalized Eigenproblem

Legacy facts:

- `src/apps/molresponse/ExcitedResponse.cpp :: rotate_excited_space` builds
  `S = <X|X> - <Y|Y>` and `A = inner(chi, lambda)`, then calls `excited_eig`.
- `src/apps/molresponse/ExcitedResponse.cpp :: reduce_subspace` performs an SVD
  of `S`, counts small singular values, and if needed transforms both `S` and
  `A` into the reduced active singular subspace before diagonalization.
- `src/apps/molresponse/ExcitedResponse.cpp :: excited_eig` then performs:
  - `sygvp`
  - diagonal-dominance swap heuristic
  - phase fixing
  - degenerate-cluster cleanup by SVD/polar-style rotation
  - final ascending energy sort

Current facts:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: build_rotation_matrices`
  reproduces the same broad structure for `S` and `A`.
- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: diagonalize_excited_bundle`
  reproduces most of legacy `excited_eig`:
  - `sygvp`
  - diagonal-dominance swap heuristic
  - phase fixing
  - near-degenerate cluster cleanup
  - final energy sort
- Phase 4 added `condition_excited_overlap_matrix`, which retries `sygvp` after
  SVD-conditioning the overlap matrix.

Comparison:

- The current code now mirrors most of legacy `excited_eig` behavior.
- The remaining gap is legacy `reduce_subspace`.
- Current code floors singular values and retries in the full space. Legacy
  actually removes near-null overlap modes and diagonalizes in a reduced
  subspace.
- This explains why the Phase 4 `condition_overlap_svd` retry helped runtime
  robustness, but has not fully restored legacy rotation behavior.

Implementation target implied by the code:

- Replace or extend overlap conditioning so the current solver can emulate
  legacy `reduce_subspace` semantics instead of only flooring singular values.

## 5. Main Update Loop

Legacy facts:

- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate`
  computes:
  - `conv_den = max(100 * thresh, dconv)`
  - `relative_max_target = max(50 * thresh, 0.5 * dconv)`
  - `max_rotation` from the protocol threshold
- The loop then:
  - rotates and updates through `update_response`
  - computes old and new densities from `old_chi` and `new_chi`
  - sets `bsh_residualsX = new_res.residual_norms`
  - updates `omega` and `Chi`
  - sets `density_residuals = norm2s_T(world, rho_omega - rho_omega_old)`
- Convergence is checked at the top of the next iteration, not immediately
  after the update.

Current facts:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: make_iteration_contract`
  reproduces the same threshold-derived density target, relative target, and
  `max_rotation`.
- `iterate(...)` now logs and persists this contract in metadata.
- For restricted variants the loop calls
  `iterate_typed_bundle_legacy_sequence`, then computes
  `max_residual`, `max_density_change`, and `max_relative_residual`, and checks
  the dual gate at the bottom of the loop.

Comparison:

- The convergence contract is now intentionally aligned with legacy.
- The control-flow is still not identical:
  - current convergence is evaluated at the end of the same iteration
  - current code has a non-legacy fallback when bundle rotation fails:
    `omega_ += 0.45 * (target_omega - omega_)`
- That fallback should be treated as scaffolding, not parity behavior.

## 6. Per-Root Update Kernel

Legacy facts:

- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::update_response`
  computes bundle potentials, rotates the bundle, builds
  `theta_X = rotated_v_x - rotated_EOX + rotated_gamma_x`, runs
  `bsh_update_excited`, computes residuals with `update_residual`, optionally
  applies KAIN, and leaves step restriction disabled with `if (false)`.
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::bsh_update_excited`
  then:
  - builds orbital shifts with `create_shift`
  - applies `theta_X.x = apply_shift(...)`
  - multiplies by `-2`
  - applies BSH operators
  - projects out the ground state

Current facts:

- `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_typed_bundle_legacy_sequence`
  computes bundle potentials, rotates the bundle, and then updates each root by
  calling `iterate_state_from_potentials`.
- `iterate_state_from_potentials` builds
  `theta = v0_state.flat - rotated_e0_flat + gamma_state.flat`,
  calls `bsh_update_from_theta`, optionally uses the per-root accelerator,
  then projects, normalizes, and evaluates density and residual metrics.
- `bsh_update_from_theta` preserves the same basic BSH update structure as
  legacy: shifts, `-2`, apply BSH, project, truncate, normalize.
- Phase 4 now leaves explicit step restriction off by default for restricted
  variants, matching the legacy `if (false)` behavior.

Comparison:

- This is structurally close, but not yet numerically equivalent.
- The largest direct mismatch in this section is upstream: the current
  `v0_state` / `gamma_state` / `rotated_e0_flat` are built from a different
  `lambda` decomposition than legacy.
- The current per-root accelerator/damping path is also still more scaffold-like
  than legacy KAIN history management.

## 7. Convergence Metrics

Legacy facts:

- Legacy convergence on the excited path is effectively:
  `max_density_change <= conv_den` and
  `max_relative_bsh <= relative_max_target`,
  unless `conv_only_dens` overrides the second gate.

Current facts:

- Restricted Stage 2c now uses the same dual-gate structure through
  `IterationContract`.
- Unrestricted variants still use scaffold/fallback convergence semantics.

Comparison:

- The convergence rule itself is no longer the primary parity gap on the
  restricted path.
- The current failure mode is that the translated kernel is not driving those
  metrics down far enough, not that the wrong gate is being checked.

## Priority Gaps For The Next Kernel Pass

1. Rebuild the legacy `Lambda_X` quantity on the current restricted path.
   Evidence:
   `src/apps/molresponse/ResponseBase.cpp :: compute_response_potentials`
   vs
   `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: compute_state_response_potentials`

2. Replace overlap-flooring-only logic with a true legacy-style reduced-subspace
   diagonalization path.
   Evidence:
   `src/apps/molresponse/ExcitedResponse.cpp :: reduce_subspace`
   vs
   `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: condition_excited_overlap_matrix`

3. After fixing 1 and 2, reassess whether the remaining Phase 4 instability is
   mostly in trial refinement or in the per-root accelerator/update sequence.
   Evidence:
   `src/apps/molresponse/ExcitedResponse.cpp :: iterate_trial`
   vs
   `src/apps/molresponse_v2/ExcitedStateBundleSolver.cpp :: iterate_trial`

## Explicit Facts vs Inference

Facts directly visible in code:

- Legacy `Lambda_X` includes `T0X`; current `lambda_flat` does not.
- Legacy guess refinement applies BSH updates; current guess refinement does not.
- Legacy diagonalization has an explicit subspace-reduction helper;
  current code does not.
- Legacy main excited update leaves step restriction disabled;
  current restricted path now does the same.

Inference:

- The missing `T0X` term is likely a major cause of poor dynamic-restricted
  convergence because it changes the bundle-rotation generalized eigenproblem.
- The weaker current guess refinement likely worsens early-iteration behavior,
  but it appears secondary to the `Lambda_X` mismatch.
