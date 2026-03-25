# Legacy Excited-State Implementation Report

Line numbers in this report refer to the current checkout in `/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next` as inspected on 2026-03-13.

## 1. Purpose and scope

This report documents the legacy excited-state implementation in the old `src/apps/molresponse` codebase, centered on the compiled path:

- `src/apps/molresponse/molresponse.cc :: main`
- `src/apps/molresponse/ExcitedResponse.hpp/.cpp :: ExcitedResponse`
- `src/apps/molresponse/ResponseBase.hpp/.cpp :: ResponseBase`
- `src/apps/molresponse/x_space.h/.cc :: X_space`
- `src/apps/molresponse/response_functions.h :: response_space`

What this report covers:

- The mathematical and algorithmic role of the old excited-state solver.
- The actual control flow from input parsing through initialization, iteration, convergence, and output.
- The class model and data containers used to represent excited states.
- How response vectors are created, stored, rotated, updated, and persisted.
- Legacy assumptions that will matter when reintegrating this behavior into `molresponse_v2`.

What this report does not cover:

- `molresponse_v2`.
- Frequency-response implementation details except where they expose shared legacy infrastructure that the excited-state code depends on.
- A redesign or migration plan.

Scope note:

- `src/apps/molresponse/calc_runner.cc` contains an older `TDDFT::solve_excited_states` implementation, but it is not part of the current build for `molresponse`; `src/apps/molresponse/CMakeLists.txt` includes `ExcitedResponse.cpp` and not `calc_runner.cc` (`src/apps/molresponse/CMakeLists.txt` lines 12-23, 42).

## 2. High-level summary

Fact:

- The current old `molresponse` excited-state solver is a bundle-based iterative solver. All requested excited states live together in one `X_space` object (`src/apps/molresponse/ResponseBase.hpp` line 715, `src/apps/molresponse/x_space.h` lines 33-458).
- Initialization creates an oversized trial space, improves it through a separate guess iteration, sorts by estimated excitation energies, and then selects the lowest `num_states` rows into the main state bundle `Chi` (`src/apps/molresponse/ExcitedResponse.cpp` lines 7-115, 479-677).
- The main iteration repeatedly:
  - builds operator-applied vectors `Lambda_X`, `V0X`, and `gamma`,
  - diagonalizes the current state bundle in a generalized subspace eigenproblem to obtain updated `omega` and a rotated basis,
  - forms a BSH right-hand side `theta_X`,
  - applies state-specific BSH Green's operators,
  - computes residuals and density changes,
  - optionally applies KAIN updates,
  - tests convergence (`src/apps/molresponse/ExcitedResponse.cpp` lines 2032-2477).

Inference supported by code comments and operator split:

- The code is solving a TDDFT/TDHF excited-state response fixed-point problem where the kinetic-energy part is handled through the BSH Green's operator, while the non-kinetic remainder is assembled into `theta_X`. This is consistent with the explicit comment in `ExcitedResponse.hpp`:
  - `chi^m = -2 G * Theta chi` (`src/apps/molresponse/ExcitedResponse.hpp` lines 106-120)
  - and with the implementation split
    - `Lambda_X = (T0X + V0X - E0X) + gamma` for subspace diagonalization (`src/apps/molresponse/ResponseBase.cpp` lines 1116-1171)
    - `theta_X = V0X - E0X + gamma` for the BSH update (`src/apps/molresponse/ExcitedResponse.cpp` lines 2420-2453).

Main architectural ideas:

- States are not standalone objects. One excited state is a row slice of a shared `X_space`.
- State identities are positional, not object-based. They are repeatedly reordered and rotated by eigenvector matrices.
- Orthogonalization, diagonalization, normalization, and persistence all operate on the whole bundle, not on independent state objects.
- The solver uses a mix of:
  - generalized eigenvalue solves in the current subspace,
  - direct Green's-function / BSH updates,
  - optional KAIN acceleration on the response vectors.

Important legacy observations:

- The excited-state path does not contain a `response_context.set_strategy(...)` call, even though `compute_gamma_full()` and `compute_gamma_static()` depend on `response_context` strategies (`src/apps/molresponse/ExcitedResponse.hpp` line 19, `src/apps/molresponse/ResponseBase.cpp` lines 687-799, `src/apps/molresponse/FrequencyResponse.hpp` lines 564-583). This is a concrete dependency gap in the old design.
- Several shared `ResponseBase` helpers decide whether `Y` is present by checking `r_params.omega() != 0.0` instead of `calc_type` or `tda` (`src/apps/molresponse/ResponseBase.cpp` lines 408-419, 1466, 1515, 1582, 1619). In an excited-state run, `ResponseParameters::omega` is a frequency-response input and is normally `0.0` (`src/apps/molresponse/response_parameters.h` line 97). This makes the full excited-state `Y` path internally inconsistent.

## 3. Relevant files and components

| File | Role in excited-state path | Key symbols |
| --- | --- | --- |
| `src/apps/molresponse/CMakeLists.txt` | Defines the compiled old `molresponse` target; confirms `ExcitedResponse.cpp` is built and `calc_runner.cc` is not. | `MOLRESPONSE_SOURCES`, `add_mad_executable(molresponse ...)` |
| `src/apps/molresponse/molresponse.cc` | User-facing driver. Creates `ExcitedResponse` when `response_parameters.excited_state()` is true and calls `solve()`. | `main` |
| `src/apps/molresponse/global_functions.h` | Defines `CalcParams` and declares shared operator helpers. | `CalcParams`, `initialize_calc_params`, `T`, `newK` |
| `src/apps/molresponse/global_functions.cc` | Parses input into `CalcParams`; provides kinetic operator helper `T`; provides `newK`. | `initialize_calc_params`, `T`, `newK` |
| `src/apps/molresponse/response_parameters.h` | Defines all response/excited-state input parameters and derives `calc_type`. | `ResponseParameters`, `set_derived_values` |
| `src/apps/molresponse/response_functions.h` | Defines `response_space`, the per-channel state matrix used inside `X_space`. | `response_space`, `response_space_inner` |
| `src/apps/molresponse/x_space.h` | Defines `X_space`, the main excited-state bundle container, plus flattening and transformation utilities. | `X_space`, `X_vector`, `to_vector`, `from_vector`, operators |
| `src/apps/molresponse/x_space.cc` | Implements conversions between `X_space`, flattened vectors, and concatenated per-state response matrices. | `to_response_matrix`, `to_X_space`, `inner` |
| `src/apps/molresponse/ResponseBase.hpp` | Declares the common response solver base class and strategy interfaces used by both excited and frequency solvers. | `ResponseBase`, `Context`, `residuals`, `save_x_space`, `load_x_space` |
| `src/apps/molresponse/ResponseBase.cpp` | Implements shared setup, Hamiltonian construction, ground-state operators, density construction, gamma/V0/Lambda assembly, residuals, KAIN, transforms, normalization, sorting, and protocol loop. | `ResponseBase::solve`, `compute_response_potentials`, `compute_V0X`, `compute_gamma_*`, `update_residual`, `kain_x_space_update`, `make_density`, `sort`, `transform` |
| `src/apps/molresponse/ExcitedResponse.hpp` | Declares the excited-state solver class and all solver-specific helpers. | `ExcitedResponse`, `omega`, `initialize`, `iterate`, `update_response`, `bsh_update_excited` |
| `src/apps/molresponse/ExcitedResponse.cpp` | Implements guess generation, guess refinement, main excited-state iteration, subspace rotation/eigenvalue extraction, BSH updates, save/load, and optional analysis. | `ExcitedResponse::initialize`, `iterate_trial`, `iterate`, `update_response`, `rotate_excited_space`, `excited_eig`, `save`, `load` |
| `src/apps/molresponse/FrequencyResponse.hpp` | Not part of the excited-state control flow, but important for understanding shared infrastructure because it is the only visible place that initializes `response_context` strategies. | `FrequencyResponse` constructor |
| `src/apps/molresponse/calc_runner.cc` | Historical precursor containing a separate `TDDFT::solve_excited_states` implementation. Useful for lineage, but not part of the current build. | `TDDFT::solve_excited_states`, `TDDFT::create_random_guess` |

## 4. Top-level execution flow

### 4.1 Entry point and parameter setup

1. `src/apps/molresponse/molresponse.cc :: main` starts MADNESS, parses command-line flags, and reads `response.in` through `initialize_calc_params(world, "response.in")` (`molresponse.cc` lines 60-101, `global_functions.cc` lines 14-25).
2. `initialize_calc_params()`:
   - constructs `ResponseParameters`,
   - reads the response input block,
   - opens the referenced ground-state archive via `GroundStateCalculation`,
   - extracts the molecule,
   - derives response-specific values such as `num_orbitals`, `spinrestricted`, `L`, `lo`, `xc`, and `calc_type` (`global_functions.cc` lines 14-25, `response_parameters.h` lines 161-188).
3. If `response_parameters.excited_state()` is true, `main` constructs `ExcitedResponse calc(world, calc_params)` and calls `calc.solve(world)` (`molresponse.cc` lines 84-102).

### 4.2 Object construction

`ExcitedResponse` does not add construction logic of its own; its constructor forwards to `ResponseBase` (`src/apps/molresponse/ExcitedResponse.hpp` line 19).

`ResponseBase::ResponseBase` does the common setup (`src/apps/molresponse/ResponseBase.cpp` lines 12-29):

- copies `ResponseParameters`, `Molecule`, and `GroundStateCalculation`,
- stores ground orbitals and ground energies,
- allocates `Chi` as `X_space(world, num_states, num_orbitals)`,
- initializes the XC functional object,
- sets global function defaults for the simulation box and truncation mode.

### 4.3 Protocol loop

`ResponseBase::solve` is the top-level solver loop (`src/apps/molresponse/ResponseBase.cpp` lines 1741-1780):

1. Read the protocol list `r_params.protocol()`.
2. For each truncation threshold:
   - call `set_protocol(world, iter_thresh)` to choose `k`, set MADNESS thresholds, create Coulomb/gradient operators, build nuclear potential, mask, and ground density (`ResponseBase.hpp` lines 722-772),
   - on the first protocol:
     - if `restart` is enabled, call `load(world, restart_file)`,
     - otherwise call `ExcitedResponse::initialize(world)`,
   - call `check_k(world, iter_thresh, k)` to reproject stored functions if the current polynomial order changed (`ResponseBase.cpp` lines 38-132),
   - create the JSON node for this protocol via `protocol_to_json(...)`,
   - call `ExcitedResponse::iterate(world)`.
3. After all protocols complete, write `converged` into the JSON tree.

Compact call flow:

```text
main
  -> initialize_calc_params
  -> ExcitedResponse::ExcitedResponse
     -> ResponseBase::ResponseBase
  -> ResponseBase::solve
     -> ResponseBase::set_protocol
     -> ExcitedResponse::initialize   [first protocol if not restart]
     -> ExcitedResponse::load         [first protocol if restart]
     -> ResponseBase::check_k
     -> ExcitedResponse::iterate
  -> ResponseBase::output_json
```

### 4.4 Initialization flow

`ExcitedResponse::initialize` (`src/apps/molresponse/ExcitedResponse.cpp` lines 7-115) performs the initial excited-state guess construction:

1. Allocate `trial` with `2 * num_states` rows and `num_orbitals` columns.
2. Choose one guess generator:
   - random guesses: `make_random_trial()` (`ExcitedResponse.cpp` lines 154-184)
   - NWChem-based guesses: `make_nwchem_trial()` (`ExcitedResponse.cpp` lines 187-312)
   - Cartesian `x/y/z` guesses: `create_trial_functions2()` (`ExcitedResponse.cpp` lines 386-476)
   - otherwise a virtual-AO-based guess via `create_virtual_ao_guess()` (`ExcitedResponse.cpp` lines 2832-2984)
3. Optionally redistribute the trial trees and ground orbitals for load balance.
4. Project each trial state's `x` row out of the occupied ground-state space using `QProjector`.
5. Apply two rounds of `gram_schmidt(world, trial.x)` and `normalize(world, trial.x)` (`ResponseBase.cpp` lines 2366-2396, 1844-1856).
6. Improve the trial subspace by calling `iterate_trial(world, trial)`.
7. Sort the refined trial states by `omega` with `sort(world, omega, trial.x)` (`ResponseBase.cpp` lines 2302-2330).
8. Select the lowest `num_states` trial rows into `Chi.x` via `select_functions(...)`; initialize `Chi.y` to zeros (`ExcitedResponse.cpp` lines 105-111, `ResponseBase.cpp` lines 2267-2300).
9. Save a hard-coded restart snapshot `guess_restart` regardless of the configured `save_file`.

### 4.5 Trial-space refinement flow

`ExcitedResponse::iterate_trial` (`src/apps/molresponse/ExcitedResponse.cpp` lines 479-677) is a simplified pre-solver used only during initialization:

1. Start from the current `guesses.x`.
2. In each guess iteration:
   - optionally load balance the guesses,
   - build TDA transition densities `rho_omega = transition_densityTDA(world, ground_orbitals, guesses.x)` (`ResponseBase.cpp` lines 2153-2176),
   - reproject out the occupied space and truncate,
   - normalize `guesses.x`,
   - call `compute_response_potentials(world, guesses, xc, "tda")` to obtain `Lambda_X`, `V0X`, and `gamma` in TDA mode,
   - call `rotate_excited_space(...)` to diagonalize the current trial bundle and obtain `omega` plus rotated vectors,
   - if this is not the last guess iteration:
     - build `theta_X = rotated_v_x - E0X + rotated_gamma_x`,
     - build energy shifts with `create_shift(...)`,
     - build per-state/per-orbital BSH operators with `create_bsh_operators(...)`,
     - apply them to get new trial vectors,
     - project them and multiply by the boundary mask.
   - run Gram-Schmidt and normalization again.
3. Exit with improved trial states and trial-space `omega`.

### 4.6 Main excited-state iteration flow

`ExcitedResponse::iterate` (`src/apps/molresponse/ExcitedResponse.cpp` lines 2032-2373) is the production excited-state solve:

1. Create convergence targets from the current threshold and `dconv`.
2. Create storage:
   - `rho_omega` for per-state densities,
   - `residuals`, `old_Chi`, `old_Lambda_X`,
   - KAIN solvers `kain_x_space`,
   - arrays for residual norms and density changes.
3. For each iteration:
   - optionally load balance `Chi`,
   - if `iter > 0`, test convergence using:
     - density residual max `< conv_den`
     - relative BSH residual max `< relative_max_target`
     - or `conv_only_dens`
   - if converged or at max iteration, optionally save and break.
   - otherwise call `update_response(...)`.
4. `update_response(...)` returns:
   - `new_omega`
   - `old_chi` (actually the rotated previous bundle)
   - `new_chi`
   - `new_res` containing both residual vectors and residual norms
5. Recompute densities from `old_chi` and `new_chi`.
6. Update:
   - `bsh_residualsX`, `bsh_residualsY`
   - `omega`
   - `Chi`
   - `density_residuals`
7. Accumulate timing and JSON iteration data.
8. After exit, print final energies and residual summaries.

### 4.7 Finalization and output

Finalization/output is split across multiple functions:

- `ExcitedResponse::save` writes a restart bundle (`ExcitedResponse.cpp` lines 2746-2782).
- `ExcitedResponse::load` reads a restart bundle (`ExcitedResponse.cpp` lines 2785-2825).
- `ExcitedResponse::excited_to_json` appends per-iteration `omega` to JSON, but stores it under `protocol_data[*]["property_data"]` (`ExcitedResponse.cpp` lines 2020-2030).
- `ResponseBase::function_data_to_json` writes norms and residuals per iteration (`ResponseBase.cpp` lines 1680-1697).
- `ResponseBase::output_json` writes the final JSON document to a fixed file `response_base.json` (`ResponseBase.cpp` lines 2125-2147).
- `ExcitedResponse::analysis` can compute dipoles, oscillator strengths, quadrupoles, and dominant orbital contributions, but the call is commented out in `iterate` (`ExcitedResponse.cpp` lines 2357-2371, 2557-2743).

## 5. Core class and data model

### 5.1 `response_space`: one channel of response vectors

`response_space` is defined in `src/apps/molresponse/response_functions.h` lines 26-636.

Fact:

- `response_space` stores a 2D matrix `x[state][orbital]` of real-space functions.
- `num_states` is the number of response states.
- `num_orbitals` is the number of occupied ground-state orbitals represented inside each state vector.
- `active` is a list of row indices used to mark which states are currently active.

Interpretation:

- A single row `response_space::x[b]` is one vector over occupied orbitals.
- In TDA mode, that row is the actual excited-state vector for state `b`.

Important operations:

- `to_vector()` and `from_vector()` flatten/unflatten the active rows into a contiguous `vector_real_function_3d` in state-major order (`response_functions.h` lines 598-635).
- `response_space_inner(a, b)` forms the matrix of row-wise orbital inner products (`response_functions.h` lines 559-596).
- `transform(world, response_space, U)` rotates the state index by a dense coefficient matrix `U` (`ResponseBase.cpp` lines 2178-2203).

### 5.2 `X_space`: paired `X`/`Y` state bundle

`X_space` is defined in `src/apps/molresponse/x_space.h` lines 33-458.

Fact:

- `X_space` owns:
  - `response_space x`
  - `response_space y`
  - `n_states`
  - `n_orbitals`
  - `active`
- `num_states()` and `num_orbitals()` describe the bundle dimensions.
- `reset_active()` and `set_active()` synchronize `active`, `x.active`, and `y.active` (`x_space.h` lines 48-69).

Representation of one excited state:

- One excited state is row `b` of the shared bundle:
  - `chi.x[b]` is the `X`-channel occupied-orbital vector.
  - `chi.y[b]` is the `Y`-channel occupied-orbital vector.
- There is no standalone "ExcitedState" class in the compiled path.
- The only dedicated `ExcitedSpace` struct in `ExcitedResponse.hpp` just wraps `chi` and `l_chi`, and no call site uses it (`ExcitedResponse.hpp` lines 9-13; `rg` search found no references beyond the definition).

Flattening and matrix views:

- `X_space::to_vector()` flattens active states in this order:
  - all `x` orbital components of active state 0,
  - then all `y` orbital components of active state 0,
  - then the next active state, and so on (`x_space.h` lines 422-457).
- `to_response_matrix(x)` converts each state into one concatenated row `[x_orbitals..., y_orbitals...]` (`x_space.cc` lines 56-69).
- `to_conjugate_response_matrix(x)` swaps the ordering to `[y..., x...]` (`x_space.cc` lines 71-85).

This is important because the solver alternates between:

- `X_space` for physically meaningful grouped storage,
- `response_space` for channel-local operations,
- `response_matrix` and flattened vectors for KAIN and state-space linear algebra.

### 5.3 State ownership and lifecycle

Ownership in the active solver:

- `ResponseBase::Chi` is the authoritative bundle of current excited states (`ResponseBase.hpp` line 715).
- `ExcitedResponse::omega` is the authoritative tensor of current excitation energies (`ExcitedResponse.hpp` line 25).
- Temporary bundles such as `trial`, `rotated_chi`, `temp_Lambda_X`, `new_chi`, `old_Chi`, and `residuals` are created as whole-bundle `X_space` objects.

Lifecycle:

1. `Chi` is allocated in `ResponseBase` construction.
2. `initialize()` fills `Chi.x` from the refined guess subspace and zeros `Chi.y`.
3. `iterate()` repeatedly replaces the entire `Chi` bundle with a new rotated/updated bundle.
4. `save()` and `load()` persist and restore the whole bundle, not individual states.

### 5.4 Role of `X_space` in solver coupling

`X_space` is the central coupling point in the legacy design:

- state reordering and state mixing happen by transforming the entire bundle with dense matrices (`ResponseBase.cpp` lines 2205-2213),
- sorting is done on the whole bundle (`ResponseBase.cpp` lines 2334-2364),
- density construction operates on the whole bundle (`ResponseBase.cpp` lines 368-395),
- convergence is decided from bundle-wide maxima (`ExcitedResponse.cpp` lines 2149-2238),
- save/load snapshots operate on the bundle (`ExcitedResponse.cpp` lines 2746-2825).

This is the main architectural reason the old code is bundle-centric rather than state-object-centric.

### 5.5 `response_context` dependency

`ResponseBase` owns `Context response_context` (`ResponseBase.hpp` line 634), which dispatches:

- inner-product strategy,
- Coulomb/exchange/XC response-kernel strategies,
- density strategy,
- load strategy (`ResponseBase.hpp` lines 429-541).

Fact:

- `FrequencyResponse` explicitly initializes these strategies in its constructor (`FrequencyResponse.hpp` lines 564-583).
- `ExcitedResponse` does not.

Inference:

- Any excited-state path that reaches `compute_gamma_full()` or `compute_gamma_static()` depends on infrastructure that is visibly initialized only for the frequency-response solver. That is a concrete legacy coupling/incompleteness issue.

## 6. Function-by-function analysis

### `main`

- Location: `src/apps/molresponse/molresponse.cc :: main` (line 60)
- Purpose: top-level executable entry point.
- Inputs: `argc`, `argv`.
- Outputs: process exit code.
- Side effects: initializes MADNESS, reads input, constructs solver objects, triggers calculation, writes JSON.
- Called by: program entry.
- Calls:
  - `initialize_calc_params`
  - `ExcitedResponse` constructor
  - `ResponseBase::solve`
  - `ResponseBase::output_json`
- Algorithmic role: selects the excited-state solver based on input flags.
- Notes / assumptions:
  - Uses `response_parameters.excited_state()` as the mode switch (`molresponse.cc` line 85).

### `initialize_calc_params`

- Location: `src/apps/molresponse/global_functions.cc :: initialize_calc_params` (line 14)
- Purpose: build the combined `CalcParams` object used by `ExcitedResponse`.
- Inputs: `World&`, input filename.
- Outputs: `CalcParams { ground_calculation, molecule, response_parameters }`.
- Side effects: reads the response input block and the ground-state archive.
- Called by: `main`.
- Calls:
  - `ResponseParameters::read_input_and_commandline_options`
  - `GroundStateCalculation` constructor
  - `ResponseParameters::set_ground_state_calculation_data`
  - `ResponseParameters::set_derived_values`
- Algorithmic role: converts user input into a solver-ready parameter bundle.
- Notes / assumptions:
  - Sets `calc_type` to `"tda"` if `excited_state && tda`, otherwise `"full"` for excited-state runs (`response_parameters.h` lines 175-180).

### `ResponseBase::solve`

- Location: `src/apps/molresponse/ResponseBase.cpp` line 1741
- Purpose: protocol-level driver shared by all legacy response solvers.
- Inputs: `World&`.
- Outputs: none.
- Side effects:
  - mutates global function defaults,
  - mutates `Chi`,
  - mutates JSON/timing state,
  - may load restart data.
- Called by: `main`.
- Calls:
  - `set_protocol`
  - `initialize` or `load`
  - `check_k`
  - `protocol_to_json`
  - virtual `iterate`
- Algorithmic role: orchestrates multi-threshold protocol execution.
- Notes / assumptions:
  - Only the first protocol performs `initialize()` or restart `load()`.
  - Later protocols reuse the existing `Chi` bundle and just reproject via `check_k`.

### `ExcitedResponse::initialize`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 7
- Purpose: build the initial excited-state bundle `Chi`.
- Inputs: `World&`.
- Outputs: none directly; writes into `Chi` and `omega`.
- Side effects:
  - allocates and mutates a temporary `trial` bundle,
  - mutates `Chi.x`, `Chi.y`, and `omega`,
  - writes `guess_restart`.
- Called by: `ResponseBase::solve`.
- Calls:
  - `make_random_trial`
  - `make_nwchem_trial`
  - `create_trial_functions2`
  - `create_virtual_ao_guess`
  - `iterate_trial`
  - `sort`
  - `select_functions`
  - `save`
- Algorithmic role: create an oversized trial subspace, refine it, and pick the lowest roots as the starting bundle.
- Notes / assumptions:
  - Always starts with `2 * num_states` trial rows unless the AO path yields fewer.
  - `Chi.y` is initialized to zeros even for full mode; the `Y` channel is generated later in the main iteration.
  - The restart snapshot name is hard-coded to `guess_restart`.

### `ExcitedResponse::make_random_trial`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 154
- Purpose: build random, localized, occupied-space-projected `X` guesses.
- Inputs: `World&`, number of trial states `m`.
- Outputs: `X_space` with populated `x` and zero `y`.
- Side effects: none outside the returned object.
- Called by: `initialize`, `make_nwchem_trial`.
- Calls:
  - `add_randomness`
  - `QProjector`
  - `normalize`
- Algorithmic role: generic guess fallback.
- Notes / assumptions:
  - Adds high-magnitude random noise, masks it, then localizes it around atoms using atom-centered Gaussians.

### `ExcitedResponse::make_nwchem_trial`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 187
- Purpose: build guesses from external NWChem virtual orbitals.
- Inputs: `World&`, number of trial states `m`.
- Outputs: `X_space` with populated `x`.
- Side effects: reads NWChem files from `r_params.nwchem_dir()`.
- Called by: `initialize`.
- Calls:
  - `slymer::NWChem_Interface`
  - `transform`
  - `make_random_trial` if NWChem does not provide enough virtuals
  - `QProjector`, `normalize`
- Algorithmic role: more structured initial subspace than purely random guesses.
- Notes / assumptions:
  - Each constructed trial row contains one virtual orbital inserted into one occupied-orbital slot.

### `ExcitedResponse::create_trial_functions2`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 386
- Purpose: build guesses from Cartesian coordinates times ground orbitals.
- Inputs: `World&`.
- Outputs: `X_space` with `3 * n * n` trial rows.
- Side effects: none outside returned object.
- Called by: `initialize`.
- Calls:
  - `make_xyz_functions`
  - `truncate`
- Algorithmic role: deterministic `x/y/z`-based trial space.
- Notes / assumptions:
  - The code creates `3 * n * n` states and populates exactly one occupied-orbital component in each row.
  - Fact: the outer loop index `i` is not used in the assignment (`ExcitedResponse.cpp` lines 414-423). This means the generated rows contain repeated blocks, and the implementation does not exactly match the explanatory comment block that follows.

### `ExcitedResponse::create_virtual_ao_guess`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 2832
- Purpose: generate an occupied-to-virtual AO-based guess subspace.
- Inputs: `World&`.
- Outputs: `X_space` whose rows correspond to virtual/occupied placements.
- Side effects: builds an AO basis (`aug-cc-pvdz`) in-memory and diagonalizes an approximate virtual-space Fock-like matrix.
- Called by: `initialize`.
- Calls:
  - `AtomicBasisSet`
  - `QProjector`
  - `svd`
  - local lambda `F(...)` for approximate Fock action
  - `syev`
- Algorithmic role: legacy attempt to construct a virtual-space-informed initial subspace without external NWChem files.
- Notes / assumptions:
  - Not used in the default configuration because `guess_xyz` defaults to `true`.
  - The guess bundle is still bundle-oriented: one row per `(virtual, occupied slot)` placement.

### `ExcitedResponse::iterate_trial`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 479
- Purpose: improve the initial guess subspace before selecting the working bundle.
- Inputs: `World&`, mutable `X_space& guesses`.
- Outputs: none directly; mutates `guesses` and writes `omega`.
- Side effects:
  - mutates `guesses.x`,
  - mutates solver member `omega`.
- Called by: `initialize`.
- Calls:
  - `transition_densityTDA`
  - `compute_response_potentials(..., "tda")`
  - `rotate_excited_space`
  - `create_shift`
  - `create_bsh_operators`
  - `apply(world, bsh_x_operators, ...)`
  - `gram_schmidt`, `normalize`
- Algorithmic role: precondition the initial subspace with repeated TDA-like diagonalize-plus-BSH steps.
- Notes / assumptions:
  - Operates only on `x`, not on `y`.
  - Uses `guess_max_iter`, not `maxiter`.

### `ExcitedResponse::iterate`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 2032
- Purpose: main excited-state solve loop.
- Inputs: `World&`.
- Outputs: none directly; mutates `Chi`, `omega`, JSON/timing state, and `all_done`.
- Side effects:
  - mutates `Chi`,
  - mutates `omega`,
  - mutates convergence arrays and JSON,
  - may save restart data.
- Called by: `ResponseBase::solve`.
- Calls:
  - `load_balance_chi`
  - `make_density`
  - `update_response`
  - `function_data_to_json`
  - `excited_to_json`
  - `save`
- Algorithmic role: repeated solve/update/convergence loop for the full excited-state bundle.
- Notes / assumptions:
  - Convergence is bundle-wide and based on maxima over all states.
  - No state is removed from the active set after converging; all states remain in the bundle until the whole solve stops.
  - `kain_x_space` uses a hard-coded `set_maxsub(10)` when KAIN is enabled (`ExcitedResponse.cpp` lines 2086-2097), not `r_params.maxsub()`.
  - Step restriction is effectively disabled because the call is wrapped in `if (false)` in `update_response` (`ExcitedResponse.cpp` lines 2463-2466).

### `ExcitedResponse::update_response`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 2375
- Purpose: compute one nonlinear response update from the current bundle.
- Inputs:
  - current `Chi`,
  - XC operator,
  - ground-state projector,
  - KAIN storage,
  - iteration index and step limit,
  - previous residual metadata.
- Outputs:
  - `new_omega`
  - `rotated_chi`
  - `new_chi`
  - `residuals { residual, residual_norms }`
- Side effects: none on members except through the returned values.
- Called by: `iterate`.
- Calls:
  - `compute_response_potentials`
  - `rotate_excited_space`
  - `bsh_update_excited`
  - `update_residual`
  - `kain_x_space_update`
  - `normalize`
- Algorithmic role:
  - build the operator-applied bundle,
  - diagonalize and rotate the current state bundle,
  - assemble the BSH right-hand side,
  - apply the BSH update,
  - measure the residual,
  - optionally apply KAIN.
- Notes / assumptions:
  - `theta_X` is formed as `rotated_v_x - rotated_EOX + rotated_gamma_x`; `T0X` is not included because the Green's operator update handles the kinetic part.
  - `Xvector` and `Xresidual` parameters are currently unused in the function body.
  - The caller passes `Tensor<double>()` for `old_residuals` (`ExcitedResponse.cpp` lines 2247-2249), which makes the intended residual history behavior unclear.

### `ResponseBase::compute_response_potentials`

- Location: `src/apps/molresponse/ResponseBase.cpp` line 1116
- Purpose: assemble the operator-applied vectors needed for excited-state diagonalization and BSH updates.
- Inputs: `World&`, `const X_space& chi`, XC operator, `calc_type`.
- Outputs: `(Lambda_X, V0X, gamma)`.
- Side effects: none outside returned objects.
- Called by: `iterate_trial`, `update_response`.
- Calls:
  - `T(world, chi_copy.x)` / `T(world, chi_copy.y)`
  - `compute_V0X`
  - `compute_gamma_full`, `compute_gamma_static`, or `compute_gamma_tda`
- Algorithmic role:
  - `T0X`: kinetic operator applied to the current bundle,
  - `E0X`: occupied-space Hamiltonian contribution,
  - `V0X`: ground-state Fock/local potential contribution,
  - `gamma`: first-order response kernel contribution,
  - `Lambda_X = T0X + V0X - E0X + gamma`.
- Notes / assumptions:
  - `compute_Y` is `true` only when `calc_type == "full"`.
  - This is the main assembly point for the excited-state linearized operator.

### `ResponseBase::compute_V0X`

- Location: `src/apps/molresponse/ResponseBase.cpp` line 1179
- Purpose: apply the ground-state potential/Fock-like operator to the response bundle.
- Inputs: `World&`, `const X_space& X`, XC operator, `compute_Y`.
- Outputs: `X_space V0`.
- Side effects:
  - may reuse stored potentials if `store_potential` is true,
  - uses the shared Coulomb operator and response macro task for exchange.
- Called by: `compute_response_potentials`, `compute_theta_X`, `compute_F0X`.
- Calls:
  - Coulomb of the ground density,
  - XC potential,
  - `ResponseComputeGroundExchange` via `MacroTask`
- Algorithmic role: evaluate the ground-state local + Coulomb + exchange + XC potential acting on the current response vectors.
- Notes / assumptions:
  - Uses `stored_v_nuc` and `stored_v_coul` when `store_potential` is enabled.
  - In full mode, acts on the concatenated `[x, y]` bundle through `X_space::to_vector()` / `from_vector()`.

### `ResponseBase::compute_gamma_tda`

- Location: `src/apps/molresponse/ResponseBase.cpp` line 924
- Purpose: compute the TDA first-order response-kernel term `gamma`.
- Inputs: `World&`, `(chi, phi0, rho1)` bundle, XC operator.
- Outputs: `X_space gamma` with populated `x`.
- Side effects: temporarily load balances orbitals and changes the MADNESS process map.
- Called by: `compute_response_potentials`, `compute_lambda_X`, `compute_theta_X`.
- Calls:
  - `transition_densityTDA`
  - Coulomb operator on each transition density
  - `newK`
  - XC kernel application
  - `QProjector`
- Algorithmic role: build the response-kernel correction `2J - c_xc K + W` in TDA form.
- Notes / assumptions:
  - Projects the result back out of the occupied subspace after construction.

### `ResponseBase::compute_gamma_full`

- Location: `src/apps/molresponse/ResponseBase.cpp` line 687
- Purpose: compute the full `X/Y` first-order response-kernel term `gamma`.
- Inputs: `World&`, `(chi, phi0, rho1)` bundle, XC operator.
- Outputs: `X_space gamma` with both `x` and `y`.
- Side effects: temporarily load balances orbitals and changes the MADNESS process map.
- Called by: `compute_response_potentials`, `compute_lambda_X`.
- Calls:
  - `make_density`
  - `response_context.compute_j1`
  - `response_context.compute_k1`
  - `response_context.compute_VXC1`
  - `QProjector`
- Algorithmic role: build the full response-kernel correction `2J - c_xc K + (1-c_xc)W`.
- Notes / assumptions:
  - Fact: this path depends on `response_context` strategies.
  - Fact: the excited-state path contains no visible `response_context.set_strategy(...)` call.
  - Inference: the full/non-TDA excited-state path is incomplete or depends on hidden setup not present in this directory.

### `ExcitedResponse::rotate_excited_space`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 813
- Purpose: diagonalize the current state bundle in the current subspace and rotate all associated vectors into the resulting eigenbasis.
- Inputs: `chi`, `lchi`, `v_chi`, `gamma_chi`.
- Outputs:
  - `new_omega`
  - rotated `chi`
  - rotated `lchi`
  - rotated `v_chi`
  - rotated `gamma_chi`
- Side effects: none outside returned values.
- Called by: `iterate_trial`, `update_response`.
- Calls:
  - `excited_eig`
  - `rotate_excited_vectors`
- Algorithmic role:
  - build overlap matrix `S = <x|x> - <y|y>`
  - build response matrix `A = <chi|lchi>`
  - symmetrize both
  - solve the generalized eigenproblem
  - rotate the whole bundle into that basis.
- Notes / assumptions:
  - The state order is not preserved across this step. State identity is whatever row index remains after the rotation and subsequent sorting.

### `ExcitedResponse::excited_eig`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 963
- Purpose: solve the generalized subspace eigenproblem used by `rotate_excited_space`.
- Inputs: overlap matrix `S`, response matrix `A`, degeneracy threshold.
- Outputs: `(new_omega, U)`.
- Side effects: mutates local copies of `S` and `A`.
- Called by: `rotate_excited_space`.
- Calls:
  - `sygvp`
  - `sort_eigenvalues`
  - `svd` for degenerate-cluster phase cleanup
- Algorithmic role:
  - solve `A c = S c omega` in the current state subspace,
  - postprocess eigenvectors to stabilize ordering and phases,
  - sort roots into ascending order.
- Notes / assumptions:
  - Uses a symmetric generalized solver (`sygvp`) after explicit symmetrization of `S` and `A`.
  - This is the active root-extraction routine in the main control path.

### `ExcitedResponse::bsh_update_excited`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 2480
- Purpose: apply the Green's-function / BSH update to the rotated excited-state bundle.
- Inputs: `omega`, mutable `theta_X`, projector.
- Outputs: new `X_space bsh_X`.
- Side effects: mutates `theta_X` in place before applying the BSH operators.
- Called by: `update_response`.
- Calls:
  - `create_shift`
  - `apply_shift`
  - `create_bsh_operators`
  - `apply(world, bsh_x_ops, ...)`
  - `QProjector`
  - `normalize`
- Algorithmic role:
  - shift the RHS so the BSH operators are well-defined,
  - build one BSH operator per `(state, occupied orbital)` pair,
  - apply them to obtain updated `x` and optionally `y`.
- Notes / assumptions:
  - For the `y` channel in full mode it uses `omega_minus = -omega` (`ExcitedResponse.cpp` lines 2510-2515).
  - The projector is applied after the BSH operator.

### `ExcitedResponse::create_shift`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 1888
- Purpose: compute additive shifts that keep the BSH denominator negative.
- Inputs: ground orbital energies, current excitation energies, channel label.
- Outputs: tensor `shift[state, orbital]`.
- Side effects: debug printing.
- Called by: `iterate_trial`, `bsh_update_excited`.
- Calls: none beyond printing.
- Algorithmic role: if `ground(p) + omega(k) > 0`, shift it to `-(ground + omega + 0.05)` so `mu = sqrt(-2*(...))` is real.
- Notes / assumptions:
  - The target `-0.05` is hard-coded and described as arbitrary in the comment (`ExcitedResponse.hpp` lines 84-95).

### `ExcitedResponse::apply_shift`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 1942
- Purpose: add the computed shift times the current function to the BSH RHS.
- Inputs: shift tensor, potential-like `response_space V`, current functions `response_space f`.
- Outputs: shifted `response_space`.
- Side effects: truncates the returned response space.
- Called by: `iterate_trial`, `bsh_update_excited`.
- Calls: none beyond arithmetic/truncation.
- Algorithmic role: implement `V + shift * f` for each `(state, orbital)` pair.

### `ExcitedResponse::create_bsh_operators`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 1976
- Purpose: create the per-state/per-orbital BSH convolution operators.
- Inputs: shift tensor, ground energies, current `omega`, `lo`, `thresh`.
- Outputs: nested vector of BSH operators indexed by `[state][orbital]`.
- Side effects: debug printing.
- Called by: `iterate_trial`, `bsh_update_excited`.
- Calls: `BSHOperatorPtr3D`.
- Algorithmic role: build the Green's operators used for the fixed-point update.
- Notes / assumptions:
  - Uses `mu = sqrt(-2.0 * (ground(p) + omega(k) + shift(k, p)))` (`ExcitedResponse.cpp` lines 2002-2005).

### `ResponseBase::update_residual`

- Location: `src/apps/molresponse/ResponseBase.cpp` line 1440
- Purpose: compute the difference between the previous and newly BSH-updated response bundle and derive per-state residual norms.
- Inputs: old `chi`, updated `g_chi`, `calc_type`, previous residual tensor, previous residual bundle.
- Outputs: `residuals { residual, residual_norms }`.
- Side effects: timing/logging.
- Called by: `ExcitedResponse::update_response`, `FrequencyResponse::update_response`.
- Calls:
  - `to_response_matrix`
  - `norm2s`
- Algorithmic role: measure nonlinear update size state-by-state.
- Notes / assumptions:
  - Fact: it ignores its `calc_type` argument and instead decides whether `Y` exists using `r_params.omega() != 0.0` (`ResponseBase.cpp` line 1466).
  - Fact: for an excited-state run, `ResponseParameters::omega` is generally `0.0`, so this function treats the problem as `x`-only even when `ExcitedResponse` is running in full mode.

### `ResponseBase::kain_x_space_update`

- Location: `src/apps/molresponse/ResponseBase.cpp` line 1570
- Purpose: apply KAIN acceleration to the state bundle.
- Inputs: current `chi`, residual bundle, per-state KAIN solver storage.
- Outputs: updated `X_space`.
- Side effects: mutates KAIN solver histories.
- Called by: `ExcitedResponse::update_response`, `FrequencyResponse::update_response`.
- Calls:
  - `to_response_matrix`
  - `XNonlinearSolver::update`
- Algorithmic role: nonlinear acceleration of the state bundle update.
- Notes / assumptions:
  - Fact: it also decides whether `Y` exists using `r_params.omega() != 0.0` (`ResponseBase.cpp` line 1582).
  - In the excited-state solver this means KAIN is effectively `x`-only unless the unrelated frequency-response `omega` input is nonzero.

### `ExcitedResponse::save`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 2746
- Purpose: write a restart snapshot of the excited-state bundle.
- Inputs: `World&`, archive filename.
- Outputs: none.
- Side effects: writes a parallel archive.
- Called by: `initialize`, `iterate`, `ResponseBase::solve` via restart/save flow.
- Calls: archive serialization only.
- Algorithmic role: bundle-level persistence.
- Notes / assumptions:
  - Stored fields are:
    - ground archive name
    - `tda` flag
    - `num_orbitals`
    - `num_states`
    - `omega`
    - `Chi.x`
    - optional `Chi.y`
  - It does not store convergence history, protocol position, KAIN history, density history, or active-state metadata.

### `ExcitedResponse::load`

- Location: `src/apps/molresponse/ExcitedResponse.cpp` line 2785
- Purpose: restore a restart snapshot.
- Inputs: `World&`, archive filename.
- Outputs: none directly; mutates `Chi` and `omega`.
- Side effects: reads a parallel archive.
- Called by: `ResponseBase::solve` on restart.
- Calls: archive deserialization only.
- Algorithmic role: restore the bundle before entering `iterate()`.
- Notes / assumptions:
  - Fact: `Chi` is reconstructed as `X_space(world, num_states, num_orbitals)` and then populated from the archive.
  - Ambiguity: the function also reads archive metadata into expressions like `r_params.archive()` and `r_params.num_states()` (`ExcitedResponse.cpp` lines 2805-2808), but the visible `ResponseParameters` accessors return by value (`response_parameters.h` lines 107-160). It is not clear from this file alone whether that metadata actually mutates `r_params`.

### Dormant diagonalization/augmentation helper stack

The following functions are defined in `ExcitedResponse`, but no call site in the current compiled excited-state path references them:

- `deflateGuesses` (`ExcitedResponse.cpp` line 696)
- `deflateTDA` (`ExcitedResponse.cpp` line 718)
- `deflateFull` (`ExcitedResponse.cpp` line 758)
- `augment` / `augment_full` (`ExcitedResponse.cpp` lines 1104, 1197)
- `unaugment` / `unaugment_full` (`ExcitedResponse.cpp` lines 1244, 1282)
- `diagonalizeFullResponseMatrix` (`ExcitedResponse.cpp` line 1365)
- `GetFullResponseTransformation` (`ExcitedResponse.cpp` line 1420)
- `diagonalizeFockMatrix` (`ExcitedResponse.cpp` line 1641)
- `get_fock_transformation` (`ExcitedResponse.cpp` line 1686)
- `reduce_subspace` (`ExcitedResponse.cpp` line 866)

Fact:

- A search in `src/apps/molresponse` found definitions for these functions but no active call sites in the current `ExcitedResponse::initialize -> iterate_trial -> iterate -> update_response` path.

Interpretation:

- These appear to be remnants of an older or alternate subspace-deflation design that was not removed after the current `rotate_excited_space`-based control path was adopted.

## 7. Solver algorithm reconstruction

### 7.1 Initialization algorithm

Reconstructed algorithm from `ExcitedResponse::initialize` and `iterate_trial`:

```text
allocate trial bundle with ~2*num_states rows
build trial.x using one of several guess generators
load-balance trial.x and ground orbitals
project trial.x out of occupied space
repeat twice:
  Gram-Schmidt trial.x
  normalize trial.x
for guess iteration = 0 .. guess_max_iter-1:
  rho <- transition_densityTDA(trial.x)
  project and normalize trial.x
  (Lambda, V0, gamma) <- compute_response_potentials(trial, calc_type="tda")
  (omega, chi_rot, lambda_rot, v_rot, gamma_rot) <- rotate_excited_space(...)
  if not last guess iteration:
    theta <- v_rot - E0X + gamma_rot
    shift <- create_shift(ground_energies, omega)
    trial.x <- BSH(theta, shift, omega)
    project and mask trial.x
  Gram-Schmidt trial.x
  normalize trial.x
sort trial.x by omega
select the lowest num_states rows into Chi.x
zero Chi.y
save guess_restart
```

### 7.2 Main excited-state algorithm

Reconstructed algorithm from `ExcitedResponse::iterate` and `update_response`:

```text
rho <- make_density(Chi)
for iter = 0 .. maxiter-1:
  optionally load-balance Chi

  if iter > 0:
    evaluate convergence from:
      density_residuals = ||rho_new - rho_old||
      relative_bsh = bsh_residual / ||Chi||
    if converged:
      optionally save
      break

  (Lambda, V0, gamma) <- compute_response_potentials(Chi, calc_type)
  (omega_new, Chi_rot, Lambda_rot, V_rot, Gamma_rot) <- rotate_excited_space(...)
  E0_rot <- 0 or Chi_rot * ham_no_diag depending on localization mode
  theta <- V_rot - E0_rot + Gamma_rot
  Chi_new <- bsh_update_excited(omega_new, theta)
  residual <- Chi_rot - Chi_new   [shared helper currently x-only if omega==0]
  if KAIN enabled and iter > 0:
    Chi_new <- KAIN(Chi_rot, residual)
  normalize/truncate Chi_new

  rho_old <- make_density(Chi_rot)
  rho_new <- make_density(Chi_new)
  density_residuals <- ||rho_new - rho_old||
  bsh_residuals <- residual norms
  Chi <- Chi_new
  omega <- omega_new
```

### 7.3 Meaning of the major assembled quantities

Fact from naming and implementation:

- `T0X`: kinetic operator applied to the response bundle (`ResponseBase.cpp` lines 1129-1139).
- `V0X`: ground-state Fock/local potential contribution (`ResponseBase.cpp` lines 1149, 1179-1332).
- `E0X`: occupied-space Hamiltonian contribution built from `hamiltonian` or `ham_no_diag` (`ResponseBase.cpp` lines 1142-1147; `ExcitedResponse.cpp` lines 2420-2428).
- `gamma`: response-kernel correction from Coulomb, exchange, and XC response (`ResponseBase.cpp` lines 1152-1163, 687-1054).
- `Lambda_X`: `T0X + V0X - E0X + gamma` (`ResponseBase.cpp` lines 1165-1170).
- `theta_X`: `V0X - E0X + gamma` in the rotated basis (`ExcitedResponse.cpp` lines 2441-2453).

Inference:

- `Lambda_X` is the operator used to form the generalized eigenproblem in the current bundle subspace.
- `theta_X` is the non-kinetic part of the fixed-point equation that the BSH Green's operator acts on.

### 7.4 Convergence logic

Convergence in `ExcitedResponse::iterate` is based on two bundle-level maxima (`ExcitedResponse.cpp` lines 2149-2238):

- density convergence:
  - `density_residuals.max() < conv_den`
  - where `conv_den = max(100 * thresh, dconv)`
- BSH update convergence:
  - `relative_max_bsh < relative_max_target`
  - where `relative_max_bsh` is `bsh_residual / ||Chi||`
  - and `relative_max_target = max(50 * thresh, 0.5 * dconv)`

Fact:

- `conv_only_dens` can bypass the BSH criterion (`ExcitedResponse.cpp` lines 2203-2205, `response_parameters.h` line 50).

Important behavior:

- The solver does not deactivate already-converged excited states.
- Convergence is all-or-nothing at the bundle level.

## 8. Data flow and state management

### 8.1 How response vectors are created

Creation sources:

- random Gaussian-masked noise: `make_random_trial`
- NWChem virtual orbitals: `make_nwchem_trial`
- Cartesian coordinate times occupied orbitals: `create_trial_functions2`
- AO/virtual-space construction: `create_virtual_ao_guess`

All of these return an `X_space` where:

- `x[state][orbital]` contains the trial occupied-orbital response amplitudes,
- `y` is initially zero or left default-zero.

### 8.2 How the main bundle is stored

The working bundle is `ResponseBase::Chi`:

- `Chi.x[b][p]` means state `b`, occupied orbital slot `p`, `X` channel.
- `Chi.y[b][p]` means state `b`, occupied orbital slot `p`, `Y` channel.

State identity is positional:

- state `b` is "whatever row `b` currently means after the latest rotation/sort".
- The code does not preserve a stable state object or UUID across rotations.

### 8.3 How states are distinguished

Fact:

- The only durable per-state identifier in the running solver is the row index into `X_space`.
- Energies are kept in `ExcitedResponse::omega`, one scalar per row (`ExcitedResponse.hpp` line 25).

Implication:

- When `rotate_excited_space` applies `U`, every state becomes a linear combination of the previous states.
- After `sort_eigenvalues` or `sort`, row ordering changes again.
- Any later code that refers to "state 0" or "state 1" is referring to the current sorted/rotated ordering, not to a persistent original state object.

### 8.4 Where frequencies / excitation energies enter

`omega` enters the computation in two places:

1. Root extraction:
   - `rotate_excited_space -> excited_eig` computes `omega` from the generalized subspace eigenproblem (`ExcitedResponse.cpp` lines 813-863, 963-1102).
2. BSH update:
   - `create_shift` uses `ground_energies` and `omega`,
   - `create_bsh_operators` uses `ground_energies`, `omega`, and the shifts to build the Green's operators (`ExcitedResponse.cpp` lines 1888-2018, 2487-2516).

In full mode:

- `bsh_update_excited` uses `omega` for `x` and `-omega` for `y` (`ExcitedResponse.cpp` lines 2510-2515).

### 8.5 How all states are updated together

Fact:

- All requested states are solved together inside a single bundle.

Coupling points:

- shared orthogonalization (`gram_schmidt`)
- shared subspace diagonalization (`rotate_excited_space`)
- shared sorting and selection (`sort`, `select_functions`)
- shared normalization (`normalize(world, X_space&)`)
- shared convergence stop condition
- shared save/load snapshot

Per-state updates inside the shared bundle:

- `create_shift` computes one shift per `(state, occupied orbital)` pair.
- `create_bsh_operators` creates one BSH operator per `(state, occupied orbital)` pair.
- `make_density` and `gamma` computations return one density / one response-kernel vector per state.

### 8.6 Data-flow summary

The core state data path is:

```text
guess generator
  -> trial.x
  -> iterate_trial
  -> sort/select
  -> Chi
  -> compute_response_potentials(Chi)
     -> Lambda_X, V0X, gamma
  -> rotate_excited_space(Chi, Lambda_X, V0X, gamma)
     -> omega, rotated bundles
  -> theta_X = V0X - E0X + gamma
  -> bsh_update_excited(theta_X, omega)
     -> new_chi
  -> update_residual(rotated_chi, new_chi)
  -> Chi := new_chi
  -> save/load as a bundle
```

## 9. Legacy design constraints relevant to reintegration

### 9.1 Bundle-centric storage is fundamental

Fact:

- The old design assumes all excited states live together in one `X_space`.

Porting consequence:

- Any new architecture that wants one object per state will need adapters for:
  - dense state-space rotations,
  - bundle-level diagonalization,
  - bundle-level save/load,
  - bundle-level convergence metrics.

### 9.2 State identity is mutable and positional

Fact:

- The old solver repeatedly rotates and sorts the entire state bundle.

Porting consequence:

- A later reintegration cannot assume "state 3" is a persistent object whose contents are updated in place without reidentification.

### 9.3 The solver logic is written around shared containers

Examples:

- `rotate_excited_space` takes the whole bundle and returns a whole rotated bundle.
- `compute_response_potentials` works on the whole `X_space`.
- `save/load` persist the whole bundle.
- `make_density` and `update_residual` consume whole bundles.

Porting consequence:

- The old code assumes tightly coupled bundle-wide operations rather than per-state methods.

### 9.4 Save/load is snapshot-based, not workflow-based

Fact:

- `ExcitedResponse::save` stores only:
  - response bundle,
  - `omega`,
  - a few scalar metadata fields.

Missing from the snapshot:

- protocol index,
- current iteration,
- residual history,
- density history,
- KAIN history,
- JSON/timing state.

Porting consequence:

- Restart reuse in a newer architecture will need a richer persistence model if it wants true iterative restarts.

### 9.5 Full-mode dependencies are entangled with unrelated shared infrastructure

Fact:

- Full excited-state gamma assembly depends on `response_context` strategies.
- `FrequencyResponse` initializes them; `ExcitedResponse` does not.

Fact:

- Several shared helpers decide whether `Y` exists from `r_params.omega()`, not from excited-state mode.

Porting consequence:

- Reintegration should explicitly separate:
  - excited-state mode,
  - TDA/full channel handling,
  - frequency-response `omega`,
  - strategy initialization.

### 9.6 Several control knobs are not actually honored

Fact:

- `maxsub` is not used in excited-state KAIN setup; the code hard-codes `10` (`ExcitedResponse.cpp` lines 2091-2097).
- `step_restrict` is effectively bypassed in `ExcitedResponse::update_response` because the branch is `if (false)` (`ExcitedResponse.cpp` lines 2463-2466).
- `analysis()` exists but is not called (`ExcitedResponse.cpp` lines 2357-2371, 2557-2743).

Porting consequence:

- Later reintegration should not assume that input parameters in the legacy file were all behaviorally active.

### 9.7 The file contains dormant alternative solver machinery

Fact:

- There is a large unused helper stack for subspace augmentation/deflation and alternate diagonalization.

Porting consequence:

- A reintegration effort should not blindly re-implement every function in `ExcitedResponse.cpp`; some are historical remnants, not part of the active path.

## 10. Open questions / ambiguities

### 10.1 Was the full excited-state path ever fully wired?

Fact:

- No `response_context.set_strategy(...)` call appears on the `ExcitedResponse` path.
- `compute_gamma_full()` and `compute_gamma_static()` depend on `response_context`.

Inference:

- The non-TDA/full path appears incomplete or stale in this checkout.

Code locations:

- `src/apps/molresponse/ExcitedResponse.hpp` line 19
- `src/apps/molresponse/ResponseBase.cpp` lines 687-799
- `src/apps/molresponse/FrequencyResponse.hpp` lines 564-583

### 10.2 Is `Y` intentionally ignored by shared residual/KAIN helpers in excited-state runs?

Fact:

- `update_residual`, `kain_rf_space_update`, `kain_x_space_update`, `x_space_step_restriction`, and part of `load_balance_chi` all key off `r_params.omega() != 0.0`.

Inference:

- For excited-state runs, where `omega` is typically not the excitation energy but the frequency-response input, the shared helpers behave as though `Y` is absent even when `calc_type == "full"`.

Code locations:

- `src/apps/molresponse/ResponseBase.cpp` lines 408-419, 1466, 1515, 1582, 1619
- `src/apps/molresponse/response_parameters.h` line 97

### 10.3 Is the residual history API in `ExcitedResponse::update_response` unfinished?

Fact:

- `ExcitedResponse::iterate` passes `Tensor<double>()` as `old_residuals` into `update_response`.
- `ResponseBase::update_residual` starts by copying `old_residuals`.

Ambiguity:

- It is unclear whether the residual tensor is supposed to be pre-sized elsewhere, or whether this is unfinished/stale code.

Code locations:

- `src/apps/molresponse/ExcitedResponse.cpp` lines 2247-2249
- `src/apps/molresponse/ResponseBase.cpp` lines 1468-1499

### 10.4 Does restart metadata actually mutate `ResponseParameters` during load?

Fact:

- `ExcitedResponse::load` and `LoadExcitedXSpace::load_x_space` deserialize into expressions like `r_params.archive()` and `r_params.num_states()`.
- The visible `ResponseParameters` accessors return by value.

Ambiguity:

- From this file set alone, it is not clear whether the archive metadata is being meaningfully restored or just read through a non-mutating expression.

Code locations:

- `src/apps/molresponse/ExcitedResponse.cpp` lines 2805-2808
- `src/apps/molresponse/ResponseBase.hpp` lines 110-127
- `src/apps/molresponse/response_parameters.h` lines 107-160

### 10.5 Is `create_trial_functions2` doing what the comment says?

Fact:

- The explanatory comment describes one `x/y/z`-modulated orbital placed into each occupied slot.
- The actual code does not use the outer loop variable `i` in the assignment.

Ambiguity:

- The implementation appears to generate repeated blocks rather than the exact pattern described in the comment.

Code locations:

- `src/apps/molresponse/ExcitedResponse.cpp` lines 413-455

### 10.6 Which diagonalization helper stack is the intended one?

Fact:

- The active control path uses `rotate_excited_space -> excited_eig`.
- A second large helper stack (`deflateTDA`, `deflateFull`, `diagonalizeFullResponseMatrix`, `get_fock_transformation`, etc.) is defined but not called from the active path.

Ambiguity:

- It is unclear whether those functions are abandoned experiments, unfinished alternates, or code paths that were active in an earlier revision.

### 10.7 Historical duplicate implementation

Fact:

- `calc_runner.cc` contains another excited-state driver (`TDDFT::solve_excited_states`) that is not in the current build.

Ambiguity:

- If later reintegration work needs the "oldest" behavior, a human review may be needed to decide whether `calc_runner.cc` captures an earlier intended algorithm better than the compiled `ExcitedResponse` path.

Code locations:

- `src/apps/molresponse/calc_runner.cc` lines 406-561
- `src/apps/molresponse/CMakeLists.txt` lines 12-23, 42

## 11. Appendix: code reference index

### Build and entry

- `src/apps/molresponse/CMakeLists.txt` lines 12-23, 42
- `src/apps/molresponse/molresponse.cc :: main` (line 60)
- `src/apps/molresponse/global_functions.cc :: initialize_calc_params` (line 14)

### Parameters and mode selection

- `src/apps/molresponse/response_parameters.h :: ResponseParameters` (line 24)
- `src/apps/molresponse/response_parameters.h :: ResponseParameters::set_derived_values` (line 170)

### Core containers

- `src/apps/molresponse/response_functions.h :: response_space` (line 26)
- `src/apps/molresponse/response_functions.h :: response_space_inner` (line 559)
- `src/apps/molresponse/x_space.h :: X_space` (line 33)
- `src/apps/molresponse/x_space.h :: X_space::to_vector` (line 422)
- `src/apps/molresponse/x_space.h :: X_space::from_vector` (line 442)
- `src/apps/molresponse/x_space.cc :: to_response_matrix` (line 56)
- `src/apps/molresponse/x_space.cc :: inner(const X_space&, const X_space&)` (line 185)

### Shared solver base

- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::ResponseBase` (line 12)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::check_k` (line 38)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::ComputeHamiltonianPair` (line 142)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::solve` (line 1741)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::make_density` (line 368)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::load_balance_chi` (line 397)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::compute_response_potentials` (line 1116)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::compute_V0X` (line 1179)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::compute_gamma_full` (line 687)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::compute_gamma_static` (line 801)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::compute_gamma_tda` (line 924)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::update_residual` (line 1440)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::kain_x_space_update` (line 1570)
- `src/apps/molresponse/ResponseBase.cpp :: normalize(World&, X_space&)` (line 1858)
- `src/apps/molresponse/ResponseBase.cpp :: transition_densityTDA` (line 2153)
- `src/apps/molresponse/ResponseBase.cpp :: transform(World&, const response_space&, const Tensor<double>&)` (line 2180)
- `src/apps/molresponse/ResponseBase.cpp :: transform(World&, const X_space&, const Tensor<double>&)` (line 2205)
- `src/apps/molresponse/ResponseBase.cpp :: select_functions` (line 2267)
- `src/apps/molresponse/ResponseBase.cpp :: sort(World&, Tensor<double>&, X_space&)` (line 2334)
- `src/apps/molresponse/ResponseBase.cpp :: ResponseBase::output_json` (line 2125)

### Excited-state solver

- `src/apps/molresponse/ExcitedResponse.hpp :: ExcitedResponse` (line 16)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::initialize` (line 7)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::make_random_trial` (line 154)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::make_nwchem_trial` (line 187)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::create_trial_functions` (line 316)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::create_trial_functions2` (line 386)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate_trial` (line 479)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::rotate_excited_space` (line 813)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::excited_eig` (line 963)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::sort_eigenvalues` (line 1845)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::create_shift` (line 1888)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::apply_shift` (line 1942)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::create_bsh_operators` (line 1976)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::excited_to_json` (line 2020)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::iterate` (line 2032)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::update_response` (line 2375)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::bsh_update_excited` (line 2480)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::analysis` (line 2557)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::save` (line 2746)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::load` (line 2785)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::create_virtual_ao_guess` (line 2832)
- `src/apps/molresponse/ExcitedResponse.cpp :: ExcitedResponse::create_response_guess` (line 2991)

### Full-path infrastructure dependency contrast

- `src/apps/molresponse/FrequencyResponse.hpp :: FrequencyResponse::FrequencyResponse` (lines 564-583)

### Historical, not in current build

- `src/apps/molresponse/calc_runner.cc :: TDDFT::solve_excited_states` (line 406)

