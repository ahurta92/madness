# Convergence metric definitions and ownership

This inventory records the implemented named moldft -> molresponse workflow as
of 2026-10-20. A common label does not establish that two measurements are
numerically comparable. Formula, state identity, units, normalization, reduction,
derivative discretization and sampling phase must also agree.

## Shared kinetic definition

`KineticEnergy.h` owns the real three-dimensional gradient quadratic:

    T_ij = 1/2 sum_(d=x,y,z) <D_d f_i | D_d f_j>

`SCF::kinetic_energy_matrix` prepares its configured derivatives, including halo
staging and compression, then calls `kinetic_energy_matrix_from_gradients`.
The result remains a distributed matrix. The response `rs::kinetic_diagonal`
wrapper calls `kinetic_energy_diagonal`, which computes only T_ii, on deep copies,
with the default derivative operators used by the response kernels. Both
contractions use the same sum and factor in `detail::contract_kinetic_gradients`.
The overload accepting explicit derivative operators supports like-for-like
comparisons with SCF's `ble` and `bspline` choices. Sharing the mathematical
formula does not silently select a different derivative discretization.

The helper applies no occupation or spin factor, amplitude normalization,
reduction across functions, iteration difference or convergence threshold.
For a normalized orbital the diagonal is in hartree. A driven response function
retains its perturbation-dependent amplitude scale; its scalar therefore carries
that scale squared, as recorded in the response metadata.

SCF forms a separate occupation-weighted per-spin total in `make_fock_matrix`;
the total-energy evaluation combines spins (doubling the restricted result).
The response reporter instead records every stored X/Y function separately.
Neither operation changes the definition of the underlying T_ii.

## Shared per-function BSH definition

`FunctionResidual.h` defines `r_i = ||f_i - F_i(f)||_L2`, where `F` is the
application's raw fixed-point update before acceleration or damping. Every
stored function carries one scalar. No occupation/spin weight, amplitude
normalization or reduction is applied. `FunctionResidualReport.h` supplies
common serialization and validation; the parent aggregate carries units and
sample provenance. The application still owns its fixed-point equation.

SCF uses its already computed residual functions and norms. The observer's
components follow its historical aggregate from the previous update. Named
response measures the raw input/update pair before KAIN, on independent copies,
and records `function_bsh_history` for each solve/rung. Restricted ground state
has `alpha`; static response has `x_alpha`; full response also has `y_alpha`.
The helper supports explicit beta blocks without inventing restricted partners.

The existing aggregates remain maximum over SCF functions/spins and quadrature
over response functions/blocks. Existing response relative gates retain their
output-state normalization; the new raw components are not normalized by it.
Saved-state validation recomputes the fixed-point difference on the archived
state independently. Historical function indices refer to the sampled working
basis; they do not establish correspondence across orbital rotations or rungs.
Older checkpoints without these optional fields remain usable, with no inferred
components. Invalid individual scalars are rejected before reuse/publication.

## Current measurements

| Measurement | moldft | Named driven response |
|---|---|---|
| Density change | Maximum of alpha/beta L2 norms of successive input density differences; units a0^(-3/2). | L2 norm of successive kernel-induced density differences, per response channel. Static closed-shell density is 4 sum(phi_i X_i); full is 2 sum(phi_i (X_i+Y_i)). No cross-channel reduction. |
| BSH residual | Maximum orbital/spin L2 norm of the raw fixed-point difference; normalized orbital convention, dimensionless. | Quadrature over stored orbital functions and response blocks of the raw fixed-point difference, per channel; response-amplitude L2 units. |
| Energy change | Absolute change of successive total-energy scalars, hartree; currently reported without a new energy gate. | No local driven-state energy-change criterion; ground-reference evidence is carried separately. |
| Kinetic value/change | Occupation/spin-weighted total and absolute change of successive total scalars, hartree. | Individual unnormalized T_ii values initially and after every iteration, plus independent saved-state evaluation. Matched absolute component changes and their maximum are also reported; no kinetic gate. |

The two density and BSH reductions above are existing solver conventions, not
interchangeable definitions. This refactor does not relabel them as equivalent
or change their thresholds. See `MoldftCheckpoint.hpp`, `shared_response_checkpoint.hpp` and
`ConvergenceReport.hpp` for the current report fields; see the native SCF and
static/full kernels for the numerical construction.

## Sampling and stopping policy

Historical solver samples and independently reevaluated saved-state samples
remain distinct. The report records iteration, representation, phase, state and
reference identity, availability and whether a measurement is an active gate.
Function indices refer to the current reference's orbital ordering. Histories
restart at a new solve/rung/reference; unrelated function identities must not be
compared implicitly. Old missing measurements remain unavailable.

The existing SCF gates and the response absolute/relative gate policy remain in
force. Reported targets alone do not activate a new criterion. The approved
kinetic diagnostic is the maximum of matched absolute per-function changes,
without occupation/spin weights. `FunctionMetricSummary.h` owns the pure scalar
maximum and component-wise difference operations; the response adapter owns
comparison eligibility, units and sampling. Thresholds, SCF orbital alignment,
and activation of a new stopping criterion remain separate scientific decisions.

Each response BSH history entry now retains its own existing absolute/relative
limits and output normalization. `maximum_comparison` replaces only the residual
reduction for a counterfactual on that sample; it does not alter solver gates.
Incomplete older histories are not backfilled from final-iteration evidence.

## Extension boundary

Use shared numerical functions for common mathematics, application adapters for
state extraction and physical weights, and separate convergence policies for
stopping decisions. An application must not override a formula while keeping
its metric identity unchanged. A future common descriptor should make differing
normalizations and reductions explicit; this slice adds no descriptor/schema
migration or base-class hierarchy.

Coverage here is SCF and the named static/full response diagnostic path. Existing
`Kinetic<T,NDIM>` operators, response subspace Gram matrices, coupled excited-state
solvers and other applications have not all been migrated to this helper. Their
representations and contractions require an explicit audit before claiming
application-wide unification.

## Direct verification

`test_function_kinetics` checks analytic Gaussian values (including amplitude
scaling), source-function nonmutation, empty input, response versus actual SCF
matrix diagonals, an analytic off-diagonal, and SCF's actual occupation-weighted
Fock-builder total. It also compares the explicit-operator helper with SCF for
`ble` and `bspline` derivatives. Tests run in serial and on two MPI ranks.

`test_function_residuals` checks analytic differences [3,0] and [0,4], input
preservation, empty/mismatched inputs, and invalid component rejection, in
serial and on two MPI ranks. Named convergence acceptance reconstructs both
existing aggregates from their components, compares independent saved-state
checks with the next actual solver evaluation, exercises continuation and
rejects corrupted evidence while retaining the previous publication.

## SCF charge and spin density comparisons

`DensityChanges.h` defines two diagnostic L2 norms from application-supplied
physical spin-density fields, including their existing occupation factors:

    charge = ||(rho_alpha(n)-rho_alpha(p)) + (rho_beta(n)-rho_beta(p))||_L2
    spin   = ||(rho_alpha(n)-rho_alpha(p)) - (rho_beta(n)-rho_beta(p))||_L2

The helper preserves input functions and applies no normalization or tolerance.
These are field norms; sums/differences of the separate spin norms lose the
cross term and cannot implement them. Restricted closed-shell charge change is
twice the single-spin change and spin change is zero. Opposing spin changes can
cancel in charge while remaining visible in spin. With no beta electrons, charge
and spin changes are the same alpha-density norm.

SCF computes these only when its collective checkpoint observer is installed.
Named ground-state reports add `charge_density_change` and `spin_density_change`
beside the existing spin-maximum `density_change`. They retain historical input
and independently reevaluated saved-state phases; missing previous samples are
unavailable. Saved-state evaluation in this frontend remains restricted HF.
The shared field algebra itself does not impose a restricted-shell assumption.

No new default convergence criteria are supplied: both diagnostics carry
`active_gate=false` and null target, effective_limit and satisfied. Existing SCF
spin gates remain unchanged. Tests specify protocol, dconv, econv and iteration
budget explicitly in their reference decks; these existing controls must not be
silently reinterpreted as total-charge or spin-density targets. Future policies
need explicit case-appropriate criteria, not universal inferred defaults.


## Kinetic trace contract (2026-10-20)

The unweighted collection trace is `K = sum_i (1/2) sum_a <d_a f_i|d_a f_i>`.
The change criterion is `abs(K_n - K_previous)`, distinct from both the raw
trace and the maximum of per-function changes. Unitary rotation within a fixed
collection preserves this trace; no orbital matching is needed. The collection,
state/root, equation, protocol, representation, and sampling phase must still
match. Different collection dimensions and missing previous samples are
unavailable, never zero convergence evidence.

Named moldft counts its stored alpha orbitals once, plus stored beta orbitals
only when unrestricted. Restricted spin multiplicity and occupations do not
enter this trace. The occupation-weighted physical kinetic energy remains
separate. Each driven response channel counts its stored X/Y and alpha/beta
blocks positively, without combining independent channels or states and without
using the signed RPA norm. Response units retain the response amplitude scale.
All individual diagonals remain diagnostics. No trace threshold or gate is
implied by reporting. Derived trace evidence is checked against its components.
