# Input options: current behavior and proposed mapping

**Source review — 2026-11-11.** Current checkout:
`76255552e16e542f6b5ae2d0d7e2b6b034c7ed48`. This includes changes after drop 6.
No programs were built or run. The target behavior comes from
[the design](design_draft.md). **Approved on 2026-11-11:** require explicit
convergence names, with no `dconv`/`econv` aliases. Read-only inspection and reuse-or-iterate solver restarts were approved on
2026-11-12, as was requiring response `maxiter` without an `excited.maxiter`
alias. Explicit `excited.maxiter`, `save` and `checkpoint.keep_all` all use
the normal fatal unknown-input error, without custom migration messages.
The 2026-11-12 consolidation in [design section 1](design_draft.md#input-routing-and-retired-controls)
fixes the remaining shared routing and retired-policy controls for review.
Current-behavior columns remain source findings, not claims of implemented target support.

SCF, Nemo and molresponse are traced through their parameter definitions and
execution paths. The final table inventories other release engines' declared
controls; their inner algorithms still need the planned per-engine adoption work.

## Release interface — 2026-11-19

Madqc is the sole executable for supported engines. Named-stage input and
supported-engine standalone applications are retirement targets. The historical
columns below still describe their existing flags; those flags are not a promise
to retain standalone CLIs. Preserve madqc options/group ownership and migrate
scientific capabilities through its workflow interface. Separate executables
remain only for unsupported engines.

## SCF implementation update — 2026-11-19

For ordinary `madqc --wf=scf`, the SCF prerequisite of `--wf=response`, and
the temporary standalone `moldft` test entry, `dft convergence.criteria`, the four
`convergence.<metric>` absolute/relative controls and
`checkpoint.every_iterations` are registered. Omitted criteria select BSH,
density and energy changes for ordinary SCF; frozen virtuals select BSH only.
Defaults are absolute, with zero relative tolerance. `maxiter` is a fresh update
allowance, and endpoint checkpoints are always written. `nv_its` belongs to
explicit virtual-count preparation and does not consume that allowance.

SCF rejects explicit `dconv`, `econv`, `save` and retired policy switches.
Unadopted engines retain their own inputs and reject the new shared SCF keys.
File and CLI precedence is unchanged; every original CLI occurrence is validated.
JSON input is broadcast and validation failures are collective.

Solver `auto`/`iterate` reuses qualifying saved metrics or performs fresh checks.
Existing native/AO/NWChem seed routes supply guesses, not convergence evidence.
`--continue` requires an existing shared SCF selection. `read_only` inspects
saved evidence and can plot the saved native functions without evaluating a new
request, creating an attempt or replacing results. `--restart_info=<directory>`
reads shared SCF evidence from that application directory.

The historical source findings below still describe the pre-adoption interfaces;
response/Nemo convergence input migration and final named-driver retirement have
not occurred in this increment.

## What the current interfaces actually expose

The workflow path accepts existing `dft`, `response` and other application
parameter groups, including CLI overrides such as `--dft="maxiter=30"` and
`--response="maxiter=30"`. File settings are read before group CLI overrides.
Keep that interface and the workflow selectors. The separate named-stage path
has additional criteria, continuation and checkpoint controls that must be moved
into the retained workflow implementation; their existence does not establish
workflow support today. In particular, `--continue` is consumed by the named
runner, not a shared run-level workflow restart controller.

**Approved 2026-11-13:** validate every original command-line occurrence before
applying overrides. Keep the ordered user arguments as well as the parser's final
maps; use owning registrations/readers for syntax, key, type and per-value checks.
A retired, unknown or malformed earlier occurrence cannot be hidden by a later
valid one. Validation does not apply discarded values. File-then-CLI precedence
and the last-occurrence replacement of an entire CLI group remain unchanged;
do not merge earlier group payloads. Check cross-setting constraints only after
resolving the accepted inputs. Preserve path case and winning-input provenance.

For example, `--response="excited.maxiter=3" --response="maxiter=3"` fails with
the normal unknown-input error. Two valid `maxiter` group overrides still use
the last one. The same rule covers retained madqc CLI forms and malformed values.

**Approved 2026-11-16:** original-occurrence validation and final CLI reading
share the file reader's collective error handling. Rank zero distributes the
error type and message before any rank leaves the input path; every rank throws
consistently on failure. Accept resolved settings only on success. Invalid input
exits nonzero before numerical work or checkpoint publication, with the existing
diagnostics and precedence. Two-rank failure tests cover retired and malformed
inputs; serial/non-MPI handling remains supported.

**Implemented prerequisite, 2026-11-18:** the parser now retains ordered human
arguments. Registered group readers validate each original payload on a private
copy, then apply only the winning whole group. File and CLI reads share collective
error propagation, preserve case-sensitive values and leave settings unchanged
when parsing fails. The combined file/CLI entry is one transaction: rank zero
selects and checks the input file, and a rejected override also rolls back the
file values. Per-value hooks are available to owning readers; cross-setting
checks still follow resolution. Existing permissive/strict group policies remain
in force until adoption. Standalone option owners still need to consume the
retained occurrences; global retired-key cleanup is not complete.

The approved target names `convergence.criteria`, `convergence.<metric>` and
`convergence.<metric>.relative` inputs are not registered in the inspected
SCF/response parameter classes. They still require implementation.

Sources: [CLI dispatch](src/apps/madqc/madqc.cpp),
[named-path dispatch](src/apps/madqc/RunNamedStages.hpp),
[CLI occurrence parsing](src/madness/mra/commandlineparser.h),
[owning parameter readers](src/madness/mra/QCCalculationParametersBase.cc),
[input precedence](src/madness/mra/QCCalculationParametersBase.h),
[maintained SCF input](src/examples/qc/scf_he_hf/scf_he_hf.in),
[maintained response input](src/examples/qc/response_he_alpha/response_he_alpha.in).

## Convergence

| Existing option | What it does now | Mapping to the new design |
| --- | --- | --- |
| `dft dconv`; standalone `molresponse --dconv`; `response dconv` | More than a density tolerance. SCF uses `max(threshold,dconv)`, multiplies the density gate by `max(5,natom)` and the BSH gate by 5. Nemo uses its selected criteria and protocol-adjusted tolerances with a different residual reduction. Response derives density, BSH, excitation-energy and FD relative targets from this one value. | **Approved: no alias.** Require explicit `convergence.density_change` and `convergence.bsh_residual`, plus applicable energy and relative controls. Explicit `dconv` is a migration error in the adopted interface; explain the former coupled roles without copying its number automatically. |
| `dft econv`; named molresponse `econv` | Ordinary SCF's engine gate uses density/BSH, not an energy test merely because `econv` is set. Nemo's selected energy criteria use it. Named SCF can treat it as a requested energy-change target; named molresponse uses it for excited-root energy change and rejects it for driven-only work. `ResponseParameters` does not register ordinary response `econv`. | **Approved: no alias.** Require `convergence.energy_change`: total energy for SCF/Nemo, excitation energy for roots, inapplicable for driven states. Explicit `econv` is a migration error in the adopted interface; preserve the existing CLI group forms. |
| `convergence_criteria [bsh_residual,total_energy,each_energy,density]` | Registered in shared DFT parameters and consumed by Nemo; it does not select the ordinary SCF gate described above. `each_energy` checks the energy-component vector, not the same quantity as total-energy change. | Use `convergence.criteria`. Document conceptual replacements `total_energy` → `energy_change`, `density` → `density_change`; remove the old key with a generic unknown-input error, without parser aliases. `each_energy` is not the same metric as total-energy change and is not silently translated. |
| `conv_only_dens` | Deprecated SCF switch that removes BSH from its engine stopping test. | Express the intended selected criteria through the common list; no second solver-policy switch. Remove `conv_only_dens` with a generic unknown-input error. The new default metric list is not the old density-only selection. |
| `accuracy.bsh_residual`, `accuracy.kinetic_trace_change` | Registered for named targets. Ordinary SCF only builds those target checks when a checkpoint observer is attached; Nemo explicitly refuses these controls. Named response has corresponding target machinery. | Canonical `convergence.bsh_residual` and `convergence.kinetic_change`; remove the old keys with generic unknown-input errors. In particular, the old kinetic trace and residual conventions are not automatically the approved physical metrics. |
| `convergence.use_default_criteria` | Selects an engine check versus written targets in the named infrastructure. Nemo refuses written-target stopping. | Remove the two-policy choice. One resolved criteria list controls stopping; no explicit list means the agreed defaults. Remove the old key with the normal fatal unknown-input error for every value. |
| Relative convergence | FD response already uses `atol + rtol * norm`, but `rtol` is derived from protocol-adjusted `dconv` and the policy factor. It is not the proposed independent per-metric input. | Optional `convergence.<metric>.relative`, default zero, with the approved metric-specific scale. Do not silently carry the old hidden relative term into an absolute-only request. |
| `progress.window`, `progress.minimum_relative_improvement`; response `stall.window`, `stall.ratio` | Experimental explicit-target progress controls and an ES plateau detector exist. Stall detection is off by default; some named intermediate paths now suppress early exits to perform fixed work. | Stall-based stopping is deferred. The production path uses intermediate threshold-dependent absolute tolerances and an update cap. Remove all listed stall/progress keys; every explicit value uses the normal fatal unknown-input error. |
| `diagnostic.iterations_per_rung` | Can force a fixed number of iterations for diagnostics; zero disables it. | It must not override final convergence or silently restore fixed-work production policy. Remove the key for every explicit value with the normal fatal unknown-input error; experimental fixed-work studies remain outside this production contract. |

**Why aliases were rejected:** for a small SCF system with threshold `1e-6`
and `dconv=1e-4`, the current density and BSH limits are both `5e-4`.
Treating the same written number as a literal density-change tolerance gives
`1e-4`, and no longer sets the BSH tolerance. Response also loses its old derived
relative term. This is a scientific meaning change, not just a spelling change.

The approved intermediate heuristic remains: adjust applicable absolute limits
at intermediate protocols, then use the actual requested final tolerances.
Current SCF/response gates floor limits by the representation threshold even on
the final rung. That final relaxation must not survive in the new controller.
Nemo's protocol helper also has a distinct intermediate `0.1*sqrt(threshold)`
relaxation. Reusing legacy heuristics therefore requires reviewed mappings by
metric and model, not one factor copied across every application.

The [November 19 source audit](docs/design/madqc-legacy-convergence.md) records
the exact SCF/driven/excited engine formulas, the derived response default,
measurement timing, and shared-input examples. It distinguishes reproducing
tolerance magnitudes from reproducing a stopping test or numerical trajectory.

Sources: [DFT parameters](src/madness/chem/CalculationParameters.h),
[SCF gate limits](src/madness/chem/SCF.h), [SCF stopping](src/madness/chem/SCF.cc),
[Nemo checks](src/madness/chem/nemo.h), [Nemo support restrictions](src/madness/chem/nemo.cc),
[Nemo protocol limits](src/madness/chem/SCFProtocol.h),
[response policy](src/apps/molresponse/solvers/convergence_policy.hpp),
[named input handling](src/apps/madqc/NamedStages.hpp).

## Iteration budgets and protocol settings

| Existing option | What it does now | Mapping to the new design |
| --- | --- | --- |
| `dft maxiter`, `response maxiter`; standalone `--maxiter` | SCF/Nemo use the DFT cap; driven response uses its own cap. In madqc, an unwritten response `maxiter` inherits the DFT value. An unwritten response `dconv` instead derives from `100 * finest_protocol_threshold`. | Keep the cap and existing madqc group/CLI ownership; standalone-only flags retire with their executables. Each executed protocol gets a fresh main-update allowance per invocation; retain historical counts without subtracting them. Resolve inherited defaults explicitly in the settings report. Numerical tolerance defaults still require validation. |
| `excited.maxiter` | An explicitly supplied value overrides the response cap for the ES bundle. Otherwise the executor inherits response `maxiter`; the registered default of 20 is not blindly applied. | **Approved 2026-11-12:** require response `maxiter` only. Reject explicit `excited.maxiter` through the normal fatal unknown-input path, without a custom migration message; it is neither an alias nor an override. Reject mixed inputs even when the values agree. Both paths use the resolved response cap without a separate excited-state default. |
| Named SCF `maxiter [N1,N2,...]`, internal `maxiter_per_rung` | Per-protocol caps supported for named moldft; named response rejects the list. `restart.start_rung` is supplied by the adapter. | Keep scalar input. Per-protocol vectors may remain internal; dropped named-stage list syntax is not introduced into workflow groups. Protocol entry stays the planner's responsibility, not a user-set rung index. |
| `response seed.start_rung`; test-driver `--seed-start-rung` | `coarse` keeps the ladder; `fine` with a DALTON seed collapses it to its finest threshold before response planning. The test driver shares `apply_seed_start_rung`; this is not a new public standalone capability. | **Approved 2026-11-16:** remove both controls and ladder rewriting with response adoption; reject every explicit value through the generic unknown-input path. Keep seed construction. The user chooses a single-entry protocol for direct fine refinement or an explicit ladder; every new seed starts at the first requested entry. Saved checkpoint progress remains the restart planner's responsibility. |
| `excited.guess_max_iter`, `seed.es_warmup` | Explicit deck settings affect a separate TDA warmup. The executor's cold-start warmup default is 10, different from the registered guess default of 5; seeded warmup defaults off. | Warmup is disabled in the initial design. Retain guess generation/seeding; remove both listed molresponse keys entirely, with generic unknown-input errors for every explicit value. Any future guess/warmup iteration allowance must remain separate from main `maxiter`. |
| DFT `protocol`, `k`; standalone response `--protocol` | Workflow response inherits the DFT protocol; an ordinary response-level `protocol` is unknown. Standalone response has its own protocol flag. | Preserve madqc ownership; retire the standalone flag with its executable. Final threshold and k are mandatory; keep per-state saved progress. Do not add a response-block protocol key merely because the shared controller accepts a protocol. |
| `nvalpha`, `nvbeta`, `nv_extra`, `nv_step`, `nv_its`, `freeze_occupied`; `aocc`, `bocc` | PR #798 adds coupled/frozen virtual solves, explicit occupations, atomic-guess padding and a count-reduction schedule. `nv_its` supplies the intermediate cap; the last membership uses `maxiter`. The current virtual-only loop counts evaluations and has a separate `5*dconv` gate. | Preserve these SCF inputs and numerical algorithms. Bring both loops into shared measurement/restart/checkpoint handling; use `nv_its` as a scoped preparation-update cap and normal `maxiter` for the final membership, with terminal evaluation and fresh restart allowances. Record both spin counts, preparation phase, occupations and fixed-reference identity. `nv_factor` remains retired. Physical fixed-reference change metrics are inapplicable to virtual convergence; no new metric aliases or automatic warmup. |
| `optimization maxiter`; inner correlated-solver iteration controls | These count geometry steps or distinct inner/outer numerical solves. | Preserve their scopes. They are not aliases for the electronic-state main-update budget. |
| `optimization gtol`, `gradient_precision` | The driver derives `gtol=max(1e-4,dconv)` and `gradient_precision=max(1e-6,0.1*dconv)`, unless explicitly overridden. At the default `dconv=1e-4`, these resolve to `1e-4` and `1e-5`. | **Approved 2026-11-13:** independent defaults `gtol=1e-4`, `gradient_precision=1e-5`. Explicit optimization-group file/CLI overrides retain their precedence. Electronic criteria/protocol changes do not alter these two defaults. Preserve other optimizer defaults and its existing precision clamp, Hessian-update and stopping behavior; replacement is deferred. |

Sources: [virtual-orbital PR #798](https://github.com/m-a-d-n-e-s-s/madness/pull/798),
[SCF scheduling](src/madness/chem/SCF.h), [SCF solvers](src/madness/chem/SCF.cc),
[optimizer driver](src/madness/chem/Drivers.hpp),
[optimizer precision handling](src/madness/chem/molopt.h),
[parameter defaults](src/madness/chem/ParameterManager.hpp),
[response parameter definitions](src/madness/chem/ResponseParameters.hpp),
[explicit ES overrides](src/apps/molresponse/deck_settings.hpp),
[executor budgets/warmup](src/apps/molresponse/calc/calc_executor.hpp),
[standalone response CLI](src/apps/molresponse/main.cpp).

For a well-converged basis-set seed, recommend a single desired final protocol,
for example DFT-group `protocol [1e-6]`. An explicit `protocol [1e-4, 1e-6]`
requests both stages even with a seed. These are illustrative thresholds, not
new defaults. Seed selection does not change the requested ladder or supply MRA
convergence evidence; matching checkpoint restart retains its own saved progress.

## Restart and guessing

| Existing option | What it does now | Mapping to the new design |
| --- | --- | --- |
| DFT `restart auto` | Default. Chooses load, iterate or a fresh/seeded start from archive metadata. Some bad/incompatible archives warn and fall back to a guess. | Keep automatic reuse for compatible saved work. An expected bad/missing checkpoint is a strict failure; a genuinely new task may start fresh. Saved measurements decide reuse versus further work under the current request. |
| DFT `restart none` | Ignores saved ground orbitals and starts from a guess. It is not a response-wide fresh-start flag. | Preserve explicit fresh-start intent for its application; fresh initialization still writes endpoint checkpoints. Register the same scoped mode in the response group during adoption; ground-only settings do not override it. New ground reference generations invalidate downstream evidence normally. |
| `--continue` | Reopens a named-stage output root and runs per-state checkpoint planning. The workflow path returns from named dispatch before this handling. | Preserve the CLI entry and connect it to the retained workflow's shared restart planner. Require a recognized run at the selected prefix; reject explicit fresh/seed/inspection conflicts. Auto can also reuse recognized runs. New tasks start normally; no named-stage syntax or separate scientific policy returns. |
| DFT `restart iterate` | Requires a readable archive and requests continued work even when its old header claims convergence. | The shared solver policy is reuse-or-iterate: reuse qualifying saved work or resume an unmet request, with fresh KAIN when numerical work is needed. Do not retain a separate forced-iteration policy through this old spelling. Require an existing checkpoint for the existing solve, then apply the common policy; newly added tasks start normally. |
| DFT `restart read_only` | Loads without iteration after basic compatibility checks, even if the saved precision is insufficient; it can return a saved energy for changed Hamiltonian settings with a warning. | **Approved 2026-11-12:** keep read-only for inspecting results and plotting saved functions, including unfinished work. Check readability and operation requirements, not new convergence tolerances. Preserve saved context/status and perform no solver iterations or checkpoint rewrites. Normal solver requests reuse qualifying work or resume solving. Route the old spelling to inspection before solver dispatch. Reject explicit requests for new scientific computation in inspection mode; default solver selection alone does not count as such a request. |
| DFT `restart ao`, `restart nwchem`; existing DALTON/basis/frequency seed controls | Choose initial guesses. Supported sources vary by engine. Ground `restartao` and `no_compute` are already retired with migration errors. | Preserve supported guess sources and explicit source failures. Seeds do not inherit convergence or protocol progress; use the agreed first-requested-protocol fallback for new/unmatched independent states. Do not revive retired ground flags or imply that another engine's `no_compute` is already retired. |
| Response saved-state reuse | Driven states and ES bundles have their own loaders and scheduler decisions; explicit per-state plans are currently passed by named continuation. The ordinary response parameter group has no matching general `restart` key. | Carry those plans into the retained workflow. Reuse final independent states whose valid saved measurements satisfy the current request; advance completed intermediate work, resume unfinished work at its own protocol, and keep coupled roots together. Register response `restart auto/none/iterate` and route them to the same task planner, with auto as default. This is new wiring, not a renamed DFT-only option. |

Sources: [ground restart modes and decisions](src/madness/chem/RestartPlan.h),
[named continuation](src/apps/madqc/RunNamedStages.hpp),
[response handoff](src/apps/molresponse/madqc_adapter.hpp),
[response state plans](src/apps/molresponse/calc/calc_manager.hpp),
[ground retired-key errors](src/madness/chem/CalculationParameters.h).

## Nemo nuclear-response controls

| Existing option or behavior | Mapping to the shared design |
| --- | --- |
| `hessian`, `purify_hessian` | Preserve Nemo-owned property selection and supported purification. CPHF tasks depend on the eligible final ground reference; direct Hessian/vibrational assembly waits for every required response. |
| `read_cphf=true`, `restart_cphf=true` | Normalize either (or both) to shared `restart iterate` semantics for CPHF tasks. Loaded functions never imply convergence: reuse qualifying evidence or iterate. Otherwise inherit the resolved application restart policy. Reject explicit fresh/seed or inspection conflicts; keep the common rule for newly added tasks. Neither flag retains a private file-selection/convergence policy or requires a legacy-file importer. |
| Hardcoded CPHF cap of ten evaluations | Use the Nemo application's resolved `maxiter` separately for each task/protocol, with main-update counting, terminal evaluation and a fresh restart allowance. No new CPHF iteration knob is needed initially. |
| Automatic preparatory pass labelled LDA | Disable it under the agreed warmup policy. Preserve direct guesses; ordinary stages use the requested model and shared protocol/control rules. |
| `RMS/norm(xi)` and Hessian-row change stopping | Replace implicit residual scaling with the shared maximum and optional relative policy after representation confirmation. Total-energy change is inapplicable. Florian must confirm response density/kinetic definitions and whether the Hessian-row change remains a scientific acceptance check or a diagnostic; no silent relabeling or removal. |

Sources: [Nemo controls](src/madness/chem/nemo.h),
[response solves and property assembly](src/madness/chem/nemo.cc),
[madqc Nemo binding](src/madness/chem/Applications.hpp).

## Checkpointing and archive controls

| Existing option | What it does now | Mapping to the new design |
| --- | --- | --- |
| DFT `save` | Boolean orbital/archive output control. Nemo saves during iterations when enabled; ordinary application paths use it for archive publication. It is not a shared checkpoint interval. | **Approved 2026-11-12:** remove the `save` input entirely, for both true and false. Use the normal fatal unknown-input error only; no alias, deprecated-key registration, migration hint or warning-and-ignore path. Endpoint/valid controlled-stop checkpoints remain required. |
| `checkpoint.every_iterations N` | Named moldft only: optional interior publication every N iterations, plus endpoints. Omission defaults to endpoints; explicitly supplied N must be positive. Named molresponse refuses the key because it only publishes at solve boundaries. | Reuse this control in the existing application groups, with shared semantics: omission means endpoints; N=1 means every update, larger N means periodic checkpoints. Response needs an evaluated-state callback inside its driver. This is adoption work, not existing response support. |
| `checkpoint.keep_all` | Named debugging escape that retains every generation. Default retention removes older unneeded generations. | **Approved 2026-11-12:** remove this input entirely, for both true and false, with the normal fatal unknown-input error only. Retain latest usable checkpoints plus necessary temporary/referenced generations. No keep-all input or special migration handling remains in the adopted interface. |
| `io backend native/hdf5`, `response hdf5`, standalone response `--hdf5` | Backend selection exists. HDF5 currently affects supported response I/O and requires the build feature; non-response workflow requests are presently reported as ignored. | **Approved 2026-11-13:** native is the sole shared restart backend. Adopted unified execution paths reject requests selecting HDF5 as unsupported, even with HDF5 build support; no silent fallback, conversion or HDF5 mirror is required. Preserve existing separate HDF5 utilities/build support. Shared HDF5 and new analysis/export integration are future work outside parameter/I/O unification. |
| Prefix/output controls and `--restart_info` | Select locations or inspect an archive; standalone response has `--calc-dir`. | Preserve madqc entry options, route paths through the new run layout, and adapt inspection to the shared records; standalone-only path flags retire with their executables. Reading a result JSON file cannot substitute for checkpoint validation. |

Sources: [named checkpoint inputs](src/apps/madqc/NamedStages.hpp),
[retention helper](src/madness/chem/CheckpointRetention.h),
[SCF publication observer](src/apps/madqc/MoldftCheckpoint.hpp),
[response publication settings](src/apps/madqc/ResponseCheckpoint.hpp),
[Nemo saving](src/madness/chem/nemo.cc),
[workflow archive publication](src/madness/chem/Applications.hpp),
[I/O configuration](src/madness/chem/ParameterManager.hpp),
[backend routing](src/apps/madqc/madqc.cpp).

## Other engines: preserve the meaning of inner controls

These are declared inputs, not a completed audit of every solver branch.
Shared lifecycle adoption does not make all similarly named parameters aliases.

| Engine | Existing controls needing adapter-specific mapping |
| --- | --- |
| Standalone MP2 | `econv`, `dconv`, `thresh`, `maxiter`, boolean `restart`, `no_compute`. Unwritten `dconv` derives from `0.1*sqrt(econv)` and `thresh` from `econv`. [Definition](src/madness/chem/mp2.h). |
| CC2 family | `econv`, `econv_pairs`, `dconv_3d`, `dconv_6d`, `iter_max`, `iter_max_3d`, `iter_max_6d`, boolean `restart`. Keep singles/pairs and inner/outer scopes distinct. [Definition](src/madness/chem/CCParameters.h). |
| TDHF/CIS | `econv`, `dconv`, final `maxiter`, `guess_maxiter`, guess tolerances, and `restart no_restart/iterate/no_compute`. [Definition](src/madness/chem/TDHF.h). |
| OEP | Own `maxiter`, boolean `restart`, `no_compute`, `conv_threshold`, and `save_iter_*` diagnostics. Diagnostic dumps are not automatically complete restart checkpoints. [Definition](src/madness/chem/oep.h). |
| PNO | `maxiter_t`, `maxiter_micro`, `maxiter_macro`; `econv_micro/pairs/macro`; `dconv`; `restart none/mp2/cispd/all`; `save_pnos`. Preserve amplitude/adaptive-function distinctions. [Definition](src/madness/chem/PNOParameters.h). |
| Znemo / Zcis | Znemo shares the underlying ground-state machinery and adds magnetic controls; Zcis declares `maxiter`, `guess_maxiter`, `dconv`. Physical metric definitions still require review. [Znemo](src/madness/chem/znemo.h), [Zcis](src/madness/chem/zcis.h). |
| Dirac | Separate parameter structure with `max_iter`, `dconv`, restart and per-iteration save controls (`do_save` / documented `no_save`, `savefile`). These are structure declarations, not a verified unified input grammar. [Definition](src/apps/dirac/DFParameters.h). |

## Consolidated target and remaining validation

[Design section 1](design_draft.md#input-routing-and-retired-controls) is the
single target contract for routing and removed inputs. The tables above retain
the old behavior only to guide source changes. There is no second legacy policy.

The consolidated revision removes warmup/stall/fixed-work switches and old
criterion-selector switches with generic unknown-input errors for all explicit
values. Migration messages remain for `dconv` and `econv`; `excited.maxiter`,
`save` and `checkpoint.keep_all` use only generic unknown-input errors.
Inspection is separate from reuse-or-iterate solver execution.
`--continue` requires an existing run, while normal auto mode can also reuse it.

Scoped inner-solver mappings are an adapter implementation prerequisite. Exact
numerical defaults and model definitions need scientific validation; this input
contract does not certify them. The [complete example](output_record_examples.md)
uses the new target syntax, which is not yet implemented.
