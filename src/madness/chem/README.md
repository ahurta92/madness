# MADchem blueprint

## System Overview

Reusable MADchem electronic-structure engines, molecular and operator primitives, convergence/checkpoint contracts, and sequential workflow adapters over MADmra.

The application layer owns command-line execution and workflow publication. This library provides the numerical methods, molecular data, operators, and shared planning/record primitives they call. Current engines retain distinct iteration and restart behavior.

## Subdirectories & System Boundaries

Chemistry source is kept directly in this directory. The principal boundaries are:

| Boundary | Responsibility |
| --- | --- |
| [Applications](../../apps/README.md) | Executable startup, workflow composition and per-application I/O. |
| [madqc](../../apps/madqc/README.md) | Named-stage adapters, requested properties, checkpoint publication and downstream release. |
| [Molecular response](../../apps/molresponse/README.md) | Driven response, coupled TDA/RPA bundles, state scheduling and response properties. |
| [MRA vector operations](../mra/vmra.h) | Distributed adaptive-function operations used by the chemistry engines. |

## Local Interface & Invariants

### Interface

- MADchem is the reusable chemistry component over mra. CMake also builds plotxc outside libraries-only builds and a benchmark_exchange_operator utility only in testing/application builds. It selects exactly one XC implementation according to Libxc::xc presence. [CMakeLists.txt:add_mad_library](./CMakeLists.txt) [CMakeLists.txt:MADCHEM_SOURCES](./CMakeLists.txt) [CMakeLists.txt:plotxc](./CMakeLists.txt) [CMakeLists.txt:benchmark_exchange_operator](./CMakeLists.txt)
- SCF exposes conventional molecular HF/DFT orbital optimization, density/Fock/operator evaluation, analytic gradients, dipoles, population analysis and orbital restart I/O. Nemo wraps SCF data with regularized orbitals F=psi/R and R-squared metric operations. [SCF.h:SCF](./SCF.h) [nemo.h:Nemo](./nemo.h) [nemo.h:NemoBase](./nemo.h)
- CalculationParameters and ParameterManager supply typed input groups and derived values. Prefix precedence is command-line then explicit dft deck prefix then input stem then mad. The dft prefix is the source of truth. [CalculationParameters.h:CalculationParameters](./CalculationParameters.h) [ParameterManager.hpp:ParameterManager](./ParameterManager.hpp) [ParameterManager.hpp:IOParameters](./ParameterManager.hpp)
- Application run/results and Driver execute compose sequential workflows through StepContext: molecule, distinct SCF and Nemo live references, absolute archive paths and JSON results. CC2/TDHF/OEP consumers validate their already-constructed reference before using context. ResponseApplication is a generic bridge whose implementation is supplied in the application layer. [Applications.hpp:Application](./Applications.hpp) [Applications.hpp:StepContext](./Applications.hpp) [Applications.hpp:check_context_matches_reference](./Applications.hpp) [Applications.hpp:ResponseApplication](./Applications.hpp) [Drivers.hpp:Workflow](./Drivers.hpp)
- WorkflowBuilders is an internal convenience factory, not installed among MADCHEM_HEADERS; response construction redirects to the application layer. CCLib has only a cc_lib label stub, MoldftLib.hpp is empty, and QCPropertyInterface has private placeholder virtual methods. [WorkflowBuilders.hpp:add_workflow_drivers](./WorkflowBuilders.hpp) [CMakeLists.txt:MADCHEM_HEADERS](./CMakeLists.txt) [CCLib.hpp:cc_lib](./CCLib.hpp) [QCPropertyInterface.h:QCPropertyInterface](./QCPropertyInterface.h)
- Raw metric primitives are separate from stopping policy: fixed-point function L2 residual components, physical charge/spin density-change norms, kinetic gradient matrices/diagonals, and maximum/trace/trace-change reductions. [FunctionResidual.h:function_residual_norms](./FunctionResidual.h) [DensityChanges.h:density_change_norms](./DensityChanges.h) [KineticEnergy.h:kinetic_energy_diagonal](./KineticEnergy.h) [FunctionMetricSummary.h:function_metric_maximum](./FunctionMetricSummary.h) [FunctionMetricSummary.h:function_metric_trace_change](./FunctionMetricSummary.h)
- ExplicitCriteria compares only user-written requested targets; ConvergenceSettings is the shared source for tables/settings/outcome JSON; LadderPlanner is a pure continuation/seed planner over saved evidence and current requested rungs. [ExplicitCriteria.h:evaluate_explicit_criteria](./ExplicitCriteria.h) [ExplicitCriteria.h:requested_accuracy_targets](./ExplicitCriteria.h) [ConvergenceSettings.h:rung_exit](./ConvergenceSettings.h) [LadderPlanner.h:plan_ladder](./LadderPlanner.h)
- Saved states live in node folders: NodeArchive replaces one folder's archive and description together, with no lock and no generation. RootCoordination holds the collective rank-zero helpers and the checked archive close; RunFiles writes results.json and timings.json at the top of the run folder and refuses a run folder from an older version (one with checkpoints/). [NodeArchive.h:Node](./NodeArchive.h) [RootCoordination.h:close_archive](./RootCoordination.h) [RunFiles.h:write_results_entry](./RunFiles.h) [RunFiles.h:write_timing_entry](./RunFiles.h)
- RestartMetadata and RestartPlan are the legacy orbital header and source-selection policy, separate from NodeArchive and LadderPlanner. Source choices include initial guess, restartdata, AO projections and NWChem. [Restart.h:RestartMetadata](./Restart.h) [RestartPlan.h:plan_restart](./RestartPlan.h) [SCF.cc:get_initial_orbitals](./SCF.cc)
- Molecule and AtomicBasisSet own geometry/atomic metadata and Gaussian basis data; slymer ESInterface/NWChem_Interface read external electronic-structure results; MolecularOrbitals packages orbital functions with energies, occupations, symmetry and localization metadata. [molecule.h:Molecule](./molecule.h) [molecularbasis.h:AtomicBasisSet](./molecularbasis.h) [ESInterface.h:ESInterface](./ESInterface.h) [NWChem.h:NWChem_Interface](./NWChem.h) [MolecularOrbitals.h:MolecularOrbitals](./MolecularOrbitals.h)
- SCFOperatorBase supplies application/matrix interfaces for local, nuclear, Coulomb, exchange, kinetic and XC operators; Kinetic and Fock have intentionally unsupported direct-function paths and should not be assumed to support every base overload. [SCFOperators.h:SCFOperatorBase](./SCFOperators.h) [SCFOperators.h:Kinetic](./SCFOperators.h) [SCFOperators.h:Fock](./SCFOperators.h)
- DispersionCorrection supports none, d3bj and d3zero with optional ATM through simple-dftd3. It adds molecular coordinate-only energy, gradient and finite-difference Hessian. PCMParameters exists without PCMSolver, but PCM construction/evaluation throws when that dependency is absent. [dispersion.cc:DispersionCorrection](./dispersion.cc) [dispersion.cc:hessian](./dispersion.cc) [dispersion.cc:reject](./dispersion.cc) [pcm.cc:PCMParameters](./pcm.cc) [pcm.cc:PCM](./pcm.cc)
- PopulationAnalysis supplies Mulliken/Lowdin/IAO analysis and Result structs serialize solver/property output; ResultsSummary presents existing JSON, while VizManifest inventories plot files rather than generating field data. [population.h:PopulationAnalysis](./population.h) [Results.h:SCFResults](./Results.h) [ResultsSummary.hpp:write_results_summary](./ResultsSummary.hpp) [VizManifest.hpp:write_viz_manifest](./VizManifest.hpp)
- OEP approximates exact exchange with a local potential over a converged, closed-shell Nemo reference. Its supported model keywords are oaep, ocep, dcep and mrks; it is not an independent general HF/DFT starting-point solver. [oep.h:OEP_Parameters](./oep.h) [oep.h:value](./oep.h) [oep.h:check_converged](./oep.h)
- Localizer provides pm, boys, new/new_sys and cholesky transformations, with optional R-weighted metric and core/valence separation. Cholesky uses deterministic density pivots and optional pivot hysteresis for stable continuation; it does not optimize the PM/Boys localization objective. [localizer.h:Localizer](./localizer.h) [localizer.cc:localize](./localizer.cc) [localizer.cc:localize_cholesky](./localizer.cc)
- GaussianFunction in the slymer import layer represents Cartesian or spherical Gaussian combinations; AtomicBasisSet is the separate native molecular basis interface. BasisFunction supplies an abstract importer-side function boundary. [gaussian.h:GaussianFunction](./gaussian.h) [gaussian.cc:GaussianFunction](./gaussian.cc) [molecularbasis.h:AtomicBasisSet](./molecularbasis.h) [basis.h:BasisFunction](./basis.h)
- NuclearCorrelationFactor projects regularization potentials for Nemo, while Diamagnetic_potential_factor supplies the complex magnetic path with a Gaussian factor and boxed diamagnetic potential. Their formulas and supported derivatives differ by factor family. [correlationfactor.h:NuclearCorrelationFactor](./correlationfactor.h) [correlationfactor.cc:create_nuclear_correlation_factor](./correlationfactor.cc) [diamagneticpotentialfactor.h:Diamagnetic_potential_factor](./diamagneticpotentialfactor.h)
- Geometry optimization helpers and vibrational analysis consume energy/gradient/Hessian boundaries: MolOpt performs quasi-Newton optimization, and compute_frequencies mass-weights Hessians with optional translation/rotation projection and negative signed frequencies for imaginary modes. [molopt.h:MolOpt](./molopt.h) [molecular_optimizer.h:MolecularOptimizer](./molecular_optimizer.h) [vibanal.cc:compute_frequencies](./vibanal.cc)
- GridBuilder wraps IntegratorXX molecular quadrature when compiled with that dependency and throws when unavailable. Compiled asymptotic-correction helpers exist, but the SCF constructor sets ac none and Nemo do_ac returns false. [IntegratorXX.h:GridBuilder](./IntegratorXX.h) [AC.h:AC](./AC.h) [AC.cc:AC](./AC.cc) [nemo.h:do_ac](./nemo.h)
- TDHF::solve_cis implements the CIS/TDA path; TDHF::solve_tdhf and make_tdhf_potentials explicitly throw NOT IMPLEMENTED. The application-layer response workflow is a separate implementation. CC2, PNO and MP3 engine declarations live here; their implementations are separate from the CIS/TDA response boundary. Znemo/Zcis provide the complex magnetic HF/CIS path. [TDHF.cc:solve_cis](./TDHF.cc) [TDHF.cc:solve_tdhf](./TDHF.cc) [TDHF.cc:make_tdhf_potentials](./TDHF.cc) [Applications.hpp:ResponseApplication](./Applications.hpp) [CC2.h:CC2](./CC2.h) [PNO.h:PNO](./PNO.h) [mp3.h:MP3](./mp3.h) [znemo.h:Znemo](./znemo.h) [zcis.h:Zcis](./zcis.h)
- CC2::solve dispatches mp2, mp3, cc2, cis, cispd, adc2 and lrcc2 calculations from CCParameters and returns a tasks JSON array; its explicit result map currently accumulates MP2, MP3, CC2 and LRCC2 excitation energies. [CC2.cc:CC2::solve](./CC2.cc) [CCParameters.h:CCParameters::calc_type](./CCParameters.h)
- CCPotentials supplies ground/excited pair construction, regularized Fock/commutator potentials, singles diagrams, pair and correlation-energy contractions, orbital bra/ket views and intermediate caching. [CCPotentials.h:CCPotentials](./CCPotentials.h)
- CCFunction and CC_vecfunction annotate 3D functions with occupied index/type, response energy and errors; CCPair contains pair labels, cusp-free 6D function, regularization components and constant part. [ccpairfunction.h:CCFunction](./ccpairfunction.h) [CCStructures.h:CC_vecfunction](./CCStructures.h) [CCStructures.h:CCPair](./CCStructures.h)
- CCPairFunction represents explicit high-dimensional or separated sum-of-products components, optionally multiplied by a two-particle operator, with projections, partial inner products, particle exchange and consolidation. [ccpairfunction.h:CCPairFunction](./ccpairfunction.h)
- LowRankFunction and LowRankFunctionFactory represent and fit sums g_i(1)h_i(2) using sampled functors, orthonormalization, contractions and alternating optimization. [lowrankfunction.h:LowRankFunction](./lowrankfunction.h) [lowrankfunction.h:LowRankFunctionFactory](./lowrankfunction.h)
- PNOPairs groups pair-specific PNOs, amplitudes, Fock/fluctuation matrices, overlap caches, frozen flags and PairEnergies; ElectronPairIterator enumerates packed lower-triangle pairs in active occupied indices. [PNOStructures.h:PNOPairs](./PNOStructures.h) [PNOStructures.cpp:ElectronPairIterator::operator ++](./PNOStructures.cpp)
- PNOParameters and F12Parameters configure PNO rank/adaptivity, thresholds, pair masks, restart/compute controls, exchange approximations and projected/Hylleraas F12 energy paths. [PNOParameters.h:PNOParameters](./PNOParameters.h) [PNOParameters.h:F12Parameters](./PNOParameters.h)
- F12Potentials supplies regularized fluctuation potentials, f12/Coulomb/exchange integrals and projected/Hylleraas energies; BasisFunctions creates Gaussian, partial-wave, predefined, file-based and excitation-operator guess orbitals. [PNOF12Potentials.h:F12Potentials](./PNOF12Potentials.h) [PNOGuessFunctions.h:BasisFunctions](./PNOGuessFunctions.h)
- MP2 is an OptimizationTargetInterface/QCPropertyInterface for correlation energies built on a HartreeFock wrapper around Nemo; ElectronPair holds first-order functions, singlet/triplet energies and restart iteration/convergence metadata. [mp2.h:MP2](./mp2.h) [mp2.h:HartreeFock](./mp2.h) [mp2.h:ElectronPair](./mp2.h)
- CorrelationFactor constructs Slater F12 and transformed electronic potentials; CorrelationFactor2 supplies an alternative factor fixed to gamma=0.5. [electronic_correlation_factor.h:CorrelationFactor](./electronic_correlation_factor.h) [electronic_correlation_factor.h:CorrelationFactor2](./electronic_correlation_factor.h)
- PNO::solve operates on a caller-owned vector of PNOPairs and dispatches MP2 or CIS(D). The solver supports fixed or adaptively enlarged PNO subspaces, orbital optimization, canonicalization and density-eigenvalue truncation. [PNO.cpp:PNO::solve](./PNO.cpp) [PNO.cpp:PNO::iterate_pairs](./PNO.cpp) [PNO.cpp:PNO::pno_compress](./PNO.cpp)
- MP3 computes the third-order energy contribution from MP2 pair functions and Info, with serial and MacroTask drivers combining CD, EF, GHIJ and KLMN terms. The macro driver rebuilds MP2 low-rank regularization components before distributing full square-index pair data. [mp3.cc:MP3::mp3_energy_contribution](./mp3.cc) [mp3.cc:MP3::mp3_energy_contribution_macrotask_driver](./mp3.cc)
- Standalone MP2::value returns only correlation energy, supports a selected occupied pair or all active pairs, and selects coupled equations for localized orbitals versus refined independent equations for canonical orbitals. [mp2.cc:MP2::value](./mp2.cc)

### Invariants

- Function handles are shallow. The input/raw-update residual overload and density-change helper copy operands to preserve caller-owned trees. The existing-residual overload instead calls norm2s directly; it preserves values but can reconstruct shared functions and change their tree state. Residual components are raw unweighted L2 fixed-point differences with no normalization; invalid/empty/nonfinite components remain invalid, not zero. [FunctionResidual.h:function_residual_norms](./FunctionResidual.h) [FunctionResidualReport.h:function_residual_record](./FunctionResidualReport.h) [DensityChanges.h:density_change_norms](./DensityChanges.h) [../mra/vmra.h:norm2s](../mra/vmra.h)
- Trace metrics are unweighted sums. Trace changes compare sums over matching collection layouts; the caller must establish the same physical state/reference, representation and sampling convention. Unitary rotations within the compared collection need no orbital-by-orbital alignment. The separate elementwise absolute-change helper does require physical correspondence between functions; equal dimensions alone do not establish it. [FunctionMetricSummary.h:function_metric_trace](./FunctionMetricSummary.h) [FunctionMetricSummary.h:function_metric_trace_change](./FunctionMetricSummary.h) [FunctionMetricSummary.h:function_metric_absolute_change](./FunctionMetricSummary.h)
- Explicit target tests require finite nonnegative measurements, matching units and quantity, and use inclusive <=. Engine checks and requested criteria are evaluated separately; convergence.use_default_criteria selects one decider, not an AND. When callers set final_protocol=false, intermediate rungs consume their supplied iteration cap; final rungs order convergence, cap, then stall. Named SCF and the response paths supply this stage role. Unnamed/standalone SCF without a checkpoint observer treats every rung as final for this decision and can converge early at intermediate protocols. A resumed run receives a fresh cap. [ExplicitCriteria.h:evaluate_explicit_criteria](./ExplicitCriteria.h) [ConvergenceSettings.h:use_default_criteria_key](./ConvergenceSettings.h) [ConvergenceSettings.h:rung_exit](./ConvergenceSettings.h) [SCF.cc:solve](./SCF.cc) [ScfContinuation.h:SavedRung](./ScfContinuation.h)
- Progress detection requires supplied explicit targets and a valid window/improvement setting, cannot invent rungs or manufacture engine success, and history must reset across representation/reference changes. Saved evidence is rejudged against current requests; diverged records cannot resume. [ExplicitCriteria.h:validate_progress_configuration](./ExplicitCriteria.h) [ExplicitCriteria.h:explicit_progress](./ExplicitCriteria.h) [LadderPlanner.h:plan_ladder](./LadderPlanner.h)
- A node save stages the new archive and description under temporary names; the complete new description is the commit point, and the next inspect finishes or discards an interrupted save. Archives are closed collectively; every rank closes, including nonwriters. [NodeArchive.h:Node](./NodeArchive.h) [RootCoordination.h:close_archive](./RootCoordination.h)
- Legacy orbital format writes version 5 and reads versions 4/5; version-4 missing evidence remains unknown. Representation numeric tags mo=1,nemo=2,znemo=3 cannot be renumbered. MO and Nemo archives are not interchangeable. Requested geometry and Hamiltonian settings are retained on load; a representation/cell mismatch throws and reprojection invalidates saved convergence evidence. [Restart.h:RestartMetadata](./Restart.h) [Restart.h:Representation](./Restart.h) [SCF.cc:load_mos](./SCF.cc)
- Explicit missing restart sources fail rather than silently start fresh; automatic load failures are agreed collectively and fall back loudly to rung zero. NWChem import is a distinct alignment path that can translate/rotate the current molecule to the imported frame. [RestartPlan.h:plan_restart](./RestartPlan.h) [SCF.cc:get_initial_orbitals](./SCF.cc) [SCF.cc:initial_guess_from_nwchem](./SCF.cc)
- SCF protocol changes mutate FunctionDefaults and rebuild parameter-bound operators. Nemo regularized orbitals normalize in the R-squared metric and must rebuild/reproject nuclear factors when threshold or k changes. [SCF.h:set_protocol](./SCF.h) [nemo.h:NemoBase](./nemo.h) [nemo.h:set_protocol](./nemo.h)
- SCF explicit target measurement is enabled through the named-stage checkpoint observer. Its checkpoint event energy/density describe an input state and may lag saved updated orbitals; the observer receives phase/comparison information to evaluate the saved state separately. A criteria-selected stop cannot create an engine convergence claim. [SCF.cc:solve](./SCF.cc) [SCF.h:IterationCheckpoint](./SCF.h) [SCF.cc:evaluate_restricted_hf](./SCF.cc)
- Occupations are validated integral 0/1 occupied-span choices; holes require canonical orbitals. Localization and symmetry are mutually exclusive. Retired gopt/restartao/no_compute/nv_factor user keys throw migration errors. Population analysis requires supported all-electron configurations. [CalculationParameters.h:set_derived_values](./CalculationParameters.h)
- Meta-GGA XC requires occupation-weighted tau before potential construction. Weak-form GGA needs an explicit allow_weak_form caller opt-in and weak_gga; multiplicative and flux terms travel together. The lda-only backend supports HF/LDA but no fxc_apply response kernel. [SCFOperators.h:XCOperator](./SCFOperators.h) [SCF.cc:apply_potential](./SCF.cc) [xcfunctional_ldaonly.cc:initialize](./xcfunctional_ldaonly.cc) [xcfunctional_ldaonly.cc:fxc_apply](./xcfunctional_ldaonly.cc)
- XC argument slot numbering is a fixed indexed internal contract with reserved holes; nuclear cusp factors are evaluated pointwise at quadrature rather than projected into cusp-product trees. [xcfunctional.h:xc_arg](./xcfunctional.h) [xcfunctional.h:xc_kernel_apply](./xcfunctional.h)
- Workflow execution changes process cwd through ScopedCWD, so drivers run sequentially. OptimizeDriver needs analytic gradients, rejects initial_hessian, and does not checkpoint optimizer history/Hessian. [PathManager.hpp:ScopedCWD](./PathManager.hpp) [Drivers.hpp:Workflow](./Drivers.hpp) [Drivers.hpp:OptimizeDriver](./Drivers.hpp)
- Nemo::value rejects user-supplied progress-window/improvement controls, accuracy.bsh_residual, accuracy.kinetic_trace_change, and convergence.use_default_criteria=false. The shared criteria vocabulary therefore does not imply every engine accepts those settings. Nemo retains its own energy/density/BSH convergence rules. [nemo.cc:value](./nemo.cc) [nemo.cc:check_convergence](./nemo.cc) [nemo.h:Nemo](./nemo.h)
- SCF engine convergence tests per-spin density change against a molecule-size-scaled effective density tolerance and the previous orbital-update BSH residual against five times that tolerance unless density-only is requested. A saved-state raw L2 metric is not interchangeable with that lagged engine check. [SCF.cc:solve](./SCF.cc) [SCF.h:IterationCheckpoint](./SCF.h)
- Nemo normalizes regularized orbitals in the R-squared metric. Its ground-state path supports meta-GGA tau, but the CPHF operator construction rejects meta-GGA; the older CPHF implementation is a separate closed-shell LDA/GGA path. [nemo.cc:solve](./nemo.cc) [nemo.cc:make_fock_operator](./nemo.cc) [nemo.cc:solve_cphf](./nemo.cc)
- The Libxc backend rejects Laplacian-dependent and range-separated hybrid functionals. fxc_apply supports the restricted LDA/GGA response path and rejects meta-GGA, despite broader ground-state energy/potential support. [xcfunctional_libxc.cc:initialize](./xcfunctional_libxc.cc) [xcfunctional_libxc.cc:fxc_apply](./xcfunctional_libxc.cc)
- The SCF Exchange constructor removes zero-occupation orbitals and rejects fractional occupations. Explicit set_bra_and_ket expects an already-conjugated bra and copies both inputs; the Nemo metric is carried by its bra. Owner-pinned exchange expects one subworld per rank, collective finalize participation even on idle ranks, and cleanup of cached Functions/reducers before their Worlds die. Final truncation follows accumulation of all tiles. [exchangeoperator.cc:Exchange](./exchangeoperator.cc) [exchangeoperator.cc:set_bra_and_ket](./exchangeoperator.cc) [exchangeoperator.h:multiworld_efficient](./exchangeoperator.h) [exchangeoperator.h:finalize](./exchangeoperator.h) [exchangeoperator.h:truncate](./exchangeoperator.h)
- Localizer core/valence enforcement requires Fock information and rejects inconsistent block separation. Its localize result marks the returned orbital span occupied, so it should not be treated as a generic occupation-preserving transform. [localizer.cc:localize](./localizer.cc) [localizer.cc:enforce_core_valence_separation](./localizer.cc)
- Molecule coordinates are stored in atomic units. Geometry groups default to atomic units while XYZ input is converted from angstroms; orientation/symmetrization may alter coordinates unless no_orient is set. NWChem import handles its own frame and units. [molecule.cc:read_file](./molecule.cc) [molecule.cc:read_xyz](./molecule.cc) [molecule.cc:orient](./molecule.cc) [NWChem.cc:read](./NWChem.cc)
- AtomicBasisSet::read_file resolves basis files through MRA_CHEMDATA_DIR, including its environment override. CorePotentialManager::read_file uses MRA_DATA_DIR/coredata. GTHPseudopotential::make_pseudo_potential opens gth.xml in the working directory. These are separate lookup contracts. [molecularbasis.cc:read_file](./molecularbasis.cc) [corepotential.cc:read_file](./corepotential.cc) [gth_pseudopotential.h:make_pseudo_potential](./gth_pseudopotential.h)
- NuclearCorrelationFactor holds a reference to Molecule; that molecule must outlive it. initialize rebuilds the projected factor potentials, and some factor families deliberately throw for unsupported higher derivatives rather than providing general Hessian support. [correlationfactor.h:NuclearCorrelationFactor](./correlationfactor.h) [correlationfactor.h:initialize](./correlationfactor.h) [correlationfactor.h:poly4erfc](./correlationfactor.h)
- Znemo requires no_orient and disallows simultaneously using an explicit magnetic Gaussian factor and a nontrivial nuclear correlation factor. It rejects dispersion; gradients require zero explicit field. Zcis uses its own CIS iteration/checkpoint conventions rather than the stage-checkpoint protocol. [znemo.cc:recompute_factors_and_potentials](./znemo.cc) [znemo.cc:value](./znemo.cc) [znemo.cc:gradient](./znemo.cc) [zcis.cc:iterate](./zcis.cc) [zcis.cc:save_guess](./zcis.cc)
- OEP rejects dispersion and requires a converged reference; its local potential and orbital iteration use the Nemo framework and its own convergence rules. The legacy no_compute branch loads orbitals but returns zero instead of a recovered stored energy. [oep.h:value](./oep.h) [oep.cc:iterate](./oep.cc) [oep.cc:load_restartdata](./oep.cc)
- projector_irrep implements a finite supported Abelian point-group set; StrongOrthogonalityProjector temporarily tightens the higher-dimensional FunctionDefaults threshold. Callers cannot treat projector operations as independent of global MRA defaults. [pointgroupsymmetry.h:projector_irrep](./pointgroupsymmetry.h) [pointgroupsymmetry.cc:projector_irrep](./pointgroupsymmetry.cc) [projector.h:StrongOrthogonalityProjector](./projector.h)
- CCPair functions[0] must be a pure 6D component; complete_pair_with_low_rank_parts requires exactly that one component and supplies reconstructible tails for MP2/CC2/LRCC2. [CCStructures.cc:CCPairBuilder::make_pair](./CCStructures.cc) [CCStructures.cc:CCPairBuilder::complete_pair_with_low_rank_parts](./CCStructures.cc)
- CCFunction, CC_vecfunction, CCPairFunction copy construction and LowRankFunction copies share underlying functions; explicit copy helpers clone. CCPairFunction constructors from separated vectors deep-copy inputs, while constructors from pure Function handles are shallow. [ccpairfunction.h:CCFunction](./ccpairfunction.h) [ccpairfunction.h:CCPairFunction](./ccpairfunction.h) [CCStructures.h:CC_vecfunction](./CCStructures.h) [lowrankfunction.h:LowRankFunction](./lowrankfunction.h)
- Convolution bra states must already carry any nuclear-correlation metric; CC potential bra orbitals are R² times Nemo kets, and exchange includes frozen occupied orbitals. [CCPotentials.h:CCPotentials::get_active_mo_bra](./CCPotentials.h) [CCPotentials.h:CCPotentials::make_bra](./CCPotentials.h) [CCPotentials.h:CCPotentials::K](./CCPotentials.h) [CCStructures.h:CCConvolutionOperator](./CCStructures.h)
- CC2 initializes frozen-core count from Fock/localization data, verifies core/valence separation and can relocalize the Nemo orbitals; response calculations check agreement with TDHF frozen count. [CC2.cc:CC2::solve](./CC2.cc) [CC2.cc:CC2::enforce_core_valence_separation](./CC2.cc)
- PNOPairs::verify rejects inconsistent pair-array sizes, unknown pair type and inconsistent CIS data. CISData requires nonnegative excitation number, positive omega and nonempty x. [PNOStructures.cpp:PNOPairs::is_consistent](./PNOStructures.cpp) [PNOStructures.h:CISData::initialized](./PNOStructures.h)
- PNO packed overlap tensors require i>=j and in-range i,j,k, and implement transpose symmetry when addressing equivalent blocks. [PNOTensors.h:Tensor_IJ_IK](./PNOTensors.h) [PNOTensors.h:Tensor_IJ_KJ](./PNOTensors.h)
- CorrelationFactor::apply_U requires a modified nonstandard screening operator; U2 explicitly throws. CorrelationFactor2 is constrained to gamma=0.5. [electronic_correlation_factor.h:CorrelationFactor::apply_U](./electronic_correlation_factor.h) [electronic_correlation_factor.h:CorrelationFactor::U2](./electronic_correlation_factor.h) [electronic_correlation_factor.h:CorrelationFactor2](./electronic_correlation_factor.h)
- Low-rank molecular-grid selection implements random and 3D DFT grids; cartesian is accepted by the parameter enum but not molecular_grid. LowRankFunction orthonormalize/reorthonormalize implement canonical/cholesky and throw for other names. [lowrankfunction.h:molecular_grid](./lowrankfunction.h) [lowrankfunction.h:LowRankFunction::orthonormalize](./lowrankfunction.h)
- BasisFunctions solid-harmonic Gaussian evaluation requires MADNESS_HAS_BOOST; predefined nonempty basis sets reject atoms beyond atomic number 18. guess_with_psi4 throws. [PNOGuessFunctions.cpp:BasisFunctions::SolidHarmonicGaussian::operator ()](./PNOGuessFunctions.cpp) [PNOGuessFunctions.h:BasisFunctions::predefined_guess](./PNOGuessFunctions.h) [PNOGuessFunctions.h:BasisFunctions::guess_with_psi4](./PNOGuessFunctions.h)
- CC pair storage uses i<=j, recovering reversed order by particle exchange; PNO ElectronPairIterator uses i>=j with local active indices and adds freeze to output names. These conventions are distinct. [CC2.cc:CC2::initialize_pairs](./CC2.cc) [CCPotentials.cc:CCPotentials::get_pair_function](./CCPotentials.cc) [PNOStructures.cpp:ElectronPairIterator::operator ++](./PNOStructures.cpp) [PNOStructures.cpp:PNOPairs::name](./PNOStructures.cpp)
- CCPairFunction dimension must split evenly into two particles. Its particle-facing operations use 1/2 while projector variable selection uses 0/1 conventions. Several operator-attached pure projection paths are explicitly unsupported, so representation affects supported operations. [ccpairfunction.h:CCPairFunction](./ccpairfunction.h) [ccpairfunction.cc:CCPairFunction::project_out](./ccpairfunction.cc) [ccpairfunction.cc:CCPairFunction::partial_inner](./ccpairfunction.cc)
- Consolidation first collects matching representation/operator types; supported transformations are op_pure_to_pure, op_dec_to_dec, op_dec_to_pure, dec_to_pure and remove_lindep. Explicit pure-component summation uses sequential fenced reconstructed gaxpy to preserve TT_2D behavior. [ccpairfunction.cc:CCPairFunction::consolidate](./ccpairfunction.cc) [ccpairfunction.cc:CCPairFunction::collect_same_types](./ccpairfunction.cc)
- PNO amplitudes and Fock/fluctuation matrices must match each pair basis rank; rotations transform PNOs, amplitudes, potentials and initialized overlap blocks together. Frozen PNO pairs retain their data and stop optimization, separately from frozen-core orbital exclusion. [PNO.cpp:PNO::t_solve](./PNO.cpp) [PNO.cpp:PNO::transform_pairs](./PNO.cpp) [PNOStructures.cpp:PNOPairs::clear_intermediates](./PNOStructures.cpp)
- PNO adaptive rank growth requires EXOP guesses. Each pair is projected out of occupied and existing PNO spaces, orthonormalized with rank-revealing Cholesky and capped by maxrank; PNO compression retains density eigenvalues above tpno (negative tpno disables truncation). [PNO.cpp:PNO::iterate_pairs](./PNO.cpp) [PNO.cpp:PNO::grow_rank](./PNO.cpp) [PNO.cpp:PNO::orthonormalize_cholesky](./PNO.cpp) [PNO.cpp:PNO::pno_compress](./PNO.cpp)
- PNO orbital update BSH parameters must have nonpositive energy; positive eigenvalues throw. CIS(D) needs active-orbital-matched response functions and nonzero excitation energy. [PNO.cpp:PNO::make_bsh_operators](./PNO.cpp) [PNO.cpp:PNO::iterate_pairs_internal](./PNO.cpp)
- F12 pair-energy arrays distinguish direct and exchange halves, singlet and triplet contributions, and regularized versus explicit-F12 pieces. Projected F12 uses fQg integrals; Hylleraas uses PNOs plus optional complementary auxiliary functions. [PNOF12Potentials.cpp:F12Potentials::compute_fQg_integrals](./PNOF12Potentials.cpp) [PNOF12Potentials.cpp:F12Potentials::compute_projected_f12_energies](./PNOF12Potentials.cpp) [PNOF12Potentials.cpp:F12Potentials::compute_hylleraas_f12_energies](./PNOF12Potentials.cpp) [PNOStructures.h:PairEnergies](./PNOStructures.h)
- Standalone MP2 nuclear-correlated projectors use R²-weighted bra orbitals and Nemo kets and require a restricted reference. Coordinate cache checks use the sum of squared coordinates, and MP2::value(x) requires matching HF geometry checksum. [mp2.cc:MP2::value](./mp2.cc) [mp2.h:HartreeFock](./mp2.h)
- MP3 consumes a complete active-space triangular set of MP2 pairs, reconstructs reversed pairs by exchange and uses R²-weighted contractions; several serial routines obtain World from the first pair and therefore assume a nonempty set. [mp3.cc:MP3::compute_mp3_cd](./mp3.cc) [mp3.cc:MP3::compute_mp3_ef_with_permutational_symmetry](./mp3.cc) [mp3.cc:MP3::compute_mp3_ghij_fast](./mp3.cc)

## Global Side Effects

- Workflow writes prefix/task_i directories and aggregate prefix.calc_info.json after every task, including a failed task record before rethrow. Root writes a temporary file and renames it, broadcasting failure. SCFApplication regenerates input in its isolated library workdir and writes calculation metadata. [Drivers.hpp:Workflow](./Drivers.hpp) [Applications.hpp:SCFApplication](./Applications.hpp) [Applications.hpp:moldft_lib](./Applications.hpp)
- SCF save writes prefix.restartdata parallel shards and, except NWChem-orbital mode, rank-zero prefix.restartaodata AO projection data. Its solve also merges prefix.fock.json (only when not bound to a node; SCFCheckpoint writes the ground-state node's archive/scf.fock.json in the node save instead) and can write orbital energies files; plotting generates density/orbital/coulomb DX and optional cube files. [SCF.cc:save_mos](./SCF.cc) [SCF.cc:solve](./SCF.cc) [SCF.cc:do_plots](./SCF.cc)
- results.json and timings.json are replaced by writing a temporary name and renaming it. Only rank zero of the main World writes them; process groups write their own timing files, which the main session merges and removes. [RunFiles.h:write_results_entry](./RunFiles.h) [RunFiles.h:write_timing_entry](./RunFiles.h)
- Dispersion evaluates on rank zero then broadcasts/cache results, and its Hessian finite-differences analytic gradients; PCM holds PCMSolver cavity/context state and compute_pcm_energy uses the previously computed surface potential/charge labels. [dispersion.cc:compute](./dispersion.cc) [dispersion.cc:hessian](./dispersion.cc) [pcm.cc:compute_pcm_potential](./pcm.cc) [pcm.cc:compute_pcm_energy](./pcm.cc)
- Build installs named Gaussian basis/coredata/structure library data; MRA_CHEMDATA_DIR points in-tree builds at chemistry data. SAP.cc contains static tabulated charges for interpolation, not an executable. [CMakeLists.txt:MRA_CHEMDATA_DIR](./CMakeLists.txt) [CMakeLists.txt:install](./CMakeLists.txt) [SAP.h:SAPCharges](./SAP.h) [SAP.cc:SAPCharges](./SAP.cc)
- Nemo writes orbital checkpoints during iterations when saving is enabled and emits individual nemo, rho and rhonemo fields on the value path. Its CPHF implementation also writes xi guess/solution and density files. These legacy outputs are distinct from node folders. [nemo.cc:solve](./nemo.cc) [nemo.cc:value](./nemo.cc) [nemo.cc:solve_cphf](./nemo.cc)
- TDHF CIS/TDA iterations save roots and transition-density fields. OEP writes a distinct unprefixed restartdata_OEP archive containing orbitals, Fock and final potential, a Slater-potential field, and prefix.oep_calc_info.json. [TDHF.cc:iterate_vectors](./TDHF.cc) [TDHF.cc:save](./TDHF.cc) [oep.cc:save_restartdata](./oep.cc) [oep.cc:solve](./oep.cc)
- Znemo writes its custom reference/orbital/AO archives and individual magnetic orbital/density/current fields; Zcis saves its own unprefixed cis_guess archive. These formats do not use node folders. [znemo.h:save_orbitals](./znemo.h) [znemo.cc:save_orbitals](./znemo.cc) [znemo.cc:analyze](./znemo.cc) [zcis.cc:save_guess](./zcis.cc)
- MAD_EXCH_TASK_PROFILE enables exchange JSONL files per rank and rank-zero CSV timing tables. MAD_EXCH_FORCE_GENERAL routes eligible symmetric exchange through the general path. The benchmark_exchange_operator executable emits benchmark JSON, but its reference-error computation is commented out and should not be cited as a numerical correctness test. [exchangeoperator.h:MAD_EXCH_TASK_PROFILE](./exchangeoperator.h) [exchangeoperator.h:MAD_EXCH_FORCE_GENERAL](./exchangeoperator.h) [benchmark_exchange_operator.cc:main](./benchmark_exchange_operator.cc)
- plotxc is a hard-coded XC sampling utility that writes a data table; it is not a chemistry workflow runner or a registered CTest. SAP.cc is static numeric data, and lda.cc contains built-in Slater/VWN numerical routines used by the fallback XC backend. [plotxc.cc:main](./plotxc.cc) [SAP.cc:SAPCharges](./SAP.cc) [lda.cc:x_rks_s__](./lda.cc) [lda.cc:c_rks_vwn5__](./lda.cc)
- CC2 solve loops write named pair amplitude files and _const files; singles use parallel binary restart archives. CCPair cloud/archive transport preserves pure 6D amplitude and optional constant part rather than all regularization tails. [CC2.cc:CC2::solve_mp2_coupled](./CC2.cc) [CC2.cc:CC2::solve_cc2](./CC2.cc) [CCStructures.h:CC_vecfunction::save_restartdata](./CCStructures.h) [CCStructures.h:CCPair::cloud_store](./CCStructures.h) [CCStructures.h:CCPair::serialize](./CCStructures.h)
- CC2::output_calc_info_schema updates `<prefix>.calc_info` with model, driver=energy, return_energy and a model-key energy. [CC2.cc:CC2::output_calc_info_schema](./CC2.cc)
- CCPairBuilder::load_function and CCPotentials::load_function may reconstruct and resave compressed restart files, then reset thresholds and truncate. [CCStructures.h:CCPairBuilder::load_function](./CCStructures.h) [CCPotentials.h:CCPotentials::load_function](./CCPotentials.h)
- CC2 and potential tasks run collective MRA operations, rebuild intermediates, print RMS/max residuals, energy changes, timings, sizes and memory use; macro tasks use serialized Info and reconstructed functions. [CC2.cc:CC2::solve_mp2_coupled](./CC2.cc) [CC2.cc:CC2::iterate_singles](./CC2.cc) [CCStructures.h:Info](./CCStructures.h) [CCStructures.cc:MacroTaskIteratePair::operator()](./CCStructures.cc)
- BasisFunctions reads basis files (default bas), and internal partial-wave guesses write virtuals* plots. gridbase::visualize writes coordinate files; randomgrid samples std::random_device-backed normal distributions. [PNOGuessFunctions.cpp:BasisFunctions::guess_virtuals_from_file](./PNOGuessFunctions.cpp) [PNOGuessFunctions.cpp:BasisFunctions::guess_virtuals_internal](./PNOGuessFunctions.cpp) [lowrankfunction.h:gridbase::visualize](./lowrankfunction.h) [lowrankfunction.h:randomgrid::gaussian_random_distribution](./lowrankfunction.h)
- PNO exchange can evaluate full Exchange, return zero for neglect, or form an XC potential from the supplied string. [PNOStructures.cpp:ParametrizedExchange::operator ()](./PNOStructures.cpp)
- CCPotentials::make_mo_ket takes shallow Nemo orbital handles and changes thresholds, truncates and reconstructs them; potential calculation updates shared intermediate caches, and operator application can change the selected particle. [CCPotentials.cc:CCPotentials::make_mo_ket](./CCPotentials.cc) [CCPotentials.cc:CCPotentials::get_CC2_singles_potential_gs](./CCPotentials.cc) [CCPotentials.cc:CCPotentials::get_CC2_singles_potential_ex](./CCPotentials.cc) [CCStructures.h:CCConvolutionOperator](./CCStructures.h)
- PNO::solve sets the global 3D threshold. Iteration mutates pair bases, amplitudes, energy arrays and caches; saves per-pair basis functions during iterations and prints pair energies, average/max ranks, residuals, eigenvalues, timings and memory. [PNO.cpp:PNO::solve](./PNO.cpp) [PNO.cpp:PNO::iterate_pairs_internal](./PNO.cpp) [PNO.cpp:PNO::save_pnos](./PNO.cpp) [PNO.cpp:PNO::print_ranks](./PNO.cpp)
- PNO restart loads basis files by PNOPairs::name, catches load failures and leaves those bases empty. CIS(D) reads singles named `<excitation>_x<occupied+freeze>`; F12 auxiliary bases can be read from auxbas_file. [PNO.cpp:PNO::load_pnos](./PNO.cpp) [PNO.cpp:PNO::solve_cispd](./PNO.cpp) [PNOF12Potentials.cpp:F12Potentials::read_cabs_from_file](./PNOF12Potentials.cpp)
- Standalone MP2 construction sets global 3D/6D truncate modes to 3, selects protocols and derives HF econv=0.01*MP2 threshold. Preparation temporarily tightens the global 6D threshold; load balancing redistributes its global process map. [mp2.cc:MP2::MP2](./mp2.cc) [mp2.cc:MP2::guess_mp1_3](./mp2.cc) [mp2.cc:load_balance](./mp2.cc)
- Standalone MP2 reads/writes pair_ij binary archives with iteration/convergence metadata, writes GVpair1 and GVpair preparation files and pair_ij_psi1_converged function files, and optionally loads mRKS_potential_final for OEP. The final amplitude filename is written even when the iteration cap is reached. [mp2.h:ElectronPair::store_pair](./mp2.h) [mp2.h:ElectronPair::load_pair](./mp2.h) [mp2.cc:MP2::guess_mp1_3](./mp2.cc) [mp2.cc:MP2::solve_residual_equations](./mp2.cc) [mp2.cc:MP2::multiply_with_0th_order_Hamiltonian](./mp2.cc)
- MP3 has no new amplitude iteration; drivers evaluate and print contribution terms and timings, and the macro driver schedules collective MacroTaskQ work and reads futures after run_all. [mp3.cc:MP3::mp3_energy_contribution_macrotask_driver](./mp3.cc)

## Source files

| File | Role |
| --- | --- |
| [AC.cc](./AC.cc) | Asymptotic-correction implementation, inactive in present SCF/Nemo constructors |
| [AC.h](./AC.h) | Asymptotic-correction parameters and functors |
| [Applications.hpp](./Applications.hpp) | Application boundaries, live reference context and workflow library adapters |
| [BSHApply.h](./BSHApply.h) | Batched BSH application and fixed-point update support |
| [CC2.cc](./CC2.cc) | Correlated workflow dispatch, pair/singles iterations, restart and energy result assembly |
| [CC2.h](./CC2.h) | Correlated CC2/CIS/MP2 engine declarations |
| [CCLib.hpp](./CCLib.hpp) | Label-only correlated library stub |
| [CCParameters.h](./CCParameters.h) | CC calculation modes, derived thresholds, convergence and restart parameters |
| [CCPotentials.cc](./CCPotentials.cc) | Regularized pair construction, singles diagrams, Fock/commutator potentials and energy kernels |
| [CCPotentials.h](./CCPotentials.h) | Public potential/operator API and orbital/metric accessors |
| [CCStructures.cc](./CCStructures.cc) | Pair completion, cached operators and distributed pair/singles work implementations |
| [CCStructures.h](./CCStructures.h) | CC pair/vector state, index maps, convolution caches, Info serialization and macro task contracts |
| [CMakeLists.txt](./CMakeLists.txt) | Chemistry component source/header membership, optional dependencies, executable and CTest registration |
| [CalculationParameters.h](./CalculationParameters.h) | SCF input defaults, validation and derived occupations/geometry |
| [MetricRecords.h](./MetricRecords.h) | JSON records of iteration measurements, for node descriptions and history rows |
| [ConvergenceSettings.h](./ConvergenceSettings.h) | Shared stopping decider/settings/outcome JSON rules |
| [DensityChanges.h](./DensityChanges.h) | Charge/spin density-field change norms |
| [Drivers.hpp](./Drivers.hpp) | Sequential workflow execution and analytic geometry optimization |
| [ESInterface.h](./ESInterface.h) | External electronic-structure reader interface |
| [ExplicitCriteria.h](./ExplicitCriteria.h) | Explicit target evaluation and progress-window validation |
| [FunctionMetricSummary.h](./FunctionMetricSummary.h) | Metric maxima, unweighted traces and changes |
| [FunctionResidual.h](./FunctionResidual.h) | Raw fixed-point function residual L2 components |
| [FunctionResidualReport.h](./FunctionResidualReport.h) | JSON residual record validity and provenance |
| [GuessFactory.cc](./GuessFactory.cc) | Polynomial/Gaussian/plane-wave guess construction; unsupported trigonometric helper |
| [GuessFactory.h](./GuessFactory.h) | Projected excitation guess functions and polynomial helpers |
| [InputWriter.cpp](./InputWriter.cpp) | Required-group JSON-to-deck serialization |
| [InputWriter.hpp](./InputWriter.hpp) | JSON-to-input serialization declarations |
| [IntegratorXX.h](./IntegratorXX.h) | Optional molecular integration grid wrapper |
| [KineticEnergy.h](./KineticEnergy.h) | Shared gradient-form kinetic diagonal/matrix primitives |
| [LadderPlanner.h](./LadderPlanner.h) | Pure restart/seed planning over supplied rung ladder |
| [MoldftLib.hpp](./MoldftLib.hpp) | Empty unlisted header; no implementation or contract |
| [MolecularOrbitals.cc](./MolecularOrbitals.cc) | Orbital transformations, native/AO archive and cube I/O |
| [MolecularOrbitals.h](./MolecularOrbitals.h) | Orbital functions with energy/occupation/symmetry/localization metadata |
| [NWChem.cc](./NWChem.cc) | External output/movecs parsing and Gaussian conversion |
| [NWChem.h](./NWChem.h) | NWChem output/movecs reader boundary |
| [PNO.cpp](./PNO.cpp) | MP2/CIS(D) PNO solving, adaptive ranks, basis optimization/compression and restart I/O |
| [PNO.h](./PNO.h) | Pair-natural-orbital engine declarations |
| [PNOF12Potentials.cpp](./PNOF12Potentials.cpp) | F12 operator setup, regularized kernels, projected/Hylleraas corrections and auxiliary bases |
| [PNOF12Potentials.h](./PNOF12Potentials.h) | F12 contraction/energy API and inline helper formulas |
| [PNOGuessFunctions.cpp](./PNOGuessFunctions.cpp) | File/internal Gaussian guess construction, basis parsing and solid harmonics |
| [PNOGuessFunctions.h](./PNOGuessFunctions.h) | Guess basis/excitation-operator API, predefined exponents and Gaussian functors |
| [PNOParameters.cpp](./PNOParameters.cpp) | Case-normalized stream conversions for pair, energy and guess enums |
| [PNOParameters.h](./PNOParameters.h) | PNO/F12 settings, type enums, rank/threshold/restart controls and derived defaults |
| [PNOStructures.cpp](./PNOStructures.cpp) | PNO state validation, packed pair iteration, cache extraction and configurable exchange |
| [PNOStructures.h](./PNOStructures.h) | Pair energies, CIS metadata, PNO basis/matrix containers and iteration/exchange declarations |
| [PNOTensors.h](./PNOTensors.h) | Packed triangular pair-overlap storage with transpose symmetry |
| [ParameterManager.hpp](./ParameterManager.hpp) | Typed input groups, prefix precedence and derived settings |
| [PathManager.hpp](./PathManager.hpp) | Filesystem paths and process-cwd scope helper |
| [QCPropertyInterface.h](./QCPropertyInterface.h) | Private placeholder property virtual interface; no public section |
| [ResponseParameters.hpp](./ResponseParameters.hpp) | Compatibility response parameter group |
| [Restart.h](./Restart.h) | Versioned legacy orbital metadata and representation tags |
| [RestartPlan.h](./RestartPlan.h) | Legacy initial/restart/AO/NWChem source policy |
| [Results.h](./Results.h) | Solver/property result records and JSON serialization |
| [ResultsSummary.hpp](./ResultsSummary.hpp) | Text presentation of existing workflow/results JSON |
| [SAP.cc](./SAP.cc) | Static numerical charge-potential tables |
| [SAP.h](./SAP.h) | Static potential-table data declarations |
| [SCF.cc](./SCF.cc) | HF/DFT iterations, orbitals/checkpoints, operators and molecular properties |
| [SCF.h](./SCF.h) | Conventional molecular HF/DFT solver and observer interface |
| [SCFOperators.cc](./SCFOperators.cc) | Operator construction, potential application and matrices |
| [SCFOperators.h](./SCFOperators.h) | Local, nuclear, Coulomb, exchange, kinetic, XC and Fock operator boundaries |
| [SCFProtocol.h](./SCFProtocol.h) | Legacy Nemo threshold protocol and derived density tolerance |
| [SCFTargetAdapter.hpp](./SCFTargetAdapter.hpp) | Geometry energy/gradient target adapters |
| [ScfContinuation.h](./ScfContinuation.h) | Rung identity and saved convergence evidence types |
| [RootCoordination.h](./RootCoordination.h) | Collective rank-zero helpers and checked native-archive close |
| [RunFiles.h](./RunFiles.h) | results.json, timings.json and the older-folder check at the top of a run folder |
| [TDHF.cc](./TDHF.cc) | CIS/TDA iteration and roots; full TDHF deliberately unimplemented |
| [TDHF.h](./TDHF.h) | CIS/TDA response engine declarations and unsupported full-TDHF entry |
| [VizManifest.hpp](./VizManifest.hpp) | Best-effort recursive plot-file manifest |
| [WorkflowBuilders.hpp](./WorkflowBuilders.hpp) | Uninstalled workflow/driver convenience construction |
| [atomutil.cc](./atomutil.cc) | Static atomic data and nuclear smoothing implementation |
| [atomutil.h](./atomutil.h) | Atomic symbols/data and smoothed nuclear potential declarations |
| [basis.h](./basis.h) | Importer-side abstract basis-function boundary |
| [benchmark_exchange_operator.cc](./benchmark_exchange_operator.cc) | Timing executable; reference-error computation disabled |
| [ccpairfunction.cc](./ccpairfunction.cc) | Pure/separated pair algebra, consolidation, operator conversion and partial contractions |
| [ccpairfunction.h](./ccpairfunction.h) | CCFunction metadata and shared pair-component representation/API |
| [corepotential.cc](./corepotential.cc) | Core-data parsing and potential/orbital evaluation |
| [corepotential.h](./corepotential.h) | Core potential/orbital data and radial functions |
| [correlationfactor.cc](./correlationfactor.cc) | Nuclear regularization factor factory |
| [correlationfactor.h](./correlationfactor.h) | Nuclear regularization factor families and potential/derivative formulas |
| [diamagneticpotentialfactor.cc](./diamagneticpotentialfactor.cc) | Magnetic factor/potential projection and regularized derivatives |
| [diamagneticpotentialfactor.h](./diamagneticpotentialfactor.h) | Magnetic Gaussian factor and boxed potential boundary |
| [dispersion.cc](./dispersion.cc) | Rank-zero D3 evaluation/broadcast/cache and coordinate derivatives |
| [dispersion.h](./dispersion.h) | Optional molecular D3 correction boundary |
| [distpm.cc](./distpm.cc) | Distributed systolic localization implementation |
| [distpm.h](./distpm.h) | Distributed orbital localization entry points |
| [electronic_correlation_factor.h](./electronic_correlation_factor.h) | Electronic F12 factors and screened transformed-potential functors |
| [exchangeoperator.cc](./exchangeoperator.cc) | Public exchange construction and algorithm dispatch |
| [exchangeoperator.h](./exchangeoperator.h) | Exchange task scheduling, cache/reduction internals and pure owner helpers |
| [gaussian.cc](./gaussian.cc) | Gaussian construction, normalization, evaluation and algebra |
| [gaussian.h](./gaussian.h) | Importer-side primitive/contracted Cartesian/spherical Gaussians |
| [gth_pseudopotential.cc](./gth_pseudopotential.cc) | Element pseudocharge lookup in supplied XML |
| [gth_pseudopotential.h](./gth_pseudopotential.h) | GTH local/nonlocal projector construction and XML loading |
| [lda.cc](./lda.cc) | Internal built-in Slater exchange and VWN5 correlation routines |
| [localizer.cc](./localizer.cc) | Localization algorithms and deterministic Cholesky pivots |
| [localizer.h](./localizer.h) | Localization method and metric/core-valence interface |
| [lowrankfunction.h](./lowrankfunction.h) | Sampling grids, low-rank functors, product-function fitting and orthonormalization |
| [masks_and_boxes.h](./masks_and_boxes.h) | Smooth scalar mask and radial box templates |
| [molecular_functors.h](./molecular_functors.h) | Molecular density/basis/nuclear derivative projection functors |
| [molecular_optimizer.h](./molecular_optimizer.h) | Legacy molecular optimizer helper |
| [molecularbasis.cc](./molecularbasis.cc) | Native Gaussian basis reading/evaluation |
| [molecularbasis.h](./molecularbasis.h) | Native atomic Gaussian basis and contraction types |
| [molecule.cc](./molecule.cc) | Geometry parsing/orientation, coordinates and nuclear data |
| [molecule.h](./molecule.h) | Molecular geometry, atomic metadata and I/O |
| [molopt.h](./molopt.h) | Active quasi-Newton molecular optimizer |
| [mp2.cc](./mp2.cc) | Standalone regularized MP2 reference setup, independent/local pair equations and checkpointing |
| [mp2.h](./mp2.h) | Standalone MP2 API, Nemo/HF adapter, parameters and serializable ElectronPair |
| [mp3.cc](./mp3.cc) | Third-order energy contractions and serial/distributed contribution drivers |
| [mp3.h](./mp3.h) | Third-order perturbation engine declarations |
| [nemo.cc](./nemo.cc) | Regularized solver, legacy checkpoint/CPHF and nuclear derivatives |
| [nemo.h](./nemo.h) | Regularized-orbital solver and R-squared metric interface |
| [oep.cc](./oep.cc) | Closed-shell OEP solve, own restart format and fields |
| [oep.h](./oep.h) | Local exact-exchange potential model/iteration parameters |
| [pcm.cc](./pcm.cc) | PCMSolver initialization, cavity and surface-charge evaluation |
| [pcm.h](./pcm.h) | PCMSolver cavity/response wrapper and parameter interface |
| [plotxc.cc](./plotxc.cc) | Hard-coded XC sampling executable |
| [pointgroupoperator.h](./pointgroupoperator.h) | Spatial point-group operations on functions and coordinates |
| [pointgroupsymmetry.cc](./pointgroupsymmetry.cc) | Symmetry validation and projected orbital basis construction |
| [pointgroupsymmetry.h](./pointgroupsymmetry.h) | Supported Abelian point groups and irreducible projection |
| [polynomial.cc](./polynomial.cc) | Polynomial coefficient implementation |
| [polynomial.h](./polynomial.h) | Polynomial coefficient operations used by guess functions |
| [population.cc](./population.cc) | AO projections and Mulliken/Lowdin/IAO computations |
| [population.h](./population.h) | Mulliken/Lowdin/IAO population analysis boundary |
| [potentialmanager.cc](./potentialmanager.cc) | Nuclear density, Gaussian/Fermi nuclear, Wigner-Seitz and SAP functor implementations |
| [potentialmanager.h](./potentialmanager.h) | Molecular/core/nuclear/SAP potential functors and manager |
| [projector.h](./projector.h) | Single-particle and strong-orthogonality projection templates |
| [vibanal.cc](./vibanal.cc) | Mass weighting, translation/rotation projection and frequencies |
| [vibanal.h](./vibanal.h) | Mass-weighted vibrational frequency and mode analysis |
| [write_test_input.h](./write_test_input.h) | Installed test-input helper, not a solver engine |
| [xcfunctional.h](./xcfunctional.h) | XC backend interface and fixed-index argument layout |
| [xcfunctional_ldaonly.cc](./xcfunctional_ldaonly.cc) | Fallback HF/LDA backend with unsupported response kernel |
| [xcfunctional_libxc.cc](./xcfunctional_libxc.cc) | Libxc functional parsing and pointwise energy/potential/kernel contractions |
| [zcis.cc](./zcis.cc) | Complex magnetic CIS iteration and own cis_guess archive |
| [zcis.h](./zcis.h) | Complex magnetic CIS declarations |
| [znemo.cc](./znemo.cc) | Complex magnetic HF iteration and legacy field/archive output |
| [znemo.h](./znemo.h) | Complex magnetic HF declarations |

## System Verification

Use an existing configured build and run CTest serially. Optional backends and feature flags control which methods and tests are available.

```sh
cmake --build <build-dir> --target MADchem
```

```sh
ctest --test-dir <build-dir> -R "^madness/test/chem/" -L short --output-on-failure
```

```sh
ctest --test-dir <build-dir> -R "^madness/test/chem/(test_node_archive|test_plan_ladder|test_convergence_settings)/run$" --output-on-failure
```

## Current limitations

Header presence is not a guarantee every route is implemented: several F12 CIS(D) off-diagonal helper branches throw, and CC pair completion currently handles MP2/CC2/LRCC2 only. [PNOF12Potentials.h:F12Potentials::compute_cispd_f12_energies](./PNOF12Potentials.h) [CCStructures.cc:CCPairBuilder::complete_pair_with_low_rank_parts](./CCStructures.cc)

Solver returns are not an unconditional convergence certificate. solve_mp2_coupled stops for absolute energy change plus max pair residual or iteration limit; solve_cc2 uses max residual, singles-change return and a signed delta<econv test. iterate_singles returns whether overall per-orbital change stayed below dconv_3D after its inner loop. [CC2.cc:CC2::solve_mp2_coupled](./CC2.cc) [CC2.cc:CC2::solve_cc2](./CC2.cc) [CC2.cc:CC2::iterate_singles](./CC2.cc)

CC2 solve_lrcc2 currently overwrites loaded singles with a CIS copy before iterations. This is current source behavior; avoid promising complete LRCC2 singles restart reuse. [CC2.cc:CC2::solve_lrcc2](./CC2.cc)

PNO amplitude convergence uses the combined Frobenius residual norm against its threshold. Failure can warn and proceed; an especially large residual at maxiter_t throws. Orbital convergence uses the maximum residual scaled by sqrt(pair-density eigenvalue), while micro convergence combines this with absolute energy change and then reenters to recompute amplitudes/energies. [PNO.cpp:PNO::t_solve](./PNO.cpp) [PNO.cpp:PNO::update_pno](./PNO.cpp) [PNO.cpp:PNO::iterate_pairs_internal](./PNO.cpp)

Adaptive PNO macro convergence uses absolute total-energy change or all per-pair changes below thresholds; converged/insignificant pairs can be frozen and maxrank or the iteration limit can also stop work. It does not promise a uniformly converged result on every return. [PNO.cpp:PNO::adaptive_solver](./PNO.cpp) [PNO.cpp:PNO::freeze_insignificant_pairs](./PNO.cpp)

PNO save/load handles bases, not a full solver checkpoint. no_compute skips the inner solve and sizes zero amplitude matrices; restarting bases rebuilds solver intermediates. [PNO.cpp:PNO::save_pnos](./PNO.cpp) [PNO.cpp:PNO::load_pnos](./PNO.cpp) [PNO.cpp:PNO::iterate_pairs_internal](./PNO.cpp)

Standalone MP2 independent-pair convergence requires absolute energy change and residual norm thresholds; the localized coupled solve uses energy change plus the sum of pair residual norms. Both can return at the iteration limit. [mp2.cc:MP2::solve_residual_equations](./mp2.cc) [mp2.cc:MP2::solve_coupled_equations](./mp2.cc)

Pair energies in standalone MP2 and PNO include explicit direct/exchange and diagonal/off-diagonal factors; report whether a value is correlation energy, F12 correction or total electronic energy instead of treating all printed energy fields interchangeably. [mp2.cc:MP2::compute_energy](./mp2.cc) [PNO.cpp:PNO::compute_projected_mp2_energies](./PNO.cpp) [PNOF12Potentials.cpp:F12Potentials::print_pair_energies](./PNOF12Potentials.cpp) [CC2.cc:CC2::solve](./CC2.cc)
