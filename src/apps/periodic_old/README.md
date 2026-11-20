# periodic_old component

## Overview

Dormant historical periodic molecular/solid HF and LDA solver experiments, local basis/atom support, and Ewald diagnostics.

## Interface

- This entire subtree is absent from ../CMakeLists.txt traversal. Its local manifest describes HartreeFock from eigsolver.cc/dft.cc, several test executables, and esolver from solver_driver.cc; none are configured by the active application build. hartreefock.cc and ewald.cc are not even included in this local manifest. [../CMakeLists.txt:add_subdirectory](../CMakeLists.txt) [CMakeLists.txt:HARTREEFOCK_SOURCES](./CMakeLists.txt)
- solver_driver.cc main initializes runtime/MRA and constructs Solver<double,3>(world,"input"). ElectronicStructureParams scans a dft block for cubic box, functional 1=LDA/2=HF, waveorder/thresh,maxits,rcriterion,nempty,basis XML,periodic/kmesh/offset,spin,restart,solver/maxsub and damping controls. [solver_driver.cc:main](./solver_driver.cc) [electronicstructureparams.h:ElectronicStructureParams::read_file](./electronicstructureparams.h)
- MolecularEntity scans geometry...end with units first (atomic default or angstrom), optional fractional coordinates in [0,1], atom symbol and coordinates, with ghost charge field. AtomicBasisSet loads TinyXML basis/atomicguess records; Cartesian Gaussian evaluation supports S/P/D/F (G is recognized by XML but evaluator rejects unsupported angular momentum). [mentity.cc:MolecularEntity::read_file](./mentity.cc) [molecularbasis.h:AtomicBasisSet::read_file](./molecularbasis.h) [molecularbasis.h:ContractedGaussianShell::eval](./molecularbasis.h)
- Solver represents k-point ranges as [begin,end), separate alpha/beta orbitals and occupations; builds initial Gaussian basis eigensystem, density-weighted Coulomb/LDA or HF exchange, BSH updates, optional KAIN then damping and per-k-point Gram-Schmidt. [esolver.h:KPoint](./esolver.h) [solver.h:Solver::initial_guess](./solver.h) [solver.h:Solver::apply_potential](./solver.h) [solver.h:Solver::update_orbitals](./solver.h)
- Older EigSolver/EigSolverOp/IEigSolverObserver form a separate operator/observer abstraction with density- or orbital-dependent potentials and serial/vector solve entrypoints; DFT and HartreeFock wrappers report energies every third iteration. EigSolver owns/deletes operator pointers, not observer pointers. [eigsolver.h:EigSolverOp](./eigsolver.h) [eigsolver.cc:EigSolver::~EigSolver](./eigsolver.cc) [dft.cc:DFT::iterateOutput](./dft.cc) [hartreefock.cc:HartreeFock::iterateOutput](./hartreefock.cc)
- Ewald experiment implements real/reciprocal cubic lattice sums, Gaussian Coulomb fits and nuclear-charge convergence comparisons; its main selects test_nuclear_potential3 with hardcoded geometry,L=6.5,thresholds and k=8. It has no configured target. [ewald.cc:EwaldNuclearPotentialFunctor](./ewald.cc) [ewald.cc:main](./ewald.cc) [CMakeLists.txt:HARTREEFOCK_SOURCES](./CMakeLists.txt)

## Invariants

- Treat these as historical experimental interfaces, not supported production paths: spin-polarized LDA throws, periodic kinetic energy reporting returns zero, and old DFT references libxc_ldaop whose definition is commented out. HartreeFock implementation signatures differ from its header; Ewald references get_entity absent from local MolecularEntity. No compile/run correctness is asserted. [solver.h:Solver::apply_potential](./solver.h) [solver.h:Solver::calculate_kinetic_energy](./solver.h) [dft.cc:XCFunctionalLDA::op_r](./dft.cc) [libxc.h:libxc_ldaop](./libxc.h) [hartreefock.cc:HartreeFockExchangeOp::op_o](./hartreefock.cc) [hartreefock.h:HartreeFockExchangeOp](./hartreefock.h) [ewald.cc:make_nuclear_charge_density_individual](./ewald.cc) [mentity.h:MolecularEntity](./mentity.h)
- Solver::init overwrites requested nelec with nuclear charge and rounds odd electron counts up; occupied states are filled per k-point without smearing by set_occs2. Explicit kpoints with nonperiodic input is rejected. [solver.h:Solver::init](./solver.h) [solver.h:Solver::set_occs2](./solver.h) [electronicstructureparams.h:ElectronicStructureParams::read_file](./electronicstructureparams.h)
- Solver stopping is mean orbital-update residual <= rcriterion or maxits. Reprojection adds two to waveorder and divides thresh by100 when residual<20*thresh or after20 iterations; the stored maxthresh/maxwaveorder are not enforced caps there. Generic EigSolver solves a fixed maxits and can _exit(0) if corrected eigenvalue remains nonnegative, so process success is not convergence evidence. [solver.h:Solver::solve](./solver.h) [solver.h:Solver::compute_residual](./solver.h) [solver.h:Solver::reproject](./solver.h) [eigsolver.cc:EigSolver::solve](./eigsolver.cc)
- Restart is fixed archive orbitals; metadata contains spin flag,kpoint ranges,orbital counts/functions, with eigenvalues reinitialized and k mismatch projected. It is a historical format without modern identity/schema/checksum guarantees; do not present compatibility with madqc checkpoints. [solver.h:Solver::save_orbitals](./solver.h) [solver.h:Solver::load_orbitals](./solver.h)

## Side Effects

- Changes global 3D defaults for cell,k,thresh,periodic boundary and process mapping; initializes/replaces Coulomb operators after reprojection. [solver_driver.cc:main](./solver_driver.cc) [solver.h:Solver::init](./solver.h) [solver.h:Solver::reproject](./solver.h)
- Historical run overwrites output.txt,matrix.txt,eigs.txt, root-only kpoints.txt and distributed orbitals archive each iteration/end; emits extensive matrices/energies/residuals/stdout diagnostics. Initial AO DX/line plots are unconditional, while final/iteration orbital plots use plotorbs. Initial density plot uses an uninitialized filename buffer, another reason not to offer an unvalidated run command. [solver.h:Solver::init](./solver.h) [solver.h:Solver::save_orbitals](./solver.h) [solver.h:Solver::initial_guess](./solver.h) [solver.h:Solver::solve](./solver.h)
- Embedded LDA utilities apply density floors and mutate tensor handles for functional/potential evaluation; local libxc.h wraps included Slater/VWN routines rather than an active external LibXC integration. Ewald helpers can enumerate a fixed 1001^3 lattice cube and print lengthy diagnostics. [lda.h:xc_rks_generic_lda](./lda.h) [libxc.h:ldaop](./libxc.h) [ewald.cc:generate_R_vectors](./ewald.cc) [ewald.cc:generate_G_vectors](./ewald.cc)

## Component boundaries

- [matlab](./matlab/README.md)

## Source files

| File | Role |
| --- | --- |
| [complexfun.h](./complexfun.h) | Legacy scalar/tensor/function complex conversion and component unary operators. |
| [dft.h](./dft.h) | Generic DFT operator and observer interfaces. |
| [dft.cc](./dft.cc) | DFT Coulomb/nuclear/XC operator composition and energy reporting. |
| [eigsolver.h](./eigsolver.h) | Generic eigenvalue solver, owned operators and nonowned observers. |
| [eigsolver.cc](./eigsolver.cc) | Fixed-iteration serial/vector BSH eigensolver and density construction. |
| [electronicstructureapp.h](./electronicstructureapp.h) | Legacy distribution map and atomic/nuclear projection functors. |
| [electronicstructureparams.h](./electronicstructureparams.h) | dft-block historical input schema/defaults and serialization. |
| [esolver.h](./esolver.h) | Complex 3D types,k-point ranges and Bloch kinetic matrix helpers. |
| [ewald.cc](./ewald.cc) | Unbuilt exploratory cubic Ewald and fitted nuclear-potential comparisons. |
| [hartreefock.h](./hartreefock.h) | Historical HF operator/wrapper declarations. |
| [hartreefock.cc](./hartreefock.cc) | Unbuilt legacy HF/exchange implementation with interface drift. |
| [lda.h](./lda.h) | Embedded generated Slater exchange/VWN5 correlation kernels and tensor adapters. |
| [libxc.h](./libxc.h) | Local LDA unary wrappers; external LibXC experiments commented out. |
| [mentity.h](./mentity.h) | Atom table/molecular coordinates, charge/potential and serialization contracts. |
| [mentity.cc](./mentity.cc) | Historical atom data,geometry parser,smoothed potential/density. |
| [molecularbasis.h](./molecularbasis.h) | Cartesian Gaussian shells,basis XML,atomic-density guesses. |
| [outputwriter.h](./outputwriter.h) | Singleton owning optional debug/log/eigenvalue ofstreams. |
| [poperator.h](./poperator.h) | Gaussian-expanded periodic Coulomb/BSH constructors. |
| [solver.h](./solver.h) | Complex k-point SCF/KAIN research implementation and archive/output lifecycle. |
| [solver_driver.cc](./solver_driver.cc) | Unbuilt esolver launcher using fixed input filename. |
| [util.h](./util.h) | Commented historical diagnostics, no active utility implementations. |
| [util.cc](./util.cc) | License-only source, no active implementation. |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

No active build target or CTest entry verifies this subtree. Inspect ../CMakeLists.txt and configured target/test inventories to confirm exclusion; local target declarations do not constitute supported build commands. Modern periodic coverage is elsewhere (e.g. moldft testperiodicdft). [../CMakeLists.txt:add_subdirectory](../CMakeLists.txt) [CMakeLists.txt:HARTREEFOCK_SOURCES](./CMakeLists.txt) [../moldft/CMakeLists.txt:testperiodicdft](../moldft/CMakeLists.txt)
