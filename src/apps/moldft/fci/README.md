# fci component

## Overview

Legacy standalone Fortran full-configuration-interaction source with C interoperability helpers; outside the active CMake graph.

## Interface

- GNUmakefile declares fci executable, libfci.a and timer; uses gfortran, gcc and bundled BLAS objects. [GNUmakefile:fci](./GNUmakefile)
- wopen_,wclose_,getwa_,putwa_,walltm_ implement synchronous word-addressed Fortran I/O; ieor_ supplies XOR and linux_cputime_ CPU seconds. [crayio.c:wopen_](./crayio.c) [ieor.c:ieor_](./ieor.c) [linux_cputime.c:linux_cputime_](./linux_cputime.c)
- The standalone driver reads orbital/electron counts, spin multiplicity, state symmetry, requested roots, tolerances, orbital symmetries and integrals from unit 5. It calls fci_setup, fci_load and fci_solve, then prints density diagnostics. [fci.F:program fci](./fci.F) [fci_setup.F:fci_setup](./fci_setup.F) [fci_load.F:fci_load](./fci_load.F)
- fci_solve uses block Davidson; fci_hv applies the Hamiltonian and fci_hami constructs its packed triangular matrix. Coupling and spin-adaptation routines share configuration through COMMON blocks. [fci_solve.F:fci_solve](./fci_solve.F) [fci_hv.F:fci_hv](./fci_hv.F) [fci_hami.F:fci_hami](./fci_hami.F) [fciinfo.fh:fciinfo](./fciinfo.fh)
- fci_onepdm computes a weighted average one-particle density matrix, normalizing the supplied weights in place; fci_1pdm/fci_2pdm provide driver diagnostics. [fci_onepdm.F:fci_onepdm](./fci_onepdm.F) [fci_1pdm.F:fci_1pdm](./fci_1pdm.F) [fci_2pdm.F:fci_2pdm](./fci_2pdm.F)

## Invariants

- Units are 0..98; word addresses start at 1, count>=0; words are 8 bytes; reads cannot exceed tracked file length. [crayio.c:CheckUnit](./crayio.c) [crayio.c:CheckAddr](./crayio.c) [crayio.c:getwa_](./crayio.c)
- fci_setup initializes the shared orbital/spin configuration before Hamiltonian, coupling or density routines run. Capacity is fixed by fciparam.fh; the driver separately limits active orbitals to 14 and clips requested roots to three. [fci_setup.F:fci_setup](./fci_setup.F) [fciparam.fh:nactmx](./fciparam.fh) [fci.F:maxroot](./fci.F)
- Integral symmetry is checked and forbidden small entries are set to zero; inconsistent entries above the hard-coded tolerances call fci_err. [fci_intchk.F:fci_intchk](./fci_intchk.F)
- fci_solve caps Davidson work at 200 iterations. fci_davids checks residual norms for each requested root while reducing its internal working tolerance toward the requested tolerance. [fci_solve.F:fci_solve](./fci_solve.F) [fci_davids.F:fci_davids](./fci_davids.F)

## Side Effects

- Opens/creates files mode 0660, seeks/reads/writes/closes; shared process-global file table and optional I/O statistics. [crayio.c:wopen_](./crayio.c) [crayio.c:putwa_](./crayio.c) [crayio.c:wclose_](./crayio.c)
- The driver writes evecs.dat and diags.dat; Davidson uses workfile for disk-backed subspace storage. fci_solve constructs fresh spin-projected guesses rather than resuming an iteration from those files. [fci.F:evecfilename](./fci.F) [fci_solve.F:fci_solve](./fci_solve.F) [fci_davids.F:fci_davids](./fci_davids.F)
- FCI routines print energies, residuals, density matrices and CPU timings. fci_err reports the error and stops with code 99. [fci_tinit.F:fci_tprint](./fci_tinit.F) [fci_err.F:fci_err](./fci_err.F)

## Source files

| File | Role |
| --- | --- |
| [crayio.c](./crayio.c) | Fortran ABI word I/O and wall time. |
| [ieor.c](./ieor.c) | Fortran ABI integer XOR. |
| [linux_cputime.c](./linux_cputime.c) | Fortran ABI CPU time. |
| [bitops.fh](./bitops.fh) | Preprocessor mapping of integer bit operations. |
| [daxpy.F](./daxpy.F) | Bundled BLAS implementation. |
| [dcopy.F](./dcopy.F) | Bundled BLAS implementation. |
| [ddot.F](./ddot.F) | Bundled BLAS implementation. |
| [dgemm.F](./dgemm.F) | Bundled BLAS implementation. |
| [dscal.F](./dscal.F) | Bundled BLAS implementation. |
| [fci.F](./fci.F) | Standalone input, solve and density-diagnostic driver. |
| [fci_1pdm.F](./fci_1pdm.F) | One-particle density diagnostics. |
| [fci_2pdm.F](./fci_2pdm.F) | Two-particle density diagnostics. |
| [fci_3pdm.F](./fci_3pdm.F) | Three-particle density diagnostic helper, not called by the driver. |
| [fci_axb.F](./fci_axb.F) | FCI helper: fci_axb. |
| [fci_cntone.F](./fci_cntone.F) | FCI helper: fci_cntone. |
| [fci_couple1.F](./fci_couple1.F) | FCI helper: fci_couple1. |
| [fci_couple2.F](./fci_couple2.F) | FCI helper: fci_couple2. |
| [fci_couple3.F](./fci_couple3.F) | FCI helper: fci_couple3. |
| [fci_davids.F](./fci_davids.F) | Disk-backed block-Davidson iteration and residual reporting. |
| [fci_detwalk.F](./fci_detwalk.F) | FCI helper: fci_det_to_walk. |
| [fci_dfill.F](./fci_dfill.F) | FCI helper: dfill. |
| [fci_diags.F](./fci_diags.F) | FCI helper: fci_diags. |
| [fci_err.F](./fci_err.F) | Diagnostic output and STOP 99. |
| [fci_hami.F](./fci_hami.F) | Explicit triangular Hamiltonian construction. |
| [fci_hmod.F](./fci_hmod.F) | Effective one-electron integral adjustment. |
| [fci_hsort.F](./fci_hsort.F) | FCI helper: fci_hsort, fci_heapify. |
| [fci_hv.F](./fci_hv.F) | Hamiltonian-vector product wrapper and scratch sizing. |
| [fci_ifill.F](./fci_ifill.F) | FCI helper: ifill. |
| [fci_intchk.F](./fci_intchk.F) | Integral symmetry checks and zeroing. |
| [fci_load.F](./fci_load.F) | Integral input and packing. |
| [fci_make_eaj.F](./fci_make_eaj.F) | Spin permutations, phases and Eaj application. |
| [fci_mfndoc.F](./fci_mfndoc.F) | FCI helper: fci_mfind_occ. |
| [fci_mfndwt.F](./fci_mfndwt.F) | FCI helper: fci_mfind_weight. |
| [fci_mgraph.F](./fci_mgraph.F) | Primitive spin graph construction. |
| [fci_num.F](./fci_num.F) | FCI helper: fci_num_dets, fci_num_csf, fci_num_owalks, fci_ndets, fci_ncsf, fci_max_num_dets. |
| [fci_ofndoc.F](./fci_ofndoc.F) | FCI helper: fci_ofind_occ. |
| [fci_ofndwt.F](./fci_ofndwt.F) | FCI helper: fci_ofind_weight. |
| [fci_ograph.F](./fci_ograph.F) | Orbital occupation graph construction. |
| [fci_onelcc.F](./fci_onelcc.F) | FCI helper: fci_onelcc. |
| [fci_onepdm.F](./fci_onepdm.F) | Weighted density matrix; normalizes weights in place. |
| [fci_output.F](./fci_output.F) | FCI helper: fci_output. |
| [fci_owalk_info.F](./fci_owalk_info.F) | FCI helper: fci_owalk_info. |
| [fci_prntri.F](./fci_prntri.F) | FCI helper: fci_prntri. |
| [fci_prspnad.F](./fci_prspnad.F) | FCI helper: fci_prspnad. |
| [fci_prteop.F](./fci_prteop.F) | FCI helper: fci_prteop. |
| [fci_prvec.F](./fci_prvec.F) | FCI helper: fci_print_vector. |
| [fci_setup.F](./fci_setup.F) | Configuration validation and graph initialization. |
| [fci_sigma.F](./fci_sigma.F) | Blocked Hamiltonian-vector contraction. |
| [fci_solve.F](./fci_solve.F) | Workspace allocation, guesses and block-Davidson solve. |
| [fci_spnad.F](./fci_spnad.F) | Determinant/CSF transformations. |
| [fci_spnfrc.F](./fci_spnfrc.F) | Spin projection and normalization. |
| [fci_srt1cc.F](./fci_srt1cc.F) | FCI helper: fci_srt1cc. |
| [fci_tinit.F](./fci_tinit.F) | Routine timing and CPU clock dispatch. |
| [fciinfo.fh](./fciinfo.fh) | Shared orbital, spin-graph and determinant COMMON state. |
| [fciparam.fh](./fciparam.fh) | Fixed FCI storage capacities. |
| [fcitime.fh](./fcitime.fh) | Shared timing counters and routine names. |
| [idamax.F](./idamax.F) | Bundled BLAS implementation. |
| [renorm.F](./renorm.F) | In-place Euclidean vector normalization. |
| [rsg.f](./rsg.f) | Generalized symmetric eigenproblem and bundled EISPACK helpers. |
| [screen.F](./screen.F) | Threshold-based vector zeroing. |
| [timer.f](./timer.f) | Standalone DAXPY timing driver. |
| [yacobi.F](./yacobi.F) | Packed symmetric Jacobi diagonalization and eigenpair ordering. |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

Legacy build recipe (historical flags/toolchain; unverified):

```sh
make -C src/apps/moldft/fci -f GNUmakefile fci
```

No CTest registration.

## Current limitations

This historical standalone code has no registered CTest coverage. The logical function fci_davids does not assign its return value, although fci_solve tests that value; its printed convergence diagnostics do not establish a reliable programmatic success result. [fci_davids.F:fci_davids](./fci_davids.F) [fci_solve.F:fci_solve](./fci_solve.F)
