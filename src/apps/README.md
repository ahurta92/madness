# MADNESS applications

`madqc` is the public executable for supported quantum-chemistry engines. The
SCF, Nemo, molecular response, CIS, OEP and CC2 algorithms remain in their
libraries; their duplicate application mains and executable targets are retired.
`madqc --wf=<name>` selects the workflow. Named-stage input is removed.

| Directory | Build purpose |
|---|---|
| [madqc](madqc/README.md) | Unified driver and workflow integration tests |
| [molresponse](molresponse/README.md) | Response implementation, private numerical tests, import and analysis tools |
| [moldft](moldft/README.md) | Private point-group tests, separate periodic DFT and historical data/utilities |
| [mp2](mp2/README.md) | Separate MP2 implementation, distinct from madqc's CC2 MP2 route |
| [pno](pno/README.md) | Pair-natural-orbital MP2/F12 engine, not yet in madqc |
| [znemo](znemo/README.md), [zcis](zcis/README.md) | Complex-orbital engines, not yet in madqc |
| [dirac](dirac/README.md) | Relativistic DFdriver, not yet in madqc |
| [tdse](tdse/README.md) | Separate time-dependent 3D/4D applications |
| [plot](plot/README.md) | Plotting and conversion tools |

The `nemo`, `cis`, `oep` and `cc2` directories retain historical decks and notes;
they no longer declare executable targets. `periodic_old` and `interior_bc`
remain outside the active CMake tree. Numerical chemistry is in
[src/madness/chem](../madness/chem/README.md).

Use [qctest decks](../examples/qc/README.md) for supported workflow examples.
Executable consolidation and shared convergence/checkpoint adoption are distinct:
SCF uses the shared session; other engines continue their current numerical
implementations until their adapters are migrated. See the
[test migration index](../../docs/design/madqc-driver-test-migration.md).
