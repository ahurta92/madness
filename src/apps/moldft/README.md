# SCF data and retained utilities

Molecular HF/DFT is run through `madqc --wf=scf`. Its numerical implementation
lives in [SCF.cc](../../madness/chem/SCF.cc) and its shared iteration/restart
session in [SCFRun.cc](../../madness/chem/SCFRun.cc). The duplicate `moldft` main
and target are retired.

This directory remains in CMake for `testpg` and `testperiodicdft`. Periodic DFT
is not yet a madqc engine and retains its separate executable. Historical
molecular decks, basis-conversion scripts, point-group/subspace utilities and
[fci data](fci/README.md) remain here. Old Python utilities may require Python 2;
they are not readers for the shared result format.

Use [curated QC decks](../../examples/qc/README.md) for molecular runs:

```sh
cmake --build <build-dir> --target madqc
<build-dir>/src/apps/madqc/madqc --wf=scf --geometry=he --dft="xc=hf"
```
