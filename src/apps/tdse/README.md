# tdse component

## Overview

Time-dependent one-electron Schrödinger applications in 3D and coupled electron/nuclear 4D.

## Interface

- Active tdse and tdse4 executables link MADmra; installation is commented and no tests registered here. tdse.confused.cc is a dormant earlier 3D variant, absent from targets. [CMakeLists.txt:add_mad_executable](./CMakeLists.txt)
- tdse reads fixed cwd files input and restart; tdse4 reads input4 and restart4. InputParameters parses whitespace tags (unknown tags throw), rank zero reads then broadcasts; controls include L,F,omega,ncycle,k,thresh,safety,cut,prefix,dump/print/load-balance cadence,nio,tScale,target_time. 3D additionally accepts nuclei, near/far plot extents,nplot and iState 1s/2s. [tdse.cc:InputParameters::read](./tdse.cc) [tdse4.cc:InputParameters::read](./tdse4.cc) [tdse.cc:doit](./tdse.cc)
- 3D constructs smoothed nuclear potential, normalized initial 1s or orthogonalized 2s, then propagates in a z-polarized sin-squared pulse with Chin-Chen splitting; dormant Trotter branch is hardcoded false. 4D models H2+ with three electronic and one mass-scaled bond coordinate, fixed R0/proton mass, x-polarized pulse and Trotter splitting. [tdse.cc:converge](./tdse.cc) [tdse.cc:propagate](./tdse.cc) [tdse4.cc:Vn](./tdse4.cc) [tdse4.cc:Ve](./tdse4.cc) [tdse4.cc:propagate](./tdse4.cc)
- Ground-state fixed-point iterations stop at residual < max(1e-5,thresh), capped at 30 in 3D and 205 with 0.3 mixing in 4D. Time step derives from cut bandlimit and tScale; stop step derives from target_time, with 10-step (3D) or 20-step (4D) field-free delay. [tdse.cc:converge](./tdse.cc) [tdse4.cc:converge](./tdse4.cc) [tdse.cc:propagate](./tdse.cc) [tdse4.cc:propagate](./tdse4.cc)

## Invariants

- Runtime startup precedes Function construction; dimension defaults set cell/order/threshold/safety and process map. Distributed operations and archive access are collective; APPLY shallow-copies its handle and reconstructs/broadens the shared source. [tdse.cc:main](./tdse.cc) [tdse.cc:doit](./tdse.cc) [tdse.cc:APPLY](./tdse.cc) [tdse4.cc:APPLY](./tdse4.cc)
- Provide complete valid input and restart integer explicitly: parser has no completeness/range validation, restart extraction does not default step0, and cadence values are divisors. 3D nucleus arrays have fixed MAXNATOM=99 without a parser bounds check. Missing nonzero-step archive aborts; missing step-zero archive triggers initial-state solve. [tdse.cc:InputParameters](./tdse.cc) [tdse.cc:doit](./tdse.cc) [tdse4.cc:doit](./tdse4.cc)
- Ground-state iteration cap alone is not a convergence guarantee; code proceeds to save/propagate after cap. Reuse requires compatible MRA defaults and original step-zero archive for overlap diagnostics. [tdse.cc:converge](./tdse.cc) [tdse.cc:propagate](./tdse.cc) [tdse4.cc:converge](./tdse4.cc) [tdse4.cc:propagate](./tdse4.cc)

## Side Effects

- Writes distributed wavefunction archives `<prefix>-<five-digit-step>`, overwrites restart/restart4 at dump boundaries, and prints time/field/energy/norm/initial overlap/dipoles; 3D also reports acceleration,4D mean bond length and energy components. [tdse.cc:wave_function_store](./tdse.cc) [tdse.cc:print_stats](./tdse.cc) [tdse4.cc:print_stats](./tdse4.cc)
- 3D writes 10001-point z-axis `<prefix>`-`<step>`.lineplot on print boundaries; doplot constructs grids but its plotdx call is commented, so nplot does not produce advertised S.dx/L.dx volumes. 4D overwrites plot.dat with potential/guess samples. Load-balancing changes FunctionDefaults maps and redistributes trees. [tdse.cc:line_plot](./tdse.cc) [tdse.cc:doplot](./tdse.cc) [tdse.cc:loadbal](./tdse.cc) [tdse4.cc:doit](./tdse4.cc)

## Source files

| File | Role |
| --- | --- |
| [tdse.cc](./tdse.cc) | Active 3D fixed-nuclei one-electron solver and Chin-Chen real-time propagation. |
| [tdse4.cc](./tdse4.cc) | Active 4D H2+ electron/bond-coordinate solver and Trotter propagation. |
| [tdse.confused.cc](./tdse.confused.cc) | Dormant 3D separated-convolution alternative with older runtime APIs. |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

```sh
cmake --build <build-dir> --target tdse tdse4
```

No registered CTest case. A numerical run needs prepared input/restart files. [CMakeLists.txt:add_mad_executable](./CMakeLists.txt) [tdse.cc:doit](./tdse.cc) [tdse4.cc:doit](./tdse4.cc)
