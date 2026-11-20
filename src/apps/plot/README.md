# Plot Functions from Chemistry Libraries

```bash
plot2cube file=mra_orbital_0
plot2plane file=mra_orbital_1
```

Requires these files:

- mra_orbital_0.00000
- input

<!-- EGM architecture contract: begin -->

## Architecture overview

Archive-to-cube and archive-to-plane conversion utilities.

## Interface

- plot2cube file=name1[,name2] writes each real 3D Function archive as cube. [plot2cube.cpp:main](./plot2cube.cpp)
- plot2plane file=name1[,name2] tries combined 3D and 6D plane output. [plot2plane.cpp:do_plot](./plot2plane.cpp)
- Both active installed targets link MADchem. [CMakeLists.txt:add_mad_executable](./CMakeLists.txt)

## Invariants

- Input archive names omit .00000; each file must hold one real Function; geometry and calculation parameters supply cell and cube atom header. [plot2cube.cpp:main](./plot2cube.cpp)
- Plane loader verifies tree state and converts loaded tensor type. [plot2plane.cpp:do_plot](./plot2plane.cpp)

## Side Effects

- Reads archives and writes `<filename>`.cube or `<firstfilename>`.plane; updates 3D/6D FunctionDefaults cell. [plot2cube.cpp:main](./plot2cube.cpp) [plot2plane.cpp:main](./plot2plane.cpp)
- Missing files print help and return zero; plane catches and reports 3D/6D failures. [plot2plane.cpp:main](./plot2plane.cpp)

## Source files

| File | Role |
| --- | --- |
| [plot2cube.cpp](./plot2cube.cpp) | Cube conversion main. |
| [plot2plane.cpp](./plot2plane.cpp) | Template do_plot and plane conversion main. |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

```sh
cmake --build <build-dir> --target plot2cube plot2plane
```

No local CTest registration.

<!-- EGM architecture contract: end -->
