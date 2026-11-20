# tools component

## Overview

Offline response/ground-state seed conversion, Dalton contraction diagnostics, bundle overlap and tree/coefficient exports.

## Interface

- dump_mra_trees loads ground MOs and optional FD/ES states and exports reconstructed/compressed JSON trees, geometry VTK, optional dense cube, HTG, adaptive AMR and native coefficients [dump_mra_trees.cpp:main](./dump_mra_trees.cpp) [dump_mra_trees.cpp:dump_function_trees](./dump_mra_trees.cpp).
- dalton_gto/read_molden parse spherical s through g Molden basis/MO data; dalton_rspvec parses Fortran sequential records and splits occupied-virtual TDA/RPA blocks; dalton_mra projects AO combinations and rotates the occupied gauge [dalton_gto.hpp:read_molden](./dalton_gto.hpp) [dalton_rspvec.hpp:read_rspvec](./dalton_rspvec.hpp) [dalton_rspvec.hpp:split_ov](./dalton_rspvec.hpp) [dalton_mra.hpp:occupied_gauge_rotation](./dalton_mra.hpp).
- seed_from_dalton converts selected raw RSPVEC entries into one TDA/Full seed bundle; dalton_es_seed_bundle calls the production ES import routine; seed_es_from_hdf5 creates a one-root TDA seed from structured per-orbital functions; seed_moldft_from_dalton writes ground-state restartdata [seed_from_dalton.cpp:main](./seed_from_dalton.cpp) [dalton_es_seed_bundle.cpp:main](./dalton_es_seed_bundle.cpp) [seed_es_from_hdf5.cpp:main](./seed_es_from_hdf5.cpp) [seed_moldft_from_dalton.cpp:main](./seed_moldft_from_dalton.cpp).
- es_overlap prints normalized response-metric overlaps and best matches; tpa_from_dalton projects external vectors and evaluates legacy/residue diagnostic contractions without iterative solves [es_overlap.cpp:run_overlap](./es_overlap.cpp) [tpa_from_dalton.cpp:main](./tpa_from_dalton.cpp).

## Invariants

- RSPVEC reader assumes little-endian 4-byte Fortran record markers; split_ov accepts exactly nocc*nvir or twice that size. Molden parser assumes pure spherical s,p,d,f,g and ignores spin headers [dalton_rspvec.hpp:read_fortran_record](./dalton_rspvec.hpp) [dalton_rspvec.hpp:split_ov](./dalton_rspvec.hpp) [dalton_gto.hpp:read_molden](./dalton_gto.hpp).
- Occupied gauge rotation uses M(M^T M)^(-1/2) and rejects minimum overlap eigenvalue below 0.5 [dalton_mra.hpp:occupied_gauge_rotation](./dalton_mra.hpp).
- seed_from_dalton, seed_moldft_from_dalton, seed_es_from_hdf5 and tpa_from_dalton enforce one rank. Seed states have converged=false and finite residual sentinels; scale defaults sqrt(2) for seed tools [seed_from_dalton.cpp:main](./seed_from_dalton.cpp) [seed_moldft_from_dalton.cpp:main](./seed_moldft_from_dalton.cpp) [seed_es_from_hdf5.cpp:main](./seed_es_from_hdf5.cpp) [tpa_from_dalton.cpp:main](./tpa_from_dalton.cpp).
- seed_from_dalton with archive rotates to its occupied gauge but does not Q-project in that tool; without archive it uses [-200,200]^3 and warns that seed is unrotated. Production dalton_es_seed_bundle uses production Q-projected seeds, but directly calls locate/seed and does not run the top-level import geometry/method gate [seed_from_dalton.cpp:main](./seed_from_dalton.cpp) [dalton_es_seed_bundle.cpp:main](./dalton_es_seed_bundle.cpp).
- Exporter collective tree operations run on all ranks with rank-zero file emission. FD export skips nonclosed shell and unknown types, and uses the restart loader; no convergence flag filter is imposed at enumeration. ES loader validates shell/type [dump_mra_trees.cpp:write_tree_json](./dump_mra_trees.cpp) [dump_mra_trees.cpp:dump_fd_states](./dump_mra_trees.cpp) [dump_mra_trees.cpp:dump_es_bundle](./dump_mra_trees.cpp).
- HTG/AMR need VTK; coefficients need HDF5 and one rank. AMR samples cell centers, rebases to populated minimum level and warns about empty interior levels; XML writers use ASCII. The documented --maxlevel is not parsed/used by main [dump_mra_trees.cpp:write_amr](./dump_mra_trees.cpp) [dump_mra_trees.cpp:write_htg](./dump_mra_trees.cpp) [dump_mra_trees.cpp:main](./dump_mra_trees.cpp).

## Side Effects

- Export creates out directory (default mra_trees), geometry.vtk and per-field files, changes active FunctionDefaults across saved protocols and reconstructs/compresses shared in-memory handles [dump_mra_trees.cpp:main](./dump_mra_trees.cpp) [dump_mra_trees.cpp:dump_function_trees](./dump_mra_trees.cpp).
- Seed tools write response_metadata.json/root archives or ground-state restartdata and may replace existing bundle paths; production bundle helper honors usable-bundle skip. Dalton discovery can extract tar members through the imported helper [seed_from_dalton.cpp:main](./seed_from_dalton.cpp) [seed_es_from_hdf5.cpp:main](./seed_es_from_hdf5.cpp) [seed_moldft_from_dalton.cpp:main](./seed_moldft_from_dalton.cpp) [dalton_es_seed_bundle.cpp:main](./dalton_es_seed_bundle.cpp).
- Overlap/TPA tools print diagnostics and modify active in-memory representations, without numerical solver iterations or property metadata writes [es_overlap.cpp:run_overlap](./es_overlap.cpp) [tpa_from_dalton.cpp:main](./tpa_from_dalton.cpp).

## Source files

| File | Role |
| --- | --- |
| [dalton_gto.hpp](./dalton_gto.hpp) | Standalone Molden parser and normalized spherical Gaussian evaluator |
| [dalton_mra.hpp](./dalton_mra.hpp) | AO/occupied-virtual projection and occupied gauge rotation |
| [dalton_rspvec.hpp](./dalton_rspvec.hpp) | Binary Fortran-record parser and block splitting |
| [dalton_es_seed_bundle.cpp](./dalton_es_seed_bundle.cpp) | Production ES seed writer standalone |
| [dump_mra_trees.cpp](./dump_mra_trees.cpp) | Offline field/tree exporter |
| [es_overlap.cpp](./es_overlap.cpp) | Normalized root-overlap diagnostic |
| [seed_es_from_hdf5.cpp](./seed_es_from_hdf5.cpp) | Structured per-orbital HDF5 to one-root TDA seed |
| [seed_from_dalton.cpp](./seed_from_dalton.cpp) | Selected RSPVEC records to TDA/Full seed bundle |
| [seed_moldft_from_dalton.cpp](./seed_moldft_from_dalton.cpp) | Molden occupied orbitals to GS seed |
| [tpa_from_dalton.cpp](./tpa_from_dalton.cpp) | Imported-vector-only TPA diagnostic |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

Core tools are targets in the parent CMake manifest; no tool-specific CTest registrations are defined. Parser/import testing is registered as test_dalton_import [../CMakeLists.txt:test_dalton_import](../CMakeLists.txt).

```sh
cmake --build <build-dir> --target dump_mra_trees seed_from_dalton dalton_es_seed_bundle seed_moldft_from_dalton tpa_from_dalton test_dalton_import
```

```sh
ctest --test-dir <build-dir> -R "^madness/test/molresponse/test_dalton_import/run$"
```

seed_es_from_hdf5 and es_overlap targets exist only with MADNESS_ENABLE_HDF5; HTG/AMR capability on dump_mra_trees additionally requires MADNESS_ENABLE_VTK [../CMakeLists.txt:MADNESS_ENABLE_HDF5](../CMakeLists.txt) [../CMakeLists.txt:MADNESS_ENABLE_VTK](../CMakeLists.txt).
