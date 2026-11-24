# tools component

## Overview

Offline response/ground-state seed conversion, Dalton contraction diagnostics, excited-state root overlap and tree/coefficient exports.

## Interface

- dump_mra_trees loads the ground MOs of a madqc run folder (--run-dir, default the current directory, older name --fd-calc-dir) from its node states/ground_state/ at the level the description records, with its Fock file from the node archive when one was saved and the molecule from analysis/ground_state/, and, with --fd/--es, its saved driven states and excited-state sets, and exports reconstructed/compressed JSON trees, geometry VTK, optional dense cube, HTG, adaptive AMR and native coefficients [dump_mra_trees.cpp:main](./dump_mra_trees.cpp) [dump_mra_trees.cpp:dump_function_trees](./dump_mra_trees.cpp).
- dalton_gto/read_molden parse spherical s through g Molden basis/MO data; dalton_rspvec parses Fortran sequential records and splits occupied-virtual TDA/RPA blocks; dalton_mra projects AO combinations and rotates the occupied gauge [dalton_gto.hpp:read_molden](./dalton_gto.hpp) [dalton_rspvec.hpp:read_rspvec](./dalton_rspvec.hpp) [dalton_rspvec.hpp:split_ov](./dalton_rspvec.hpp) [dalton_mra.hpp:occupied_gauge_rotation](./dalton_mra.hpp).
- seed_from_dalton converts selected raw RSPVEC entries into a TDA/Full seed of the excited-state set, and seed_es_from_hdf5 a one-root TDA seed from structured per-orbital functions; both write <run>/seeds/response/excited-<model>/ with the run's ground-state model as its settings. seed_moldft_from_dalton writes ground-state restartdata [seed_from_dalton.cpp:main](./seed_from_dalton.cpp) [seed_es_from_hdf5.cpp:main](./seed_es_from_hdf5.cpp) [seed_moldft_from_dalton.cpp:main](./seed_moldft_from_dalton.cpp).
- es_overlap prints the response-metric overlap of every root of an excited-state set or seed folder with every root of another; tpa_from_dalton projects external vectors and evaluates legacy/residue diagnostic contractions without iterative solves [es_overlap.cpp:main](./es_overlap.cpp) [tpa_from_dalton.cpp:main](./tpa_from_dalton.cpp).

## Invariants

- RSPVEC reader assumes little-endian 4-byte Fortran record markers; split_ov accepts exactly nocc*nvir or twice that size. Molden parser assumes pure spherical s,p,d,f,g and ignores spin headers [dalton_rspvec.hpp:read_fortran_record](./dalton_rspvec.hpp) [dalton_rspvec.hpp:split_ov](./dalton_rspvec.hpp) [dalton_gto.hpp:read_molden](./dalton_gto.hpp).
- Occupied gauge rotation uses M(M^T M)^(-1/2) and rejects minimum overlap eigenvalue below 0.5 [dalton_mra.hpp:occupied_gauge_rotation](./dalton_mra.hpp).
- seed_from_dalton, seed_moldft_from_dalton, seed_es_from_hdf5 and tpa_from_dalton enforce one rank. Seeds are starting guesses only (kind "seed", no metrics); scale defaults sqrt(2) for seed tools [seed_from_dalton.cpp:main](./seed_from_dalton.cpp) [seed_moldft_from_dalton.cpp:main](./seed_moldft_from_dalton.cpp) [seed_es_from_hdf5.cpp:main](./seed_es_from_hdf5.cpp) [tpa_from_dalton.cpp:main](./tpa_from_dalton.cpp).
- seed_from_dalton with archive rotates to its occupied gauge but does not Q-project in that tool; without archive it uses [-200,200]^3 and warns that seed is unrotated [seed_from_dalton.cpp:main](./seed_from_dalton.cpp).
- Exporter collective tree operations run on all ranks with rank-zero file emission [dump_mra_trees.cpp:write_tree_json](./dump_mra_trees.cpp).
- dump_mra_trees, es_overlap and the kernel tests' --calc-dir (which also gives them the ground state) find saved states by listing the node folders under <run>/states/response/ and reading their descriptions; they never rebuild node names, never read response_metadata.json, and never write into the run folder (a save in progress is read past, not finished). A half-saved or unreadable node stops the tool. Each loaded node is printed with its folder, kind, protocol and save time [saved_states.hpp:response_nodes](../solvers/saved_states.hpp) [saved_states.hpp:read_saved_node](../solvers/saved_states.hpp).
- HTG/AMR need VTK; coefficients need HDF5 and one rank. AMR samples cell centers, rebases to populated minimum level and warns about empty interior levels; XML writers use ASCII. The documented --maxlevel is not parsed/used by main [dump_mra_trees.cpp:write_amr](./dump_mra_trees.cpp) [dump_mra_trees.cpp:write_htg](./dump_mra_trees.cpp) [dump_mra_trees.cpp:main](./dump_mra_trees.cpp).

## Side Effects

- Export creates out directory (default mra_trees), geometry.vtk and per-field files, changes active FunctionDefaults across saved protocols and reconstructs/compresses shared in-memory handles [dump_mra_trees.cpp:main](./dump_mra_trees.cpp) [dump_mra_trees.cpp:dump_function_trees](./dump_mra_trees.cpp).
- Seed tools replace the seed folder they write, or write ground-state restartdata. Dalton discovery can extract tar members through the imported helper [seed_from_dalton.cpp:main](./seed_from_dalton.cpp) [seed_es_from_hdf5.cpp:main](./seed_es_from_hdf5.cpp) [seed_moldft_from_dalton.cpp:main](./seed_moldft_from_dalton.cpp).
- The TPA tool prints diagnostics and modifies active in-memory representations, without numerical solver iterations or property metadata writes [tpa_from_dalton.cpp:main](./tpa_from_dalton.cpp).

## Source files

| File | Role |
| --- | --- |
| [dalton_gto.hpp](./dalton_gto.hpp) | Standalone Molden parser and normalized spherical Gaussian evaluator |
| [dalton_mra.hpp](./dalton_mra.hpp) | AO/occupied-virtual projection and occupied gauge rotation |
| [dalton_rspvec.hpp](./dalton_rspvec.hpp) | Binary Fortran-record parser and block splitting |
| [dump_mra_trees.cpp](./dump_mra_trees.cpp) | Offline field/tree exporter |
| [es_overlap.cpp](./es_overlap.cpp) | Root-overlap diagnostic between two excited-state set or seed folders |
| [seed_es_from_hdf5.cpp](./seed_es_from_hdf5.cpp) | Structured per-orbital HDF5 to a one-root TDA seed in seeds/ |
| [seed_from_dalton.cpp](./seed_from_dalton.cpp) | Selected RSPVEC records to a TDA/Full seed in seeds/ |
| [seed_moldft_from_dalton.cpp](./seed_moldft_from_dalton.cpp) | Molden occupied orbitals to GS seed |
| [tpa_from_dalton.cpp](./tpa_from_dalton.cpp) | Imported-vector-only TPA diagnostic |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

Core tools are targets in the parent CMake manifest; no tool-specific CTest registrations are defined. Parser/import testing is registered as test_dalton_import [../CMakeLists.txt:test_dalton_import](../CMakeLists.txt).

```sh
cmake --build <build-dir> --target dump_mra_trees es_overlap seed_from_dalton seed_moldft_from_dalton tpa_from_dalton test_dalton_import
```

```sh
ctest --test-dir <build-dir> -R "^madness/test/molresponse/test_dalton_import/run$"
```

The seed_es_from_hdf5 target exists only with MADNESS_ENABLE_HDF5; HTG/AMR capability on dump_mra_trees additionally requires MADNESS_ENABLE_VTK [../CMakeLists.txt:MADNESS_ENABLE_HDF5](../CMakeLists.txt) [../CMakeLists.txt:MADNESS_ENABLE_VTK](../CMakeLists.txt).
