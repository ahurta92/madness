# src/apps/molresponse/CMakeLists.txt

if(NOT TARGET MADchem)
  message(STATUS "molresponse: MADchem not available, skipping")
  return()
endif()

# Optional VTK for native HyperTreeGrid (.htg) export from dump_mra_trees (--htg).
# OFF by default so a stock build needs no VTK; turn on with -DMADNESS_ENABLE_VTK=ON
# (point CMake at the cluster VTK via -DVTK_DIR=...). Only dump_mra_trees links it.
option(MADNESS_ENABLE_VTK "Build dump_mra_trees with native VTK HyperTreeGrid (.htg) export" OFF)
if(MADNESS_ENABLE_VTK)
  find_package(VTK 9.0 REQUIRED COMPONENTS CommonCore CommonDataModel IOXML)
  if(VTK_FOUND)
    message(STATUS "molresponse: VTK ${VTK_VERSION} found -> dump_mra_trees --htg enabled")
  endif()
endif()

# The solver core is covered by the shared-path tests and the madqc qctests.
add_mad_executable(test_protocol_key
                   "tests/test_protocol_key.cpp"
                   "MADchem")
add_mad_executable(test_bsh_shift
                   "tests/test_bsh_shift.cpp"
                   "MADchem")
# Per-root ES plateau detection: a stuck root is caught even when its |dw| has
# met the target (the solve-level max tracks miss it). Pure logic.
add_mad_executable(test_es_plateau
                   "tests/test_es_plateau.cpp"
                   "MADchem")
# The relative FD convergence gate (review C1): ||dx|| <= atol + rtol ||x||.
add_mad_executable(test_fd_gate
                   "tests/test_fd_gate.cpp"
                   "MADchem")
# Deck excited.* iteration keys -> ExecutorSettings (review finding C2).
add_mad_executable(test_deck_settings
                   "tests/test_deck_settings.cpp;GroundState.cpp"
                   "MADchem")
add_mad_executable(test_response_metadata
                   "tests/test_response_metadata.cpp"
                   "MADchem")
# dalton.dir import contract (showcase W3): locate-by-convention (+ tarball
# extraction), geometry fingerprint, .out method/basis parse, exact 1-to-1
# frequency matching. Pure filesystem/parsing — no World, no fixtures.
add_mad_executable(test_dalton_import
                   "tests/test_dalton_import.cpp;GroundState.cpp"
                   "MADchem")
add_mad_executable(test_response_property_planner
                   "tests/test_response_property_planner.cpp"
                   "MADchem")
add_mad_executable(test_state_metrics
                   "tests/test_state_metrics.cpp"
                   "MADchem")
# Full Raman Tier-B formulas (invariants, mode projection) on synthetic tensors.
add_mad_executable(test_raman_spectrum
                   "tests/test_raman_spectrum.cpp"
                   "MADtensor")
add_mad_executable(test_calc_manager
                   "tests/test_calc_manager.cpp"
                   "MADchem")
add_mad_executable(test_kernel_equivalence
                   "tests/test_kernel_equivalence.cpp;GroundState.cpp"
                   "MADchem")
# Pins the three algebraic identities (M1/M2/M3) behind the state-grouped
# reformulation of the 2PA two-electron residue.
add_mad_executable(test_tpa_regroup_identities
                   "tests/test_tpa_regroup_identities.cpp;GroundState.cpp"
                   "MADchem")
# Gate for the promoted state-free (P,Q) residue source (source-unification
# step 3, kernels/tpa_source_spec.hpp): equality vs tpa_e3_residue (2e) and
# tpa_moment_residue_1e (1e), plus the (P,Q)-vs-V^{bc} comparison report.
add_mad_executable(test_tpa_pq_vs_vbc
                   "tests/test_tpa_pq_vs_vbc.cpp;GroundState.cpp"
                   "MADchem")
# Term-RESOLVED comparison of the two quadratic sources (V^BC vs (P,Q)):
# every SourceEntry evaluated separately, aligned slots, per-term norms and
# contractions (2026-09-04 review: catch compensation the contracted scalar
# gates cannot see). Informational — the equality gates stay separate.
add_mad_executable(test_vbc_pq_terms
                   "tests/test_vbc_pq_terms.cpp;GroundState.cpp"
                   "MADchem")
add_mad_executable(test_step_restriction
                   "tests/test_step_restriction.cpp"
                   "MADchem")
# Slot-identity contract of rs::diagonalize (dominance order, phase fix,
# rank-deficient dim padding) — guards against re-introducing the
# ascending-omega re-sort that caused the coarse-rung slot churn.
add_mad_executable(test_diagonalize_slots
                   "tests/test_diagonalize_slots.cpp"
                   "MADchem")
# Inc 1 of the ES node-φ state-parallel design (docs/21): node-aligned subworld
# creation (one subworld per physical node via Split_type SHARED) + map test.
add_mad_executable(test_node_subworlds
                   "tests/test_node_subworlds.cpp"
                   "MADchem")
# Inc 2 (docs/21): prove ground state shipped into node-aligned subworlds is
# ONE distributed copy per node (not per-rank replicated); de-risk cross-world copy.
add_mad_executable(test_node_phi
                   "tests/test_node_phi.cpp"
                   "MADchem")
# Keystone of the persistent-subworld recompute design (docs/22 §8): assemble the
# subspace matrix from per-node partial columns + universe allreduce (no function
# movement), == direct universe matrix.
add_mad_executable(test_subspace_allreduce
                   "tests/test_subspace_allreduce.cpp"
                   "MADchem")
# Inc 1 of the FD exchange tensor layer (docs 26/27/28): the tensor-layer θ
# (build_ctx + assemble_theta; shared Tx + cached g0, dot-contracted) matches the
# reference compute_V0x − compute_E0x + compute_gamma (Static ClosedShell).
add_mad_executable(test_exchange_ctx
                   "tests/test_exchange_ctx.cpp"
                   "MADchem")
# S1 of the persistent-subworld ES design (docs/23): prove the GS-dependent build
# kernels (V0x/γ/E0x/T0x → Λ) reproduce bit-for-bit in a node-aligned subworld,
# the one piece the φ-count + A/S-scalar proofs did not cover.
add_mad_executable(test_es_build_subworld
                   "tests/test_es_build_subworld.cpp"
                   "MADchem")
# Offline MRA tree-skeleton dumper (R2/L4 structural export for viz/partitioning).
# Loads converged orbitals and dumps the adaptive octree per orbital; never solves.
add_mad_executable(dump_mra_trees
                   "tools/dump_mra_trees.cpp;GroundState.cpp"
                   "MADchem")
# tpa_from_dalton: contraction-only 2PA harness — reads DALTON's own converged
# vectors (RSPVEC: EXCITLAB X_f + DIPLEN N^A(w_f/2); molden.inp: occupied MOs),
# projects them to MRA, and passes them through tpa_sources/tpa_moment with ZERO
# solves. Equation-vs-vector discriminator against DALTON's printed S tensors
# (TPA_SCOPING.md §5k). No HDF5 needed (text molden + Fortran-binary RSPVEC).
add_mad_executable(tpa_from_dalton
                   "tools/tpa_from_dalton.cpp;GroundState.cpp"
                   "MADchem")
# Native HTG export is optional: only when -DMADNESS_ENABLE_VTK=ON found VTK.
# Guarded by MADNESS_HAS_VTK in the source; a no-VTK build is unchanged.
if(MADNESS_ENABLE_VTK AND VTK_FOUND)
  target_compile_definitions(dump_mra_trees PRIVATE MADNESS_HAS_VTK)
  target_link_libraries(dump_mra_trees PRIVATE ${VTK_LIBRARIES})
  vtk_module_autoinit(TARGETS dump_mra_trees MODULES ${VTK_LIBRARIES})
endif()
# --- io-hdf5 thread: HDF5-backed I/O (opt-in; docs/30) -----------------------
# Native MADNESS archive stays the default. No-op unless configured with
# -DMADNESS_ENABLE_HDF5=ON, so a normal build is byte-for-byte unchanged.
option(MADNESS_ENABLE_HDF5 "Enable HDF5-backed function/checkpoint I/O in molresponse (opt-in; native archive remains the default)" OFF)
if(MADNESS_ENABLE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C)
  message(STATUS "molresponse: HDF5 enabled — ${HDF5_LIBRARIES}")
  # P1 build-wiring smoke: pure HDF5 C-API round-trip, no MADNESS function code.
  add_mad_executable(test_hdf5_smoke "tests/test_hdf5_smoke.cpp" "MADworld")
  target_include_directories(test_hdf5_smoke PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_hdf5_smoke PRIVATE ${HDF5_LIBRARIES})
  target_compile_definitions(test_hdf5_smoke PRIVATE MADNESS_HAS_HDF5)
  # P2 Layer A: Function raw-coeff round-trip vs the native archive.
  add_mad_executable(test_function_hdf5 "tests/test_function_hdf5.cpp;GroundState.cpp" "MADchem")
  target_include_directories(test_function_hdf5 PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_function_hdf5 PRIVATE ${HDF5_LIBRARIES})
  target_compile_definitions(test_function_hdf5 PRIVATE MADNESS_HAS_HDF5)
  # Make the real solver binaries HDF5-aware so cm_run / madqc can use the opt-in
  # restart path in response_state.hpp. MADNESS_HAS_HDF5 only enables the path;
  # default (env MADRESPONSE_IO_HDF5 unset) still writes the legacy archive.
  # dump_mra_trees joins the list so its --coeffs native-coefficient export can call
  # save_function_hdf5 (the structured /meta+/keys+/coeffs writer in
  # solvers/function_hdf5_io.hpp); guarded by MADNESS_HAS_HDF5 in the source.
  foreach(_hdf5_tgt dump_mra_trees)
    if(TARGET ${_hdf5_tgt})
      target_include_directories(${_hdf5_tgt} PRIVATE ${HDF5_INCLUDE_DIRS})
      target_link_libraries(${_hdf5_tgt} PRIVATE ${HDF5_LIBRARIES})
      target_compile_definitions(${_hdf5_tgt} PRIVATE MADNESS_HAS_HDF5)
    endif()
  endforeach()
  # dalton-cis-seed thread: converter that turns pymra per-occ-orbital .mad.h5
  # response functions into a TDA excited-state seed, <run>/seeds/response/
  # excited-tda/, that the set starts from. HDF5-only (uses load_function_hdf5).
  add_mad_executable(seed_es_from_hdf5 "tools/seed_es_from_hdf5.cpp;GroundState.cpp" "MADchem")
  target_include_directories(seed_es_from_hdf5 PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(seed_es_from_hdf5 PRIVATE ${HDF5_LIBRARIES})
  target_compile_definitions(seed_es_from_hdf5 PRIVATE MADNESS_HAS_HDF5)
  # dalton-cis-seed (b): overlap-based pinning — identify which converged ES root
  # matches the operator-selected seed by |<seed|root>| (character, not energy).
  # Its reader of the older bundle layout is gone; node folders arrive in step 4c part 5d.
  add_mad_executable(es_overlap "tools/es_overlap.cpp;GroundState.cpp" "MADchem")
  target_include_directories(es_overlap PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(es_overlap PRIVATE ${HDF5_LIBRARIES})
  target_compile_definitions(es_overlap PRIVATE MADNESS_HAS_HDF5)
endif()

# dalton-cis-seed (c): native DALTON reader — reads RSPVEC + molden.inp directly,
# projects x_i onto MRA, and emits the ES seed without a Python/HDF5
# intermediate. Writes <run>/seeds/response/excited-<model>/ (native archives)
# that the set starts from — no HDF5 at runtime, and every HDF5 path in
# response_state.hpp is #ifdef MADNESS_HAS_HDF5 guarded, so this compiles in a
# plain (non-HDF5) build. Kept OUT of the MADNESS_ENABLE_HDF5 block so it is
# always available (enabling HDF5 here drags parallel-HDF5's Intel-MPI include
# path in and breaks the xeonmax toolchain — see es_bench MPI notes).
add_mad_executable(seed_from_dalton "tools/seed_from_dalton.cpp;GroundState.cpp" "MADchem")

# dalton-cis-seed (GS): project the occupied DALTON/molden MOs and write a
# moldft version-4 .restartdata (native path only, no HDF5). Ground-state analogue
# of seed_from_dalton; moldft --restart resumes from it (SCF::load_mos).
add_mad_executable(seed_moldft_from_dalton "tools/seed_moldft_from_dalton.cpp;GroundState.cpp" "MADchem")



# --- CTest registration (M3): tiered so a bare `ctest` is fast and green on ---
# --- any machine, and heavier tiers are opt-in by label / cache variable. ----
if (BUILD_TESTING)
  # Tier 1 — "unit": pure logic, no fixtures, no MPI launcher required.
  foreach(_t test_calc_manager test_protocol_key test_bsh_shift test_deck_settings test_es_plateau test_fd_gate
             test_raman_spectrum
             test_response_metadata
             test_response_property_planner test_state_metrics
             test_diagonalize_slots test_dalton_import)
    add_test(NAME madness/test/molresponse/${_t}/run COMMAND ${_t})
    set_property(TEST madness/test/molresponse/${_t}/run
                 PROPERTY LABELS "unit;response;short" APPEND)
  endforeach()

  # Tier 2 — "response-small": self-contained solver tests (SYNTHETIC ground
  # state; no fixture archives). Registered SERIAL on purpose: launcher choice
  # is cluster-specific (mpirun hangs in batch on some systems) and the
  # multi-rank coverage of these paths lives in the es_bench state-parallel
  # gates. Serial runs still validate the numerics end-to-end.
  foreach(_t test_exchange_ctx test_node_subworlds test_node_phi
             test_subspace_allreduce test_es_build_subworld test_step_restriction)
    add_test(NAME madness/test/molresponse/${_t}/run COMMAND ${_t})
    set_property(TEST madness/test/molresponse/${_t}/run
                 PROPERTY LABELS "response-small;response;medium" APPEND)
  endforeach()

  # Tier 3 — "response-fixture": need a moldft ground-state archive; registered
  # only when the fixture tree is available so CI without fixtures stays green.
  set(MADNESS_RESPONSE_FIXTURE_DIR "" CACHE PATH
      "fixture tree for response regression tests (<dir>/h2o/mad.restartdata.*)")
  if (MADNESS_RESPONSE_FIXTURE_DIR)
    add_test(NAME madness/test/molresponse/test_kernel_equivalence/run
             COMMAND test_kernel_equivalence
                     --archive=${MADNESS_RESPONSE_FIXTURE_DIR}/h2o/mad.restartdata)
    # Source-unification gates: the (P,Q) source and the M1/M3 regrouping
    # identities the (P,Q) promotion rests on.
    add_test(NAME madness/test/molresponse/test_tpa_pq_vs_vbc/run
             COMMAND test_tpa_pq_vs_vbc
                     --archive=${MADNESS_RESPONSE_FIXTURE_DIR}/h2o/mad.restartdata)
    add_test(NAME madness/test/molresponse/test_vbc_pq_terms/run
             COMMAND test_vbc_pq_terms
                     --archive=${MADNESS_RESPONSE_FIXTURE_DIR}/h2o/mad.restartdata)
    add_test(NAME madness/test/molresponse/test_tpa_regroup_identities/run
             COMMAND test_tpa_regroup_identities
                     --archive=${MADNESS_RESPONSE_FIXTURE_DIR}/h2o/mad.restartdata)
    foreach(_t test_kernel_equivalence test_tpa_pq_vs_vbc test_vbc_pq_terms
               test_tpa_regroup_identities)
      set_property(TEST madness/test/molresponse/${_t}/run
                   PROPERTY LABELS "response-fixture;response;long" APPEND)
    endforeach()
  endif()
endif (BUILD_TESTING)

# Per-function kinetic diagonals: analytic functions, no external fixtures.
if(BUILD_TESTING)
  add_mad_executable(test_es_kinetic_evidence "tests/test_es_kinetic_evidence.cpp" "MADchem")
  add_test(NAME madness/test/molresponse/test_es_kinetic_evidence/run COMMAND test_es_kinetic_evidence)
  set_tests_properties(madness/test/molresponse/test_es_kinetic_evidence/run
    PROPERTIES LABELS "unit;response;short" TIMEOUT 30)

  add_mad_executable(test_fd_explicit_accuracy "tests/test_fd_explicit_accuracy.cpp" "MADchem")
  add_test(NAME madness/test/molresponse/test_fd_explicit_accuracy/run COMMAND test_fd_explicit_accuracy)
  set_tests_properties(madness/test/molresponse/test_fd_explicit_accuracy/run
    PROPERTIES LABELS "unit;response;short" TIMEOUT 30)

  add_mad_executable(test_shared_iteration "tests/test_shared_iteration.cpp" "MADworld;nlohmann_json::nlohmann_json")
  add_mad_executable(test_shared_es_iteration "tests/test_shared_es_iteration.cpp" "MADchem")
  add_test(NAME madness/test/molresponse/test_shared_es_iteration/run COMMAND test_shared_es_iteration)
  set_tests_properties(madness/test/molresponse/test_shared_es_iteration/run PROPERTIES LABELS "unit;response;long" TIMEOUT 240)
  if(ENABLE_MPI AND MPIEXEC_EXECUTABLE)
    add_test(NAME madness/test/molresponse/test_shared_es_iteration/mpi2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_shared_es_iteration> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(madness/test/molresponse/test_shared_es_iteration/mpi2 PROPERTIES LABELS "unit;response;long" TIMEOUT 240
      ENVIRONMENT "OMPI_MCA_hwloc_base_binding_policy=none")
  endif()
  add_mad_executable(test_es_comparison "tests/test_es_comparison.cpp" "MADchem")
  add_test(NAME madness/test/molresponse/test_es_comparison/run COMMAND test_es_comparison)
  set_tests_properties(madness/test/molresponse/test_es_comparison/run PROPERTIES LABELS "unit;response;short" TIMEOUT 30)
  add_mad_executable(test_shared_response_review "tests/test_shared_response_review.cpp" "MADchem")
  add_test(NAME madness/test/molresponse/test_shared_response_review/run COMMAND test_shared_response_review)
  set_tests_properties(madness/test/molresponse/test_shared_response_review/run PROPERTIES LABELS "unit;response;short" TIMEOUT 60)
  if(ENABLE_MPI AND MPIEXEC_EXECUTABLE)
    add_test(NAME madness/test/molresponse/test_shared_response_review/mpi2
      COMMAND ${MPIEXEC_EXECUTABLE} -n 2 $<TARGET_FILE:test_shared_response_review>)
    set_tests_properties(madness/test/molresponse/test_shared_response_review/mpi2 PROPERTIES
      LABELS "unit;response;short" TIMEOUT 60 ENVIRONMENT "MAD_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;OMPI_MCA_hwloc_base_binding_policy=none")
  endif()
  add_mad_executable(test_shared_response_records "tests/test_shared_response_records.cpp" "MADchem")
  add_test(NAME madness/test/molresponse/test_shared_response_records/run COMMAND test_shared_response_records)
  set_tests_properties(madness/test/molresponse/test_shared_response_records/run PROPERTIES LABELS "unit;response;short" TIMEOUT 60)
  if(ENABLE_MPI AND MPIEXEC_EXECUTABLE)
    add_test(NAME madness/test/molresponse/test_shared_response_records/mpi2
      COMMAND ${MPIEXEC_EXECUTABLE} -n 2 $<TARGET_FILE:test_shared_response_records>)
    set_tests_properties(madness/test/molresponse/test_shared_response_records/mpi2 PROPERTIES
      LABELS "unit;response;short" TIMEOUT 60 ENVIRONMENT "MAD_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;OMPI_MCA_hwloc_base_binding_policy=none")
  endif()
  add_mad_executable(test_shared_fd_planning "tests/test_shared_fd_planning.cpp" "MADchem")
  add_test(NAME madness/test/molresponse/test_shared_fd_planning/run COMMAND test_shared_fd_planning)
  set_tests_properties(madness/test/molresponse/test_shared_fd_planning/run PROPERTIES LABELS "unit;response;short" TIMEOUT 30)
  add_test(NAME madness/test/molresponse/test_shared_iteration/run COMMAND test_shared_iteration)
  set_tests_properties(madness/test/molresponse/test_shared_iteration/run PROPERTIES LABELS "unit;molresponse;short" TIMEOUT 30)
  add_mad_executable(test_response_property_metadata "tests/test_response_property_metadata.cpp" "MADworld;nlohmann_json::nlohmann_json")
  if(ENABLE_MPI AND MPIEXEC_EXECUTABLE)
    foreach(mode read write)
      add_test(NAME madness/test/molresponse/test_response_property_metadata/${mode}_mpi2
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_response_property_metadata> ${mode} ${MPIEXEC_POSTFLAGS})
      set_tests_properties(madness/test/molresponse/test_response_property_metadata/${mode}_mpi2 PROPERTIES
        LABELS "unit;response;short" TIMEOUT 15
        ENVIRONMENT "MAD_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;OMPI_MCA_hwloc_base_binding_policy=none")
    endforeach()
  endif()

  add_mad_executable(test_iteration_protocol "tests/test_iteration_protocol.cpp" "MADworld;nlohmann_json::nlohmann_json")
  add_test(NAME madness/test/molresponse/test_iteration_protocol/run COMMAND test_iteration_protocol)
  set_tests_properties(madness/test/molresponse/test_iteration_protocol/run
                       PROPERTIES LABELS "unit;response;short")
  add_mad_executable(test_function_metric_summary "tests/test_function_metric_summary.cpp" "MADchem")
  add_test(NAME madness/test/molresponse/test_function_metric_summary/run COMMAND test_function_metric_summary)
  set_tests_properties(madness/test/molresponse/test_function_metric_summary/run
    PROPERTIES LABELS "unit;response;short" TIMEOUT 30)

  add_mad_executable(test_function_kinetics "tests/test_function_kinetics.cpp" "MADchem")
  add_test(NAME madness/test/molresponse/test_function_kinetics/run COMMAND test_function_kinetics)
  set_tests_properties(madness/test/molresponse/test_function_kinetics/run
    PROPERTIES LABELS "response;medium" TIMEOUT 120)
  if(ENABLE_MPI)
    add_test(NAME madness/test/molresponse/test_function_kinetics/mpi2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_function_kinetics> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(madness/test/molresponse/test_function_kinetics/mpi2
      PROPERTIES LABELS "response;medium" TIMEOUT 120)
  endif()

  add_mad_executable(test_function_residuals "tests/test_function_residuals.cpp" "MADchem")
  add_test(NAME madness/test/molresponse/test_function_residuals/run COMMAND test_function_residuals)
  set_tests_properties(madness/test/molresponse/test_function_residuals/run
    PROPERTIES LABELS "response;medium" TIMEOUT 120)
  if(ENABLE_MPI)
    add_test(NAME madness/test/molresponse/test_function_residuals/mpi2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_function_residuals> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(madness/test/molresponse/test_function_residuals/mpi2
      PROPERTIES LABELS "response;medium" TIMEOUT 120)
  endif()
endif()
