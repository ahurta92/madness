# src/apps/molresponse_v3/CMakeLists.txt

if(NOT TARGET MADchem)
  message(STATUS "molresponse_v3: MADchem not available, skipping")
  return()
endif()

# Optional VTK for native HyperTreeGrid (.htg) export from dump_mra_trees (--htg).
# OFF by default so a stock build needs no VTK; turn on with -DMADNESS_ENABLE_VTK=ON
# (point CMake at the cluster VTK via -DVTK_DIR=...). Only dump_mra_trees links it.
option(MADNESS_ENABLE_VTK "Build dump_mra_trees with native VTK HyperTreeGrid (.htg) export" OFF)
if(MADNESS_ENABLE_VTK)
  find_package(VTK 9.0 REQUIRED COMPONENTS CommonCore CommonDataModel IOXML)
  if(VTK_FOUND)
    message(STATUS "molresponse_v3: VTK ${VTK_VERSION} found -> dump_mra_trees --htg enabled")
  endif()
endif()

add_mad_executable(molresponse_v3 "main.cpp;GroundState.cpp" "MADresponse2;MADchem")
# Legacy direct-FD test binaries (test_v3_solver / test_es_solver) retired in R0b
# along with the top-level FDSolver.hpp / ESSolver.hpp they exercised; the solver
# core is covered by test_calc_manager_run / test_run_response / test_v3_*_skeleton.
add_mad_executable(test_v3_es_skeleton
                   "tests/test_v3_es_skeleton.cpp;GroundState.cpp"
                   "MADresponse2;MADchem")
add_mad_executable(test_v3_fd_skeleton
                   "tests/test_v3_fd_skeleton.cpp;GroundState.cpp"
                   "MADresponse2;MADchem")
add_mad_executable(test_protocol_key
                   "tests/test_protocol_key.cpp"
                   "MADresponse2;MADchem")
add_mad_executable(test_response_metadata
                   "tests/test_response_metadata.cpp"
                   "MADresponse2;MADchem")
add_mad_executable(test_response_property_planner
                   "tests/test_response_property_planner.cpp"
                   "MADresponse2;MADchem")
add_mad_executable(test_state_metrics
                   "tests/test_state_metrics.cpp"
                   "MADresponse2;MADchem")
add_mad_executable(test_calc_manager
                   "tests/test_calc_manager.cpp"
                   "MADresponse2;MADchem")
add_mad_executable(test_calc_manager_run
                   "tests/test_calc_manager_run.cpp;GroundState.cpp"
                   "MADresponse2;MADchem")
add_mad_executable(test_run_response
                   "tests/test_run_response.cpp;GroundState.cpp"
                   "MADresponse2;MADchem")
add_mad_executable(test_kernel_equivalence
                   "tests/test_kernel_equivalence.cpp;GroundState.cpp"
                   "MADresponse2;MADchem")
add_mad_executable(test_step_restriction
                   "tests/test_step_restriction.cpp"
                   "MADresponse2;MADchem")
# Inc 1 of the ES node-φ state-parallel design (docs/21): node-aligned subworld
# creation (one subworld per physical node via Split_type SHARED) + map test.
add_mad_executable(test_node_subworlds
                   "tests/test_node_subworlds.cpp"
                   "MADresponse2;MADchem")
# Inc 2 (docs/21): prove ground state shipped into node-aligned subworlds is
# ONE distributed copy per node (not per-rank replicated); de-risk cross-world copy.
add_mad_executable(test_node_phi
                   "tests/test_node_phi.cpp"
                   "MADresponse2;MADchem")
# Keystone of the persistent-subworld recompute design (docs/22 §8): assemble the
# subspace matrix from per-node partial columns + universe allreduce (no function
# movement), == direct universe matrix.
add_mad_executable(test_subspace_allreduce
                   "tests/test_subspace_allreduce.cpp"
                   "MADresponse2;MADchem")
# Offline MRA tree-skeleton dumper (R2/L4 structural export for viz/partitioning).
# Loads converged orbitals and dumps the adaptive octree per orbital; never solves.
add_mad_executable(dump_mra_trees
                   "tools/dump_mra_trees.cpp;GroundState.cpp"
                   "MADresponse2;MADchem")
# Native HTG export is optional: only when -DMADNESS_ENABLE_VTK=ON found VTK.
# Guarded by MADNESS_HAS_VTK in the source; a no-VTK build is unchanged.
if(MADNESS_ENABLE_VTK AND VTK_FOUND)
  target_compile_definitions(dump_mra_trees PRIVATE MADNESS_HAS_VTK)
  target_link_libraries(dump_mra_trees PRIVATE ${VTK_LIBRARIES})
  vtk_module_autoinit(TARGETS dump_mra_trees MODULES ${VTK_LIBRARIES})
endif()
add_dependencies(applications-madness molresponse_v3)

install(TARGETS molresponse_v3
        DESTINATION "${MADNESS_INSTALL_BINDIR}")
