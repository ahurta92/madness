# src/apps/molresponse_v3/CMakeLists.txt

if(NOT TARGET MADchem)
  message(STATUS "molresponse_v3: MADchem not available, skipping")
  return()
endif()

add_mad_executable(molresponse_v3 "main.cpp;GroundState.cpp" "MADresponse2;MADchem")
add_mad_executable(test_v3_solver "test_solver.cpp;GroundState.cpp" "MADresponse2;MADchem")
add_mad_executable(test_es_solver "test_es_solver.cpp;GroundState.cpp" "MADresponse2;MADchem")
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
add_mad_executable(test_kernel_equivalence
                   "tests/test_kernel_equivalence.cpp;GroundState.cpp"
                   "MADresponse2;MADchem")
add_mad_executable(test_step_restriction
                   "tests/test_step_restriction.cpp"
                   "MADresponse2;MADchem")
add_dependencies(applications-madness molresponse_v3)

install(TARGETS molresponse_v3
        DESTINATION "${MADNESS_INSTALL_BINDIR}")
