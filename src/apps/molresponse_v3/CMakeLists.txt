# src/apps/molresponse_v3/CMakeLists.txt

if(NOT TARGET MADchem)
  message(STATUS "molresponse_v3: MADchem not available, skipping")
  return()
endif()

add_mad_executable(molresponse_v3 "main.cpp;GroundState.cpp" "MADresponse2;MADchem")
add_mad_executable(test_v3_solver "test_solver.cpp;GroundState.cpp" "MADresponse2;MADchem")
add_dependencies(applications-madness molresponse_v3)

install(TARGETS molresponse_v3
        DESTINATION "${MADNESS_INSTALL_BINDIR}")
