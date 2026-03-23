# src/apps/molresponse_v2/CMakeLists.txt

if(CMAKE_CXX_STANDARD GREATER_EQUAL 17)
    message(STATUS "Using C++17 or higher standard")
else()
    message(FATAL_ERROR "C++17 or higher is required.")
endif()

set(MOLRESPONSE_SOURCES2
    ResponseManager.cpp
    GroundStateData.cpp
    ResponseSolver.cpp
    ResponseState.cpp
    ResponseInitializer.cpp
    FrequencyLoop.cpp
    ExcitedStateBundleSolver.cpp
    ExcitedResponse.cpp
)

set(MOLRESPONSE_HEADERS2
    GroundStateData.hpp
    MolecularProperty.hpp
    ResponseBundle.hpp
    ResponseManager.hpp
    ResponseState.hpp
    ResponseIO.hpp
        ResponseRecord.hpp
    DerivedStatePlanner.hpp
    StateParallelPlanner.hpp
    FrequencyLoop.hpp
    PropertyManager.hpp
    ResponseSolver.hpp
    ResponseSolverUtils.hpp
    ResponseKernels.hpp
    ResponseVectorKernels.hpp
    excited_ops/RestrictedBundleCommon.hpp
    excited_ops/RestrictedFullBundleOps.hpp
    excited_ops/RestrictedTDABundleOps.hpp
    ResponseInitializer.hpp
    ResponseDebugLogger.hpp
    ResponseDebugLoggerMacros.hpp
    broadcast_json.hpp
#    ResponseParameters.hpp
)

# Create the MADresponse2 library
add_mad_library(response2 MOLRESPONSE_SOURCES2 MOLRESPONSE_HEADERS2 "chem" "madness/chem/molresponse2/")
set(targetname MADresponse2)
option(MADNESS_ENABLE_LEGACY_EXCITED_BUNDLE_ADAPTER
       "Enable legacy molresponse-based excited-bundle adapter"
       OFF)
if (MADNESS_ENABLE_LEGACY_EXCITED_BUNDLE_ADAPTER)
    if (NOT TARGET MADresponse)
        message(FATAL_ERROR "MADNESS_ENABLE_LEGACY_EXCITED_BUNDLE_ADAPTER requires MADresponse target")
    endif()
    target_link_libraries(MADresponse2 PUBLIC MADresponse)
    target_compile_definitions(MADresponse2 PUBLIC MADNESS_ENABLE_LEGACY_EXCITED_BUNDLE_ADAPTER=1)
endif()
#add_dependencies(applications-madness test_preliminaries)
add_dependencies(applications-madness test_parameter_manager)
add_dependencies(applications-madness molresponse2)
add_dependencies(applications-madness test_tda_h2)
#add_mad_executable(test_preliminaries "test_preliminaries.cpp" "MADresponse2;MADchem")
add_mad_executable(test_parameter_manager "test_parameter_manager.cpp" "MADresponse2;MADchem")
add_mad_executable(molresponse2 "molresponse2.cpp" "MADresponse2;MADchem")
add_mad_executable(test_tda_h2 "test_tda_h2.cpp" "MADresponse2;MADchem")

# Add install targets for executables
install(TARGETS molresponse2 test_parameter_manager test_tda_h2
        DESTINATION "${MADNESS_INSTALL_BINDIR}")
