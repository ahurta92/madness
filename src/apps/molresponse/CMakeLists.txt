# src/apps/molresponse

if(CMAKE_CXX_STANDARD GREATER_EQUAL 17)
    # C++ standard is 17 or higher
    # Add your code or configurations here
else()
    # C++ standard is below 17
    # Handle the case where the required standard is not met
endif()


set(MOLRESPONSE_SOURCES
        ResponseBase.cpp
        ExcitedResponse.cpp
        FrequencyResponse.cpp
        property.cc
        global_functions.cc
        timer.cc
        Plot_VTK.cc
        basic_operators.cc
        response_parameters.cpp
        x_space.cc
        )
set(MOLRESPONSE_HEADERS
        ResponseBase.hpp
        ExcitedResponse.hpp
        FrequencyResponse.hpp
        response_macrotask.hpp
        property.h
        x_space.h
        global_functions.h
        timer.h
        Plot_VTK.h
        basic_operators.h
        response_parameters.h
        )
# Create the MADchem library
add_mad_library(response MOLRESPONSE_SOURCES MOLRESPONSE_HEADERS "chem" "madness/chem/molresponse/")
set(targetname MADresponse)
add_dependencies(applications-madness molresponse)

add_mad_executable(molresponse "molresponse.cc" "MADresponse;MADchem")

install(TARGETS molresponse DESTINATION "${MADNESS_INSTALL_BINDIR}")


set(R_BASE_SOURCES
        ResponseBase.cpp
        ExcitedResponse.cpp
        FrequencyResponse.cpp
        property.cc
        global_functions.cc
        timer.cc
        basic_operators.cc
        Plot_VTK.cc
        response_parameters.cpp
        x_space.cc
        )

set(R_BASE_HEADERS
        ResponseBase.hpp
        ExcitedResponse.hpp
        FrequencyResponse.hpp
        response_macrotask.hpp
        property.h
        global_functions.h
        x_space.h
        timer.h
        basic_operators.h
        Plot_VTK.h
        response_parameters.h
        )
# Create the MADchem library
add_mad_library(response_base R_BASE_SOURCES R_BASE_HEADERS "chem" "molresponse/")

set(targetname MADresponse_base)

set(RALL_SOURCES
        ResponseBase.cpp
        ExcitedResponse.cpp
        FrequencyResponse.cpp
        property.cc
        global_functions.cc
        timer.cc
        Plot_VTK.cc
        basic_operators.cc
        response_parameters.cpp
        x_space.cc
        )

set(RALL_HEADERS
        ResponseBase.hpp
        ExcitedResponse.hpp
        FrequencyResponse.hpp
        response_macrotask.hpp
        property.h
        global_functions.h
        timer.h
        basic_operators.h
        Plot_VTK.h
        response_parameters.h
        x_space.h
        )
add_mad_library(all_response RALL_SOURCES RALL_HEADERS "chem" "molresponse/")
set(targetname MADall_response)
