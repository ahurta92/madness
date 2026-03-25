# CMakeLists.txt for molresponse_legacy
#
# Frozen compatibility baseline imported from bsundahl/madness
# (src/apps/molresponse).
# Intended as a validation reference for molresponse_v2 closed-shell
# excited-state reintegration.
# Do not add new features here; keep local changes minimal and mark them
# // LEGACY_PATCH.

cmake_minimum_required(VERSION 3.12)

if(NOT TARGET MADchem)
  message(STATUS "molresponse_legacy: MADchem not available, skipping")
  return()
endif()

set(MOLRESPONSE_LEGACY_SOURCES
  TDDFT.cc
  density.cc
  property.cc
  # property_functions.cc excluded: depends on property_operators.h (missing from bsundahl repo)
  global_functions.cc
  basic_operators.cc
  calc_runner.cc
  timer.cc
  Plot_VTK.cc
  iter_freq2.cc
  # iterate_freq.cc excluded: calls CreateCoulombDerivativeRF etc. (not defined anywhere)
  iterate_excited.cc
  iterate_gamma.cc
  # iterate_xy.cc excluded: depends on missing CreateCoulombDerivativeRF/CreateExchangeDerivativeRF
  #   which are not in any bsundahl source; functions defined here (CreateGamma, ComputeHf,
  #   ComputeGf) are also not called from any other file.
  Lambda_X.cc
  Theta_X.cc
  molresponse.cc
)

add_executable(molresponse_legacy ${MOLRESPONSE_LEGACY_SOURCES})

target_link_libraries(molresponse_legacy MADchem)

target_include_directories(molresponse_legacy PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}   # also provides molresponse/ symlink -> .
  ${CMAKE_SOURCE_DIR}/src       # LEGACY_PATCH: <chem/...> -> src/chem/ symlink -> madness/chem/
  ${CMAKE_SOURCE_DIR}/src/madness  # LEGACY_PATCH: <chem/...> angle-bracket includes
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_BINARY_DIR}/src/madness
)

# LEGACY_PATCH: force-include preamble before all TUs to inject 'using namespace madness'
# without editing every source file (bsundahl code assumes this everywhere).
target_compile_options(molresponse_legacy PRIVATE
  -w
  -include "${CMAKE_CURRENT_SOURCE_DIR}/legacy_preamble.h"
)

add_dependencies(applications-madness molresponse_legacy)
