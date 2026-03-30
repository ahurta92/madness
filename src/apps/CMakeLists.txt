# src/apps

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(plot)
add_subdirectory(tdse)
add_subdirectory(moldft)
#add_subdirectory(molresponse)
# Frozen compatibility baseline for excited-state parity/reference work.
add_subdirectory(molresponse_legacy)
add_subdirectory(molresponse_v2)
add_subdirectory(molresponse_v3)
add_subdirectory(madqc_v2)
add_subdirectory(pno)
add_subdirectory(nemo)
add_subdirectory(cis)
add_subdirectory(dirac)
add_subdirectory(mp2)
add_subdirectory(cc2)
add_subdirectory(oep)
add_subdirectory(znemo)
add_subdirectory(zcis)
#add_subdirectory(hf)
#add_subdirectory(interior_bc)
