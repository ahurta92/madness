#examples/myexamples

set(EXAMPLE_SOURCES
  mhehf
)

#if(LIBXC_FOUND)
#  list(APPEND EXAMPLE_SOURCES hefxc)
#endif()

# Create executables for example applications
foreach(example ${EXAMPLE_SOURCES})
  add_executable(${example} ${example}.cc)
  target_link_libraries(${example} MADchem MADgit ${DISABLEPIE_LINKER_FLAG})

endforeach()
