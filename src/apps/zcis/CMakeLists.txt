# src/apps/moldft

add_mad_executable(zcis zcis.cc "MADchem")
add_dependencies(applications-madness zcis)


install(TARGETS zcis DESTINATION "${MADNESS_INSTALL_BINDIR}")
