# src/apps/bench
#
# Microbenchmarks for MADNESS MRA operations (scaling vs number of functions
# and protocol = thresh, k).  See README.md.

add_mad_executable(bench_mra_ops bench_mra_ops.cc "MADchem")
add_dependencies(applications-madness bench_mra_ops)

install(TARGETS bench_mra_ops DESTINATION "${MADNESS_INSTALL_BINDIR}")
