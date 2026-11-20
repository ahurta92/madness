# SCF is exposed through madqc. Keep private symmetry and unsupported periodic
# DFT programs; removing the molecular executable must not drop these builds.
add_mad_executable(testpg testpg.cc MADchem)
add_mad_executable(testperiodicdft testperiodicdft.cc MADchem)
