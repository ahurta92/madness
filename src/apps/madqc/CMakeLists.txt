# src/apps/moldft

# madqc also compiles molresponse's GroundState.cpp so the v3 engine
# (molresponse_v3_lib, included header-only via madqc_adapter.hpp) links (R3).
add_mad_executable(madqc "madqc.cpp;../molresponse/GroundState.cpp" "MADchem")
add_dependencies(applications-madness madqc)

if (BUILD_TESTING)
    add_mad_executable(test_workflow_builders test_workflow_builders.cpp "MADchem")
if (BUILD_TESTING)
  add_test(NAME madness/test/madqc/test_workflow_builders/run
           COMMAND test_workflow_builders)
  set_property(TEST madness/test/madqc/test_workflow_builders/run
               PROPERTY LABELS "unit;madqc;short" APPEND)
endif()
    add_dependencies(applications-madness test_workflow_builders)
    add_test(NAME test_workflow_builders COMMAND test_workflow_builders)

    # integration tests
    add_scripted_tests(test_nemo_energy.py madqc "medium;applications")
    add_scripted_tests(test_nemo_localization.py madqc "verylong;applications")
    add_scripted_tests(test_moldft_energy.py madqc "medium;applications")
    add_scripted_tests(test_moldft_energy_mpi_parallel.py madqc "short;applications")
    add_scripted_tests(test_cc2_callable.py madqc "short;applications")
    add_scripted_tests(test_oep_energy.py madqc "long;applications")
    add_scripted_tests(test_cis_energy_he.py madqc "medium;applications")
    add_scripted_tests(test_cis_symmetry_h2o.py madqc "verylong;applications")
    add_scripted_tests(test_mp2_helium.py madqc "verylong;applications")
    add_scripted_tests(test_lrcc2_helium.py madqc "verylong;applications")
    add_scripted_tests(test_molresponse_h2o_alpha_beta_z.py madqc "verylong;applications")
    # RE-ENABLED (R-1 close-out, 2026-07-10 vbc gate 2083856): rescoped to
    # molresponse_v3 — alpha xyz + static beta + SINGLE-COMPONENT raman
    # (atom 0, z) via requested_properties; reference regenerated from a v3 run
    # validated against the archived v2 values (beta_zzz 315.595 vs v2 315.559
    # = 0.011%; raman z_z_nuc0z -13.549 vs v2 -13.502 = 0.35%). v2's full
    # per-atom tensor remains archived at
    # madness_studies/refs/raman_v2_h2o_full_tensor.json (state-parallel scope).
    add_scripted_tests(test_molresponse_lih_alpha_raman_beta_xyz.py madqc "verylong;applications")
    # excited-state tests are WIP — disabled until Phase 4 is complete
    # add_scripted_tests(test_molresponse_excited_metadata_smoke.py madqc "short;applications")
    # add_scripted_tests(test_molresponse_excited_restart_reuse.py madqc "short;applications")
    # add_scripted_tests(test_molresponse_excited_protocol_projection.py madqc "short;applications")
    # add_scripted_tests(test_molresponse_excited_guess_fallback.py madqc "short;applications")
    # add_scripted_tests(test_molresponse_excited_restart_mode_compat.py madqc "short;applications")
    # add_scripted_tests(test_znemo_energy.py madqc "short;applications")
endif(BUILD_TESTING)

install(TARGETS madqc DESTINATION "${MADNESS_INSTALL_BINDIR}")
