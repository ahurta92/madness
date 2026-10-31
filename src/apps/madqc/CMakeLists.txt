# src/apps/moldft

# madqc also compiles molresponse's GroundState.cpp so the v3 engine
# (molresponse_v3_lib, included header-only via madqc_adapter.hpp) links (R3).
add_mad_executable(madqc "madqc.cpp;../molresponse/GroundState.cpp" "MADchem")
add_dependencies(applications-madness madqc)

# io-hdf5 thread: make madqc HDF5-aware so the deck opt-in (`io backend hdf5`
# / the `response.hdf5` alias) works through the madqc driver. The option is
# declared in ../molresponse/CMakeLists.txt (processed before this dir), and
# this mirrors its per-target wiring — without the define, madqc's TU compiles
# the no-HDF5 branch and throws on the opt-in even in an HDF5 build.
if(MADNESS_ENABLE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C)
  target_include_directories(madqc PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(madqc PRIVATE ${HDF5_LIBRARIES})
  target_compile_definitions(madqc PRIVATE MADNESS_HAS_HDF5)
endif()

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
    #
    # The "run a deck, compare result keys" tests below moved to qctest cases in
    # src/examples/qc (label `qctest`, see that README). Each one is now an input
    # deck a user can copy, run and diff, instead of CLI flags embedded in a
    # Python file, and the assertions are declared in check.json:
    #
    #   test_moldft_energy.py            -> scf_he_hf
    #   test_moldft_energy_mpi_parallel  -> scf_he_hf_mpi
    #   test_nemo_energy.py              -> nemo_he_hf
    #   test_nemo_localization.py        -> nemo_h2o_{canon,boys,new}
    #   test_cis_energy_he.py            -> cis_he_singlets
    #   test_oep_energy.py               -> oep_be_oaep
    #   test_mp2_helium.py               -> mp2_he_corr
    #   test_lrcc2_helium.py             -> lrcc2_he_excited
    #   test_cc2_callable.py             -> dropped, not ported: every numeric
    #       value in its reference was 0.0 (maxiter=1 with freeze 1 freezes
    #       helium's only orbital), so all four of its numeric assertions
    #       compared 0.0 against 0.0 and could not fail. A `calc_type cc2`
    #       ground state costs 1833 s, and lrcc2_he_excited covers it anyway.
    #
    # Retiring these also retires their references: test_mp2_helium's
    # mp2_correlation_energy was stale by 1e-4 against a current run.
    #
    # What stays here are the tests of *behavior* rather than of numbers -- they
    # chain runs or compare structures that a declarative check.json cannot
    # express.
    add_scripted_tests(test_cis_symmetry_h2o.py madqc "verylong;applications")
    add_scripted_tests(test_molresponse_h2o_alpha_beta_z.py madqc "verylong;applications")
    # Driver-hardening gate: SCF restart, refusal to answer a changed
    # Hamiltonian out of the checkpoint, and corrupt-checkpoint recovery. Four
    # sequential madqc runs on helium, so it is heavier than a single-shot test.
    # Registered commented out when it was written ("pending on-alloc build")
    # and never enabled, which is how two wrong paths in it went unnoticed --
    # it looked for the archive and the checkpoint in the cwd rather than in
    # <prefix>/task_0/moldft/. Both fixed; enabled here.
    add_scripted_tests(test_moldft_restart_e2e.py madqc "long;applications")
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

if(BUILD_TESTING)
  add_mad_executable(test_named_stages "test_named_stages.cc" "MADchem")
  add_test(NAME madness/test/madqc/test_named_stages/run COMMAND test_named_stages)
  set_tests_properties(madness/test/madqc/test_named_stages/run
    PROPERTIES LABELS "unittests;short" TIMEOUT 60)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_named_stages/mpi2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_named_stages> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(madness/test/madqc/test_named_stages/mpi2
      PROPERTIES LABELS "unittests;short" TIMEOUT 60)
  endif()
  # D-078, D-079: verdicts replayed from recorded files (the from-store tier).
  add_mad_executable(test_named_replay "test_named_replay.cc" "MADchem")
  add_test(NAME madness/test/madqc/test_named_replay/run
    COMMAND test_named_replay ${CMAKE_CURRENT_SOURCE_DIR}/test_named_replay_fixture.json)
  set_tests_properties(madness/test/madqc/test_named_replay/run
    PROPERTIES LABELS "unittests;short" TIMEOUT 60)
  add_mad_executable(test_response_properties "test_response_properties.cc" "MADchem")
  add_test(NAME madness/test/madqc/test_response_properties/run COMMAND test_response_properties)
  set_tests_properties(madness/test/madqc/test_response_properties/run
    PROPERTIES LABELS "unittests;short" TIMEOUT 60)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_response_properties/mpi2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_response_properties> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(madness/test/madqc/test_response_properties/mpi2
      PROPERTIES LABELS "unittests;short" TIMEOUT 60)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_named_stages_e2e/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_stages_e2e.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_named_stages_e2e/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_named_stages_e2e/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_stages_e2e.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_named_stages_e2e/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  # D-079, D-080: continued from the stored generations (make_test_store.py);
  # skipped unless MADQC_TEST_STORE is set.
  add_test(NAME madness/test/madqc/test_named_from_store/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_from_store.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
      --mpiexec ${MPIEXEC_EXECUTABLE})
  set_tests_properties(madness/test/madqc/test_named_from_store/run
    PROPERTIES LABELS "applications;medium" TIMEOUT 900 SKIP_RETURN_CODE 77)
  add_test(NAME madness/test/madqc/test_named_checkpoint_e2e/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_checkpoint_e2e.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_named_checkpoint_e2e/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 1200)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_named_checkpoint_e2e/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_checkpoint_e2e.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_named_checkpoint_e2e/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 1200)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_named_econv/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_econv.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_named_econv/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 1200)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_named_econv/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_econv.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_named_econv/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 1200)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_named_es_econv/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_es_econv.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_named_es_econv/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 1200)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_named_es_econv/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_es_econv.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_named_es_econv/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 1200)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_named_es_restart/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_es_restart.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_named_es_restart/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 1200)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_named_es_restart/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_es_restart.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_named_es_restart/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 1200)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_named_checkpoint_cadence/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_checkpoint_cadence.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_named_checkpoint_cadence/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 1200)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_named_checkpoint_cadence/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_checkpoint_cadence.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_named_checkpoint_cadence/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 1200)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_convergence_modes/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_convergence_modes.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_convergence_modes/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 1500)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_convergence_modes/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_convergence_modes.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_convergence_modes/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 1500)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_response_convergence_modes/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_response_convergence_modes.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_response_convergence_modes/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 1800)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_response_convergence_modes/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_response_convergence_modes.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_response_convergence_modes/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 1800)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_named_reference_from_generation/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_reference_from_generation.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_named_reference_from_generation/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 1800)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_named_reference_from_generation/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_reference_from_generation.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_named_reference_from_generation/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 1800)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE AND MADNESS_ENABLE_HDF5)
  add_test(NAME madness/test/madqc/test_named_checkpoint_hdf5/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_checkpoint_hdf5.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_named_checkpoint_hdf5/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_named_response_checkpoint/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_response_checkpoint.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_named_response_checkpoint/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 1800)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_named_response_checkpoint/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_response_checkpoint.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_named_response_checkpoint/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 1800)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_named_convergence/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_convergence.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_named_convergence/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 1200)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_named_convergence/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_convergence.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_named_convergence/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 1200)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_named_response_eligibility/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_response_eligibility.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_named_response_eligibility/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_named_response_eligibility/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_response_eligibility.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_named_response_eligibility/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_named_density_changes/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_density_changes.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_named_density_changes/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 600)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_named_density_changes/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_named_density_changes.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_named_density_changes/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 600)
  endif()
endif()

if(BUILD_TESTING AND Python3_Interpreter_FOUND)
  add_test(NAME madness/test/madqc/test_forced_rungs/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_forced_rungs.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_forced_rungs/run PROPERTIES LABELS "madqc;applications;verylong" TIMEOUT 700)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_forced_rungs/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_forced_rungs.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_forced_rungs/mpi2
      PROPERTIES LABELS "madqc;applications;verylong" TIMEOUT 700)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_explicit_accuracy/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_explicit_accuracy.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_explicit_accuracy/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_explicit_accuracy/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_explicit_accuracy.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_explicit_accuracy/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_response_configured_progress/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_response_configured_progress.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_response_configured_progress/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_response_configured_progress/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_response_configured_progress.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_response_configured_progress/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
  endif()
endif()

if(BUILD_TESTING AND Python3_EXECUTABLE)
  add_test(NAME madness/test/madqc/test_response_explicit_accuracy/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_response_explicit_accuracy.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_response_explicit_accuracy/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_response_explicit_accuracy/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_response_explicit_accuracy.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_response_explicit_accuracy/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
  endif()
  add_test(NAME madness/test/madqc/test_configured_progress/run
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_configured_progress.py
      --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/madqc/test_configured_progress/run
    PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
  if(ENABLE_MPI)
    add_test(NAME madness/test/madqc/test_configured_progress/mpi2
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_configured_progress.py
        --binary $<TARGET_FILE:madqc> --source ${CMAKE_CURRENT_SOURCE_DIR}
        --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
    set_tests_properties(madness/test/madqc/test_configured_progress/mpi2
      PROPERTIES LABELS "applications;verylong" TIMEOUT 900)
  endif()
endif()
