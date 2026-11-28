# The sole public driver for supported chemistry engines. Response remains an
# app-layer dependency so MADchem does not depend on src/apps.
add_mad_executable(madqc "madqc.cpp;../molresponse/GroundState.cpp" "MADchem")
add_dependencies(applications-madness madqc)

# Shared restart is native; optional response analysis/import tools still use
# the existing HDF5 implementation.
if(MADNESS_ENABLE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C)
  target_include_directories(madqc PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(madqc PRIVATE ${HDF5_LIBRARIES})
  target_compile_definitions(madqc PRIVATE MADNESS_HAS_HDF5)
endif()

install(TARGETS madqc DESTINATION "${MADNESS_INSTALL_BINDIR}")

if(BUILD_TESTING)
  add_mad_executable(test_shared_reference "test_shared_reference.cc;../molresponse/GroundState.cpp" "MADchem")
  add_mad_executable(test_workflow_completion test_workflow_completion.cc "MADchem")
  add_mad_executable(test_workflow_builders test_workflow_builders.cpp "MADchem")
  foreach(test test_shared_reference test_workflow_completion test_workflow_builders)
    add_test(NAME madness/test/madqc/${test}/run COMMAND $<TARGET_FILE:${test}>)
    set_tests_properties(madness/test/madqc/${test}/run PROPERTIES
      LABELS "madqc;applications;long" TIMEOUT 300)
  endforeach()
  set_tests_properties(madness/test/madqc/test_workflow_completion/run
    madness/test/madqc/test_workflow_builders/run PROPERTIES LABELS "unit;madqc;short")
  add_test(NAME madness/test/madqc/test_workflow_nemo_repeat/run
    COMMAND test_workflow_completion --nemo-repeat)
  set_tests_properties(madness/test/madqc/test_workflow_nemo_repeat/run
    PROPERTIES LABELS "madqc;applications;medium" TIMEOUT 180)
  if(ENABLE_MPI AND MPIEXEC_EXECUTABLE)
    foreach(test test_shared_reference test_workflow_completion)
      add_test(NAME madness/test/madqc/${test}/mpi2
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:${test}> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(madness/test/madqc/${test}/mpi2 PROPERTIES
        LABELS "madqc;applications;long" TIMEOUT 300
        ENVIRONMENT "MAD_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;OMPI_MCA_hwloc_base_binding_policy=none")
    endforeach()
  endif()

  # Scientific numbers live in src/examples/qc. These tests cover chained
  # runs, checkpoint integrity and relationships between reported quantities.
  # Unregistered (files kept) until they are ported to results.json (stage 5e):
  #   test_workflow_completion.py           expects the retired checkpoints/ layout
  #   test_shared_response_review.py        expects the retired checkpoints/ layout
  #   test_workflow_response.py             reads calc_info tasks[1] 2PA rows and metadata.excited_states
  #   test_shared_response.py               reads response_metadata.json counters and calc_info tasks[*].metadata
  #   test_shared_response_ladder.py        reads/edits response_metadata.json and calc_info response_properties
  #   test_shared_response_ladder_review.py reads/edits response_metadata.json and calc_info response_properties
  #   test_shared_quadratic_response.py     reads calc_info response_properties beta/raman/tpa
  #   test_shared_derived_response.py       reads calc_info response_properties beta/raman/tpa
  #   test_shared_response_numerics.py      translates calc_info check keys and references into results.json
  #   test_molresponse_h2o_alpha_beta_z.py  compares calc_info response_properties with a calc_info reference
  #   test_molresponse_lih_alpha_raman_beta_xyz.py  compares calc_info response_properties with a calc_info reference
  if(Python3_EXECUTABLE)
    foreach(test test_stdio_review test_shared_scf test_scf_recovery test_workflow_requalification
        test_driver_contract test_workflow_seed test_shared_response_refinement)
      set(script_args)
      if(test MATCHES "^test_(workflow_seed)$")
        list(APPEND script_args --source ${PROJECT_SOURCE_DIR})
      endif()
      set(test_name ${test})
      foreach(mode run mpi2)
        if(mode STREQUAL "mpi2" AND NOT (ENABLE_MPI AND MPIEXEC_EXECUTABLE))
          continue()
        endif()
        set(mpi_args)
        if(mode STREQUAL "mpi2")
          set(mpi_args --ranks 2 --mpiexec ${MPIEXEC_EXECUTABLE})
        endif()
        add_test(NAME madness/test/madqc/${test_name}/${mode}
          COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.py
            --binary $<TARGET_FILE:madqc> ${script_args} ${mpi_args})
        set_tests_properties(madness/test/madqc/${test_name}/${mode} PROPERTIES
          LABELS "madqc;applications;verylong" TIMEOUT 1200
          ENVIRONMENT "MAD_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;OMPI_MCA_hwloc_base_binding_policy=none")
      endforeach()
    endforeach()
    add_test(NAME madness/test/madqc/test_shared_response_default_families/run
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_shared_response_defaults.py
        --binary $<TARGET_FILE:madqc> --families)
    set_tests_properties(madness/test/madqc/test_shared_response_default_families/run PROPERTIES
      LABELS "madqc;applications;verylong" TIMEOUT 660)
    foreach(check defaults seed)
      set(_shared_extra_args)
      if(check STREQUAL "seed")
        list(APPEND _shared_extra_args --source ${PROJECT_SOURCE_DIR})
      endif()
      add_test(NAME madness/test/madqc/test_shared_response_${check}/run
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_shared_response_${check}.py
          --binary $<TARGET_FILE:madqc> ${_shared_extra_args})
      set_tests_properties(madness/test/madqc/test_shared_response_${check}/run PROPERTIES
        LABELS "madqc;applications;verylong" TIMEOUT 1800
        ENVIRONMENT "MAD_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;OMPI_MCA_hwloc_base_binding_policy=none")
    endforeach()
    foreach(model tda full)
      set(_shared_es_args)
      if(model STREQUAL "full")
        list(APPEND _shared_es_args --full)
      endif()
      add_test(NAME madness/test/madqc/test_shared_excited_response_${model}/run
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_shared_excited_response.py
          --binary $<TARGET_FILE:madqc> ${_shared_es_args})
      set_tests_properties(madness/test/madqc/test_shared_excited_response_${model}/run PROPERTIES
        LABELS "madqc;applications;verylong" TIMEOUT 1200
        ENVIRONMENT "MAD_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;OMPI_MCA_hwloc_base_binding_policy=none")
    endforeach()
    add_test(NAME madness/test/madqc/test_shared_response_nuclear/run
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_shared_response_variants.py
        --binary $<TARGET_FILE:madqc> --nuclear)
    set_tests_properties(madness/test/madqc/test_shared_response_nuclear/run PROPERTIES
      LABELS "madqc;applications;long" TIMEOUT 600
      ENVIRONMENT "MAD_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;OMPI_MCA_hwloc_base_binding_policy=none")
    foreach(atom H F C)
      add_test(NAME madness/test/madqc/test_shared_response_variants_${atom}/run
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_shared_response_variants.py
          --binary $<TARGET_FILE:madqc> --atom ${atom})
      set(_variant_tier verylong)
      if(atom STREQUAL "H")
        set(_variant_tier long)
      endif()
      set_tests_properties(madness/test/madqc/test_shared_response_variants_${atom}/run PROPERTIES
        LABELS "madqc;applications;${_variant_tier}" TIMEOUT 900
        ENVIRONMENT "MAD_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;OMPI_MCA_hwloc_base_binding_policy=none")
    endforeach()
    add_test(NAME madness/test/madqc/test_scf_final_evidence/run
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_scf_final_evidence.py
        --binary $<TARGET_FILE:madqc>)
    set_tests_properties(madness/test/madqc/test_scf_final_evidence/run
      PROPERTIES LABELS "madqc;applications;long" TIMEOUT 240)
    add_test(NAME madness/test/madqc/test_workflow_seed_accuracy/run
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_workflow_seed.py
        --binary $<TARGET_FILE:madqc> --source ${PROJECT_SOURCE_DIR} --numerical)
    set_tests_properties(madness/test/madqc/test_workflow_seed_accuracy/run
      PROPERTIES LABELS "madqc;applications;verylong" TIMEOUT 900)
  endif()

  # These larger scientific comparisons retain their original references.
  add_scripted_tests(test_cis_symmetry_h2o.py madqc "verylong;applications")
endif()
