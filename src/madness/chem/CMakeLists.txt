# src/apps/chem

add_definitions(-DMRA_CHEMDATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

# Set the CHEM sources and header files
set(MADCHEM_HEADERS
    AC.h
    atomutil.h
    basis.h
    BSHApply.h
    CalculationParameters.h
    CalculationOutcome.h
    CC2.h
    mp3.h
    ccpairfunction.h
    CCPotentials.h
    CCStructures.h
    CCParameters.h
    corepotential.h
    correlationfactor.h
    diamagneticpotentialfactor.h
    dispersion.h
    distpm.h
    electronic_correlation_factor.h
    ESInterface.h
    exchangeoperator.h
    gaussian.h
    gth_pseudopotential.h
    GuessFactory.h
    IntegratorXX.h
    KineticEnergy.h
    DensityChanges.h
    FunctionResidual.h
    FunctionResidualReport.h
    ExplicitCriteria.h
    FunctionMetricSummary.h
    localizer.h
    lowrankfunction.h
    masks_and_boxes.h
    molecularbasis.h
    molecular_functors.h
    molecular_optimizer.h
    MolecularOrbitals.h
    molecule.h
    molopt.h
    mp2.h
    nemo.h
    NWChem.h
    oep.h
    pcm.h
    PNOF12Potentials.h
    PNOGuessFunctions.h
    PNO.h
    PNOParameters.h
    PNOStructures.h
    PNOTensors.h
    pointgroupoperator.h
    pointgroupsymmetry.h
    polynomial.h
    population.h
    potentialmanager.h
    ParameterManager.hpp
    ResponseParameters.hpp
    Restart.h
    MetricRecords.h
    RootCoordination.h
    RunFiles.h
    IterationMetrics.h
    IterationHistory.h
    IterationOutput.h
    NodeArchive.h
    TimingRecords.h
    MetricHelpers.h
    ConvergenceSettings.h
    LadderPlanner.h
    ScfContinuation.h
    RestartPlan.h
    PathManager.hpp
    Applications.hpp
    Drivers.hpp
    ResultsSummary.hpp
    VizManifest.hpp
    SCFTargetAdapter.hpp
    InputWriter.hpp
    CCLib.hpp
    projector.h
    Results.h
    QCPropertyInterface.h
    SAP.h
    SCF.h
    SCFCheckpoint.h
    SCFScientificSettings.h
    SCFInput.h
    SCFRun.h
    SCFEvaluationPolicy.h
    SCFOperators.h
    SCFProtocol.h
    TDHF.h
    vibanal.h
    write_test_input.h
    xcfunctional.h
    zcis.h
    znemo.h
)
set(MADCHEM_SOURCES
    SCFCheckpoint.cc
    SCFScientificSettings.cc
    SCFRun.cc
    AC.cc
    atomutil.cc
    CC2.cc
    mp3.cc
    ccpairfunction.cc
    CCPotentials.cc
    CCStructures.cc
    corepotential.cc
    correlationfactor.cc
    diamagneticpotentialfactor.cc
    dispersion.cc
    distpm.cc
    exchangeoperator.cc
    gaussian.cc
    gth_pseudopotential.cc
    GuessFactory.cc
    lda.cc
    localizer.cc
    molecularbasis.cc
    MolecularOrbitals.cc
    molecule.cc
    mp2.cc
    nemo.cc
    NWChem.cc
    InputWriter.cpp
    oep.cc
    pcm.cc
    pointgroupsymmetry.cc
    potentialmanager.cc
    polynomial.cc
    population.cc
    SAP.cc
    SCF.cc
    SCFOperators.cc
    TDHF.cc
    vibanal.cc
    zcis.cc
    znemo.cc
    PNO.cpp  PNOF12Potentials.cpp  PNOGuessFunctions.cpp  PNOParameters.cpp  PNOStructures.cpp
)
if(TARGET Libxc::xc)
  list(APPEND MADCHEM_SOURCES xcfunctional_libxc.cc)
else()
  list(APPEND MADCHEM_SOURCES xcfunctional_ldaonly.cc)
endif()

# Create the MADchem library
add_mad_library(chem MADCHEM_SOURCES MADCHEM_HEADERS "mra" "madness/chem/")

  # N.B. MADchem too flaky for unity builds
  set_target_properties(MADchem-obj PROPERTIES UNITY_BUILD FALSE)

  set(targetname MADchem)

  # Add dependencies for MADchem
  target_include_directories(${targetname} PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/apps/>)
  if (TARGET Libxc::xc)
    target_link_libraries(${targetname} PUBLIC Libxc::xc)
  endif ()

  if (PCM_FOUND)
    target_include_directories(${targetname} PUBLIC ${PCM_INCLUDE_DIRS})
    target_link_libraries(${targetname} PUBLIC ${PCM_LIBRARIES})
  endif ()

  if (INTEGRATORXX_FOUND)
    target_include_directories(${targetname} PUBLIC ${INTEGRATORXX_INCLUDE_DIRS})
    target_link_libraries(${targetname} PUBLIC ${INTEGRATORXX_LIBRARIES})
  endif ()

  # simple-dftd3 contributes a library to link and one header, used by exactly
  # one translation unit (dispersion.h deliberately does not include it).
  #
  # The header is reached by absolute path rather than by adding -I, the same
  # way MADNESS_MPI_HEADER reaches mpi.h. Its usual install prefix is a conda
  # environment whose include/ also carries a copy of nlohmann/json.hpp, and an
  # -I -- target-wide or source-scoped, the latter being placed ahead of the
  # target's -- would shadow the nlohmann_json MADNESS actually builds against.
  # One TU disagreeing with the rest about a header-only library's version is an
  # ODR violation, so the -I is avoided entirely rather than ordered around.
  if (TARGET s-dftd3::s-dftd3)
    target_link_libraries(${targetname} PUBLIC s-dftd3::s-dftd3)
  elseif (DFTD3_FOUND)
    target_link_libraries(${targetname} PUBLIC ${DFTD3_LIBRARIES})
  endif ()
  if (DFTD3_HEADER)
    set_source_files_properties(dispersion.cc PROPERTIES COMPILE_DEFINITIONS
        "MADNESS_DFTD3_HEADER=\"${DFTD3_HEADER}\"")
  endif ()

  if (BOOST_FOUND)
    target_link_libraries(${targetname} PUBLIC Boost::headers)
  endif ()

# Add executables
if (NOT MADNESS_BUILD_LIBRARIES_ONLY)
  add_mad_executable(plotxc "plotxc.cc;xcfunctional.h" "MADchem")
endif()

# Install files
install(FILES sto-3g sto-6g 6-31g coredata/mcp coredata/mcp2 coredata/mcp_guess structure_library
    DESTINATION "${MADNESS_INSTALL_DATADIR}"
    COMPONENT chem)

# Add unit tests    
if(BUILD_TESTING)
  
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  # The list of unit test source files
  set(CHEM_TEST_SOURCES_SHORT test_pointgroupsymmetry.cc test_masks_and_boxes.cc
      test_exchange_asymmetric.cc
  test_qc.cc test_MolecularOrbitals.cc test_BSHApply.cc test_projector.cc
  test_results_summary.cc
  test_exchange_owner_assign.cc test_exchange_batch_cache.cc test_restart.cc
  test_population.cc)
  set(CHEM_TEST_SOURCES_LONG test_localizer.cc test_ccpairfunction.cc test_low_rank_function.cc test_scf_evaluation.cc test_scf_iteration.cc test_scf_defaults.cc test_scf_frozen_iteration.cc test_scf_checkpoint.cc)
  if (TARGET Libxc::xc)
    list(APPEND CHEM_TEST_SOURCES_SHORT test_dft.cc )
    list(APPEND CHEM_TEST_SOURCES_LONG test_SCFOperators.cc)
  endif(TARGET Libxc::xc)
  if (INTEGRATORXX_FOUND)
    list(APPEND CHEM_TEST_SOURCES_SHORT test_IntegratorXX.cc )
  endif(INTEGRATORXX_FOUND)
  if (MADNESS_HAS_DFTD3)
    list(APPEND CHEM_TEST_SOURCES_SHORT test_dispersion.cc )
  endif(MADNESS_HAS_DFTD3)

  add_unittests(chem "${CHEM_TEST_SOURCES_SHORT}" "MADchem;MADgtest" "unittests;short")
  add_unittests(chem "test_density_changes.cc" "MADchem" "unittests;short")
  if(ENABLE_MPI)
    add_mpi_tests(chem test_density_changes "2" "MADchem" "unittests;short")
  endif()
  # The planner, plan_ladder (DESIGN B6, D7): pure, a table of cases.
  add_unittests(chem "test_plan_ladder.cc" "MADworld;nlohmann_json::nlohmann_json" "unittests;short")
  # Who decides convergence and the settings table (DESIGN B8, B9): pure.
  add_unittests(chem "test_convergence_settings.cc" "MADworld;nlohmann_json::nlohmann_json" "unittests;short")
  # Measurement-based shared policy; pure scalar fixtures, no numerical solve.
  add_unittests(chem "test_iteration_metrics.cc" "MADworld;nlohmann_json::nlohmann_json" "unittests;short")
  # Convergence decision module (step 1): header-only, plain values, no runtime.
  add_unittests(chem "test_convergence_decision.cc" "MADworld" "unittests;short")
  # Per-update history row (step 2): header-only, plain values, no runtime.
  add_unittests(chem "test_iteration_history.cc" "MADworld;nlohmann_json::nlohmann_json" "unittests;short")
  # Timing record of one launch, verbosity, printing rules (step 3): header-only, plain values.
  add_unittests(chem "test_iteration_output.cc" "MADworld;nlohmann_json::nlohmann_json" "unittests;short")
  # Node archive (step 4a): header-only, plain values and a temporary folder.
  add_unittests(chem "test_node_archive.cc" "MADworld;nlohmann_json::nlohmann_json" "unittests;short")
  add_unittests(chem "test_metric_helpers.cc" "MADworld;nlohmann_json::nlohmann_json" "unittests;short")
  add_unittests(chem "test_optimizer_settings.cc" "MADchem" "unittests;short")
  add_unittests(chem "test_timing_records.cc" "MADchem" "unittests;short")
  add_unittests(chem "test_scf_settings.cc" "MADchem" "unittests;short")
add_unittests(chem "test_scf_inputs.cc" "MADchem" "unittests;short")
add_unittests(chem "test_scf_run.cc" "MADchem" "unittests;long")
  add_unittests(chem "test_iteration_restart.cc" "MADworld;nlohmann_json::nlohmann_json" "unittests;short")
  add_unittests(chem "test_scf_evaluation_policy.cc" "MADworld;nlohmann_json::nlohmann_json" "unittests;short")
  # the rotating operand only moves between ranks, so the distributed path needs more than one
  add_mpi_tests(chem test_exchange_asymmetric "2" "MADchem;MADgtest" "unittests;short")
  add_unittests(chem "${CHEM_TEST_SOURCES_LONG}" "MADchem;MADgtest" "unittests;long")
  set_tests_properties(madness/test/chem/test_scf_evaluation/run PROPERTIES TIMEOUT 600)
  set_tests_properties(madness/test/chem/test_scf_iteration/run PROPERTIES TIMEOUT 600)
  set_tests_properties(madness/test/chem/test_scf_defaults/run PROPERTIES TIMEOUT 600)
  set_tests_properties(madness/test/chem/test_scf_frozen_iteration/run PROPERTIES TIMEOUT 600)
  set_tests_properties(madness/test/chem/test_scf_checkpoint/run PROPERTIES
    LABELS "unittests;medium" TIMEOUT 600)
  foreach(regression IN ITEMS frozen empty occupations lifetime)
    add_test(NAME madness/test/chem/test_scf_frozen_${regression}/run
      COMMAND $<TARGET_FILE:test_scf_frozen_iteration> --review ${regression})
    set_tests_properties(madness/test/chem/test_scf_frozen_${regression}/run PROPERTIES
      LABELS "unittests;short" TIMEOUT 120)
  endforeach()
  # Curated scientific controls use unchanged QC decks/references. Coarse Li is
  # intentionally incomplete; refining the representation must meet the request.
  foreach(variant IN ITEMS hybrid water coarse refined kinetic invalid_stage)
    if(NOT TARGET Libxc::xc AND variant MATCHES "^(hybrid|coarse|refined)$")
      continue()
    endif()
    if(variant STREQUAL "hybrid")
      set(default_args --case=scf_he_pbe0)
    elseif(variant STREQUAL "water")
      set(default_args --case=scf_h2o_hf)
    elseif(variant STREQUAL "coarse")
      set(default_args --case=scf_li_tpss --expect-limit=true)
    elseif(variant STREQUAL "refined")
      set(default_args --case=scf_li_tpss --protocol=1e-6)
    elseif(variant STREQUAL "kinetic")
      set(default_args --case=scf_he_hf --kinetic-required=true)
    else()
      set(default_args --case=scf_he_hf --invalid-stage=true)
    endif()
    add_test(NAME madness/test/chem/test_scf_defaults_${variant}/run
      COMMAND $<TARGET_FILE:test_scf_defaults> ${default_args})
    set_tests_properties(madness/test/chem/test_scf_defaults_${variant}/run PROPERTIES
      LABELS "unittests;verylong" TIMEOUT 600)
  endforeach()
  if(ENABLE_MPI)
    add_mpi_tests(chem test_scf_evaluation "2" "MADchem;MADgtest" "unittests;verylong")
    add_mpi_tests(chem test_scf_iteration "2" "MADchem;MADgtest" "unittests;verylong")
    add_mpi_tests(chem test_scf_run "2" "MADchem" "unittests;long")
    set_tests_properties(madness/test/chem/test_scf_run_mpi2/run PROPERTIES
      TIMEOUT 600 ENVIRONMENT "MAD_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;OMPI_MCA_hwloc_base_binding_policy=none")
    add_mpi_tests(chem test_scf_checkpoint "2" "MADchem;MADgtest" "unittests;long")
    set_tests_properties(madness/test/chem/test_scf_checkpoint_mpi2/run PROPERTIES
      TIMEOUT 600 ENVIRONMENT "MAD_NUM_THREADS=2;OMPI_MCA_hwloc_base_binding_policy=none")
    add_mpi_tests(chem test_scf_frozen_iteration "2" "MADchem;MADgtest" "unittests;long")
    set_tests_properties(madness/test/chem/test_scf_frozen_iteration_mpi2/run PROPERTIES
      TIMEOUT 600 ENVIRONMENT "MAD_NUM_THREADS=2;OMPI_MCA_hwloc_base_binding_policy=none")
    set_tests_properties(madness/test/chem/test_scf_iteration_mpi2/run PROPERTIES
      TIMEOUT 600 ENVIRONMENT "MAD_NUM_THREADS=2;OMPI_MCA_hwloc_base_binding_policy=none")
    set_tests_properties(madness/test/chem/test_scf_evaluation_mpi2/run PROPERTIES
      TIMEOUT 600 ENVIRONMENT "MAD_NUM_THREADS=2;OMPI_MCA_hwloc_base_binding_policy=none")
  endif()

  # Create other test executables not included in the unit tests ... consider these executables (unlike unit tests)
  if (NOT MADNESS_BUILD_LIBRARIES_ONLY)
    set(CHEM_OTHER_TESTS testxc benchmark_exchange_operator)
    foreach(_test ${CHEM_OTHER_TESTS})
      add_mad_executable(${_test} "${_test}.cc" "MADchem")
    endforeach()
  endif()

endif()
