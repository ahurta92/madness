# src/madness/misc

set(MADMISC_HEADERS misc.h ran.h cfft.h info.h gnuplot.h array_of_bools.h kahan_accumulator.h)
set(MADMISC_SOURCES
    checksum_file.cc position_stream.cc gprofexit.cc ran.cc cfft.cc info.cc numerics.cc unique_filename.cc)
# The git revision is written at build time, not configure time, so every
# rebuild names the commit that was built. Configure-time values are only the
# fallback for a tree without git. Only info.cc includes the generated header.
include(GetGitMetadata)
vgkit_cmake_git_metadata()
find_package(Git QUIET)
set(MADNESS_GIT_REVISION_HEADER ${CMAKE_CURRENT_BINARY_DIR}/madness_git_revision.h)
add_custom_target(madness-git-revision
        COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${PROJECT_SOURCE_DIR} -DOUTPUT=${MADNESS_GIT_REVISION_HEADER}
                "-DGIT_EXECUTABLE=${GIT_EXECUTABLE}" "-DFALLBACK_REVISION=${MADNESS_GIT_REVISION}"
                "-DFALLBACK_DESCRIPTION=${MADNESS_GIT_DESCRIPTION}"
                -P ${PROJECT_SOURCE_DIR}/cmake/modules/WriteGitRevision.cmake
        BYPRODUCTS ${MADNESS_GIT_REVISION_HEADER}
        COMMENT "Recording the git revision"
        VERBATIM)
set_source_files_properties(info.cc PROPERTIES OBJECT_DEPENDS ${MADNESS_GIT_REVISION_HEADER})

# Create the MADmisc library.
# MADmisc is a leaf: it must not include anything from world/tensor/mra/chem, so that
# MADworld can depend on it (see world/CMakeLists.txt).
add_mad_library(misc MADMISC_SOURCES MADMISC_HEADERS "" "madness/misc/")
target_link_libraries(MADmisc PUBLIC Threads::Threads)
add_dependencies(MADmisc-obj madness-git-revision)
target_include_directories(MADmisc-obj PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# No unit tests registered here: interp3.cc moved to src/madness/mra along with
# interpolation_1d.h, and test_gnuplot.cc stays disabled (gnuplot missing in CI).
