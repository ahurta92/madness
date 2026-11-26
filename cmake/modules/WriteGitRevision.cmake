#
# Script mode, run at every build by src/madness/misc/CMakeLists.txt:
#
#   cmake -DSOURCE_DIR=<repo> -DOUTPUT=<header> -DGIT_EXECUTABLE=<git>
#         -DFALLBACK_REVISION=<rev> -DFALLBACK_DESCRIPTION=<desc> -P WriteGitRevision.cmake
#
# Writes <header> defining MADNESS_GIT_REVISION and MADNESS_GIT_DESCRIPTION from
# the source tree as it is now, so a rebuild names the commit that was built.
# MADNESS_GIT_REVISION is the full commit hash, with "-dirty" appended when
# tracked files have uncommitted changes.
# The file is rewritten only when its text changes, so an unchanged revision
# recompiles nothing. Without git (a source tarball) the configure-time values
# are used.
#
set(revision "${FALLBACK_REVISION}")
set(description "${FALLBACK_DESCRIPTION}")
if(GIT_EXECUTABLE AND EXISTS "${SOURCE_DIR}/.git")
  execute_process(COMMAND "${GIT_EXECUTABLE}" rev-parse -q HEAD
                  WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE status
                  OUTPUT_VARIABLE output ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(status EQUAL 0)
    set(revision "${output}")
    # Uncommitted changes to tracked files (as git describe --dirty counts
    # them; untracked files do not count) mark the build: <commit>-dirty.
    execute_process(COMMAND "${GIT_EXECUTABLE}" status --porcelain --untracked-files=no
                    WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE status
                    OUTPUT_VARIABLE output ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(status EQUAL 0 AND NOT output STREQUAL "")
      string(APPEND revision "-dirty")
    endif()
  endif()
  execute_process(COMMAND "${GIT_EXECUTABLE}" describe --dirty
                  WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE status
                  OUTPUT_VARIABLE output ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(status EQUAL 0)
    set(description "${output}")
  endif()
endif()
if(NOT revision)
  set(revision "unknown")
endif()
if(NOT description)
  set(description "unknown")
endif()
set(text "// Generated at build time by cmake/modules/WriteGitRevision.cmake; do not edit.
#define MADNESS_GIT_REVISION \"${revision}\"
#define MADNESS_GIT_DESCRIPTION \"${description}\"
")
set(previous "")
if(EXISTS "${OUTPUT}")
  file(READ "${OUTPUT}" previous)
endif()
if(NOT previous STREQUAL text)
  file(WRITE "${OUTPUT}" "${text}")
endif()
