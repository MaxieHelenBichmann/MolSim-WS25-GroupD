set(OUTPUT_DIR     "${PROJECT_BINARY_DIR}/coverage")
set(PROFDATA       "${OUTPUT_DIR}/coverage.profdata")
file(GLOB_RECURSE PROFRAW_FILES "${PROJECT_BINARY_DIR}/*.profraw")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${PROJECT_BINARY_DIR}/*.profraw")

if(DEFINED CMAKE_SCRIPT_MODE_FILE)

  file(MAKE_DIRECTORY "${OUTPUT_DIR}")

  execute_process(
    COMMAND ${LLVM_PROFDATA_EXECUTABLE} merge -sparse ${PROFRAW_FILES} -o ${PROFDATA}
    RESULT_VARIABLE PROF_MERGE_RESULT
  )
  if(NOT PROF_MERGE_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to merge profraw files")
  endif()

  
  # Generate coverage report
  set(_cov_export_args export --format=lcov --instr-profile=${PROFDATA})

  # Tell llvm-cov how to rewrite any lingering absolute paths to `.`
  if(DEFINED PROJECT_SOURCE_DIR)
  list(APPEND _cov_export_args
    "--path-equivalence=${PROJECT_SOURCE_DIR},."
    "--compilation-dir=${PROJECT_SOURCE_DIR}"
  )
  endif()

  # optionally ignore external headers / build dirs
  if(NOT DEFINED COVERAGE_IGNORE_REGEX)
    set(COVERAGE_IGNORE_REGEX ".*/(tests?|build|_deps|external)/.*")
  endif()
  list(APPEND _cov_export_args "--ignore-filename-regex=${COVERAGE_IGNORE_REGEX}")

  list(APPEND _cov_export_args ${COVERAGE_TARGET_FILES})

  execute_process(
    COMMAND ${LLVM_COV_EXECUTABLE} ${_cov_export_args}
    OUTPUT_FILE "${OUTPUT_DIR}/coverage.lcov"
    ERROR_VARIABLE COV_EXPORT_ERROR
    RESULT_VARIABLE COV_EXPORT_RESULT
  )
  if(NOT COV_EXPORT_RESULT EQUAL 0)
   message(FATAL_ERROR "Failed to export coverage data: ${COV_EXPORT_ERROR}")
  endif()

  if(DEFINED PROJECT_SOURCE_DIR)
    set(_lcov_file "${OUTPUT_DIR}/coverage.lcov")
    if(EXISTS "${_lcov_file}")
      file(READ "${_lcov_file}" _lcov_contents)
      set(_needle "SF:${PROJECT_SOURCE_DIR}")
      string(FIND "${_lcov_contents}" "${_needle}" _idx)
      if(NOT _idx EQUAL -1)
        string(REPLACE "SF:${PROJECT_SOURCE_DIR}/" "SF:" _lcov_contents "${_lcov_contents}")
        string(REPLACE "SF:${PROJECT_SOURCE_DIR}" "SF:" _lcov_contents "${_lcov_contents}")
        file(WRITE "${_lcov_file}" "${_lcov_contents}")
      endif()
    endif()
  endif()
  return()
endif()

include_guard(GLOBAL)

function(enable_coverage)
  set(_targets ${ARGV})

  foreach(t IN LISTS _targets)
    target_compile_options(${t} PRIVATE
      -fprofile-instr-generate
      -fcoverage-mapping
      -fcoverage-prefix-map=${PROJECT_SOURCE_DIR}=.
      -fdebug-prefix-map=${PROJECT_SOURCE_DIR}=.
    )
    target_link_options(${t} PRIVATE
      -fprofile-instr-generate
      -fcoverage-mapping
    )
  endforeach()

  if(_targets)
    set_property(GLOBAL APPEND PROPERTY COVERAGE_ENABLED_TARGETS "${_targets}")
  endif()
endfunction()

# Create a standalone target that converts collected raw profiles into an LCOV file.
# Usage: add_coverage_report_target(<targets...>)
# - This defines a custom target named `coverage_report` with NO DEPENDS.
# - Invoke after running your instrumented binaries/tests to generate coverage.
function(add_coverage_report_target)
  # Resolve required tools when creating the conversion target
  find_program(LLVM_PROFDATA_EXECUTABLE NAMES llvm-profdata REQUIRED)
  find_program(LLVM_COV_EXECUTABLE     NAMES llvm-cov       REQUIRED)

  set(_targets ${ARGV})
  if(NOT _targets)
    # Try to pick up any targets that previously had enable_coverage() applied
    get_property(_targets GLOBAL PROPERTY COVERAGE_ENABLED_TARGETS)
  endif()

  # Build generator-expression paths for provided targets
  set(_files "")
  foreach(t IN LISTS _targets)
    list(APPEND _files "$<TARGET_FILE:${t}>")
  endforeach()

  # Deduplicate
  if(_files)
    list(REMOVE_DUPLICATES _files)
  endif()

  # Escape list for -D passing
  string(REPLACE ";" "\\;" _files_escaped "${_files}")

  if(_files_escaped STREQUAL "")
   message(FATAL_ERROR "No target files specified for coverage report generation")
  endif()

  add_custom_target(coverage_report
    COMMAND ${CMAKE_COMMAND}
      -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
      -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
      -DLLVM_PROFDATA_EXECUTABLE=${LLVM_PROFDATA_EXECUTABLE}
      -DLLVM_COV_EXECUTABLE=${LLVM_COV_EXECUTABLE}
      -DCOVERAGE_TARGET_FILES=${_files_escaped}
      -P "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/GenerateCoverage.cmake"
    COMMENT "Generating LLVM coverage"
    VERBATIM
  )
  message(STATUS "Added coverage_report target to generate coverage data for targets: ${_targets}")
endfunction()