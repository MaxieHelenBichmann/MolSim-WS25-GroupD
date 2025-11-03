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
  execute_process(
    COMMAND ${LLVM_COV_EXECUTABLE} export --format=lcov --instr-profile=${PROFDATA} ${COVERAGE_TARGET_FILES}
    OUTPUT_FILE "${OUTPUT_DIR}/coverage.lcov"
    ERROR_VARIABLE COV_EXPORT_ERROR
    RESULT_VARIABLE COV_EXPORT_RESULT
  )
  if(NOT COV_EXPORT_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to export coverage data: ${COV_EXPORT_ERROR}")
  endif()
  return()
endif()

include_guard(GLOBAL)

function(enable_coverage)
  set(_targets ${ARGV})

  foreach(t IN LISTS _targets)
    target_compile_options(${t} PRIVATE -fprofile-instr-generate -fcoverage-mapping)
    target_link_options(${t} PRIVATE -fprofile-instr-generate -fcoverage-mapping)
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
      -DLLVM_PROFDATA_EXECUTABLE=${LLVM_PROFDATA_EXECUTABLE}
      -DLLVM_COV_EXECUTABLE=${LLVM_COV_EXECUTABLE}
      -DCOVERAGE_TARGET_FILES=${_files_escaped}
      -P "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/GenerateCoverage.cmake"
    COMMENT "Generating LLVM coverage"
    VERBATIM
  )
endfunction()