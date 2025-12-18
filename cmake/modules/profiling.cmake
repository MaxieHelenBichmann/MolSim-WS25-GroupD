include_guard(GLOBAL)

option(ENABLE_PROFILING "Enable perf and valgrind profiling targets" OFF)


# Check for profiling tools
find_program(PERF_EXECUTABLE NAMES perf)
find_program(VALGRIND_EXECUTABLE NAMES valgrind)

if(PERF_EXECUTABLE)
    message(STATUS "Found perf: ${PERF_EXECUTABLE}")
else()
    message(STATUS "perf not found - profiling targets will not be available")
endif()

if(VALGRIND_EXECUTABLE)
    message(STATUS "Found valgrind: ${VALGRIND_EXECUTABLE}")
else()
    message(STATUS "valgrind not found - memory check targets will not be available")
endif()

# Function to add perf profiling targets
function(add_perf_targets TARGET_NAME)
    if(NOT PERF_EXECUTABLE)
        return()
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING " Performance profiling (perf) should be used with RelWithDebInfo")
    endif()
    set(PERF_INPUT_FILE "${PROJECT_SOURCE_DIR}/input/particles.yaml" CACHE STRING "Input file for perf profiling")
    set(PERF_ARGS "-d 0.0002 -t 5" CACHE STRING "Arguments for perf profiling runs")
    
    add_custom_target(perf-record
        COMMAND ${PERF_EXECUTABLE} record --call-graph=dwarf -F 997
                $<TARGET_FILE:${TARGET_NAME}> ${PERF_INPUT_FILE} ${PERF_ARGS}
        COMMAND ${CMAKE_COMMAND} -E echo "Performance data saved to perf.data"
        COMMAND ${CMAKE_COMMAND} -E echo "Run 'make perf-report' to analyze results"
        DEPENDS ${TARGET_NAME}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Recording performance profile with perf"
        VERBATIM
    )

    add_custom_target(perf-report
        COMMAND ${PERF_EXECUTABLE} report
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Analyzing performance profile (interactive)"
        VERBATIM
    )

    add_custom_target(perf-report-text
        COMMAND ${PERF_EXECUTABLE} report --stdio > ${PROJECT_BINARY_DIR}/perf-report.txt
        COMMAND ${CMAKE_COMMAND} -E echo "Report saved to ${PROJECT_BINARY_DIR}/perf-report.txt"
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Generating text performance report"
        VERBATIM
    )

    add_custom_target(perf-stat
        COMMAND ${PERF_EXECUTABLE} stat -d
                $<TARGET_FILE:${TARGET_NAME}> ${PERF_INPUT_FILE} ${PERF_ARGS}
        DEPENDS ${TARGET_NAME}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Running with perf stat (hardware counters)"
        VERBATIM
    )

    add_custom_target(perf-cache
        COMMAND ${PERF_EXECUTABLE} stat -e cache-references,cache-misses,L1-dcache-loads,L1-dcache-load-misses
                $<TARGET_FILE:${TARGET_NAME}> ${PERF_INPUT_FILE} ${PERF_ARGS}
        DEPENDS ${TARGET_NAME}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Analyzing cache performance"
        VERBATIM
    )

    message(STATUS "Added perf targets: perf-record, perf-report, perf-report-text, perf-stat, perf-cache")
endfunction()

function(add_valgrind_targets TARGET_NAME)
    if(NOT VALGRIND_EXECUTABLE)
        return()
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    message(WARNING "Memory Debugging depends on -g symbols so only use with Debug or RelWithDebInfo")
    endif()
    set(VALGRIND_INPUT_FILE "${PROJECT_SOURCE_DIR}/input/particles.yaml" CACHE STRING "Input file for valgrind analysis")
    set(VALGRIND_ARGS "-d 0.0002 -t 1" CACHE STRING "Arguments for valgrind runs (use short simulation)")
    set(VALGRIND_OUTPUT_DIR "${PROJECT_BINARY_DIR}/valgrind" CACHE STRING "Output directory for valgrind reports")

    add_custom_target(valgrind-memcheck
        COMMAND ${CMAKE_COMMAND} -E make_directory ${VALGRIND_OUTPUT_DIR}
        COMMAND ${VALGRIND_EXECUTABLE} 
                --leak-check=full
                --show-leak-kinds=all
                --track-origins=yes
                --verbose
                --log-file=${VALGRIND_OUTPUT_DIR}/memcheck.log
                $<TARGET_FILE:${TARGET_NAME}> ${VALGRIND_INPUT_FILE} ${VALGRIND_ARGS}
        COMMAND ${CMAKE_COMMAND} -E echo "Memory check complete. Report: ${VALGRIND_OUTPUT_DIR}/memcheck.log"
        DEPENDS ${TARGET_NAME}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Running valgrind memory check"
        VERBATIM
    )

    add_custom_target(valgrind-massif
        COMMAND ${CMAKE_COMMAND} -E make_directory ${VALGRIND_OUTPUT_DIR}
        COMMAND ${VALGRIND_EXECUTABLE}
                --tool=massif
                --massif-out-file=${VALGRIND_OUTPUT_DIR}/massif.out
                $<TARGET_FILE:${TARGET_NAME}> ${VALGRIND_INPUT_FILE} ${VALGRIND_ARGS}
        COMMAND ${CMAKE_COMMAND} -E echo "Heap profiling complete. Use 'ms_print ${VALGRIND_OUTPUT_DIR}/massif.out' to analyze"
        DEPENDS ${TARGET_NAME}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Running massif heap profiling"
        VERBATIM
    )

    message(STATUS "Added valgrind targets: valgrind-memcheck, valgrind-massif")
endfunction()
