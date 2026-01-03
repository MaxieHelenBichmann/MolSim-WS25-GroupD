include_guard(GLOBAL)

option(ENABLE_PROFILING "Enable perf and valgrind profiling targets" OFF)
option(ENABLE_INTEL_VTUNE "Enable Intel VTune profiling targets" OFF)
option(ENABLE_INTEL_ADVISOR "Enable Intel Advisor profiling targets" OFF)


# Check for profiling tools
find_program(PERF_EXECUTABLE NAMES perf)
find_program(VALGRIND_EXECUTABLE NAMES valgrind)
find_program(VTUNE_EXECUTABLE NAMES vtune)
find_program(ADVISOR_EXECUTABLE NAMES advisor advisor-cl)

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

if(VTUNE_EXECUTABLE)
    message(STATUS "Found VTune: ${VTUNE_EXECUTABLE}")
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        message(STATUS "VTune works best with Intel compilers (icpx/icx). Current: ${CMAKE_CXX_COMPILER_ID}")
    endif()
else()
    message(STATUS "VTune not found - Intel profiling targets will not be available")
endif()

if(ADVISOR_EXECUTABLE)
    message(STATUS "Found Advisor: ${ADVISOR_EXECUTABLE}")
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        message(STATUS "Advisor works best with Intel compilers (icpx/icx). Current: ${CMAKE_CXX_COMPILER_ID}")
    endif()
else()
    message(STATUS "Advisor not found - Intel analysis targets will not be available")
endif()

# Intel compiler optimization flags for profiling
if(ENABLE_INTEL_VTUNE OR ENABLE_INTEL_ADVISOR)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        # Intel-specific flags for better profiling
        # -g: debug symbols, -O2: optimization level
        # Note: -qopt-report is only useful during linking, omitted to avoid warnings
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -g -O2")
        message(STATUS "Added Intel compiler profiling flags (-g -O2)")
    endif()
endif()

# Function to add perf profiling targets
function(add_perf_targets TARGET_NAME)
    if(NOT PERF_EXECUTABLE)
        return()
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING " Performance profiling (perf) should be used with RelWithDebInfo")
    endif()
    set(PERF_INPUT_FILE "${PROJECT_SOURCE_DIR}/input/benchmark_10k.yaml" CACHE STRING "Input file for perf profiling")
    set(PERF_ARGS "" CACHE STRING "Arguments for perf profiling runs")
    
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
    set(VALGRIND_INPUT_FILE "${PROJECT_SOURCE_DIR}/input/benchmark_10k.yaml" CACHE STRING "Input file for valgrind analysis")
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

# Function to add Intel VTune profiling targets
function(add_vtune_targets TARGET_NAME)
    if(NOT VTUNE_EXECUTABLE)
        return()
    endif()
    
    set(VTUNE_INPUT_FILE "${PROJECT_SOURCE_DIR}/input/benchmark_10k.yaml" CACHE STRING "Input file for VTune profiling")
    set(VTUNE_ARGS "" CACHE STRING "Arguments for VTune profiling runs")
    set(VTUNE_OUTPUT_DIR "${PROJECT_BINARY_DIR}/vtune" CACHE STRING "Output directory for VTune results")
    
    add_custom_target(vtune-hotspots
        COMMAND ${VTUNE_EXECUTABLE} -collect hotspots -result-dir ${VTUNE_OUTPUT_DIR}/hotspots
                -- $<TARGET_FILE:${TARGET_NAME}> ${VTUNE_INPUT_FILE} ${VTUNE_ARGS}
        DEPENDS ${TARGET_NAME}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Running VTune hotspots analysis"
        VERBATIM
    )
    
    add_custom_target(vtune-memory
        COMMAND ${VTUNE_EXECUTABLE} -collect memory-access -result-dir ${VTUNE_OUTPUT_DIR}/memory
                -- $<TARGET_FILE:${TARGET_NAME}> ${VTUNE_INPUT_FILE} ${VTUNE_ARGS}
        DEPENDS ${TARGET_NAME}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Running VTune memory access analysis"
        VERBATIM
    )
    
    add_custom_target(vtune-uarch
        COMMAND ${VTUNE_EXECUTABLE} -collect uarch-exploration -result-dir ${VTUNE_OUTPUT_DIR}/uarch
                -- $<TARGET_FILE:${TARGET_NAME}> ${VTUNE_INPUT_FILE} ${VTUNE_ARGS}
        DEPENDS ${TARGET_NAME}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Running VTune microarchitecture exploration"
        VERBATIM
    )
    
    message(STATUS "Added VTune targets: vtune-hotspots, vtune-memory, vtune-uarch")
endfunction()

# Function to add Intel Advisor profiling targets
function(add_advisor_targets TARGET_NAME)
    if(NOT ADVISOR_EXECUTABLE)
        return()
    endif()
    
    set(ADVISOR_INPUT_FILE "${PROJECT_SOURCE_DIR}/input/benchmark_10k.yaml" CACHE STRING "Input file for Advisor analysis")
    set(ADVISOR_ARGS "" CACHE STRING "Arguments for Advisor runs")
    set(ADVISOR_OUTPUT_DIR "${PROJECT_BINARY_DIR}/advisor" CACHE STRING "Output directory for Advisor results")
    
    add_custom_target(advisor-survey
        COMMAND ${ADVISOR_EXECUTABLE} --collect=survey --project-dir=${ADVISOR_OUTPUT_DIR}/survey
                -- $<TARGET_FILE:${TARGET_NAME}> ${ADVISOR_INPUT_FILE} ${ADVISOR_ARGS}
        DEPENDS ${TARGET_NAME}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Running Advisor survey analysis"
        VERBATIM
    )
    
    add_custom_target(advisor-tripcounts
        COMMAND ${ADVISOR_EXECUTABLE} --collect=tripcounts --project-dir=${ADVISOR_OUTPUT_DIR}/survey
                -- $<TARGET_FILE:${TARGET_NAME}> ${ADVISOR_INPUT_FILE} ${ADVISOR_ARGS}
        DEPENDS ${TARGET_NAME}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Running Advisor trip counts analysis"
        VERBATIM
    )
    
    add_custom_target(advisor-roofline
        COMMAND ${ADVISOR_EXECUTABLE} --collect=roofline --project-dir=${ADVISOR_OUTPUT_DIR}/roofline
                -- $<TARGET_FILE:${TARGET_NAME}> ${ADVISOR_INPUT_FILE} ${ADVISOR_ARGS}
        DEPENDS ${TARGET_NAME}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Running Advisor roofline analysis"
        VERBATIM
    )
    
    message(STATUS "Added Advisor targets: advisor-survey, advisor-tripcounts, advisor-roofline")
endfunction()

# Function to add combined profiling targets
function(add_intel_batch_targets TARGET_NAME)
    if(NOT VTUNE_EXECUTABLE AND NOT ADVISOR_EXECUTABLE)
        return()
    endif()
    
    set(PROFILING_OUTPUT_DIR "${PROJECT_BINARY_DIR}/profiling-results" CACHE STRING "Output directory for all profiling results")
    
    if(VTUNE_EXECUTABLE)
        add_custom_target(vtune-all
            COMMAND ${CMAKE_COMMAND} -E make_directory ${PROFILING_OUTPUT_DIR}/vtune
            COMMAND ${VTUNE_EXECUTABLE} -collect hotspots -result-dir ${PROFILING_OUTPUT_DIR}/vtune/hotspots
                    -- $<TARGET_FILE:${TARGET_NAME}> ${VTUNE_INPUT_FILE} ${VTUNE_ARGS}
            COMMAND ${VTUNE_EXECUTABLE} -collect memory-access -result-dir ${PROFILING_OUTPUT_DIR}/vtune/memory
                    -- $<TARGET_FILE:${TARGET_NAME}> ${VTUNE_INPUT_FILE} ${VTUNE_ARGS}
            COMMAND ${VTUNE_EXECUTABLE} -collect uarch-exploration -result-dir ${PROFILING_OUTPUT_DIR}/vtune/uarch
                    -- $<TARGET_FILE:${TARGET_NAME}> ${VTUNE_INPUT_FILE} ${VTUNE_ARGS}
            DEPENDS ${TARGET_NAME}
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            COMMENT "Running all VTune analyses"
            VERBATIM
        )
    endif()
    
    if(ADVISOR_EXECUTABLE)
        add_custom_target(advisor-all
            COMMAND ${CMAKE_COMMAND} -E make_directory ${PROFILING_OUTPUT_DIR}/advisor
            COMMAND ${ADVISOR_EXECUTABLE} --collect=survey --project-dir=${PROFILING_OUTPUT_DIR}/advisor/survey
                    -- $<TARGET_FILE:${TARGET_NAME}> ${ADVISOR_INPUT_FILE} ${ADVISOR_ARGS}
            COMMAND ${ADVISOR_EXECUTABLE} --collect=tripcounts --project-dir=${PROFILING_OUTPUT_DIR}/advisor/survey
                    -- $<TARGET_FILE:${TARGET_NAME}> ${ADVISOR_INPUT_FILE} ${ADVISOR_ARGS}
            COMMAND ${ADVISOR_EXECUTABLE} --collect=roofline --project-dir=${PROFILING_OUTPUT_DIR}/advisor/roofline
                    -- $<TARGET_FILE:${TARGET_NAME}> ${ADVISOR_INPUT_FILE} ${ADVISOR_ARGS}
            DEPENDS ${TARGET_NAME}
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            COMMENT "Running all Advisor analyses"
            VERBATIM
        )
    endif()
    
    if(VTUNE_EXECUTABLE AND ADVISOR_EXECUTABLE)
        add_custom_target(intel-profile-all
            COMMAND ${CMAKE_COMMAND} -E make_directory ${PROFILING_OUTPUT_DIR}/vtune
            COMMAND ${VTUNE_EXECUTABLE} -collect hotspots -result-dir ${PROFILING_OUTPUT_DIR}/vtune/hotspots
                    -- $<TARGET_FILE:${TARGET_NAME}> ${VTUNE_INPUT_FILE} ${VTUNE_ARGS}
            COMMAND ${VTUNE_EXECUTABLE} -collect memory-access -result-dir ${PROFILING_OUTPUT_DIR}/vtune/memory
                    -- $<TARGET_FILE:${TARGET_NAME}> ${VTUNE_INPUT_FILE} ${VTUNE_ARGS}
            COMMAND ${VTUNE_EXECUTABLE} -collect uarch-exploration -result-dir ${PROFILING_OUTPUT_DIR}/vtune/uarch
                    -- $<TARGET_FILE:${TARGET_NAME}> ${VTUNE_INPUT_FILE} ${VTUNE_ARGS}
            COMMAND ${CMAKE_COMMAND} -E make_directory ${PROFILING_OUTPUT_DIR}/advisor
            COMMAND ${ADVISOR_EXECUTABLE} --collect=survey --project-dir=${PROFILING_OUTPUT_DIR}/advisor/survey
                    -- $<TARGET_FILE:${TARGET_NAME}> ${ADVISOR_INPUT_FILE} ${ADVISOR_ARGS}
            COMMAND ${ADVISOR_EXECUTABLE} --collect=tripcounts --project-dir=${PROFILING_OUTPUT_DIR}/advisor/survey
                    -- $<TARGET_FILE:${TARGET_NAME}> ${ADVISOR_INPUT_FILE} ${ADVISOR_ARGS}
            COMMAND ${ADVISOR_EXECUTABLE} --collect=roofline --project-dir=${PROFILING_OUTPUT_DIR}/advisor/roofline
                    -- $<TARGET_FILE:${TARGET_NAME}> ${ADVISOR_INPUT_FILE} ${ADVISOR_ARGS}
            COMMAND ${CMAKE_COMMAND} -E tar czf ${PROJECT_BINARY_DIR}/profiling-results.tar.gz profiling-results
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            COMMENT "Running Intel profiling and creating archive"
            VERBATIM
        )
        message(STATUS "Added combined target: intel-profile-all")
    elseif(VTUNE_EXECUTABLE)
        add_custom_target(intel-profile-all
            COMMAND ${CMAKE_COMMAND} -E make_directory ${PROFILING_OUTPUT_DIR}/vtune
            COMMAND ${VTUNE_EXECUTABLE} -collect hotspots -result-dir ${PROFILING_OUTPUT_DIR}/vtune/hotspots
                    -- $<TARGET_FILE:${TARGET_NAME}> ${VTUNE_INPUT_FILE} ${VTUNE_ARGS}
            COMMAND ${VTUNE_EXECUTABLE} -collect memory-access -result-dir ${PROFILING_OUTPUT_DIR}/vtune/memory
                    -- $<TARGET_FILE:${TARGET_NAME}> ${VTUNE_INPUT_FILE} ${VTUNE_ARGS}
            COMMAND ${VTUNE_EXECUTABLE} -collect uarch-exploration -result-dir ${PROFILING_OUTPUT_DIR}/vtune/uarch
                    -- $<TARGET_FILE:${TARGET_NAME}> ${VTUNE_INPUT_FILE} ${VTUNE_ARGS}
            COMMAND ${CMAKE_COMMAND} -E tar czf ${PROJECT_BINARY_DIR}/profiling-results.tar.gz profiling-results
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            COMMENT "Running VTune profiling and creating archive"
            VERBATIM
        )
        message(STATUS "Added combined target: intel-profile-all (VTune only)")
    elseif(ADVISOR_EXECUTABLE)
        add_custom_target(intel-profile-all
            COMMAND ${CMAKE_COMMAND} -E make_directory ${PROFILING_OUTPUT_DIR}/advisor
            COMMAND ${ADVISOR_EXECUTABLE} --collect=survey --project-dir=${PROFILING_OUTPUT_DIR}/advisor/survey
                    -- $<TARGET_FILE:${TARGET_NAME}> ${ADVISOR_INPUT_FILE} ${ADVISOR_ARGS}
            COMMAND ${ADVISOR_EXECUTABLE} --collect=tripcounts --project-dir=${PROFILING_OUTPUT_DIR}/advisor/survey
                    -- $<TARGET_FILE:${TARGET_NAME}> ${ADVISOR_INPUT_FILE} ${ADVISOR_ARGS}
            COMMAND ${ADVISOR_EXECUTABLE} --collect=roofline --project-dir=${PROFILING_OUTPUT_DIR}/advisor/roofline
                    -- $<TARGET_FILE:${TARGET_NAME}> ${ADVISOR_INPUT_FILE} ${ADVISOR_ARGS}
            COMMAND ${CMAKE_COMMAND} -E tar czf ${PROJECT_BINARY_DIR}/profiling-results.tar.gz profiling-results
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            COMMENT "Running Advisor profiling and creating archive"
            VERBATIM
        )
        message(STATUS "Added combined target: intel-profile-all (Advisor only)")
    endif()
endfunction()
