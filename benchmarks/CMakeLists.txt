if(ENABLE_BENCHMARK)
    set(BENCHMARK_EXECUTABLE_NAME MolSimBench)
    file(GLOB_RECURSE BENCHMARK_SOURCES "src/*.cpp")

    add_executable(${BENCHMARK_EXECUTABLE_NAME} ${BENCHMARK_SOURCES})

    target_link_libraries(${BENCHMARK_EXECUTABLE_NAME} PRIVATE LibMolSim benchmark::benchmark)
    message(STATUS "Benchmarking enabled and executable MolSimBench available")
endif()
