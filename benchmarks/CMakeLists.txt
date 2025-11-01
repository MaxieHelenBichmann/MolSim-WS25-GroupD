if(ENABLE_BENCHMARK)
    file(GLOB_RECURSE BENCHMARK_SOURCES "src/*.cpp")

    add_executable(MolSimBench ${BENCHMARK_SOURCES})

    target_link_libraries(MolSimBench PRIVATE
    LibMolSim
    benchmark::benchmark
    benchmark::benchmark_main
    )   
    target_compile_definitions(MolSimBench PRIVATE DISABLE_IO)
    message(STATUS "Benchmarking enabled and executable MolSimBench available")
endif()
