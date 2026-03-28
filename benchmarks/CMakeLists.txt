if(ENABLE_BENCHMARK)
    file(GLOB_RECURSE BENCHMARK_SOURCES "src/*.cpp" "code/*.cpp")

    add_executable(MolSimBench ${BENCHMARK_SOURCES})

    target_link_libraries(MolSimBench PRIVATE
    LibMolSim
    benchmark::benchmark
    benchmark::benchmark_main
    )   
    target_compile_options(MolSimBench PRIVATE -Wextra -Wall)
    add_clang_tidy_support(MolSimBench)
    message(STATUS "Benchmarking enabled and executable MolSimBench available")
endif()
