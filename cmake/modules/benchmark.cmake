option(ENABLE_BENCHMARK "Enable benchmarking" ON)
if(ENABLE_BENCHMARK)
  include(FetchContent)
  message(STATUS "Benchmarking enabled")
  FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.8.0
    SYSTEM
    )
  FetchContent_MakeAvailable(benchmark)
endif()
