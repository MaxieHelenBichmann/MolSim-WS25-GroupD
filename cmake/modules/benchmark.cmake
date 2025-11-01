option(ENABLE_BENCHMARK "Enable benchmarking" OFF)
if(ENABLE_BENCHMARK)
  include(FetchContent)

  FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.8.0)

  FetchContent_MakeAvailable(benchmark)
  message(STATUS "Benchmarking enabled")
endif()
