option(ENABLE_BENCHMARK "Enable benchmarking" ON)
if(ENABLE_BENCHMARK)
  include(FetchContent)
  message(STATUS "Benchmarking enabled")
  
  # Disable benchmark's internal googletest to avoid CMake compatibility issues
  set(BENCHMARK_DOWNLOAD_DEPENDENCIES OFF CACHE BOOL "Disable benchmark downloading googletest")
  set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "Disable benchmark gtests")
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable benchmark testing")
  
  FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.1
    SYSTEM
    )
  FetchContent_MakeAvailable(benchmark)
endif()
