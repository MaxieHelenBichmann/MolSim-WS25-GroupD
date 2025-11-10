cmake_minimum_required(VERSION 3.10)

# Enable C++20. You may raise the standard if you want
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

# define project name, version
project(MolSim_GroupD VERSION 0.0.1 LANGUAGES CXX)

include_directories(${CMAKE_SOURCE_DIR}/include)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-O0)
else()
    add_compile_options(-O3)
endif()

add_compile_definitions(
    SPDLOG_USE_STD_FORMAT
    $<$<CONFIG:Debug>:SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE>
    $<$<CONFIG:RelWithDebInfo>:SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE>
    $<$<CONFIG:Release>:SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_WARN>
    $<$<CONFIG:MinSizeRel>:SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_WARN>
)


# add additional modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(CTest)
include(doxygen)
include(gtest)
include(clang-format)
include(yaml-cpp)
include(spdlog)
include(benchmark)
include(clang-tidy)

option(COVERAGE "Enable generation of LLVM coverage artifacts" OFF)


add_subdirectory(src project)
add_subdirectory(tests)
add_subdirectory(benchmarks)

if(COVERAGE)
    include(GenerateCoverage)
    add_coverage_report_target()
endif()
