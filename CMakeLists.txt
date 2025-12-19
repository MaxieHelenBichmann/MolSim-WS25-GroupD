cmake_minimum_required(VERSION 3.10)

option(COVERAGE "Enable generation of LLVM coverage artifacts" OFF)
option(ENABLE_UBSAN "Enable undefined behavior sanitizer (debug only)" OFF)

# Enable C++20. You may raise the standard if you want
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POLICY_VERSION_MINIMUM 3.10)

# define project name, version
project(MolSim_GroupD VERSION 0.0.1 LANGUAGES CXX)

include_directories(${CMAKE_SOURCE_DIR}/include)

# let ccmake and cmake-gui offer the default build type options (single-config generators only)
if(NOT CMAKE_CONFIGURATION_TYPES)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
    if(NOT CMAKE_BUILD_TYPE)
        set(default_build_type "Release")
        message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
        set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    endif()
endif()

# -march=native and -Ofast did not yield better results in benchmarks on Rechnerhalle, try again on Cluster
add_compile_options(
    $<$<CONFIG:Debug>:-O0>
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-flto>
    $<$<CONFIG:RelWithDebInfo>:-O3>
    $<$<CONFIG:MinSizeRel>:-O3>
    $<$<CONFIG:MinSizeRel>:-flto>
)

if(ENABLE_UBSAN)
    add_compile_options($<$<CONFIG:Debug>:-fsanitize=undefined>)
    add_link_options($<$<CONFIG:Debug>:-fsanitize=undefined>)
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
include(clang-tidy)
include(clang-format)
include(gtest)
include(yaml-cpp)
include(spdlog)
include(benchmark)
include(cli11)


option(COVERAGE "Enable generation of LLVM coverage artifacts" OFF)
option(ENABLE_IO "Enable file input/output operations for result files" ON)
option(ENABLE_CHECKPOINTS "Enable creating checkpoints" ON)

add_subdirectory(src project)
add_subdirectory(tests)
add_subdirectory(benchmarks)

if(COVERAGE)
    include(GenerateCoverage)
    add_coverage_report_target()
endif()
