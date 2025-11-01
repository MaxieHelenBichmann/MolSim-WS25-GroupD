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

add_compile_options(-Ofast -march=native)

# add additional modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(CTest)
include(doxygen)
include(gtest)
include(clang-format)
include(yaml-cpp)
include(clang-tidy)

add_subdirectory(src project)
add_subdirectory(tests)

