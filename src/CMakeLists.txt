
file(GLOB_RECURSE SRCS 
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/*.h"
    )
list(FILTER SRCS EXCLUDE REGEX ".*\/MolSim\.cpp")


add_library(LibMolSim STATIC ${SRCS})
target_link_libraries(LibMolSim PUBLIC yaml-cpp::yaml-cpp spdlog::spdlog CLI11::CLI11)
target_compile_options(LibMolSim PRIVATE -Wextra -Wall -Werror)
add_clang_tidy_support(LibMolSim)
include(vtk)

add_executable(MolSim MolSim.cpp)
include(openmp)

target_link_libraries(MolSim PRIVATE LibMolSim)
add_clang_tidy_support(MolSim)


target_compile_options(MolSim
        PRIVATE
            $<$<CXX_COMPILER_ID:GNU>:-Wall -Werror>
            # disable some spam warnings for icpc...
            $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)

include(profiling)

if(ENABLE_IO)
    target_compile_definitions(MolSim PRIVATE ENABLE_IO)
endif()
if(ENABLE_CHECKPOINTS)
    target_compile_definitions(MolSim PRIVATE ENABLE_CHECKPOINTING)
endif()
if(ENABLE_STATS)
    target_compile_definitions(MolSim PRIVATE ENABLE_STATS)
endif()

if(COVERAGE)
    include(GenerateCoverage)
    enable_coverage(MolSim)
endif()


if(ENABLE_PROFILING)
    add_perf_targets(MolSim)
    add_valgrind_targets(MolSim)
endif()

if(ENABLE_INTEL_VTUNE)
    add_vtune_targets(MolSim)
endif()

if(ENABLE_INTEL_ADVISOR)
    add_advisor_targets(MolSim)
endif()

if(ENABLE_INTEL_VTUNE OR ENABLE_INTEL_ADVISOR)
    add_intel_batch_targets(MolSim)
endif()
