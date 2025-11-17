
file(GLOB_RECURSE SRCS 
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/*.h"
    )
list(FILTER SRCS EXCLUDE REGEX ".*\/MolSim\.cpp")


add_library(LibMolSim STATIC ${SRCS})
target_link_libraries(LibMolSim PUBLIC yaml-cpp::yaml-cpp spdlog::spdlog)
target_compile_options(LibMolSim PRIVATE -Wextra -Wall -Werror)
include(vtk)

add_executable(MolSim MolSim.cpp)

target_link_libraries(MolSim PRIVATE LibMolSim)


target_compile_options(MolSim
        PRIVATE
            $<$<CXX_COMPILER_ID:GNU>:-Wall -Werror>
            # disable some spam warnings for icpc...
            $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)

include(profiling)


if(COVERAGE)
    include(GenerateCoverage)
    enable_coverage(MolSim)
endif()


if(ENABLE_PROFILING)
    add_perf_targets(MolSim)
    add_valgrind_targets(MolSim)
endif()
