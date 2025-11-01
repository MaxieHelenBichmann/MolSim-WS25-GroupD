
file(GLOB_RECURSE SRCS 
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/*.h"
    )
list(FILTER SRCS EXCLUDE REGEX ".*\/MolSim\.cpp")


add_library(LibMolSim STATIC ${SRCS})
target_link_libraries(LibMolSim PUBLIC yaml-cpp::yaml-cpp spdlog::spdlog)
include(vtk)

add_executable(MolSim MolSim.cpp)

target_link_libraries(MolSim PRIVATE LibMolSim)


target_compile_options(MolSim
        PRIVATE
            $<$<CXX_COMPILER_ID:GNU>:-Wall>
            # disable some spam warnings for icpc...
            $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)
