message(STATUS "clang-tidy checks enabled")

find_program(CLANG_TIDY_EXE clang-tidy)
if(CLANG_TIDY_EXE)
    #For other options like automatic fixes add --fix at the end
    set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_EXE} -p ${CMAKE_BINARY_DIR})
    message(STATUS "clang-tidy found and enabled: ${CLANG_TIDY_EXE}")

    file(GLOB_RECURSE ALL_CXX_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_SOURCE_DIR}/src/*.cpp"
        "${CMAKE_SOURCE_DIR}/include/*.h"
        "${CMAKE_SOURCE_DIR}/tests/*.cpp"
        "${CMAKE_SOURCE_DIR}/benchmarks/*.cpp"
        "${CMAKE_SOURCE_DIR}/benchmarks/*.h)"

    )

    add_custom_target(
        fix
        COMMAND ${CLANG_TIDY_EXE}
        --fix
        --fix-errors
        --fix-notes
        -p=${CMAKE_BINARY_DIR}
        --config-file=${CMAKE_SOURCE_DIR}/.clang-tidy
        ${ALL_CXX_SOURCES}
        COMMENT "Running clang-tidy with fixes"
    )
    message(STATUS "clang-tidy fix target created. Build it with 'make fix'.")
else ()
    message(WARNING "clang-tidy not found, no static analysis")
endif ()
