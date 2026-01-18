option(ENABLE_CLANG_TIDY "Enable clang-tidy static analysis during build" OFF)
option(CLANG_TIDY_WARNINGS_AS_ERRORS "Treat clang-tidy warnings as errors (applies to clang-tidy and run-clang-tidy)" OFF)

if(ENABLE_CLANG_TIDY)
find_program(CLANG_TIDY_EXE clang-tidy REQUIRED)

    set(CLANG_TIDY_COMMAND
        ${CLANG_TIDY_EXE}
        --config-file=${CMAKE_SOURCE_DIR}/.clang-tidy
        --header-filter=${CMAKE_SOURCE_DIR}/include/.*
    )

    message(STATUS "clang-tidy enabled: ${CLANG_TIDY_EXE}")
endif()

function(add_clang_tidy_support target)
    if(ENABLE_CLANG_TIDY)
        set_target_properties(${target} PROPERTIES
            CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}"
        )
    endif()
endfunction()


find_program(RUN_CLANG_TIDY_EXE run-clang-tidy)

set(RUN_CLANG_TIDY_CMD
    ${RUN_CLANG_TIDY_EXE}
    -clang-tidy-binary=${CLANG_TIDY_EXE}
    -p=${CMAKE_BINARY_DIR}
    -quiet
    -config-file=${CMAKE_SOURCE_DIR}/.clang-tidy
    -header-filter=${CMAKE_SOURCE_DIR}/include/.*
)

if (CLANG_TIDY_WARNINGS_AS_ERRORS)
    list(APPEND CLANG_TIDY_COMMAND "--warnings-as-errors=*")
    list(APPEND RUN_CLANG_TIDY_CMD "-warnings-as-errors=*")
endif()

# add fix and lint targets that runs run-clang-tidy
add_custom_target(lint COMMAND ${RUN_CLANG_TIDY_CMD} 
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running clang-tidy analysis"
)

add_custom_target(fix COMMAND ${RUN_CLANG_TIDY_CMD} -fix -warnings-as-errors="" -extra-arg=-Wno-error
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running clang-tidy and applying fixes"
)
