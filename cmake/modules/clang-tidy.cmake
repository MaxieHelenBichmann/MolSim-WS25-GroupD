option(ENABLE_CLANG_TIDY "Enable clang-tidy static analysis during build" OFF)

if(ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE clang-tidy REQUIRED)
    find_program(RUN_CLANG_TIDY_EXE run-clang-tidy)

    set(CLANG_TIDY_COMMAND
        ${CLANG_TIDY_EXE}
        --config-file=${CMAKE_SOURCE_DIR}/.clang-tidy
        --header-filter=${CMAKE_SOURCE_DIR}/include/.*
    )

    set(RUN_CLANG_TIDY_CMD
        ${RUN_CLANG_TIDY_EXE}
        -clang-tidy-binary=${CLANG_TIDY_EXE}
        -p=${CMAKE_BINARY_DIR}
        -quiet
        -warnings-as-errors=*
        -config-file=${CMAKE_SOURCE_DIR}/.clang-tidy
        -header-filter=${CMAKE_SOURCE_DIR}/include/.*
        -source-filter=${CMAKE_SOURCE_DIR}/src/.*
    )

    # add fix and lint targets that runs run-clang-tidy
    add_custom_target(lint COMMAND ${RUN_CLANG_TIDY_CMD} 
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Running clang-tidy analysis"
    )

    add_custom_target(fix COMMAND ${RUN_CLANG_TIDY_CMD} -fix
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Running clang-tidy and applying fixes"
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

