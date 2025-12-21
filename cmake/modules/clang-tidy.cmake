option(ENABLE_CLANG_TIDY "Enable clang-tidy static analysis during build" OFF)
option(ENABLE_CLANG_TIDY_FIXES "Apply clang-tidy fixes automatically" OFF)

if(ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE clang-tidy)

    if(CLANG_TIDY_EXE)
        set(CLANG_TIDY_COMMAND
            "${CLANG_TIDY_EXE}"
            "--config-file=${CMAKE_SOURCE_DIR}/.clang-tidy"
            "--header-filter=${CMAKE_SOURCE_DIR}/include/.*"
        )

        if(ENABLE_CLANG_TIDY_FIXES)
            list(APPEND CLANG_TIDY_COMMAND "--fix" "--format-style=file")
            message(STATUS "clang-tidy fixes enabled")
        endif()

        message(STATUS "clang-tidy enabled: ${CLANG_TIDY_EXE}")
    else()
        message(FATAL_ERROR "ENABLE_CLANG_TIDY is ON but clang-tidy was not found")
    endif()
endif()

function(add_clang_tidy_support target)
    if(ENABLE_CLANG_TIDY AND CLANG_TIDY_EXE)
        set_target_properties(${target} PROPERTIES
            CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}"
        )
    endif()
endfunction()
