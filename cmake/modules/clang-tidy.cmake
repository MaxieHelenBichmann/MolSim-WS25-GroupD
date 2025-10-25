option(ENABLE_CLANG_TIDY "Enable clang-tidy checks" ON)

if(ENABLE_CLANG_TIDY)
    message(STATUS "clang-tidy checks enabled")
    
    find_program(CLANG_TIDY_EXE clang-tidy)
    if(CLANG_TIDY_EXE)
        #For other options like automatic fixes add --fix at the end
        set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_EXE} -p ${CMAKE_BINARY_DIR})
        message(STATUS "clang-tidy found and enabled: ${CLANG_TIDY_EXE}")

        file(GLOB_RECURSE ALL_CXX_SOURCES CONFIGURE_DEPENDS
            "${CMAKE_SOURCE_DIR}/src/*.cpp"
        )

        add_custom_target(
            fix
            COMMAND ${CLANG_TIDY_EXE}
            --fix
            --fix-errors
            --fix-notes
            -p=${CMAKE_BINARY_DIR}
            ${ALL_CXX_SOURCES}
            COMMENT "Running clang-tidy with fixes"
        )
    else ()
        message(WARNING "clang-tidy not found, no static analysis")
    endif ()
endif()
