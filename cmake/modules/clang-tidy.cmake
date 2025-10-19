option(ENABLE_CLANG_TIDY "Enable clang-tidy checks" ON)

if(ENABLE_CLANG_TIDY)
    message(STATUS "clang-tidy checks enabled")
    
    find_program(CLANG_TIDY_EXE clang-tidy)
    if(CLANG_TIDY_EXE)

        #For other options like automatic fixes add -fix at the end
        set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_EXE} )
        message(STATUS "clang-tidy found and enabled: ${CLANG_TIDY_EXE}")
    else ()
        message(WARNING "clang-tidy not found, no static analysis")
    endif ()
endif()
