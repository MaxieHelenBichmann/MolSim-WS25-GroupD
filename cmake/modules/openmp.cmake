option(ENABLE_OPENMP "Enable parallelization with OpenMP" OFF)
if(ENABLE_OPENMP)
    message(STATUS "Parallelization with OpenMP enabled")
    find_package(OpenMP REQUIRED)

    if(OpenMP_CXX_FOUND)
        message (STATUS "Found OpenMP Version: ${OpenMP_CXX_VERSION}")
    else ()
        message(FATAL_ERROR "OpenMP_CXX not found")
    endif ()

    target_link_libraries(LibMolSim
            PUBLIC
            OpenMP::OpenMP_CXX
    )
endif()
