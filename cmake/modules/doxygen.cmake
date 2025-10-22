option(ENABLE_DOXYGEN "Enable Doxygen generation" ON)

find_package(Doxygen)

if(ENABLE_DOXYGEN)
    message(STATUS "Doxygen generation enabled")
    if(DOXYGEN_FOUND)
        add_custom_target(doc_doxygen
            COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )
        message(STATUS "Doxygen target created. Build it with 'make doc_doxygen'.")
    else()
        message(WARNING "Doxygen not found. The 'doc_doxygen' target will not be available.")
    endif()
endif()
