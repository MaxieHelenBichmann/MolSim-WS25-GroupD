find_package(Doxygen)

if(DOXYGEN_FOUND)
    add_custom_target(doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating documentation with Doxygen"
        VERBATIM
    )
    message(STATUS "Doxygen target created. Build it with 'make doxygen'.")
else()
    message(WARNING "Doxygen not found. The 'doxygen' target will not be available.")
endif()

