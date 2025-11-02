if(TEST_ENABLED)
        include_directories(
                ${PROJECT_SOURCE_DIR}/src
                ${CMAKE_CURRENT_SOURCE_DIR})

        include("../CTestConfig.cmake")        

        file(GLOB_RECURSE TEST_SRC
                "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

        enable_testing()

        add_executable(tests ${TEST_SRC})

        target_compile_definitions(tests PRIVATE -DTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
        target_compile_options(tests PRIVATE -Wall -Wextra -Werror)

        target_link_libraries(
                tests
                LibMolSim
                GTest::gtest_main)

        include(GoogleTest)
        gtest_discover_tests(tests)
        message(STATUS "Testing enabled an tests executable created")
endif()
