if(ENABLE_TESTING)
        include_directories(
                ${PROJECT_SOURCE_DIR}/src
                ${PROJECT_SOURCE_DIR}/benchmarks/code
                ${CMAKE_CURRENT_SOURCE_DIR})

        include("../CTestConfig.cmake")

        file(GLOB_RECURSE TEST_SRC
                "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

        # Include benchmark container implementations for testing
        file(GLOB_RECURSE BENCHMARK_CONTAINER_SRC
                "${PROJECT_SOURCE_DIR}/benchmarks/code/linkedcellimpl/*.cpp")

        enable_testing()

        add_executable(tests ${TEST_SRC} ${BENCHMARK_CONTAINER_SRC})

        target_compile_definitions(tests PRIVATE -DTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
        target_compile_options(tests PRIVATE -Wall -Wextra -Werror)

        target_link_libraries(
                tests
                LibMolSim
                GTest::gtest_main
                gmock_main)

        include(GoogleTest)
        gtest_discover_tests(tests)

        if(COVERAGE)
                include(GenerateCoverage)
                enable_coverage(tests)
        endif()

        message(STATUS "Testing enabled an tests executable created")
endif()
