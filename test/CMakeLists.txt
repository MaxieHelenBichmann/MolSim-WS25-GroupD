include_directories(
        ${PROJECT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB_RECURSE TEST_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

enable_testing()

add_executable(tests ${TEST_SRC})

target_link_libraries(
        tests
        LibMolSim
        GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(tests)
