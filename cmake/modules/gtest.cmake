option(ENABLE_TESTING "Enable Testing with gtest" ON)
if(ENABLE_TESTING)
  message(STATUS "Building tests is enabled")
  include(FetchContent)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.17.0
    SYSTEM
  )
  FetchContent_MakeAvailable(googletest)
endif()
