include(FetchContent)
FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG 0.8.0
  DOWNLOAD_EXTRACT_TIMESTAMP true
  SYSTEM
)
set(YAML_CPP_CLANG_FORMAT_EXE OFF 
    CACHE BOOL "Force clang-format to not be found by yaml-cpp")
FetchContent_MakeAvailable(yaml-cpp)
set_property(DIRECTORY ${yaml-cpp_SOURCE_DIR} PROPERTY SYSTEM TRUE)
if(TARGET yaml-cpp)
    target_compile_options(yaml-cpp PRIVATE -Wno-shadow)
endif()
