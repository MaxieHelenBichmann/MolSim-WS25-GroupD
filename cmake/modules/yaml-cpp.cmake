FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG yaml-cpp-0.7.0
  DOWNLOAD_EXTRACT_TIMESTAMP true
)
set(YAML_CPP_CLANG_FORMAT_EXE OFF 
    CACHE BOOL "Force clang-format to not be found by yaml-cpp")
FetchContent_MakeAvailable(yaml-cpp)
