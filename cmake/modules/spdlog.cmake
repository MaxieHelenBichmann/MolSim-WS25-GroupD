include(FetchContent)

FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG v1.14.1
  DOWNLOAD_EXTRACT_TIMESTAMP true
  SYSTEM
)

set(SPDLOG_FMT_EXTERNAL OFF CACHE BOOL "Use external fmt library" FORCE)
FetchContent_MakeAvailable(spdlog)
