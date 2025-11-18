include(FetchContent)
  FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11
    GIT_TAG        v2.6.1
    SYSTEM
  )
  FetchContent_MakeAvailable(cli11)
