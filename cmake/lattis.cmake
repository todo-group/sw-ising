include(FetchContent)
FetchContent_Declare(
  lattis
  GIT_REPOSITORY https://github.com/todo-group/lattis.git
  GIT_TAG        main
)
list(APPEND FetchContents lattis)
