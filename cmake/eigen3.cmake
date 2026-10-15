if(TARGET Eigen3::Eigen)
  return()
endif()

find_package(Eigen3 QUIET)
if(Eigen3_FOUND)
  message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
else(Eigen3_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
    SOURCE_SUBDIR none
  )
  FetchContent_MakeAvailable(eigen3)
  if(NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
    set_target_properties(Eigen3::Eigen PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${eigen3_SOURCE_DIR}"
    )
  endif()
endif(Eigen3_FOUND)
