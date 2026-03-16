# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.15...3.27)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops")

# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

find_package( PoseLib REQUIRED )
include_directories(${CMAKE_SOURCE_DIR})


# Add a library using FindPython's tooling (pybind11 also provides a helper like
# this)
python_add_library(_core MODULE src/main.cpp WITH_SOABI)
target_link_libraries(_core PRIVATE pybind11::headers)

# Add library
add_library(HomLib STATIC IMPORTED)
set_property(TARGET HomLib PROPERTY IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/../build/libHomLib.a")

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../includes/HomLib
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/helpers
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/robust
    ${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty/RansacLib
)
target_link_libraries(_core PRIVATE HomLib)
target_link_libraries(_core PRIVATE PoseLib::PoseLib)

# The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION homlib)
