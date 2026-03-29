cmake_minimum_required(VERSION 3.15)

project(HomLib)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
    message("USING GNUCXX compiler")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-strict-aliasing")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-write-strings")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -ansi")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)
set(EIGEN_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

find_package(PoseLib REQUIRED)
include_directories(${CMAKE_SOURCE_DIR})


set(HOMLIB_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/includes/HomLib
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/helpers
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/robust
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/solvers/fitzgibbon_cvpr_2001
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/solvers/kukelova_cvpr_2015
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/solvers/nakano_icpr_2025
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/solvers/valtonenornhag_icpr_2020
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/solvers/valtonenornhag_wacv_2021
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/solvers/wadenback_3dv_2026
)

add_subdirectory(HomLib)

option(BUILD_EXAMPLES "Build examples" ON)
option(BUILD_TESTS "Build tests" ON)

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(BUILD_TESTS)
    add_subdirectory(tests)
endif()

if(SKBUILD)
    add_subdirectory(python)
endif()
