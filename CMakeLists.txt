cmake_minimum_required(VERSION 3.2)
project(HomLib)

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()

include(CTest)


IF(CMAKE_COMPILER_IS_GNUCXX)
  message("USING GNUCXX compiler")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-strict-aliasing")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-write-strings")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -ansi")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)
set(EIGEN_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

find_package( PoseLib REQUIRED )
include_directories(${CMAKE_SOURCE_DIR})

include_directories("thirdparty/RansacLib")


include_directories(
    includes/HomLib
    src/helpers
    src/robust
    src/solvers/fitzgibbon_cvpr_2001
    src/solvers/kukelova_cvpr_2015
    src/solvers/nakano_icpr_2025
    src/solvers/valtonenornhag_arxiv_2020a
    src/solvers/valtonenornhag_arxiv_2020b
    src/solvers/wadenback_2025
)

set(SRCS
    src/helpers/generate_problem_instance.cpp
    src/helpers/gj.cpp
    src/helpers/normalize2dpts.cpp
    src/helpers/radial.cpp
    src/helpers/roots.cpp
    src/robust/refinement.cpp
    src/solvers/fitzgibbon_cvpr_2001/get_fitzgibbon_cvpr_2001.cpp
    src/solvers/fitzgibbon_cvpr_2001/get_fitzgibbon_cvpr_2001_single.cpp
    src/solvers/kukelova_cvpr_2015/get_kukelova_cvpr_2015.cpp
    src/solvers/kukelova_cvpr_2015/get_kukelova_cvpr_2015_6pt.cpp
    src/solvers/kukelova_cvpr_2015/solver_kukelova_cvpr_2015.cpp
    src/solvers/nakano_icpr_2025/get_nakano_icpr_2025.cpp
    src/solvers/valtonenornhag_arxiv_2020a/get_valtonenornhag_arxiv_2020a_fHf.cpp
    src/solvers/valtonenornhag_arxiv_2020a/solver_valtonenornhag_arxiv_2020a_fHf.cpp
    src/solvers/valtonenornhag_arxiv_2020b/get_valtonenornhag_arxiv_2020b_fHf.cpp
    src/solvers/valtonenornhag_arxiv_2020b/get_valtonenornhag_arxiv_2020b_frHfr.cpp
    src/solvers/valtonenornhag_arxiv_2020b/solver_valtonenornhag_arxiv_2020b_fHf.cpp
    src/solvers/valtonenornhag_arxiv_2020b/solver_valtonenornhag_arxiv_2020b_frHfr.cpp
    src/solvers/wadenback_2025/get_wadenback_2025_double_sided.cpp
    src/solvers/wadenback_2025/get_wadenback_2025_double_sided_equal.cpp
    src/solvers/wadenback_2025/get_wadenback_2025_one_sided.cpp
    src/solvers/wadenback_2025/get_wadenback_2025_unsided.cpp
    src/solvers/wadenback_2025/wadenback_common.cpp
)

set(TEST_SRCS
    tests/test_fitzgibbon_cvpr_2001.cpp
    tests/test_helpers.cpp
    tests/test_kukelova_cvpr_2015.cpp
    tests/test_main.cpp
    tests/test_valtonenornhag_arxiv_2020a.cpp
    tests/test_valtonenornhag_arxiv_2020b.cpp
)

add_library(HomLib STATIC ${SRCS})
target_link_libraries(HomLib PRIVATE Eigen3::Eigen)
target_link_libraries(HomLib PRIVATE PoseLib::PoseLib)
target_compile_options(HomLib PRIVATE -fPIC)

add_executable(example example.cpp)
target_link_libraries(example LINK_PUBLIC HomLib Eigen3::Eigen PoseLib::PoseLib)

add_executable(example_ransac example_ransac.cpp)
target_link_libraries(example_ransac LINK_PUBLIC HomLib Eigen3::Eigen PoseLib::PoseLib)



add_executable(tests ${SRCS} ${TEST_SRCS})
target_include_directories(tests PRIVATE tests)
target_link_libraries(tests PRIVATE Eigen3::Eigen PoseLib::PoseLib)

if(CMAKE_COMPILER_IS_GNUCXX AND CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("Code Coverage ")
  LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")
  include(CodeCoverage)
  APPEND_COVERAGE_COMPILER_FLAGS()
endif()

add_test(NAME tests COMMAND tests)
enable_testing()
