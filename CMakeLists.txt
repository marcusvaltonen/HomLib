cmake_minimum_required(VERSION 3.2)
project(HomLib)

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()

include(CTest)

SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

find_package(Eigen3 REQUIRED)
set(EIGEN_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

include_directories(
    includes/HomLib
    src/helpers
    src/solvers/fitzgibbon_cvpr_2001
    src/solvers/kukelova_cvpr_2015
    src/solvers/valtonenornhag_arxiv_2020a
    src/solvers/valtonenornhag_arxiv_2020b
)

set(SRCS
    src/helpers/generate_problem_instance.cpp
    src/helpers/gj.cpp
    src/helpers/normalize2dpts.cpp
    src/helpers/radial.cpp
    src/helpers/roots.cpp
    src/solvers/fitzgibbon_cvpr_2001/get_fitzgibbon_cvpr_2001.cpp
    src/solvers/kukelova_cvpr_2015/get_kukelova_cvpr_2015.cpp
    src/solvers/kukelova_cvpr_2015/solver_kukelova_cvpr_2015.cpp
    src/solvers/valtonenornhag_arxiv_2020a/get_valtonenornhag_arxiv_2020a_fHf.cpp
    src/solvers/valtonenornhag_arxiv_2020a/solver_valtonenornhag_arxiv_2020a_fHf.cpp
    src/solvers/valtonenornhag_arxiv_2020b/get_valtonenornhag_arxiv_2020b_fHf.cpp
    src/solvers/valtonenornhag_arxiv_2020b/get_valtonenornhag_arxiv_2020b_frHfr.cpp
    src/solvers/valtonenornhag_arxiv_2020b/solver_valtonenornhag_arxiv_2020b_fHf.cpp
    src/solvers/valtonenornhag_arxiv_2020b/solver_valtonenornhag_arxiv_2020b_frHfr.cpp
)

set(EXAMPLE_SRCS example.cpp)

set(TEST_SRCS
    tests/test_fitzgibbon_cvpr_2001.cpp
    tests/test_helpers.cpp
    tests/test_kukelova_cvpr_2015.cpp
    tests/test_main.cpp
    tests/test_valtonenornhag_arxiv_2020a.cpp
    tests/test_valtonenornhag_arxiv_2020b.cpp
)

add_executable(example ${SRCS} ${EXAMPLE_SRCS})
target_link_libraries(example PRIVATE Eigen3::Eigen)

add_executable(tests ${SRCS} ${TEST_SRCS})
target_include_directories(tests PRIVATE tests)
target_link_libraries(tests PRIVATE Eigen3::Eigen)

if(CMAKE_COMPILER_IS_GNUCXX AND CMAKE_BUILD_TYPE STREQUAL "Debug")
message("debug mode")
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")
include(CodeCoverage)
APPEND_COVERAGE_COMPILER_FLAGS()
#set(COVERAGE_LCOV_EXCLUDES 'dir1/*' 'dir2/*')
#SETUP_TARGET_FOR_COVERAGE_LCOV(NAME foo_coverage
#                              EXECUTABLE foo
#                              DEPENDENCIES foo)
endif()
message("hej")

add_test(NAME tests COMMAND tests)
enable_testing()
