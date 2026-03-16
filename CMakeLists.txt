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
    src/solvers/valtonenornhag_icpr_2020
    src/solvers/valtonenornhag_wacv_2021
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
    src/solvers/valtonenornhag_icpr_2020/get_valtonenornhag_icpr_2020_fHf.cpp
    src/solvers/valtonenornhag_icpr_2020/solver_valtonenornhag_icpr_2020_fHf.cpp
    src/solvers/valtonenornhag_wacv_2021/get_valtonenornhag_wacv_2021_fHf.cpp
    src/solvers/valtonenornhag_wacv_2021/get_valtonenornhag_wacv_2021_frHfr.cpp
    src/solvers/valtonenornhag_wacv_2021/solver_valtonenornhag_wacv_2021_fHf.cpp
    src/solvers/valtonenornhag_wacv_2021/solver_valtonenornhag_wacv_2021_frHfr.cpp
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
    tests/test_valtonenornhag_icpr_2020.cpp
    tests/test_valtonenornhag_wacv_2021.cpp
)

add_library(HomLib STATIC ${SRCS})
target_link_libraries(HomLib PRIVATE Eigen3::Eigen)
target_link_libraries(HomLib PRIVATE PoseLib::PoseLib)
target_compile_options(HomLib PRIVATE -fPIC)

add_executable(example example.cpp)
target_link_libraries(example LINK_PUBLIC HomLib Eigen3::Eigen PoseLib::PoseLib)

add_executable(example_ransac example_ransac.cpp)
target_link_libraries(example_ransac LINK_PUBLIC HomLib Eigen3::Eigen PoseLib::PoseLib)

if(CMAKE_COMPILER_IS_GNUCXX AND CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("Code Coverage ")
  LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")
  include(CodeCoverage)
  APPEND_COVERAGE_COMPILER_FLAGS()
endif()

#add_test(NAME tests COMMAND tests)
#enable_testing()

# Option to disable tests (e.g., for production builds)
option(BUILD_TESTING "Build tests" ON)

if(BUILD_TESTING)
    # Enable CTest - this creates the 'test' target
    include(CTest)
    
    set(Catch2_DIR "" CACHE INTERNAL "")
    
    # Either FetchContent...
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.13.0
    )
    FetchContent_MakeAvailable(Catch2)
    
    # OR find_package if installed system-wide:
    # find_package(Catch2 3 REQUIRED)
    
    # Add your test executable
    add_executable(HomLibTests 
        tests/test_fitzgibbon_cvpr_2001.cpp
        tests/test_helpers.cpp
        tests/test_kukelova_cvpr_2015.cpp
        tests/test_valtonenornhag_icpr_2020.cpp
        tests/test_valtonenornhag_wacv_2021.cpp
    )
    
    # Link against Catch2WithMain - this provides the main() function
    target_link_libraries(HomLibTests PRIVATE HomLib Catch2::Catch2WithMain Eigen3::Eigen PoseLib::PoseLib)
    
    # Automatically discover and register all TEST_CASE macros with CTest
    # This requires the Catch2 extras directory in the module path
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    include(Catch)
    catch_discover_tests(HomLibTests)
endif()
