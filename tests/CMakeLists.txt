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
    test_bisect_sturm_scaled.cpp
    test_cai_pami_2025.cpp
    test_fitzgibbon_cvpr_2001.cpp
    test_helpers.cpp
    test_kukelova_cvpr_2015.cpp
    test_nakano_icpr_2025.cpp
    test_refinement.cpp
    test_valtonenornhag_icpr_2020.cpp
    test_valtonenornhag_icpr_2026.cpp
    test_valtonenornhag_wacv_2021.cpp
    test_wadenback_3dv_2026.cpp
)

# Link against Catch2WithMain - this provides the main() function
target_link_libraries(HomLibTests PRIVATE HomLib Catch2::Catch2WithMain Eigen3::Eigen PoseLib::PoseLib)
target_include_directories(HomLibTests PRIVATE ${HOMLIB_INCLUDE_DIRS})

# Automatically discover and register all TEST_CASE macros with CTest
# This requires the Catch2 extras directory in the module path
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(HomLibTests)

# Backend integration tests — separate executable with optional backend dependencies
add_executable(HomLibBackendTests test_backends.cpp)
target_link_libraries(HomLibBackendTests PRIVATE
    HomLib Catch2::Catch2WithMain Eigen3::Eigen PoseLib::PoseLib)
target_include_directories(HomLibBackendTests PRIVATE
    ${HOMLIB_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/thirdparty/RansacLib)

if(HOMLIB_USE_GCRANSAC)
    target_compile_definitions(HomLibBackendTests PRIVATE HOMLIB_HAS_GCRANSAC)
    target_include_directories(HomLibBackendTests PRIVATE ${GCRANSAC_INCLUDE_DIR})
    target_link_libraries(HomLibBackendTests PRIVATE GCoptimization ${OpenCV_LIBS})
endif()

if(HOMLIB_USE_MAGSAC)
    target_compile_definitions(HomLibBackendTests PRIVATE HOMLIB_HAS_MAGSAC)
    target_include_directories(HomLibBackendTests PRIVATE ${GCRANSAC_INCLUDE_DIR} ${MAGSAC_INCLUDE_DIR})
    target_link_libraries(HomLibBackendTests PRIVATE gflags OpenMP::OpenMP_CXX)
    if(NOT HOMLIB_USE_GCRANSAC)
        target_compile_definitions(HomLibBackendTests PRIVATE HOMLIB_HAS_GCRANSAC)
        target_include_directories(HomLibBackendTests PRIVATE ${GCRANSAC_INCLUDE_DIR})
        target_link_libraries(HomLibBackendTests PRIVATE GCoptimization ${OpenCV_LIBS})
    endif()
endif()

if(HOMLIB_USE_SUPERANSAC)
    target_compile_definitions(HomLibBackendTests PRIVATE HOMLIB_HAS_SUPERANSAC)
    target_include_directories(HomLibBackendTests PRIVATE ${SUPERANSAC_INCLUDE_DIR})
    target_link_libraries(HomLibBackendTests PRIVATE SupeRANSAC)
endif()

catch_discover_tests(HomLibBackendTests)
