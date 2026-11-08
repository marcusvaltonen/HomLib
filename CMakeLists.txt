cmake_minimum_required(VERSION 3.15)

# ----------------------------------------------------------------------
# Extract version from pyproject.toml
# ----------------------------------------------------------------------
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml" PYPROJECT_CONTENT)
string(REGEX MATCH "version[ \t]*=[ \t]*[\"']([^\"']+)[\"']" _ ${PYPROJECT_CONTENT})
set(PROJECT_VERSION_FULL ${CMAKE_MATCH_1})

if(NOT PROJECT_VERSION_FULL)
    message(FATAL_ERROR "Version not found in pyproject.toml")
endif()

message(STATUS "Project version: ${PROJECT_VERSION_FULL}")

# ----------------------------------------------------------------------
# Extract numeric version for CMake
# ----------------------------------------------------------------------
# Look for the pattern "major.minor.patch" at the beginning
string(REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)" _ "${PROJECT_VERSION_FULL}")
if(CMAKE_MATCH_1)
    set(VERSION_MAJOR ${CMAKE_MATCH_1})
    set(VERSION_MINOR ${CMAKE_MATCH_2})
    set(VERSION_PATCH ${CMAKE_MATCH_3})
    set(PROJECT_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
else()
    # Fallback if no numeric parts found
    set(PROJECT_VERSION "0.0.0")
endif()

project(HomLib VERSION ${PROJECT_VERSION} LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wno-strict-aliasing -Wno-write-strings -Wno-deprecated)
    add_compile_options($<$<CONFIG:Release>:-funroll-loops>)
    if(HOMLIB_NATIVE_ARCH)
        add_compile_options(-march=native)
    endif()
endif()

# C++17 by default; configure with -DCMAKE_CXX_STANDARD=20 or 23 to build in a newer standard.
if(NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
elseif(CMAKE_CXX_STANDARD LESS 17 OR CMAKE_CXX_STANDARD EQUAL 98)
    message(FATAL_ERROR "HomLib needs C++17 or later (CMAKE_CXX_STANDARD is ${CMAKE_CXX_STANDARD})")
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(HOMLIB_STRICT_STANDARD "GCC and clang: compile HomLib's own code as ISO C++ (no compiler extensions), with -pedantic-errors and deprecations as errors (MSVC ignores these flags)" OFF)

# HomLib's own targets call this; the third-party code is included as system headers, so that it
# isn't checked
function(homlib_strict_standard target)
    if(HOMLIB_STRICT_STANDARD)
        set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
        target_compile_options(${target} PRIVATE -pedantic-errors -Werror=deprecated -Werror=deprecated-declarations
            # GCC (11 and later knows these flags) keeps these C++20 deprecations out of -Wdeprecated
            $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,11>>:-Werror=volatile;-Werror=deprecated-enum-enum-conversion;-Werror=deprecated-enum-float-conversion;-Werror=deprecated-copy-dtor>)
    endif()
endfunction()

option(HOMLIB_USE_SUPERANSAC "Enable SupeRANSAC backend (requires Boost, OpenCV)" ON)
option(HOMLIB_USE_GCRANSAC "Enable GC-RANSAC backend (requires OpenCV)" ON)
option(HOMLIB_USE_MAGSAC "Enable MAGSAC backend (requires OpenCV, gflags)" ON)
option(HOMLIB_NATIVE_ARCH "Use -march=native for CPU-specific optimizations (not portable)" OFF)

find_package(Eigen3 REQUIRED)
set(EIGEN_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

find_package(PoseLib REQUIRED)
include_directories(${CMAKE_SOURCE_DIR})

if(HOMLIB_USE_SUPERANSAC)
    set(SUPERANSAC_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/superansac/include)
    add_subdirectory(thirdparty/superansac/src)
    # The submodule's src/CMakeLists.txt uses PROJECT_SOURCE_DIR/include which
    # resolves incorrectly when built as a subdirectory. Fix the include path.
    target_include_directories(SupeRANSAC PUBLIC ${SUPERANSAC_INCLUDE_DIR})
    # SupeRANSAC and GC-RANSAC (also used by MAGSAC) each bundle a different version of
    # Kolmogorov's max-flow classes Graph<> and Energy<> in the global namespace
    # (SupeRANSAC's is index-based, GC-RANSAC's pointer-based). Both end up in the
    # Python module, where the linker keeps only one definition of each inline member
    # (an ODR violation): rename SupeRANSAC's, in every target that uses its headers.
    target_compile_definitions(SupeRANSAC PUBLIC Graph=SuperansacGraph Energy=SuperansacEnergy)
endif()

if(HOMLIB_USE_GCRANSAC OR HOMLIB_USE_MAGSAC)
    find_package(OpenCV REQUIRED)
    set(GCRANSAC_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/graph-cut-ransac/src/pygcransac/include)

    # GCoptimization sources needed by GCRANSAC headers
    set(GCOPTIMIZATION_SOURCES
        ${GCRANSAC_INCLUDE_DIR}/GCoptimization.cpp
        ${GCRANSAC_INCLUDE_DIR}/LinkedBlockList.cpp
        ${GCRANSAC_INCLUDE_DIR}/graph.cpp
        ${GCRANSAC_INCLUDE_DIR}/maxflow.cpp
    )
    add_library(GCoptimization STATIC ${GCOPTIMIZATION_SOURCES})
    target_include_directories(GCoptimization PUBLIC ${GCRANSAC_INCLUDE_DIR})
    target_link_libraries(GCoptimization PUBLIC Eigen3::Eigen ${OpenCV_LIBS})
    set_target_properties(GCoptimization PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

if(HOMLIB_USE_MAGSAC)
    find_package(gflags REQUIRED)
    # Apple clang has no OpenMP runtime of its own, and Homebrew doesn't link its libomp into
    # Homebrew's prefix (keg-only): point FindOpenMP to it, unless OpenMP_ROOT is given.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND NOT DEFINED OpenMP_ROOT AND NOT DEFINED ENV{OpenMP_ROOT})
        foreach(_libomp $ENV{HOMEBREW_PREFIX}/opt/libomp /opt/homebrew/opt/libomp /usr/local/opt/libomp)
            if(EXISTS ${_libomp}/include/omp.h)
                set(OpenMP_ROOT ${_libomp})
                break()
            endif()
        endforeach()
    endif()
    find_package(OpenMP REQUIRED)
    set(MAGSAC_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/magsac/src/pymagsac/include)
    # magsac.h reads a static constexpr member function through a reference parameter
    # (constexpr size_t sample_size = estimator_.sampleSize()). P2280 (C++23, adopted as a defect
    # report) allows it; GCC applies it in every mode, Clang only in C++23 mode. If the compiler
    # rejects it, use a copy of magsac.h that calls the function through the type, in a directory
    # searched before the original one (unless magsac.h no longer has the construct).
    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles("
        struct E { static constexpr int n() { return 4; } };
        template <class T> int f(const T &e) { constexpr int n = e.n(); return n; }
        int main() { return f(E()) - 4; }" HOMLIB_CXX_CONSTEXPR_THROUGH_REFERENCE)
    if(NOT HOMLIB_CXX_CONSTEXPR_THROUGH_REFERENCE)
        set(_magsac_fixed_dir ${CMAKE_CURRENT_BINARY_DIR}/magsac_fixed)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${MAGSAC_INCLUDE_DIR}/magsac.h)
        file(READ ${MAGSAC_INCLUDE_DIR}/magsac.h _magsac)
        string(REPLACE "estimator_.sampleSize()" "ModelEstimator::sampleSize()" _magsac_fixed "${_magsac}")
        if(_magsac_fixed STREQUAL _magsac)
            file(REMOVE_RECURSE ${_magsac_fixed_dir})
            message(STATUS "MAGSAC: magsac.h has no estimator_.sampleSize() any more, using it as is "
                "(the workaround in ${CMAKE_CURRENT_LIST_FILE} is no longer needed)")
        else()
            # configure_file only writes magsac.h when its content changes, so that a new
            # configure doesn't recompile the files that include it.
            file(WRITE ${_magsac_fixed_dir}/magsac.h.in "${_magsac_fixed}")
            configure_file(${_magsac_fixed_dir}/magsac.h.in ${_magsac_fixed_dir}/magsac.h COPYONLY)
            message(STATUS "MAGSAC: the compiler rejects magsac.h, using a fixed copy in ${_magsac_fixed_dir}")
            set(MAGSAC_INCLUDE_DIR ${_magsac_fixed_dir} ${MAGSAC_INCLUDE_DIR})
        endif()
    endif()
endif()


set(HOMLIB_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/includes/HomLib
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/helpers
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/robust
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/solvers/cai_pami_2025
    ${CMAKE_CURRENT_SOURCE_DIR}/HomLib/solvers/dlt
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
