# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.15...3.27)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX
)

option(HOMLIB_NATIVE_ARCH "Use -march=native for CPU-specific optimizations (not portable)" OFF)

if(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  string(APPEND CMAKE_CXX_FLAGS_RELEASE " -funroll-loops")
  if(HOMLIB_NATIVE_ARCH)
    add_compile_options(-march=native)
  endif()
endif()

message(STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")

# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.12.0
  )
  FetchContent_MakeAvailable(pybind11)
endif()

# Common include dirs for all binding compilation units
set(BINDING_COMMON_INCLUDES
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${HOMLIB_INCLUDE_DIRS}
)
set(BINDING_THIRDPARTY_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty/RansacLib)

# Helper: configure common properties for a binding object library
function(configure_binding_object target)
    target_include_directories(${target} PRIVATE ${BINDING_COMMON_INCLUDES})
    target_include_directories(${target} SYSTEM PRIVATE ${BINDING_THIRDPARTY_INCLUDES} ${Python_INCLUDE_DIRS})
    target_link_libraries(${target} PRIVATE pybind11::headers HomLib PoseLib::PoseLib Eigen3::Eigen)
    target_compile_definitions(${target} PRIVATE HOMLIB_VERSION="${PROJECT_VERSION_FULL}")
    set_target_properties(${target} PROPERTIES
        POSITION_INDEPENDENT_CODE ON)
    if(MSVC)
        target_compile_definitions(${target} PRIVATE _USE_MATH_DEFINES)
        target_compile_options(${target} PRIVATE /bigobj)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE $<$<CONFIG:Release>:-funroll-loops>)
    endif()
    homlib_strict_standard(${target})
endfunction()

# Main module — core bindings, LOMSAC (no backend-specific headers)
python_add_library(_core MODULE src/main.cpp WITH_SOABI)
target_include_directories(_core PRIVATE ${BINDING_COMMON_INCLUDES})
target_include_directories(_core SYSTEM PRIVATE ${BINDING_THIRDPARTY_INCLUDES})
target_link_libraries(_core PRIVATE pybind11::headers)
target_link_libraries(_core PRIVATE HomLib)
target_link_libraries(_core PRIVATE PoseLib::PoseLib)
target_compile_definitions(_core PRIVATE HOMLIB_VERSION="${PROJECT_VERSION_FULL}")
if(MSVC)
    target_compile_definitions(_core PRIVATE _USE_MATH_DEFINES)
endif()
homlib_strict_standard(_core)

# SupeRANSAC backend — separate compilation unit with only SupeRANSAC includes
if(HOMLIB_USE_SUPERANSAC)
    add_library(bind_superansac_obj OBJECT src/bind_superansac.cpp)
    configure_binding_object(bind_superansac_obj)
    target_include_directories(bind_superansac_obj SYSTEM PRIVATE ${SUPERANSAC_INCLUDE_DIR})
    target_link_libraries(bind_superansac_obj PRIVATE SupeRANSAC)
    target_compile_definitions(bind_superansac_obj PRIVATE HOMLIB_HAS_SUPERANSAC)

    target_link_libraries(_core PRIVATE bind_superansac_obj SupeRANSAC)
    target_compile_definitions(_core PRIVATE HOMLIB_HAS_SUPERANSAC)
    if(APPLE)
        set_target_properties(_core PROPERTIES
            INSTALL_RPATH "@loader_path"
            BUILD_WITH_INSTALL_RPATH TRUE)
    elseif(NOT MSVC)
        set_target_properties(_core PROPERTIES
            INSTALL_RPATH "$ORIGIN"
            BUILD_WITH_INSTALL_RPATH TRUE)
    endif()
endif()

# GC-RANSAC backend — separate compilation unit with only GCRANSAC includes
if(HOMLIB_USE_GCRANSAC)
    add_library(bind_gcransac_obj OBJECT src/bind_gcransac.cpp)
    configure_binding_object(bind_gcransac_obj)
    target_include_directories(bind_gcransac_obj SYSTEM PRIVATE ${GCRANSAC_INCLUDE_DIR})
    target_link_libraries(bind_gcransac_obj PRIVATE GCoptimization ${OpenCV_LIBS})
    target_compile_definitions(bind_gcransac_obj PRIVATE HOMLIB_HAS_GCRANSAC)

    target_link_libraries(_core PRIVATE bind_gcransac_obj GCoptimization ${OpenCV_LIBS})
    target_compile_definitions(_core PRIVATE HOMLIB_HAS_GCRANSAC)
endif()

# MAGSAC backend — separate compilation unit with MAGSAC includes, and GC-RANSAC's headers
# (MAGSAC++ uses GC-RANSAC's samplers and model classes)
if(HOMLIB_USE_MAGSAC)
    add_library(bind_magsac_obj OBJECT src/bind_magsac.cpp)
    configure_binding_object(bind_magsac_obj)
    target_include_directories(bind_magsac_obj SYSTEM PRIVATE ${GCRANSAC_INCLUDE_DIR} ${MAGSAC_INCLUDE_DIR})
    target_link_libraries(bind_magsac_obj PRIVATE ${OpenCV_LIBS} gflags OpenMP::OpenMP_CXX)
    target_compile_definitions(bind_magsac_obj PRIVATE HOMLIB_HAS_MAGSAC USE_OPENMP)

    target_link_libraries(_core PRIVATE bind_magsac_obj ${OpenCV_LIBS} gflags)
    target_compile_definitions(_core PRIVATE HOMLIB_HAS_MAGSAC)
endif()

# The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION homlib)
if(HOMLIB_USE_SUPERANSAC)
    install(TARGETS SupeRANSAC DESTINATION homlib)
endif()
