cmake_policy(SET CMP0077 NEW)
cmake_minimum_required(VERSION 3.14)
project(dg_euler LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
option(BUILD_TESTS "Build unit tests" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# --- Dependencies ---
include(cmake/eigen.cmake)

# --- OpenMP ---
# Apple Clang does not bundle OpenMP; point CMake at the Homebrew libomp package.
if(APPLE)
  set(OpenMP_C_FLAGS          "-Xpreprocessor -fopenmp")
  set(OpenMP_C_LIB_NAMES      "omp")
  set(OpenMP_CXX_FLAGS        "-Xpreprocessor -fopenmp")
  set(OpenMP_CXX_LIB_NAMES    "omp")
  # Prefer the Homebrew arm64 prefix; fall back to x86_64.
  foreach(_omp_prefix /opt/homebrew/opt/libomp /usr/local/opt/libomp)
    if(EXISTS "${_omp_prefix}/lib/libomp.dylib")
      set(OpenMP_omp_LIBRARY "${_omp_prefix}/lib/libomp.dylib")
      set(OpenMP_C_INCLUDE_DIR   "${_omp_prefix}/include")
      set(OpenMP_CXX_INCLUDE_DIR "${_omp_prefix}/include")
      break()
    endif()
  endforeach()
endif()
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found – parallel build enabled.")
else()
  message(STATUS "OpenMP NOT found – building single-threaded.")
endif()

# --- Core physics / mesh libraries ---
add_library(flux_physics flux.cpp physics.cpp)
target_include_directories(flux_physics PUBLIC ${CMAKE_SOURCE_DIR} ${EIGEN_INCLUDE_DIR})

add_library(mesh_io readgri.cpp)
target_include_directories(mesh_io PUBLIC ${CMAKE_SOURCE_DIR})

add_library(problem problem.cpp)
target_include_directories(problem PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(problem PRIVATE flux_physics)

# --- DG utilities (quadrature, mass matrix, surface/volume integrals) ---
add_library(dg_utils
    geometry.cpp mm.cpp term2.cpp surf.cpp quad.cpp
    givens/quad2d.c givens/shape.c
)
target_include_directories(dg_utils PUBLIC ${CMAKE_SOURCE_DIR} ${EIGEN_INCLUDE_DIR})
target_link_libraries(dg_utils PRIVATE flux_physics mesh_io)
if(OpenMP_CXX_FOUND)
  target_link_libraries(dg_utils PUBLIC OpenMP::OpenMP_CXX)
endif()

# --- DG time advance ---
add_library(advance advance.cpp)
target_include_directories(advance PUBLIC ${CMAKE_SOURCE_DIR} ${EIGEN_INCLUDE_DIR})
target_link_libraries(advance PRIVATE flux_physics mesh_io problem dg_utils)
if(OpenMP_CXX_FOUND)
  target_link_libraries(advance PUBLIC OpenMP::OpenMP_CXX)
endif()

# --- Main executable ---
add_executable(main main.cpp write_vtu.cpp parseinput.cpp restart.cpp)
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR} ${EIGEN_INCLUDE_DIR})
target_link_libraries(main PRIVATE flux_physics mesh_io problem advance dg_utils)

add_custom_target(run COMMAND $<TARGET_FILE:main>
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                  DEPENDS main)

# --- Optional tests ---
if(BUILD_TESTS)
  add_subdirectory(test)
endif()
