cmake_minimum_required(VERSION 3.14) # FetchContent works best with 3.14+
project(flux LANGUAGES C CXX)         # Enable C for quad2d.c and shape.c

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# --- 1. Dependencies ---
include(cmake/eigen.cmake)

# --- 2. Existing Libraries ---
add_library(flux_physics flux.cpp physics.cpp)
target_include_directories(flux_physics PUBLIC ${CMAKE_SOURCE_DIR})

add_library(mesh_io readgri.cpp)
target_include_directories(mesh_io PUBLIC ${CMAKE_SOURCE_DIR})

add_library(problem problem.cpp)
target_include_directories(problem PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(problem PRIVATE flux_physics)

# --- 3. New DG Components ---
# Include the C files for quadrature and shape functions
add_library(dg_utils 
    mm.cpp 
    term2.cpp 
    quad.cpp 
    givens/quad2d.c 
    givens/shape.c
)
target_include_directories(dg_utils PUBLIC ${CMAKE_SOURCE_DIR} ${EIGEN_INCLUDE_DIR})
target_link_libraries(dg_utils PRIVATE flux_physics mesh_io)

# --- 4. Solver Logic ---
add_library(advance advance.cpp)
target_include_directories(advance PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(advance PRIVATE flux_physics mesh_io problem dg_utils)

# --- 5. Executable ---
add_executable(main main.cpp)
target_link_libraries(main PRIVATE flux_physics mesh_io problem advance dg_utils)

# Run target
add_custom_target(run COMMAND $<TARGET_FILE:main> WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} DEPENDS main)

add_subdirectory(test)
