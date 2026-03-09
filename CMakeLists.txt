cmake_minimum_required(VERSION 3.10)
project(flux LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

# Release by default for performance
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Core library: flux + physics
add_library(flux_physics flux.cpp physics.cpp)
target_include_directories(flux_physics PUBLIC ${CMAKE_SOURCE_DIR})

# Mesh I/O (read .gri)
add_library(mesh_io readgri.cpp)
target_include_directories(mesh_io PUBLIC ${CMAKE_SOURCE_DIR})

# Reconstruction (no-recon, MUSCL, etc.)
#add_library(reconstruct reconstruct.cpp)
#target_include_directories(reconstruct PUBLIC ${CMAKE_SOURCE_DIR})

# Problem parameters and initialization
add_library(problem problem.cpp)
target_include_directories(problem PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(problem PRIVATE flux_physics)

# BC type mapping
#add_library(bc bc.cpp)
#target_include_directories(bc PUBLIC ${CMAKE_SOURCE_DIR})

# VTU export for Paraview
#add_library(write_vtu write_vtu.cpp)
#target_include_directories(write_vtu PUBLIC ${CMAKE_SOURCE_DIR})
#target_link_libraries(write_vtu PRIVATE flux_physics mesh_io)

# Time advance and residual
add_library(advance advance.cpp)
target_include_directories(advance PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(advance PRIVATE flux_physics mesh_io problem)
#target_link_libraries(advance PRIVATE flux_physics mesh_io problem bc reconstruct write_vtu)

# Main solver executable
add_executable(main main.cpp)
#target_link_libraries(main PRIVATE flux_physics mesh_io problem advance write_vtu)
target_link_libraries(main PRIVATE flux_physics mesh_io problem advance)
add_custom_target(run COMMAND $<TARGET_FILE:main> WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} DEPENDS main)
