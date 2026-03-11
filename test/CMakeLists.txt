add_executable(check_eigen test_eigen.cpp)
target_include_directories(check_eigen PRIVATE ${EIGEN_INCLUDE_DIR})

add_executable(test_quad test_quad.cpp)
target_link_libraries(test_quad PRIVATE dg_utils)

add_executable(test_freestream test_freestream.cpp)
target_link_libraries(test_freestream PRIVATE dg_utils problem flux_physics mesh_io)
