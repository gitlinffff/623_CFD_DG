add_executable(check_eigen test_eigen.cpp)

# Link to Eigen (use the variable set in your eigen.cmake)
target_include_directories(check_eigen PRIVATE ${EIGEN_INCLUDE_DIR})

# add DG utils if test is needed
# target_link_libraries(check_eigen PRIVATE dg_utils)
