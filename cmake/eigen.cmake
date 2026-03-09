include(FetchContent)

set(FETCHCONTENT_QUIET TRUE)

# Define Eigen package details
set(PACKAGE_NAME eigen)
set(REPO_URL "https://gitlab.com/libeigen/eigen.git")
set(REPO_TAG "3.4.0")

# Optimization: Disable Eigen's internal builds (tests, docs, etc.)
set(EIGEN_BUILD_TESTING OFF)
set(EIGEN_MPL2_ONLY ON)
set(EIGEN_BUILD_PKGCONFIG OFF)
set(EIGEN_BUILD_DOC OFF)

# Download and make Eigen available
FetchContent_Declare(
  ${PACKAGE_NAME}
  GIT_REPOSITORY ${REPO_URL}
  GIT_TAG        ${REPO_TAG}
)

FetchContent_MakeAvailable(${PACKAGE_NAME})

# Export the source directory for inclusion
set(EIGEN_INCLUDE_DIR ${eigen_SOURCE_DIR})
