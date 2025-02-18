add_library(stp INTERFACE)

# Include stp headers
target_include_directories(stp SYSTEM INTERFACE ${PROJECT_SOURCE_DIR}/include)




# Try to find Eigen3 on the system
find_package(Eigen3 QUIET)
message(STATUS "Eigen version: ${EIGEN3_VERSION}")

# If Eigen3 is not found, use ExternalProject to clone and build Eigen3
if(NOT Eigen3_FOUND)
  message("-- Eigen3 is not found. Please clone and build from GitLab...")
  message("git clone https://gitlab.com/libeigen/eigen.git")
  message("cd eigen")
  message("mkdir build")
  message("cd build")
  message("cmake ..")
  message("make install")
  message(FATAL_ERROR "Please visit https://eigen.tuxfamily.org/ for more details.")
else()
message("found")
message(STATUS "Eigen3 found at: ${EIGEN3_INCLUDE_DIRS}")
 # Add Eigen3 include directory to stp target
  target_include_directories(stp INTERFACE ${EIGEN3_INCLUDE_DIRS})
endif()

target_link_libraries(stp INTERFACE fmt Eigen3::Eigen pthread)
