add_library(stp INTERFACE)

# Include stp headers
target_include_directories(stp SYSTEM INTERFACE ${PROJECT_SOURCE_DIR}/include)

# Try to find Eigen3 on the system
find_package(Eigen3 QUIET)

# If Eigen3 is not found, use ExternalProject to clone and build Eigen3
if(NOT Eigen3_FOUND)
    message("Eigen3 not found. Cloning and building from GitLab...")

    include(ExternalProject)
    ExternalProject_Add(
        Eigen3Project
        PREFIX ${CMAKE_BINARY_DIR}/eigen_external
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG master  # Or any specific version you want
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/eigen_install
    )

    # Set Eigen3_DIR to the installation path
    set(Eigen3_DIR ${CMAKE_BINARY_DIR}/eigen_external)

    # Include the Eigen3Config.cmake after building
    find_package(Eigen3 REQUIRED CONFIG)
    message("Eigen3 cloned, built, and configured.")
endif()


target_link_libraries(stp INTERFACE fmt Eigen3::Eigen)
