# CMakeLists.txt
cmake_minimum_required(VERSION 3.12)
# enable_language(CUDA)
# project(stp LANGUAGES CXX CUDA)
project(stp LANGUAGES CXX)

# Set C++ standard
set(STP_CXX_STANDARD "20" CACHE STRING "C++ standard")
set(CMAKE_CXX_STANDARD ${STP_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Set Eigen include directory
include_directories(/usr/local/include/eigen3)

# Options
option(STP_ENABLE_CUDA "Enable CUDA" OFF)
option(STP_EXAMPLES "Build examples" ON)
option(STP_TEST "Build tests" ON)


if(STP_ENABLE_CUDA)
  set(CUDA_FLAG "-DSTP_ENABLE_CUDA")
endif()
add_compile_options(-g -O3 -W -Wall -Wextra ${CUDA_FLAG})


# Add main include files
add_subdirectory(include)

# Add third-party libraries
add_subdirectory(lib)

if(STP_EXAMPLES)
  add_subdirectory(examples)
endif()

if(STP_TEST)
  add_subdirectory(test)
endif()
