# CMakeLists.txt
cmake_minimum_required(VERSION 3.12)
project(stp LANGUAGES CXX)

# Set C++ standard
set(STP_CXX_STANDARD "17" CACHE STRING "C++ standard")
set(CMAKE_CXX_STANDARD ${STP_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(STP_EXAMPLES "Build examples" ON)
option(STP_TEST "Build tests" ON)

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
