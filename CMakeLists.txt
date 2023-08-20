# CMakeLists.txt
cmake_minimum_required(VERSION 3.12)
project(stp)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)

# Create a target for the header-only library
add_library(stp INTERFACE)
target_include_directories(stp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Add third-party libraries
add_subdirectory(lib)

# Add an executable that uses the header-only library
add_executable(test test.cpp)
target_link_libraries(test PRIVATE stp fmt)
