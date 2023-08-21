add_library(stp INTERFACE)

# Include stp headers
target_include_directories(stp SYSTEM INTERFACE ${PROJECT_SOURCE_DIR}/include)

# Include Eigen headers
target_include_directories(stp INTERFACE ${PROJECT_SOURCE_DIR}/lib/eigen)

target_link_libraries(stp INTERFACE fmt Eigen3::Eigen)
