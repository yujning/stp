add_library(stp INTERFACE)
target_include_directories(stp SYSTEM INTERFACE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(stp INTERFACE fmt)

