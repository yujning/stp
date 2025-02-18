
#Set CUDA architecture (for RTX 3090)
set(CMAKE_CUDA_ARCHITECTURES 86)
#Set CUDA compiler options to generate position-independent code
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fPIC")
include_directories(../include)
add_library(  stp_cuda  STATIC stp_cuda.cu)
set_target_properties( stp_cuda
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON
                      )
target_compile_features(stp_cuda PUBLIC cxx_std_17)

file(GLOB FILENAMES *.cpp)

foreach(filename ${FILENAMES})
  get_filename_component(basename ${filename} NAME_WE)
  add_executable(${basename} ${filename})
  target_link_libraries(${basename} PUBLIC stp)

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${basename}.cmake")
    include(${basename}.cmake)
  endif()
endforeach()

target_link_libraries(stp_simulator PUBLIC  stp_cuda)




    

    
    




