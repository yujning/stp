#pragma once

#include <iostream>

#include <chrono>
#include "stp_utils.hpp"

#ifdef ENABLE_CUDA 


extern "C"  
std::vector<int32_t> cuda_In_KR_Matrix(int32_t dim, std::vector<int32_t>& A);

extern "C"  
void Get_Total_Thread_Num(void);

extern "C"  
std::vector<int32_t> cuda_semi_tensor_product(std::vector<int32_t>& A, std::vector<int32_t>& B);

#endif

