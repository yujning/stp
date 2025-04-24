
#ifndef STP_CUDA_H
#define STP_CUDA_H

#include <iostream>




#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Sparse>



extern "C"
std::vector<std::vector<int32_t>> Matrix_to_Vec_chain(std::vector<Eigen::MatrixXi>& A);

extern "C"                                                                              
std::vector<int32_t> my_chain_multiply_by_multi_thread( std::vector<std::vector<int32_t>>& mc, bool verbose = false);

extern "C"  
std::vector<int32_t> Matrix_to_Vec(Eigen::MatrixXi& A);

extern "C"  
std::vector<int32_t> In_KR_Matrix(int32_t dim, std::vector<int32_t>& A);

extern "C"  
void Get_Total_Thread_Num(void);

extern "C"  
std::vector<int32_t> my_semi_tensor_product(std::vector<int32_t>& A, std::vector<int32_t>& B);

#endif