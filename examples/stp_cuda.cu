#include <iostream>
#include <chrono>


#include <stp/stp_cuda.hpp>
#include <stp/stp_eigen.hpp>
#include <stp/stp_utils.hpp>

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuda_runtime_api.h>
#include <cusparse.h>

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"

uint64_t Total_Thread=0; //total supported threads


extern "C"
//get total thread number
void Get_Total_Thread_Num(void)
{
    int deviceCount;
    //get the number of CUDA devices
    cudaGetDeviceCount(&deviceCount); 

    for (int i = 0; i < deviceCount; ++i) 
    {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, i);

        std::cout << "Device " << i << ": " << deviceProp.name << std::endl;
        std::cout << "Max threads per block: " << deviceProp.maxThreadsPerBlock << std::endl;
        std::cout << "Max blocks per grid: " << (int64_t)deviceProp.maxGridSize[0] * deviceProp.maxGridSize[1] * deviceProp.maxGridSize[2] << std::endl;
        std::cout << "  X : " << deviceProp.maxGridSize[0] << std::endl;
        std::cout << "  Y : " << deviceProp.maxGridSize[1] << std::endl;
        std::cout << "  Z : " << deviceProp.maxGridSize[2] << std::endl;
        std::cout << "Total max threads: " << deviceProp.maxThreadsPerMultiProcessor * deviceProp.multiProcessorCount << std::endl;
        std::cout << std::endl;
        Total_Thread = deviceProp.maxThreadsPerMultiProcessor * deviceProp.multiProcessorCount;
    }    
}

extern "C"
//matrix encoding
std::vector<int32_t> Matrix_to_Vec(matrix& A)
{
    int32_t size = A.cols() + 1;
    std::vector<int32_t> Result(size);
    Result[0] = A.rows();
    //operate on each column
    for(int32_t i = 0; i < A.cols(); i++)
    {
        for(int32_t j=0;j<A.rows();j++)
        {
            //non-zero element
            if(A(j,i))
            {
                Result[i+1] = j;
                break;
            }
        }
    }
    return Result;
}


extern "C"
//matrix chain encoding
std::vector<std::vector<int32_t>> Matrix_to_Vec_chain(matrix_chain& A)
{
    std::vector<std::vector<int32_t>> result;
    for(int32_t i=0; i<A.size();i++)
    {
        result.push_back(Matrix_to_Vec(A[i]));
    }
    return result;

}


//In_KR_Matrix_Kernel     
__global__ void In_KR_Matrix_Kernel(int32_t sub_dim, int32_t idx_offset, int32_t *A, int32_t A_val_len, int32_t *C, int32_t C_val_len)
{
    int32_t ix = blockIdx.x * 1024 + threadIdx.x;
    //index
    int32_t idx = ix + idx_offset;

    int32_t x_code = idx / A_val_len;

    int32_t y_code = idx % A_val_len;

    //boundary check  
    if(idx < C_val_len)
    {    
        //XP+Y
        C[idx] = x_code * A[0] + A[y_code + 1]; 
    }
}


extern "C"
//In_KR_Matrix
std::vector<int32_t> In_KR_Matrix(int32_t dim, std::vector<int32_t>& A)
{
    //Get the dimensions of matrix A
    int32_t A_row = A[0];
    int32_t A_col = A.size() - 1;

    //Calculate the size of result matrix
    int32_t C_len = dim * A_col + 1;

    //Define the result matrix
    std::vector<int32_t> C(C_len);

    //Assign the number of rows of result matrix
    C[0] = A_row * dim;

    if(C_len <= 20000)
    {
        for (int32_t i = 0; i < dim; i++)
        {
            int32_t temp = i * A_row;
            int32_t idx = i * A_col + 1;
            for (int32_t j = 0; j < A_col; j++)
            {
                C[idx + j] = temp + A[j + 1];
            }
        }
        return C;
    }

    //std::cout << "In_KR_Matrix threads > 1025"<< std::endl;

    int32_t *d_A,*d_C;

    //compute space 
    size_t size_A = A.size() * sizeof(int32_t);
    size_t size_C = (C_len - 1) * sizeof(int32_t);

    //allocate memory
    cudaError_t errA, errC;
    errA = cudaMalloc((void **)&d_A, size_A);
    //Error checking
    if (errA != cudaSuccess) 
    {
        std::cerr << "Error allocating memory for In_KR_Matrix d_A: " << cudaGetErrorString(errA) << std::endl;

        std::vector<int32_t> Temp;
        //Set the size to improve the speed of computation
        Temp.resize(1);
        Temp[0] = -1;

        return Temp;
    }
    errC = cudaMalloc((void **)&d_C, size_C);
    if (errC != cudaSuccess) 
    {
        std::cerr << "Error allocating memory for In_KR_Matrix d_C: " << cudaGetErrorString(errC) << std::endl;
        std::vector<int32_t> Temp;
        //Set the size to improve the speed of computation
        Temp.resize(1);
        Temp[0] = -1;

        return Temp;
    }

    //Copy parameters
    cudaMemcpy(d_A, A.data(), size_A, cudaMemcpyHostToDevice);

    //Calculate the block size (maximum 1024)
    dim3 threadsPerBlock(1024, 1);

    //Can be done in one go (with each element of C as a thread)
    if((C_len - 1) <=Total_Thread)
    {
        //Calculate the grid size (for large scale)
        dim3 numBlocks0(( C_len - 1 + 1024 -1 ) / 1024, 1);

        //Launch GPU
        In_KR_Matrix_Kernel<<<numBlocks0, threadsPerBlock>>>(dim, 0, d_A, A_col , d_C, C_len - 1);
        cudaDeviceSynchronize(); //Wait for the kernel to complete
    }
    //Divide into blocks (by column)
    else
    {
        int32_t remain_num = C_len - 1; //remaining unassigned columns in C
        int32_t idx_offset = 0;  //Thread offset

        while(remain_num)
        {
            //the last  
            if(remain_num <= Total_Thread)
            {
                //Calculate the grid size (for large scale)
                dim3 numBlocks(( remain_num + 1024 -1 ) / 1024, 1);

                In_KR_Matrix_Kernel<<<numBlocks, threadsPerBlock>>>(dim, idx_offset, d_A, A_col, d_C, C_len - 1); 
                cudaDeviceSynchronize(); //Wait for the kernel to complete
                //Calculate the thread offset
                idx_offset += remain_num;
                remain_num = 0; //exit the loop
            }
            // Total thread
            else
            {
                //Calculate the grid size (for large scale)
                dim3 numBlocks(( Total_Thread + 1024 -1 ) / 1024, 1);

                In_KR_Matrix_Kernel<<<numBlocks, threadsPerBlock>>>(dim, idx_offset, d_A, A_col, d_C, C_len - 1);  
                cudaDeviceSynchronize(); //Wait for the kernel to complete

                //Calculate thread offset
                idx_offset += Total_Thread;
                remain_num -= Total_Thread; //exit the loop                
            }
        }     

    }

    //Segmented copy with pointer offset
    cudaMemcpy(C.data() + 1, d_C, size_C, cudaMemcpyDeviceToHost);

    //Free resources
    cudaFree(d_A);
    cudaFree(d_C);

    return C;
}



//Matrix_KR_In_Kernel
__global__ void Matrix_KR_In_Kernel(int32_t dim, int32_t idx_offset, int32_t *A, int32_t A_val_len, int32_t *C, int32_t C_val_len)
{
    int32_t ix = blockIdx.x * 1024 + threadIdx.x;
    //index
    int32_t idx = ix + idx_offset;

    //calculate x_code
    int32_t x_code = idx / dim;

    //calculate y_code
    int32_t y_code = idx % dim;

    //boundary check
    if(idx < C_val_len)
    {    
        //xp+y
        C[idx] = A[x_code] * dim + y_code;
    }
    
}

//Matrix_KR_In
std::vector<int32_t> Matrix_KR_In(int32_t dim,  std::vector<int32_t>& A)
{
    //get dimensions of matrix A
    int32_t A_row = A[0];
    int32_t A_col = A.size() - 1;

    //calculate size of result matrix
    //int32_t C_len = A_row * dim + 1;
    int32_t C_len = A_col * dim + 1;

    std::vector<int32_t> C(C_len);

    //assign number of rows of result matrix
    C[0] = A_row * dim;

    if(C_len <= 20000)
    {
        for (int32_t i = 0; i < A_col; i++)
        {
            int32_t temp = A[i + 1] * dim;
            int32_t idx = i * dim + 1;
            for (int32_t j = 0; j < dim; j++)
            {
                C[idx + j] = temp + j;
            }
        }
        return C;

    }

    //std::cout << "Matrix_KR_In threads > 1025"<< std::endl;

    int32_t *d_A,*d_C;
    //compute space
    size_t size_A = A.size() * sizeof(int32_t);
    size_t size_C = ( C_len - 1 ) * sizeof(int32_t);

    //allocate memory
    cudaError_t errA, errC;

    errA = cudaMalloc((void **)&d_A, size_A);
    //error checking
    if (errA != cudaSuccess) 
    {
        std::cerr << "Error allocating memory for Matrix_KR_In d_A: " << cudaGetErrorString(errA) << std::endl;
        std::vector<int32_t> Temp;
        //set size to improve computation speed
        Temp.resize(1);
        Temp[0] = -1;

        return Temp;
    }
    errC = cudaMalloc((void **)&d_C, size_C);
    if (errC != cudaSuccess) 
    {
        std::cerr << "Error allocating memory for Matrix_KR_In d_C: " << cudaGetErrorString(errC) << std::endl;
        std::vector<int32_t> Temp;
        //set size to improve computation speed
        Temp.resize(1);
        Temp[0] = -1;

        return Temp;
    }

    //copy parameters
    cudaMemcpy(d_A, A.data(), size_A, cudaMemcpyHostToDevice);

    //calculate block size (maximum 1024)
    dim3 threadsPerBlock(1024, 1);

    //can be done in one go (with each element of C as a thread)
    if((C_len - 1) <= Total_Thread)
    {
        //calculate grid size (for large scale)
        dim3 numBlocks(( C_len - 1 + 1024 -1 ) / 1024, 1);

        Matrix_KR_In_Kernel<<<numBlocks, threadsPerBlock>>>(dim, 0, d_A + 1, A_col, d_C, C_len - 1); 
        cudaDeviceSynchronize(); //wait for the kernel to complete
    }
    //divide into blocks (by column)
    else
    {
        int32_t remain_num = C_len - 1; //C remaining unassigned columns
        int32_t idx_offset = 0;  //thread offset

        while(remain_num)
        {
            //last time  quantity remain_num
            if(remain_num <= Total_Thread)
            {
                //calculate grid size (for large scale) 
                dim3 numBlocks(( remain_num + 1024 -1 ) / 1024, 1);

                Matrix_KR_In_Kernel<<<numBlocks, threadsPerBlock>>>(dim, idx_offset, d_A + 1, A_col, d_C, C_len - 1);  
                cudaDeviceSynchronize(); //wait for the kernel to complete
                //calculate thread offset
                idx_offset += remain_num;
                remain_num = 0; //exit the loop
            }
            //total thread
            else
            {
                //calculate grid size (for large scale)
                dim3 numBlocks(( Total_Thread + 1024 -1 ) / 1024, 1);

                Matrix_KR_In_Kernel<<<numBlocks, threadsPerBlock>>>(dim, idx_offset, d_A + 1, A_col, d_C, C_len - 1);  
                cudaDeviceSynchronize(); //wait for the kernel to complete

                //calculate thread offset
                idx_offset += Total_Thread;
                remain_num -= Total_Thread; //exit the loop                
            }
        }
    }

    //segmented copy with pointer offset
    cudaMemcpy(C.data() + 1, d_C, size_C, cudaMemcpyDeviceToHost);

    //free resources
    cudaFree(d_A);
    cudaFree(d_C);

    return C;
}         


//Matrix_Multipiy_Kernel
__global__ void Matrix_Multipiy_Kernel(int32_t idx_offset, int32_t *A, int32_t A_val_len, int32_t *B, int32_t B_val_len, int32_t *C, int32_t C_val_len, int32_t t)
{
    int32_t ix = blockIdx.x * 1024 + threadIdx.x;
    //index
    int32_t idx = ix + idx_offset;

    //calculate x_code
    int32_t x_code = idx % t ;

    //calculate y_code
    int32_t y_code = idx / t;

    //boundary check
    if(idx < C_val_len )
    {

        //result
        C[idx] = A[ B[ y_code + 1 ] * t + x_code + 1];
    }
}

extern "C"
//my_semi_tensor_product
std::vector<int32_t> my_semi_tensor_product(std::vector<int32_t>& A, std::vector<int32_t>& B)
{
    //get dimensions of matrix A and B
    int32_t A_row = A[0];
    int32_t A_col = A.size() - 1;
    int32_t B_row = B[0];
    int32_t B_col = B.size() - 1;

    if(A_col % B_row == 0)
    {
        //calculate size of result matrix
        int32_t C_len = (int64_t)A_col * B_col / B_row + 1;

        if(C_len <= 20000)
        {
            std::vector<int32_t> C(C_len);

            C[0] = A_row;
            int32_t times = A_col / B_row;

            for (int32_t i = 0; i < B_col; i++)
            {
                //t = n/p
                for(int32_t j = 0; j < times; j++)
                {
                    C[times * i + j + 1] = A[1 + B[i+1] * times + j] ;
                }
            }

            return C;
        }

        //caculate size of result matrix
        size_t size_A = A.size() * sizeof(int32_t);
        size_t size_B = B.size() * sizeof(int32_t);
        size_t size_C = (C_len - 1) * sizeof(int32_t);       

        int32_t *d_A,*d_B,*d_C;

        //allocate memory
        cudaError_t errA, errB, errC;
        errA = cudaMalloc((void **)&d_A, size_A);

        //error checking
        if (errA != cudaSuccess) 
        {
            std::cerr << "Error allocating memory for my_semi_tensor_product d_A: " << cudaGetErrorString(errA) << std::endl;
            std::vector<int32_t> Temp;
            //set size to improve computation speed
            Temp.resize(1);
            Temp[0] = -1;

            return Temp;
        }
        errB = cudaMalloc((void **)&d_B, size_B);
        //error checking
        if (errB != cudaSuccess) 
        {
            std::cerr << "Error allocating memory for my_semi_tensor_product d_B: " << cudaGetErrorString(errB) << std::endl;
            std::vector<int32_t> Temp;
            //set size to improve computation speed
            Temp.resize(1);
            Temp[0] = -1;

            return Temp;
        }
        errC = cudaMalloc((void **)&d_C, size_C);
        if (errC != cudaSuccess) 
        {
            std::cerr << "Error allocating memory for my_semi_tensor_product d_C: " << cudaGetErrorString(errC) << std::endl;
            std::vector<int32_t> Temp;
            //set size to improve computation speed
            Temp.resize(1);
            Temp[0] = -1;

            return Temp;
        }

        //copy parameters
        cudaMemcpy(d_A, A.data(), size_A, cudaMemcpyHostToDevice);
        cudaMemcpy(d_B, B.data(), size_B, cudaMemcpyHostToDevice);

        //calculate block size (maximum 1024)
        dim3 threadsPerBlock(1024, 1);

        if((C_len - 1) <= Total_Thread)
        {
            //calculate grid size (for large scale)
            dim3 numBlocks((C_len - 1 + 1024 - 2 ) / 1024, 1);
            Matrix_Multipiy_Kernel<<<numBlocks, threadsPerBlock>>>(0, d_A, A_col, d_B, B_col, d_C, C_len - 1, A_col / B_row);

            //wait for all threads to complete
            cudaDeviceSynchronize();          
        }
        else
        {
            int32_t remain_num = C_len - 1; //remaining unassigned columns in C
            int32_t idx_offset = 0;  //thread offset

            while(remain_num)
            {
                //the last  
                if(remain_num <= Total_Thread)
                {
                    //calculate grid size (for large scale) 
                    dim3 numBlocks(( remain_num + 1024 -1 ) / 1024, 1);

                    Matrix_Multipiy_Kernel<<<numBlocks, threadsPerBlock>>>(idx_offset, d_A, A_col, d_B, B_col, d_C, C_len - 1, A_col / B_row);
                    cudaDeviceSynchronize(); //wait for the kernel to complete
                    //calculate thread offset
                    idx_offset += remain_num;
                    remain_num = 0; //exit the loop
                }
                //total thread
                else
                {
                    //calculate grid size (for large scale)
                    dim3 numBlocks(( Total_Thread + 1024 -1 ) / 1024, 1);

                    Matrix_Multipiy_Kernel<<<numBlocks, threadsPerBlock>>>(idx_offset, d_A, A_col, d_B, B_col, d_C, C_len - 1, A_col / B_row);
                    cudaDeviceSynchronize(); //wait for the kernel to complete
                    //calculate thread offset
                    //idx_offset += Total_Thread;
                    cudaDeviceSynchronize(); //wait for the kernel to complete
                    //计算线程偏移
                    idx_offset += Total_Thread;
                    remain_num -= Total_Thread; //exit the loop               
                }
            }            
        }

        std::vector<int32_t> C;
        //set size to improve computation speed
        C.resize(C_len);

        //assign result matrix row number
        C[0] = A_row;

        //segmented copy with pointer offset
        cudaMemcpy(C.data() + 1, d_C, size_C, cudaMemcpyDeviceToHost);
        //free resources
        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);

        return C;
    }
    else if(B_row % A_col == 0)
    {
        std::vector<int32_t> temp = Matrix_KR_In(B_row / A_col, A);
        std::vector<int32_t> C = my_semi_tensor_product(temp, B);

        return C;
    }
    else
    {
        //error
        std::cout << "Error" << std::endl;
        std::vector<int32_t> C;
        //set size to improve computation speed
        C.resize(1);
        C[0] = -1;

        return C;
    }
    
}




extern "C"
//my_chain_multiply_by_multi_thread                                                
std::vector<int32_t> my_chain_multiply_by_multi_thread( std::vector<std::vector<int32_t>>& mc, bool verbose)
{
  std::vector<int32_t> result=mc[0];

  if(mc.size()<2)
  {
    return result;
  }

  for (size_t i = 1; i < mc.size(); i++ )
  {
    result = my_semi_tensor_product(result,mc[i]);

  }
  return result;
}




