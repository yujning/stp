#pragma once
#include <iostream>
#include <chrono>
#include <vector>
#include <cstdint>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

#include "excute_cuda.hpp"

// matrix encoding
std::vector<stp_data> Matrix_to_Vec(matrix &A)
{
    stp_data size = A.cols() + 1;
    std::vector<stp_data> Result(size);
    Result[0] = A.rows();
    // operate on each column
    for (stp_data i = 0; i < A.cols(); i++)
    {
        for (stp_data j = 0; j < A.rows(); j++)
        {
            // non-zero element
            if (A(j, i))
            {
                Result[i + 1] = j;
                break;
            }
        }
    }
    return Result;
}

// matrix chain encoding
std::vector<std::vector<int32_t>> Matrix_to_Vec_chain(matrix_chain &A)
{
    std::vector<std::vector<int32_t>> result;
    for (int32_t i = 0; i < A.size(); i++)
    {
        result.push_back(Matrix_to_Vec(A[i]));
    }
    return result;
}

#pragma region excute with vector

//generate Mn    (k - 1) - p
inline matrix generate_Mn_matrix( stp_data k)
{
    matrix Mn_matrix = matrix::Zero(k,k);
    for(size_t i = 0; i < k; i++)
    {
        Mn_matrix(k - 1 - i, i) = 1;
    }
    return Mn_matrix;
}


//generate Mc min(P,Q)
inline matrix generate_Mc_matrix( stp_data k = 2)
{
    matrix Mc_matrix = matrix::Zero(k, k*k);

    for (int32_t i = k - 1; i >= 0; i--)
    {
        for (int32_t j = k - 1; j >= 0; j--)
        {
            if(j < i)
            {
                Mc_matrix( k - 1 - j, (k - j -1) + (k - i -1) * k ) = 1;
            }
            else
            {
                Mc_matrix( k - 1 - i, (k - j - 1) + (k - i - 1) * k ) = 1;
            }
        }
    }
    return Mc_matrix;
}


// generate Md max(P,Q)
inline matrix generate_Md_matrix(stp_data k = 2)
{
    matrix Md_matrix = matrix::Zero(k, k * k);
    for (stp_data i = 0; i < k; i++)
    {
        for (stp_data j = 0; j < k; j++)
        {
            if (j > i)
            {
                Md_matrix( k - 1 - j, (k - j -1) + (k - i -1) * k ) = 1;
            }
            else
            {
                Md_matrix( k - 1 - i, (k - j -1) + (k - i -1) * k ) = 1;
            }
        }
    }
    return Md_matrix;
}




inline matrix generate_swap_matrix(const stp_data &m, const stp_data &n)
{
    matrix swap_matrixXi = matrix::Zero(m * n, m * n);
    stp_data p, q;
    for (stp_data i = 0; i < m * n / 2 + 1; i++)
    {
        p = i / m;
        q = i % m;
        stp_data j = q * n + p;
        swap_matrixXi(i, j) = 1;
        swap_matrixXi(m * n - 1 - i, m * n - 1 - j) = 1;
    }
    return swap_matrixXi;
}

inline std::vector<stp_data> generate_swap_vec(const stp_data &m, const stp_data &n)
{
    stp_data dim = m * n;
    std::vector<stp_data> swap_matrix(dim + 1);
    swap_matrix[0] = dim;
    stp_data I, J;

    for (size_t i = 0; i < dim; ++i)
    {
        I = i / n;
        J = i % n;

        swap_matrix[i + 1] = J * m + I;
    }

    return swap_matrix;
}

inline matrix generate_Mr_matrix(const stp_data &k = 2)
{
    stp_data dim = k;
    matrix Mr_matrix = matrix::Zero(dim*dim, dim);

    for (stp_data i = 0; i < dim; i++)
    {
        Mr_matrix((i+1)*(i+1), i) = 1;
    }

    return Mr_matrix;
}


inline std::vector<stp_data> generate_Mr_vec(const stp_data &k)
{
    stp_data dim = k;
    std::vector<stp_data> Mr_matrix(dim + 1);
    Mr_matrix[0] = dim * dim;

    for (size_t i = 1; i < dim + 1; ++i)
    {
        Mr_matrix[i] = (i - 1) * (k + 1);
    }

    return Mr_matrix;
}


// In_KR_Matrix
inline std::vector<stp_data> In_KR_Vec(stp_data dim, std::vector<stp_data> &A)
{
    // Get the dimensions of matrix A
    stp_data A_row = A[0];
    stp_data A_col = A.size() - 1;

    // Calculate the size of result matrix
    stp_data C_len = dim * A_col + 1;

    // Define the result matrix
    std::vector<stp_data> C(C_len);

    // Assign the number of rows of result matrix
    C[0] = A_row * dim;

    for (stp_data i = 0; i < dim; i++)
    {
        stp_data temp = i * A_row;
        stp_data idx = i * A_col + 1;
        for (stp_data j = 0; j < A_col; j++)
        {
            C[idx + j] = temp + A[j + 1];
        }
    }
    return C;
}

// Matrix_KR_In
std::vector<stp_data> Vec_KR_In(stp_data dim, std::vector<stp_data> &A)
{
    // get dimensions of matrix A
    int32_t A_row = A[0];
    int32_t A_col = A.size() - 1;

    // calculate size of result matrix
    // int32_t C_len = A_row * dim + 1;
    int32_t C_len = A_col * dim + 1;

    std::vector<int32_t> C(C_len);

    // assign number of rows of result matrix
    C[0] = A_row * dim;

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


std::vector<stp_data> Vec_semi_tensor_product(std::vector<stp_data> &A, std::vector<stp_data> &B)
{
    // get dimensions of matrix A and B
    int32_t A_row = A[0];
    int32_t A_col = A.size() - 1;
    int32_t B_row = B[0];
    int32_t B_col = B.size() - 1;

    if (A_col % B_row == 0)
    {
        // calculate size of result matrix
        int32_t C_len = (int64_t)A_col * B_col / B_row + 1;

        std::vector<int32_t> C(C_len);

        C[0] = A_row;
        int32_t times = A_col / B_row;

        for (int32_t i = 0; i < B_col; i++)
        {
            // t = n/p
            for (int32_t j = 0; j < times; j++)
            {
                C[times * i + j + 1] = A[1 + B[i + 1] * times + j];
            }
        }

        return C;
    }
    else if (B_row % A_col == 0)
    {
        std::vector<int32_t> temp = Vec_KR_In(B_row / A_col, A);
        std::vector<int32_t> C = Vec_semi_tensor_product(temp, B);

        return C;
    }
    else
    {
        // error
        std::cout << "Error" << std::endl;
        std::vector<int32_t> C;
        // set size to improve computation speed
        C.resize(1);
        C[0] = -1;

        return C;
    }
}


std::vector<stp_data> Vec_chain_multiply(std::vector<std::vector<stp_data>> &mc, bool verbose)
{
    std::vector<int32_t> result = mc[0];

    if (mc.size() < 2)
    {
        return result;
    }

    for (size_t i = 1; i < mc.size(); i++)
    {
        result = Vec_semi_tensor_product(result, mc[i]);
    }
    return result;
}

#pragma endregion

#pragma region excute with matrix

// Matrix_KR_In
matrix Matrix_KR_In(int32_t dim, matrix &A)
{
    matrix Ia = matrix::Identity(dim, dim);
    matrix KPa = Eigen::kroneckerProduct(A, Ia);

    return KPa;
}

// Matrix_KR_In
matrix In_KR_Matrix(int32_t dim, const matrix &A)
{
    matrix Ia = matrix::Identity(dim, dim);
    matrix KPa = Eigen::kroneckerProduct(Ia, A);

    return KPa;
}

inline matrix kronecker_product( const matrix& A, const matrix& B )
{
  /* trivial cases */
  auto a_dimensions = A.rows() * A.cols();
  auto b_dimensions = B.rows() * B.cols();

  if ( a_dimensions == 1u )
    {
      return B;
    }
  if ( b_dimensions == 1u )
    {
      return A;
    }

  Eigen::SparseMatrix<int> sparse_A = A.sparseView();
  Eigen::SparseMatrix<int> sparse_B = B.sparseView();

  Eigen::SparseMatrix<int> KP = Eigen::kroneckerProduct( sparse_A, sparse_B );

  matrix result( KP );


  return result;
}


matrix Matrix_semi_tensor_product(const matrix &A, const matrix &B)
{
    int32_t m = A.rows();
    int32_t n = A.cols();
    int32_t p = B.rows();
    int32_t q = B.cols();

    unsigned t = stp::get_lcm( n, p );

    matrix Ia = matrix::Identity( t / n, t / n );
    matrix Ib = matrix::Identity( t / p, t / p );

    matrix KPa = Eigen::kroneckerProduct(A, Ia);
    matrix KPb = Eigen::kroneckerProduct(B, Ib);

    matrix result = KPa * KPb;

    return result;
}

matrix Matrix_chain_multiply(std::vector<matrix> &mc, bool verbose = false)
{
    matrix result = mc[0];

    if (mc.size() < 2)
    {
        return result;
    }

    for (size_t i = 1; i < mc.size(); i++)
    {
        result = Matrix_semi_tensor_product(result, mc[i]);
    }
    return result;
}

#pragma endregion
