#pragma once

#include <numeric>

#include <Eigen/Dense>

using Eigen::MatrixXi;

namespace stp
{
  unsigned get_lcm( unsigned m, unsigned n )
  {
    return ( m * n ) / std::gcd( m , n);
  }
  
  MatrixXi generate_swap_matrix( unsigned m, unsigned n)
  {
    MatrixXi swap_matrixXi(m * n, m * n);
    swap_matrixXi = MatrixXi::Zero(m * n, m * n);
    int p, q;
    for (int i = 0; i < m * n / 2 + 1; i++)
    {
      p = i / m;
      q = i % m;
      int j = q * n + p;
      swap_matrixXi(i, j) = 1;
      swap_matrixXi(m * n - 1 - i, m * n - 1 - j) = 1;
    }
    return swap_matrixXi;
  }

  MatrixXi kronecker_product( MatrixXi A, MatrixXi B )
  {
    /* trivial cases */
    auto a_dimensions = A.rows() * A.cols();
    auto b_dimensions = B.rows() * B.cols();

    if( a_dimensions == 1u ) return B;
    if( b_dimensions == 1u ) return A;

    MatrixXi KP( A.rows() * B.rows(), A.cols() * B.cols() );

    for (int i = 0; i < A.rows(); ++i) 
    {
        for (int j = 0; j < A.cols(); ++j) 
        {
            KP.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
        }
    }

    return KP;
  }

  MatrixXi semi_tensor_product( MatrixXi A, MatrixXi B )
  {
    unsigned m = A.rows();
    unsigned n = A.cols();
    
    unsigned p = B.rows();
    unsigned q = B.cols();
  
    unsigned t = get_lcm( n, p ); 

    MatrixXi Ia = MatrixXi::Identity( t / n, t / n );
    MatrixXi Ib = MatrixXi::Identity( t / p, t / p );

    MatrixXi KPa = kronecker_product( A, Ia );
    MatrixXi KPb = kronecker_product( B, Ib );

    return KPa * KPb;
  }
}
