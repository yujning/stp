#pragma once

#include <Eigen/Dense>

#include <iostream> 
#include <numeric>

using matrix = Eigen::MatrixXi;           // Defines the type of matrix to use
using matrix_chain = std::vector<matrix>; // Defined matrix chain


namespace stp
{

  inline unsigned get_lcm( unsigned m, unsigned n )
  {
    return ( m * n ) / std::gcd( m, n );
  }

  inline matrix generate_swap_matrix( const int& m, const int& n )
  {
    matrix swap_matrixXi = matrix::Zero( m * n, m * n );
    int p, q;
    for ( int i = 0; i < m * n / 2 + 1; i++ )
    {
      p = i / m;
      q = i % m;
      int j = q * n + p;
      swap_matrixXi( i, j ) = 1;
      swap_matrixXi( m * n - 1 - i, m * n - 1 - j ) = 1;
    }
    return swap_matrixXi;
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

    matrix KP( A.rows() * B.rows(), A.cols() * B.cols() );

    for ( int i = 0; i < A.rows(); ++i )
    {
      for ( int j = 0; j < A.cols(); ++j )
      {
        KP.block( i * B.rows(), j * B.cols(), B.rows(), B.cols() ) = A( i, j ) * B;
      }
    }

    return KP;
  }

  inline matrix semi_tensor_product( const matrix& A, const matrix& B )
  {
    unsigned m = A.rows();
    unsigned n = A.cols();

    unsigned p = B.rows();
    unsigned q = B.cols();

    unsigned t = get_lcm( n, p );

    matrix Ia = matrix::Identity( t / n, t / n );
    matrix Ib = matrix::Identity( t / p, t / p );

    matrix KPa = kronecker_product( A, Ia );
    matrix KPb = kronecker_product( B, Ib );

    return KPa * KPb;
  }

  inline matrix matrix_chain_multiply( const matrix_chain& mc )
  {
    assert( mc.size() > 0 );
    matrix result_matrix;

    if ( mc.size() == 1 )
    {
      return mc[0];
    }

    result_matrix = semi_tensor_product( mc[0], mc[1] );

    for ( int i = 2; i < mc.size(); i++ )
    {
      result_matrix = semi_tensor_product( result_matrix, mc[i] );
    }
    return result_matrix;
  }

} // namespace stp

