#pragma once

#include <Eigen/Dense>

#include <iostream> 
#include <numeric>
#include <bitset>

using matrix = Eigen::MatrixXi;           // Defines the type of matrix to use
using matrix_chain = std::vector<matrix>; // Defined matrix chain


namespace stp
{

  unsigned get_lcm( unsigned m, unsigned n )
  {
    return ( m * n ) / std::gcd( m, n );
  }

  void print_binary( const matrix& mat, std::ostream& os = std::cout )
  {
    for( auto i = mat.cols() - 1; i >= 0; i-- )
    {
      os << mat( 0, i );
    }
  }

  std::string to_binary( const matrix& mat )
  {
    std::stringstream st;
    print_binary( mat, st );
    return st.str();
  }

  void print_hex( const std::string& binary_string, std::ostream& os = std::cout )
  {
    assert( binary_string.length() % 4 == 0 );   

    for ( size_t i = 0; i < binary_string.length(); i += 4 ) 
    {
      std::string block = binary_string.substr( i, 4 );

      int decimal_value = std::bitset<4>(block).to_ulong(); 
      char hex_digit;

      if (decimal_value < 10) 
      {
        hex_digit = '0' + decimal_value; 
      } 
      else 
      {
        hex_digit = 'A' + ( decimal_value - 10 );
      }

      os << hex_digit;
    }
  }

  std::string to_hex( const matrix& mat )
  {
    std::stringstream st;
    print_hex( to_binary( mat ), st );
    return st.str();
  }

  matrix generate_swap_matrix( const int& m, const int& n )
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

  matrix kronecker_product( const matrix& A, const matrix& B )
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

  matrix semi_tensor_product( const matrix& A, const matrix& B )
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

  matrix matrix_chain_multiply( const matrix_chain& mc )
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

