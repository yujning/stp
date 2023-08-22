#pragma once

#include <numeric>

#include <Eigen/Dense>

using Eigen::MatrixXd;

namespace stp
{
  unsigned get_lcm( unsigned m, unsigned n )
  {
    return ( m * n ) / std::gcd( m , n);
  }

  MatrixXd kronecker_product( MatrixXd A, MatrixXd B )
  {
    /* trivial cases */
    auto a_dimensions = A.rows() * A.cols();
    auto b_dimensions = B.rows() * B.cols();

    if( a_dimensions == 1u ) return B;
    if( b_dimensions == 1u ) return A;

    MatrixXd KP( A.rows() * B.rows(), A.cols() * B.cols() );

    for (int i = 0; i < A.rows(); ++i) 
    {
        for (int j = 0; j < A.cols(); ++j) 
        {
            KP.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
        }
    }

    return KP;
  }

  MatrixXd stp_calculation( MatrixXd A, MatrixXd B )
  {
    unsigned m = A.rows();
    unsigned n = A.cols();
    
    unsigned p = B.rows();
    unsigned q = B.cols();
  
    unsigned t = get_lcm( n, p ); 

    MatrixXd Ia = MatrixXd::Identity( t / n, t / n );
    MatrixXd Ib = MatrixXd::Identity( t / p, t / p );

    MatrixXd KPa = kronecker_product( A, Ia );
    MatrixXd KPb = kronecker_product( B, Ib );

    return KPa * KPb;
  }
}
