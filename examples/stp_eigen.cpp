#include <iostream>
#include <stp/stp_eigen.hpp>

using Eigen::MatrixXi;

int main()
{
  MatrixXi I( 2, 2 );
  I << 1, 0,
      0, 1;

  MatrixXi Theta1( 2, 1 );
  Theta1 << 1, 0;

  MatrixXi Theta2( 2, 1 );
  Theta1 << 0, 1;

  MatrixXi A = stp::kronecker_product( I, Theta1 );
  MatrixXi B = stp::kronecker_product( I, Theta2 );

  assert( A.rows() == B.rows() );

  MatrixXi W( A.rows(), A.cols() + B.cols() );
  W << A, B; // Merge

  std::cout << "W[2,4] is \n"
            << stp::generate_swap_matrix( 2, 4 ) << std::endl;
  std::cout << "W[8,2] is \n"
            << stp::generate_swap_matrix( 8, 2 ) << std::endl;
}
