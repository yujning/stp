#include <iostream>
#include <stp/stp_eigen.hpp>
#include <stp/stp_logic_expr.hpp>

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

  std::cout << "----------------------------------------------------------------------------" << std::endl;
  
  std::string expr1 = "(a & b) | (a & ~c) | (~b & ~c)";
  std::vector<std::string> inputs_order1 = { "a", "b", "c" };
  MatrixXi mat1 = stp::from_exp_to_nmx( expr1, inputs_order1, false );
  std::cout << "The normal matrix for expression " << expr1 << " is " << std::endl;
  std::cout << mat1 << std::endl;
  std::cout << "The binary truth table is: " << stp::to_binary( mat1 ) << std::endl;
  std::cout << "The hex truth table is: 0X" << stp::to_hex( mat1 ) << std::endl;
  
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  
  std::string expr2 = "(a & b & ~c) | (a & ~c) | (~b & ~c)";
  std::vector<std::string> inputs_order2 = { "a", "b", "c" };
  MatrixXi mat2 = stp::from_exp_to_nmx( expr2, inputs_order2, false );
  std::cout << "The normal matrix for expression " << expr2 << " is " << std::endl;
  std::cout << mat2 << std::endl;
  std::cout << "The binary truth table is: " << stp::to_binary( mat2 ) << std::endl;
  std::cout << "The hex truth table is: 0X" << stp::to_hex( mat2 ) << std::endl;
  
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  stp::test();
}
