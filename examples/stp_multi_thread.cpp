#include <iostream>
#include <stp/stp_logic_expr.hpp>
#include <stp/stp_multi_thread.hpp>
#include <stp/stp_eigen.hpp>

using matrix = Eigen::MatrixXi;

int main()
{
  std::string expr = "m_d m_c x_1 x_2 m_d m_c x_1 m_n x_3 m_c m_n x_2 m_n x_3";
  std::vector<std::string> input_names{ "x_1", "x_2", "x_3" };
  auto mc = stp::expr_normalize_to_chain( expr, input_names, true );

  auto mat = stp::matrix_chain_multiply_by_multi_thread( mc, 10, true);
  
  std::cout << "The normal matrix for expression " << expr << " is " << std::endl;
  std::cout << mat << std::endl;
  std::cout << "The binary truth table is: " << stp::to_binary( mat ) << std::endl;
  std::cout << "The hex truth table is: 0X" << stp::to_hex( mat ) << std::endl;
  return 0;
}
