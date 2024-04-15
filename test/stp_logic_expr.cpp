#include <catch.hpp>

#include <stp/stp_logic_expr.hpp>
#include <stp/stp_utils.hpp>

using namespace stp;
using matrix = Eigen::MatrixXi;


TEST_CASE( "expression to tt", "[expr]" )
{
  std::string expr1 = "(a & b) | (a & ~c) | (~b & ~c)";
  
  std::vector<std::string> inputs_order1 = { "a", "b", "c" };
  matrix mat1 = stp::from_exp_to_nmx( expr1, inputs_order1 );
  CHECK( stp::to_binary( mat1, true ) == "10001011" );
  CHECK( stp::to_hex( mat1, true ) == "8B" );
  
  std::vector<std::string> inputs_order2 = { "a", "c", "b" };
  matrix mat2 = stp::from_exp_to_nmx( expr1, inputs_order2 );
  CHECK( stp::to_binary( mat2, true ) == "10001101" );
  CHECK( stp::to_hex( mat2, true ) == "8D" );
  
  std::vector<std::string> inputs_order3 = { "c", "b", "a" };
  matrix mat3 = stp::from_exp_to_nmx( expr1, inputs_order3 );
  CHECK( stp::to_binary( mat3, true ) == "11010001" );
  CHECK( stp::to_hex( mat3, true ) == "D1" );
}

TEST_CASE( "STP expression to tt with different input orders", "[expr]" )
{
    std::string input = "m_d m_c x_1 x_2 m_d m_c x_1 m_n x_3 m_c m_n x_2 m_n x_3";
    
    std::vector<std::string> input_names1{ "x_1", "x_2", "x_3" };
    matrix mat1 = expr_normalize( input, input_names1 );
    CHECK( stp::to_binary( mat1, true ) == "10001011" );
    CHECK( stp::to_hex( mat1, true ) == "8B" );
    
    std::vector<std::string> input_names2{ "x_1", "x_3", "x_2" };
    matrix mat2 = expr_normalize( input, input_names2 );
    CHECK( stp::to_binary( mat2, true ) == "10001101" );
    CHECK( stp::to_hex( mat2, true ) == "8D" );
    
    std::vector<std::string> input_names3{ "x_3", "x_2", "x_1" };
    matrix mat3 = expr_normalize( input, input_names3 );
    CHECK( stp::to_binary( mat3, true ) == "11010001" );
    CHECK( stp::to_hex( mat3, true ) == "D1" );
}

TEST_CASE( "STP expression to tt with more expresssions", "[expr]" )
{
    std::string input = "m_c m_d x_1 x_2 m_d x_3 x_4";
    std::vector<std::string> input_names{ "x_1", "x_2", "x_3", "x_4" };
    matrix mat1 = expr_normalize( input, input_names, true );
    CHECK( stp::to_binary( mat1 ) == "1110111011100000" );
    CHECK( stp::to_hex( mat1 ) == "EEE0" );
}
