#include <catch.hpp>

#include <stp/stp_eigen.hpp>

using namespace stp;
using Eigen::MatrixXi;

TEST_CASE( "stp calculation using eigen", "[eigen]" )
{
    MatrixXi A(2, 4);
    A << 1, 0, 0, 0,
         0, 1, 1, 1;

    MatrixXi B(2, 4);
    B << 1, 1, 0, 1,
         0, 0, 1, 0;

    MatrixXi C = semi_tensor_product( A, B );

    CHECK( C.rows() == 2u );
    CHECK( C.cols() == 8u );
}

TEST_CASE( "swap matrix for column vector", "[eigen]" )
{
  /* A stp B = W (swap matrix) stp B stp A */
  MatrixXi A( 2, 1 );
  A << 1, 0;

  MatrixXi B( 2, 1 );
  B << 0, 1;

  MatrixXi C = generate_swap_matrix( 2, 2 );

  MatrixXi stp1 = semi_tensor_product( A, B );
  MatrixXi stp2 = semi_tensor_product( semi_tensor_product( C, B ), A );

  CHECK( stp1.rows() == stp2.rows() );
  CHECK( stp1.cols() == stp2.cols() );
  CHECK( stp1 == stp2 );
  
  MatrixXi stp3 = semi_tensor_product( B, A );
  MatrixXi stp4 = semi_tensor_product( semi_tensor_product( C, A ), B );

  CHECK( stp3.rows() == stp4.rows() );
  CHECK( stp3.cols() == stp4.cols() );
  CHECK( stp3 == stp4 );
}

TEST_CASE( "power reducing matrix", "[eigen]" )
{
  /* A stp A = PR_2 stp A, note that A is a Boolean value, cannot be a Matrix */
  MatrixXi A( 2, 1 );
  A << 1, 0;

  MatrixXi PR2( 4, 2 );
  PR2 << 1, 0,
         0, 0,
         0, 0,
         0, 1;

  MatrixXi stp1 = semi_tensor_product( A, A );
  MatrixXi stp2 = semi_tensor_product( PR2, A );

  CHECK( stp1 == stp2 );
}

/*
TEST_CASE( "expression to tt", "[eigen]" )
{
  std::string expr1 = "(a & b) | (a & ~c) | (~b & ~c)";
  
  std::vector<std::string> inputs_order1 = { "a", "b", "c" };
  MatrixXi mat1 = stp::from_exp_to_nmx( expr1, inputs_order1 );
  CHECK( stp::to_binary( mat1 ) == "10001011" );
  CHECK( stp::to_hex( mat1 ) == "8B" );
  
  std::vector<std::string> inputs_order2 = { "a", "c", "b" };
  MatrixXi mat2 = stp::from_exp_to_nmx( expr1, inputs_order2 );
  CHECK( stp::to_binary( mat2 ) == "10001101" );
  CHECK( stp::to_hex( mat2 ) == "8D" );
  
  std::vector<std::string> inputs_order3 = { "c", "b", "a" };
  MatrixXi mat3 = stp::from_exp_to_nmx( expr1, inputs_order3 );
  CHECK( stp::to_binary( mat3 ) == "11010001" );
  CHECK( stp::to_hex( mat3 ) == "D1" );
}*/
