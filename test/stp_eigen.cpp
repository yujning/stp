#include <catch.hpp>

#include <stp/stp_eigen.hpp>

using namespace stp;

TEST_CASE( "stp calculation using eigen", "[eigen]" )
{
    Eigen::MatrixXd A(2, 4);
    A << 1, 0, 0, 0,
         0, 1, 1, 1;

    Eigen::MatrixXd B(2, 4);
    B << 1, 1, 0, 1,
         0, 0, 1, 0;

    Eigen::MatrixXd C = stp_calculation( A, B );

    CHECK( C.rows() == 2u );
    CHECK( C.cols() == 8u );
}
