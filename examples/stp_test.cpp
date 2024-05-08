#include <iostream>
#include <stp/stp_eigen.hpp>

using Eigen::MatrixXi;

matrix random_generation(int row, int col)
{
  matrix result(row, col);
  for(int i = 0; i < row; i++)
    for(int j = 0; j < col; j++)
      result(i, j) = rand() % 2;
  return result;
}

int main()
{
  //n % p = 0 type matrix multiplication test
  for(int i = 0; i < 10; i++)
  {
    int m = 200;
    int n = 4 << i;
    int p = 4;
    int q = 100;
    matrix A = random_generation(m, n);
    matrix B = random_generation(p, q);
    stp::semi_tensor_product(A, B, true, stp::stp_method::definition_1);
    stp::semi_tensor_product(A, B, true, stp::stp_method::definition_2); 
    stp::semi_tensor_product(A, B, true);
  }

  std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
  
  //p % n = 0 type matrix multiplication test
  for(int i = 0; i < 10; i++)
  {
    int m = 200;
    int n = 4;
    int p = 4 << i;
    int q = 100;
    matrix A = random_generation(m, n);
    matrix B = random_generation(p, q);
    stp::semi_tensor_product(A, B, true, stp::stp_method::definition_1);
    stp::semi_tensor_product(A, B, true, stp::stp_method::definition_2); 
    stp::semi_tensor_product(A, B, true);
  }
  return 0;
}
