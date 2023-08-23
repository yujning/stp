#include <stp/stp.hpp>
#include <iostream>
#include <vector>
using namespace std;

int main()
{
  stp::MyStp stp_engine;
  stp_engine.getVersion();
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  string equ1 = "(a & b) | (a & ~c) | (~b & ~c)";
  vector<string> inputs_order1 = {"a", "b", "c"};
  matrix mat1 =  stp::from_exp_to_nmx(equ1, inputs_order1, true);
  cout << "equ1: " << endl;
  cout << mat1 << endl;
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  string equ2 = "(a & b & ~c) | (a & ~c) | (~b & ~c)";
  vector<string> inputs_order2 = {"a", "b", "c"};
  matrix mat2 =  stp::from_exp_to_nmx(equ2, inputs_order2, true);
  cout << "equ2: " << endl;
  cout << mat2 << endl;

  std::cout << "----------------------------------------------------------------------------" << std::endl;
  return 0;
}
