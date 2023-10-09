#pragma once

#include <Eigen/Dense>

#include <iostream> 
#include <numeric>

using matrix = Eigen::MatrixXi;      // Defines the type of matrix to use
using m_chain = std::vector<matrix>; // Defined matrix chain

/*
  Some static global functions for accessibility.
  Its scope is in this folder only and is used to implement encapsulation.
*/

/* Obtain the least common multiple of m and n. */
//static unsigned get_lcm( unsigned m, unsigned n );

/* Identify whether a matrix is representing a variable. */
//static bool is_variable( const matrix& mat );

/* If a matrix represents a variable, return the variable, otherwise return 0. */
//static int get_variable( const matrix& mat );

/* String segmentation function, can achieve multiple characters as delimiters to partition. */
//static std::vector<std::string> m_split( const std::string& input, const std::string& pred );

namespace stp
{
std::vector<std::string> m_split( const std::string& input, const std::string& pred )
{
  std::vector<std::string> result;
  std::string temp{ "" };
  unsigned count1 = input.size();
  unsigned count2 = pred.size();
  unsigned j;
  for ( size_t i = 0; i < count1; i++ )
  {
    for ( j = 0; j < count2; j++ )
    {
      if ( input[i] == pred[j] )
      {
        break;
      }
    }
    // if input[i] != pred中的任何一个 该字符加到temp上
    if ( j == count2 )
      temp += input[i];
    else
    {
      if ( !temp.empty() )
      {
        result.push_back( temp );
        temp.clear();
      }
    }
  }
  if ( !temp.empty() )
  {
    result.push_back( temp );
    temp.clear();
  }
  return result;
}

unsigned get_lcm( unsigned m, unsigned n )
{
  return ( m * n ) / std::gcd( m, n );
}

bool is_variable( const matrix& mat )
{
  return mat.rows() == 2 && mat.cols() == 1;
}

int get_variable( const matrix& mat )
{
  if ( mat.rows() == 2 && mat.cols() == 1 )
    return mat( 0, 0 );
  return 0;
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
    return B;
  if ( b_dimensions == 1u )
    return A;

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

matrix matrix_chain_mul( const m_chain& matrix_chain )
{
  assert( matrix_chain.size() > 0 );
  matrix result_matrix;
  if ( matrix_chain.size() == 1 )
    return matrix_chain[0];
  result_matrix = semi_tensor_product( matrix_chain[0], matrix_chain[1] );
  for ( int i = 2; i < matrix_chain.size(); i++ )
    result_matrix = semi_tensor_product( result_matrix, matrix_chain[i] );
  return result_matrix;
}

matrix normalize_matrix( m_chain matrix_chain )
{
  matrix Mr( 4, 2 ); // Reduced power matrix
  Mr << 1, 0, 0, 0, 0, 0, 0, 1;
  matrix I2( 2, 2 ); // Identity matrix
  I2 << 1, 0, 0, 1;
  matrix normal_matrix;
  int p_variable;
  int p;
  int max = 0; // the max is the number of variable
  for ( int i = 0; i < matrix_chain.size(); i++ )
  {
    if ( matrix_chain[i]( 0, 0 ) > max )
      max = matrix_chain[i]( 0, 0 );
  }
  std::vector<int> idx( max + 1 ); // id[0] is the max of idx
  p_variable = matrix_chain.size() - 1;
  while ( p_variable >= 0 )
  {
    bool find_variable = false;
    matrix& matrix = matrix_chain[p_variable];
    int var = get_variable( matrix );
    if ( var != 0 ) // 1:find a variable
    {
      if ( idx[var] == 0 ) // the variable appears for the first time ：end : not_end
      {
        idx[var] = idx[0] + 1;
        idx[0]++;
        if ( p_variable == matrix_chain.size() - 1 ) // 第一次出现的变量在矩阵链尾部
        {
          matrix_chain.pop_back();
          p_variable--;
          continue;
        }
      }
      else // the variable appears for the not first time
      {
        if ( idx[var] == idx[0] )
        {
          find_variable = true;
        }
        else
        {
          find_variable = true;
          matrix_chain.push_back( generate_swap_matrix( 2, 1 << ( idx[0] - idx[var] ) ) );
          for ( int i = 1; i <= max; i++ )
          {
            if ( idx[i] != 0 && idx[i] > idx[var] )
              idx[i]--;
          }
          idx[var] = idx[0];
        }
      }
      m_chain matrix_chain1;
      matrix_chain1.clear();
      for ( p = p_variable + 1; p < matrix_chain.size(); p++ )
      {
        matrix_chain1.push_back( matrix_chain[p] );
      }
      while ( p > p_variable + 1 )
      {
        matrix_chain.pop_back();
        p--;
      }
      if ( matrix_chain1.size() > 0 )
      {
        matrix_chain.push_back( matrix_chain_mul( matrix_chain1 ) );
      }
      if ( p_variable != matrix_chain.size() - 1 )
      {
        matrix_chain[p_variable] = kronecker_product( I2, matrix_chain[p_variable + 1] );
        matrix_chain.pop_back();
      }
      if ( find_variable )
        matrix_chain.push_back( Mr );
      continue;
    }
    else
    {
      p_variable--;
    }
  }
  // normal_matrix = matrix_chain_mul(matrix_chain);
  // matrix_chain.clear();
  // matrix_chain.push_back(normal_matrix);
  for ( int i = max; i > 0; i-- ) //!
  {
    matrix_chain.push_back( generate_swap_matrix( 2, pow( 2, idx[0] - idx[i] ) ) );
    for ( int j = 1; j <= max; j++ )
    {
      if ( idx[j] != 0 && idx[j] > idx[i] )
        idx[j]--;
    }
    idx[i] = max;
  }
  normal_matrix = matrix_chain_mul( matrix_chain );
  return normal_matrix;
}

// 表达式字符串   +   给定排序的变量名字   +    是否需要打印生成的矩阵链
// example:             (a & b) | (a & ~c) | (~b & ~c)                                 a b c
matrix from_exp_to_nmx( const std::string& expression, const std::vector<std::string>& input_names, bool print = false )
{
  // str -> vec<str> 区分变量和运算符
  std::vector<std::string> equation;
  equation.push_back( "(" );
  std::string temp = "";
  for ( int i = 0; i < expression.size(); i++ )
  {
    if ( expression[i] == ' ' )
      continue;
    else if ( expression[i] == '(' || expression[i] == ')' || expression[i] == '|' || expression[i] == '&' || expression[i] == '~' )
    {
      if ( !temp.empty() )
        equation.push_back( temp );
      std::string p1 = "";
      p1 = p1 + expression[i];
      equation.push_back( p1 );
      temp = "";
    }
    else
      temp = temp + expression[i];
  }
  equation.push_back( ")" );
  // equation = ["(",   "(", "a", "&", "b", ")",   "(", "a", "&", "~", "c", ")",  "(", "~", "a", "&", "~", "c", ")"   ")"]
  /***********/
  // deal with "~"
  for ( int i = 0; i < equation.size(); i++ )
  {
    if ( equation[i] == "~" )
    {
      equation[i] = "~(" + equation[i + 1] + ")";
      equation.erase( equation.begin() + i + 1 );
    }
    else if ( equation[i] == "&" || equation[i] == "|" || equation[i] == "(" || equation[i] == ")" )
      continue;
    else
      equation[i] = "(" + equation[i] + ")";
  }
  // equation = ["(",   "(", "(a)", "&", "(b)", ")",   "(", "(a)", "&", "~(c)", ")",  "(", "~(a)", "&", "~(c)", ")"   ")"]
  //?
  std::vector<int> left_bracket;
  for ( int i = 0; i < equation.size(); i++ )
  {
    if ( equation[i] == "(" )
    {
      left_bracket.push_back( i );
    }
    if ( equation[i] == ")" )
    {
      std::string equ = "";
      for ( int j = left_bracket[left_bracket.size() - 1] + 1; j < i; j++ ) // 遍历左右括号之间的表达式
      {
        if ( equation[j] == "|" || equation[j] == "&" )
          equ = "(" + equation[j] + ")" + equ;
        else
          equ = equ + "(" + equation[j] + ")";
      }
      //"(" -> equ
      equation[left_bracket[left_bracket.size() - 1]] = equ;
      equation.erase( equation.begin() + left_bracket[left_bracket.size() - 1] + 1, equation.begin() + i + 1 );
      i = left_bracket[left_bracket.size() - 1] - 1;
      left_bracket.pop_back();
    }
  }
  std::string equ = equation[0];
  equation.clear();
  equation = m_split( equ, "()" );
  m_chain matrix_chain;
  Eigen::MatrixXi Mc( 2, 4 );
  Mc << 1, 0, 0, 0, 0, 1, 1, 1;
  Eigen::MatrixXi Md( 2, 4 );
  Md << 1, 1, 1, 0, 0, 0, 0, 1;
  Eigen::MatrixXi Mn( 2, 2 );
  Mn << 0, 1, 1, 0;
  for ( int i = 0; i < equation.size(); i++ )
  {
    if ( equation[i] == "&" )
      matrix_chain.push_back( Mc );
    else if ( equation[i] == "|" )
      matrix_chain.push_back( Md );
    else if ( equation[i] == "~" )
      matrix_chain.push_back( Mn );
    else if ( equation[i] == "" || equation[i] == " " )
      ;
    else // variable
    {
      for ( int j = 0; j < input_names.size(); j++ )
      {
        if ( equation[i] == input_names[j] )
        {
          Eigen::MatrixXi M( 2, 1 );
          M << j + 1, j + 1;
          matrix_chain.push_back( M );
          break;
        }
      }
    }
  }
  // print m_chain
  if ( print )
  {
    auto mat_to_var = [&]( const matrix& mat ) { return input_names[mat( 0, 0 ) - 1]; };
    for ( const auto& mat : matrix_chain )
    {
      if ( is_variable( mat ) )
        std::cout << mat_to_var( mat ) << std::endl;
      else
        std::cout << mat << std::endl;
    }
  }

  matrix p = normalize_matrix( matrix_chain );
  return p;
}
} // namespace stp

