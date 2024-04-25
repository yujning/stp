/* stp: C++ semi-tensor product library for electronic design automation (EDA)
 * Copyright (C) 2023-  Ningbo University, Zhejiang, China
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file stp_eigen.hpp
  \brief header file for some basic STP operations
  \author Zhufei Chu
  \author Ruibing Zhang
*/

#pragma once

#include <Eigen/Dense>

#include <iostream> 
#include <numeric>
#include <stp/stp_timer.hpp>

using matrix = Eigen::MatrixXi;           // Defines the type of matrix to use
using matrix_chain = std::vector<matrix>; // Defined matrix chain


namespace stp
{

  inline unsigned get_lcm( unsigned m, unsigned n )
  {
    return ( m * n ) / std::gcd( m, n );
  }

  inline matrix generate_swap_matrix( const int& m, const int& n )
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

  inline matrix kronecker_product( const matrix& A, const matrix& B )
  {
    /* trivial cases */
    auto a_dimensions = A.rows() * A.cols();
    auto b_dimensions = B.rows() * B.cols();

    if ( a_dimensions == 1u )
    {
      return B;
    }
    if ( b_dimensions == 1u )
    {
      return A;
    }

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

  inline matrix semi_tensor_product( const matrix& A, const matrix& B )
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

  inline matrix semi_tensor_product2( const matrix& A, const matrix& B )
  {
    int m = A.rows(); 
    int n = A.cols(); 

    int p = B.rows(); 
    int q = B.cols(); 

    int row, col;

    matrix result_matrix;

    if (n % p == 0)
    {
      int times = n / p;
      row = m;
      col = times * q;
      result_matrix = matrix::Zero(row, col);
      for (int i = 0; i < q; ++i) 
      {
        for (int j = 0; j < p; ++j)
        {
          result_matrix.block(0, i * times, m, times) += B(j, i) * A.block(0, j * times, m, times);
        }
      }
    }
    else if (p % n == 0)
    {
      int times = p / n;
      row = m * times;
      col = q;
      result_matrix = matrix::Zero(row, col);
      for (int i = 0; i < m; ++i) //
      {
        for (int j = 0; j < n; ++j)
        {
          result_matrix.block(i * times, 0, times, q) += A(i, j) * B.block(j * times, 0, times, q);
        }
      }
    }
    else
    {
      std::cout << "matrix type error!" << std::endl;
    }
    return result_matrix;
  }

  inline matrix matrix_chain_multiply( const matrix_chain& mc )
  {
    assert( mc.size() > 0 );
    matrix result_matrix;

    if ( mc.size() == 1 )
    {
      return mc[0];
    }

    result_matrix = semi_tensor_product( mc[0], mc[1] );

    for ( int i = 2; i < mc.size(); i++ )
    {
      result_matrix = semi_tensor_product( result_matrix, mc[i] );
    }
    return result_matrix;
  }

  class stp_eigen 
  {
    public:
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
        {
          return B;
        }
        if ( b_dimensions == 1u )
        {
          return A;
        }

        stp_timer timer(true);

        matrix KP( A.rows() * B.rows(), A.cols() * B.cols() );

        for ( int i = 0; i < A.rows(); ++i )
        {
          for ( int j = 0; j < A.cols(); ++j )
          {
            KP.block( i * B.rows(), j * B.cols(), B.rows(), B.cols() ) = A( i, j ) * B;
          }
        }

        kp_time = timer.get_elapsed_ms();
        total_time = kp_time + stp_time;

        return KP;
      }
    
      matrix semi_tensor_product( const matrix& A, const matrix& B )
      {
        int n = A.cols(); 
        int p = B.rows();
        matrix result;

        stp_timer timer(true);

        if(std::gcd( n, p ) < 4) 
        {
          result = semi_tensor_product1( A, B );
        }
        else
        {
          result = semi_tensor_product2( A, B );
        }

        stp_time = timer.get_elapsed_ms();
        total_time = kp_time + stp_time;

        return result;
      }
      
      matrix matrix_chain_multiply( const matrix_chain& mc )
      {
        assert( mc.size() > 0 );

        if ( mc.size() == 1 )
        {
          return mc[0];
        }

        if ( mc.size() == 2 )
        {
          return semi_tensor_product( mc[0], mc[1] );
        }

        std::vector<int> orders = matrix_chain_order( mc );
        return parse_multiplication_order(orders, mc);
      }

      matrix normalize_matrix( matrix_chain mc )
      {
        matrix Mr( 4, 2 ); // Reduced power matrix
        Mr << 1, 0, 0, 0, 0, 0, 0, 1;

        matrix I2( 2, 2 ); // Identity matrix
        I2 << 1, 0, 0, 1;

        matrix normal_matrix;
        int p_variable;
        int p;

        int max = 0; // the max is the number of variable

        for ( int i = 0; i < mc.size(); i++ )
        {
          if ( mc[i]( 0, 0 ) > max )
          {
            max = mc[i]( 0, 0 );
          }
        }

        std::vector<int> idx( max + 1 ); // id[0] is the max of idx
        p_variable = mc.size() - 1;

        while ( p_variable >= 0 )
        {
          bool find_variable = false;
          matrix& matrix = mc[p_variable];
          int var = get_variable( matrix );

          if ( var != 0 ) // 1:find a variable
          {
            if ( idx[var] == 0 ) // the variable appears for the first time ：end : not_end
            {
              idx[var] = idx[0] + 1;
              idx[0]++;

              if ( p_variable == mc.size() - 1 ) // the variable shows in the end
              {
                mc.pop_back();
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
                mc.push_back( generate_swap_matrix( 2, 1 << ( idx[0] - idx[var] ) ) );

                for ( int i = 1; i <= max; i++ )
                {
                  if ( idx[i] != 0 && idx[i] > idx[var] )
                    idx[i]--;
                }

                idx[var] = idx[0];
              }
            }

            matrix_chain mc_temp;
            mc_temp.clear();

            for ( p = p_variable + 1; p < mc.size(); p++ )
            {
              mc_temp.push_back( mc[p] );
            }

            while ( p > p_variable + 1 )
            {
              mc.pop_back();
              p--;
            }

            if ( mc_temp.size() > 0 )
            {
              mc.push_back( matrix_chain_multiply( mc_temp ) );
            }

            if ( p_variable != mc.size() - 1 )
            {
              mc[p_variable] = kronecker_product( I2, mc[p_variable + 1] );
              mc.pop_back();
            }

            if ( find_variable )
            {
              mc.push_back( Mr );
            }
            continue;
          }
          else
          {
            p_variable--;
          }
        }

        for ( int i = max; i > 0; i-- ) 
        {
          mc.push_back( generate_swap_matrix( 2, pow( 2, idx[0] - idx[i] ) ) );

          for ( int j = 1; j <= max; j++ )
          {
            if ( ( idx[j] != 0 ) && ( idx[j] > idx[i] ) )
            {
              idx[j]--;
            }
          }

          idx[i] = max;
        }

        normal_matrix = matrix_chain_multiply( mc );
        return normal_matrix;
      }

      void print_stats()
      {
        std::cout << "kronecker   product: " << kp_time << "ms" <<  std::endl;
        std::cout << "semi_tensor product: " << stp_time << "ms" <<  std::endl;
        std::cout << "              total: " << total_time << "ms" <<  std::endl;
      }

    private:
      unsigned get_lcm( unsigned m, unsigned n )
      {
        return ( m * n ) / std::gcd( m, n );
      }

      // the method used to implement semi tensor product
      matrix semi_tensor_product1( const matrix& A, const matrix& B )
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

      matrix semi_tensor_product2( const matrix& A, const matrix& B )
      {
        int m = A.rows(); 
        int n = A.cols(); 

        int p = B.rows(); 
        int q = B.cols(); 

        int row, col;

        matrix result_matrix;

        if (n % p == 0)
        {
          int times = n / p;
          row = m;
          col = times * q;
          result_matrix = matrix::Zero(row, col);
          for (int i = 0; i < q; ++i) 
          {
            for (int j = 0; j < p; ++j)
            {
              result_matrix.block(0, i * times, m, times) += B(j, i) * A.block(0, j * times, m, times);
            }
          }
        }
        else if (p % n == 0)
        {
          int times = p / n;
          row = m * times;
          col = q;
          result_matrix = matrix::Zero(row, col);
          for (int i = 0; i < m; ++i) 
          {
            for (int j = 0; j < n; ++j)
            {
              result_matrix.block(i * times, 0, times, q) += A(i, j) * B.block(j * times, 0, times, q);
            }
          }
        }
        else
        {
          std::cout << "matrix type error!" << std::endl;
          assert(false);
        }
        return result_matrix;
      }

      // the method used to implement matrix chain multiply
      std::vector<int> matrix_chain_order(const matrix_chain& mc)
      {
        struct mc_info
        {
          uint32_t op_num = 0;  // Records the total complexity of the current operation
          uint32_t flag   = 0;  // Record the matrix multiplication sequence
          uint32_t row    = 0; 
          uint32_t col    = 0;  
        };
        
        int length = mc.size();
        std::vector<std::vector<mc_info>> dp(length, std::vector<mc_info>(length));

        for (int i = 0; i < length; i++)
        {
          dp[i][i].row = mc[i].rows();
          dp[i][i].col = mc[i].cols();
        }

        for (int l = 2; l <= length; l++) /* 矩阵链的长度 */
        {
          for (int i = 0; i < length - l + 1; i++)
          {
            int j = i + l - 1; 
            dp[i][j].op_num = INT_MAX;
            std::vector<int> temp;
            for (int k = i; k < j; k++)
            {
              temp = complexity_analysis(dp[i][k].row, dp[i][k].col, dp[k+1][j].row,  dp[k+1][j].col);
              int q = dp[i][k].op_num + dp[k+1][j].op_num + temp[0];
              if (q < dp[i][j].op_num)
              {
                dp[i][j].op_num = q;
                dp[i][j].flag = k;
                dp[i][j].row = temp[1];
                dp[i][j].col = temp[2];
              }
            }
          }
        }
        
        // Parse the Dynamic programming table
        std::vector<int> result;
        std::vector<std::vector<int>> flags(length, std::vector<int>(length));
        for(int i = 0; i < length; i++)
          for(int j = 0; j < length; j++)
            flags[i][j] = dp[i][j].flag;
        optimal_parens(flags, 0, length - 1, result);
        return result;
      }

      std::vector<int> complexity_analysis(int m, int n, int p, int q)
      {
        /*
          temp[0] operational complexity
          temp[1] rows
          temp[2] cols
        */ 
        std::vector<int> temp(3, 0);
        if (n % p == 0)
        {
          temp[0] = n * q / p;   
          temp[1] = m;
          temp[2] = n * q / p;
          return temp;
        }
        else if (p % n == 0)
        {
          temp[0] = q;
          temp[1] = m * p / n;
          temp[2] = q;
          return temp;
        }
        else
        {
          std::cout << "matrix type error!" << std::endl;
        }
        return temp;
      }

      void optimal_parens(const std::vector<std::vector<int>> &s, int i, int j, std::vector<int> &result)
      {
        if (i == j)
            result.push_back(i);
        else
        {
          result.push_back(-1);
          optimal_parens(s, i, s[i][j], result);
          optimal_parens(s, s[i][j] + 1, j, result);
          result.push_back(-2);
        }
      }

      matrix parse_multiplication_order(std::vector<int>& order, const std::vector<matrix>& mc)
      {
        std::vector<matrix> stacks;
        std::vector<int> idx;
        for (unsigned int i = 0; i < order.size(); ++i)
        {
          if (order[i] == -1)
            idx.push_back(i);
          else if (order[i] == -2)
          {
            if (stacks.size() >= 2)
            {
              stacks[stacks.size() - 2] = semi_tensor_product(stacks[stacks.size() - 2], stacks[stacks.size() - 1]);
              stacks.pop_back();
            }
            order.erase(order.begin() + idx.back(), order.begin() + i);
            i = idx.back();
            idx.pop_back();
          }
          else
            stacks.push_back(mc[order[i]]);
        }
        assert(stacks.size() == 1);
        return stacks[0];
      }

      int get_variable( const matrix& mat )
      {
        if( is_variable( mat ) )
        {
          return mat( 0, 0 );
        }
        else
        {
          return 0;
        }
      }

      bool is_variable( const matrix& mat )
      {
        return mat.rows() == 2 && mat.cols() == 1;
      }
      
    private:
      uint64_t stp_time = 0u;
      uint64_t kp_time = 0u;
      uint64_t total_time = 0;
  };

} // namespace stp

