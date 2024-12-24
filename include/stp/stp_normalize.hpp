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
  \file stp_normalize.hpp
  \brief header file for various stp expression normalization
  \author Zhufei Chu
*/

#pragma once

#include <Eigen/Dense>
#include <algorithm>

#include "stp/stp_eigen.hpp"
#include "stp/stp_utils.hpp"

using matrix = Eigen::MatrixXi;            // Defines the type of matrix to use
using matrix_chain = std::vector<matrix>;  // Defined matrix chain
using ids = std::vector<uint32_t>;

namespace stp
{
/*
 * Given a input matrix chain expressed in vector, the class is used to
 * normalize it by multiple methods.
 *
 * The input matrxi format:
 * 1000 1011 1 2
 * 0111 0100 1 2
 *
 * which corresponds to structural matrices and varaibles (with 2*1
 * dimensions, the first value is larger and equal than 1)
 * */

class stp_normalize_impl
{
 public:
  stp_normalize_impl( const stp_circuit& circuit, bool& verbose )
      : circuit( circuit ), verbose( verbose )
  {
    // initialization
    vars_order.clear();
    chain.clear();

    I2 = matrix::Zero( 2, 2 );
    I2 << 1, 0, 0, 1;

    PR2 = matrix::Zero( 4, 2 );
    PR2 << 1, 0, 0, 0, 0, 0, 0, 1;

    unsigned num_vars = 0u;
    for ( const id& pi : circuit.get_inputs() )
      {
        num_vars++;
        vars_order[ pi ] = circuit.get_inputs().size() - num_vars + 1;
        std::cout << "pi: " << pi << " order: " << vars_order[ pi ]
                  << std::endl;
      }
  }

  void run()
  {
    std::cout << "Hello" << std::endl;

    assert( circuit.get_outputs().size() == 1 );
    const id po = circuit.get_outputs()[ 0 ];
    get_chain_id( po );
    print_chain_ids();
    print_expr();
    get_chain();
    // print_chain( chain );

    auto new_chain = chain_normalize_method1();
    // print_chain( new_chain );

    std::cout << "The final results: \n"
              << matrix_chain_multiply( legal_chain( new_chain ) ) << "\n";
  }

  matrix_chain chain_normalize_method1()
  {
    matrix_chain input_chain = chain;

    while ( !is_normal_chain( input_chain ) )
      {
        auto output_chain = reorder( input_chain );
        input_chain = output_chain;
      }

    return input_chain;
  }

  matrix_chain legal_chain( const matrix_chain& mc )
  {
    matrix_chain order;

    for ( const auto& m : mc )
      {
        if ( !is_variable( m ) )
          {
            order.push_back( m );
          }
      }

    return order;
  }

  bool is_normal_chain( const matrix_chain& mc )
  {
    for ( unsigned i = 0; i < mc.size() - 1u; i++ )
      {
        if ( is_variable( mc[ i ] ) && is_matrix( mc[ i + 1 ] ) )
          {
            return false;
          }
        else if ( is_variable( mc[ i ] ) && is_variable( mc[ i + 1 ] ) )
          {
            auto a = get_variable( mc[ i ] );
            auto b = get_variable( mc[ i + 1 ] );

            if ( a >= b )
              {
                return false;
              }
          }
      }

    return true;
  }

  matrix_chain reorder( const matrix_chain& mc )
  {
    matrix_chain order = mc;

    std::vector<unsigned> insert_idx;

    for ( unsigned i = 0; i < mc.size() - 1u; i++ )
      {
        // swap variable and matrix
        if ( is_variable( mc[ i ] ) && is_matrix( mc[ i + 1 ] ) )
          {
            std::cout << "swap " << i << " and " << i + 1 << "\n";
            // swap
            matrix temp = kronecker_product( I2, order[ i + 1 ] );
            order[ i + 1 ] = order[ i ];
            order[ i ] = temp;
          }
        // variables power reducing or swapping
        else if ( is_variable( mc[ i ] ) && is_variable( mc[ i + 1 ] ) )
          {
            auto a = get_variable( mc[ i ] );
            auto b = get_variable( mc[ i + 1 ] );

            if ( a == b )
              {
                std::cout << "PR2"
                          << "\n";
                order[ i ] = PR2;
              }
            else
              {
                // order the variables
                if ( a > b )
                  {
                    std::cout << "reoder " << i << " and " << i + 1 << "\n";
                    matrix temp = order[ i + 1 ];
                    order[ i + 1 ] = order[ i ];
                    order[ i ] = temp;

                    insert_idx.push_back( i );

                    // order.insert( order.begin() + i + count,
                    // generate_swap_matrix( 2, 2 ) );
                  }
              }
          }
        else
          {
            continue;
          }
      }

    // insert the swap matrix in reverse order
    sort( insert_idx.rbegin(), insert_idx.rend() );
    matrix swap_matrix2 = generate_swap_matrix( 2, 2 );
    for ( const auto& j : insert_idx )
      {
        order.insert( order.begin() + j, swap_matrix2 );
      }

    // std::cout <<
    // "----------------------------------------------------------------------\n";
    // print_chain( order );

    return order;
  }

  bool is_variable( const matrix& m ) { return m.rows() == 2 && m.cols() == 1; }

  int get_variable( const matrix& m )
  {
    assert( is_variable( m ) );

    return m( 0, 0 );
  }

  // check whether it is a matrix
  bool is_matrix( const matrix& m ) { return m.cols() > 1; }

  void print_chain( const matrix_chain& mc )
  {
    for ( const auto& m : mc )
      {
        std::cout << m << "\n" << std::endl;
      }
  }

  void print_chain_ids()
  {
    assert( nodes.size() != 0u );
    for ( const auto n : nodes )
      {
        std::cout << " " << n;
      }
    std::cout << std::endl;
  }

  bool is_mtx( const uint32_t t )
  {
    if ( circuit.is_pi( t ) )
      {
        return false;
      }
    if ( t >= circuit.get_nodes().size() )
      {
        return false;
      }
    return true;
  }

  void print_expr()
  {
    std::cout << "[stp_expr]: ";
    for ( const uint32_t t : nodes )
      {
        if ( circuit.is_pi( t ) )
          {
            std::cout << circuit.get_node( t ).get_name() << " ";
          }
        else if ( is_mtx( t ) )
          {
            std::cout << "m" << t << " ";
          }
        else
          {
            std::cout << "UNKNOWN ";
          }
      }
    std::cout << "\n";
  }

 private:
  // get the chain expressed using node ids
  void get_chain_id( const id& n )
  {
    nodes.push_back( n );

    if ( !circuit.is_pi( n ) )
      {
        for ( const auto& input : circuit.get_node( n ).get_inputs() )
          {
            get_chain_id( input.index );
          }
      }
    else
      {
        return;
      }
  }

  matrix get_matrix( const uint32_t n )
  {
    if ( circuit.is_pi( n ) )
      {
        uint32_t order = vars_order[ n ];
        matrix result( 2, 1 );
        result << order, order;
        return result;
      }

    return circuit.get_node( n ).get_mtx();
  }

  void get_chain()
  {
    for ( const id& n : nodes )
      {
        chain.push_back( get_matrix( n ) );
      }
  }

 private:
  matrix_chain chain;
  bool verbose;
  const stp_circuit& circuit;
  ids nodes;
  std::unordered_map<uint32_t, uint32_t>
      vars_order;  // record the variables order
  matrix I2;
  matrix PR2;
};

void stp_normalize( const stp_circuit& circuit, bool verbose = false )
{
  stp_normalize_impl p( circuit, verbose );
  p.run();
}

}  // namespace stp
