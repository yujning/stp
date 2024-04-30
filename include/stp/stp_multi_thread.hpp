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
  \file stp_multi_thread.hpp
  \brief header file for stp multi thread calculation
  \author Zhufei Chu
*/

#pragma once

#include <Eigen/Dense>

#include <thread>
#include <mutex>

#include "stp/stp_eigen.hpp"

using matrix = Eigen::MatrixXi;           // Defines the type of matrix to use
using matrix_chain = std::vector<matrix>; // Defined matrix chain


namespace stp
{

  struct temp_result
  {
    matrix product; 
    int start_index;
  };

  class matrix_chain_multiply_by_multi_thread_impl
  {
    public:
    matrix_chain_multiply_by_multi_thread_impl( const matrix_chain& mc, int num_threads, bool& verbose )
      :mc( mc ), num_threads( num_threads ), verbose( verbose )
    {
    }

     matrix_chain get_sub_chain( int start, int end )
     {
       assert( start >=0 );
       assert( start < mc.size() );
       assert( end > 0 );
       assert( end <= mc.size() );
       
       std::lock_guard<std::mutex> lock( mtx );
       matrix_chain sub_chain;

       for( int i = start; i < end; i++ )
       {
         sub_chain.push_back( mc[i] );
       }

       return sub_chain;
     }

     matrix sub_chain_product( const auto& schain )
     {
       std::lock_guard<std::mutex> lock( mtx );
       return matrix_chain_multiply( schain );
     }
    
     matrix_chain run()
     {
       std::vector<std::thread> threads;

       int block_size = mc.size() / num_threads;
       int remaining_elements = mc.size() % num_threads;

       std::vector<temp_result> part( num_threads );

       if( verbose )
       {
         std::cout << "[i] The matrix chain size: " << mc.size() << " #threads: " << num_threads << std::endl; 
       }


       int start = 0;

       for( int i = 0; i < num_threads; i++ )
       {
         int num_elements = block_size + ( i < remaining_elements ? 1 : 0 );
         int end = start + num_elements;

         if( verbose )
         {
           std::cout << "[i] Thread " << i << " starts to process " << num_elements << " elements.\n";
         }
         
         threads.emplace_back( [this, i, &part, start, end]() 
             {
             part[i].product = sub_chain_product( get_sub_chain( start, end ) );
             part[i].start_index = start;
             });

         start = end;
       }

       for( auto& thread : threads )
       {
         thread.join();
       }

       std::sort( part.begin(), part.end(), []( const temp_result& a, const temp_result& b ){
            return a.start_index < b.start_index;
           } );

       matrix_chain result;

       for( const auto& pr : part )
       {
         result.push_back( pr.product );
       }
       return result;
     }

    private:
      matrix_chain mc;
      std::mutex mtx;
      int num_threads;
      bool verbose;
  };

  matrix_chain matrix_chain_multiply_by_multi_thread( const matrix_chain& mc, const int& num_threads, bool verbose = false )
  {
    matrix_chain_multiply_by_multi_thread_impl p( mc, num_threads, verbose );
    return p.run();
  }

} //end of namespace
