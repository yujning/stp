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
  \file stp_logic_expr.hpp
  \brief header file for stp based logic expression
*/

#pragma once

#include <Eigen/Dense>
#include <cmath>

#include "stp/stp_utils.hpp"

using matrix = Eigen::MatrixXi;           // Defines the type of matrix to use
using matrix_chain = std::vector<matrix>; // Defined matrix chain

namespace stp
{

  void print_strings( const std::vector<std::string>& inputs )
  {
    for( auto s : inputs )
    {
      std::cout << " " << s << " "; 
    }
    std::cout << std::endl;
  }

  class expr_normalize_impl
  {
    public:
      expr_normalize_impl( const std::string& expr, const std::vector<std::string>& input_names )
        : expr( expr ), input_names( input_names )
      {
      }

      bool is_operation( const std::string& word )
      {
        return word.substr( 0, 2 ) == "m_";
      }
      
      bool is_variable( const std::string& word )
      {
        return word.substr( 0, 2 ) == "x_";
      }

      bool is_valid_string()
      {
        std::istringstream iss( expr );
        std::string word;

        while (iss >> word) 
        {
          if ( ( !is_operation( word ) && !is_variable( word ) ) 
            || ( word.length() <= 2 )
            || ( is_operation( word ) && word.length() > 3 ) ) 
          {
            return false;
          }
        }

        return true;
      }

      void initialization()
      {
        Mr = matrix::Zero( 4, 2 );
        Mr << 1, 0, 0, 0, 0, 0, 0, 1;
        
        I2 = matrix::Zero( 2, 2 );
        I2 << 1, 0, 0, 1;
        
        Mc = matrix::Zero( 2, 4 );
        Mc << 1, 0, 0, 0, 0, 1, 1, 1;
        
        Md = matrix::Zero( 2, 4 );
        Md << 1, 1, 1, 0, 0, 0, 0, 1;
        
        Mi = matrix::Zero( 2, 4 );
        Mi << 1, 0, 1, 1, 0, 1, 0, 0;
        
        Mn = matrix::Zero( 2, 2 );
        Mn << 0, 1, 1, 0;
        
        all_tokens = parse_tokens( expr, "" );
      }

      matrix get_variable_matrix( const std::string& var )
      {
          for( int i = 0; i < input_names.size(); i++ )
          {
            if( input_names[i] == var )
            {
              matrix mv( 2, 1 );
              mv << i + 1, i + 1;
              return mv;
            }
          }

          assert( false );
      }

      void expr_to_chain()
      {

        for( const auto& token : all_tokens )
        {
          if( token == "m_c" ) { chain.push_back( Mc ); }
          else if( token == "m_d" ) { chain.push_back( Md ); }
          else if( token == "m_n" ) { chain.push_back( Mn ); }
          else if( token == "m_i" ) { chain.push_back( Mi ); }
          else if( token.substr( 0, 2 ) == "x_" ) //variable
          {
            chain.push_back( get_variable_matrix( token ) );
          }
          else
          {
            std::cerr << "Unknown token.\n";
          }
        }
      }

      void print_chain()
      {
        assert( chain.size() != 0 );

        auto mat_to_var = [&]( const matrix& mat ) { return input_names[ mat( 0, 0 ) - 1 ]; };
        
        for( const auto& m : chain )
        {
          if( m.rows() == 2 && m.cols() == 1 )
          {
            std::cout << mat_to_var( m ) << std::endl;
          }
          else
          {
            std::cout << m << std::endl;
          }
        }
      }

      int get_the_number_vars_before_operation( const std::vector<std::string>& inputs, int op_idx )
      {
        int count = 0;
        for( int i = 0; i < op_idx; i++ )
        {
          if( is_variable( inputs[i] ) )
          {
            count++;
          }
        }
        return count;
      }

      std::vector<std::string> move_vars_to_rightside( const std::vector<std::string>& inputs )
      {
        std::vector<std::string> new_expr;
        std::vector<std::string> temp_vars;

        for( int i = 0; i < inputs.size(); i++ )
        {
          if( !is_variable( inputs[i] ) )
          {
            auto count = get_the_number_vars_before_operation( inputs, i );
            if( count == 0 )
            {
              new_expr.push_back( inputs[i] );
            }
            else
            {
              int dim = std::pow( 2, count );
              new_expr.push_back( "I" + std::to_string( dim ) );
              new_expr.push_back( inputs[i] );
            }
          }
          else
          {
            if( temp_vars.size() > 0 && inputs[i] == temp_vars.back() ) //power reducing
            { 
              auto count = get_the_number_vars_before_operation( inputs, i );
                
              if( count == 0 )
              {
                temp_vars.insert( temp_vars.end() - 1, "PR2" );
              }
              else
              {
                int dim = std::pow( 2, count - 1 ); //note x_i already appear two times, so count must minus 1
                temp_vars.insert( temp_vars.end() - count, 
                                  "I" + std::to_string( dim ) + " PR2" );
              }
            }
            else
            {
              temp_vars.push_back( inputs[i] );
            }
          }
        }
        
        new_expr.insert( new_expr.end(), temp_vars.begin(), temp_vars.end() );
        return new_expr;
      }

      bool variable_compare( const std::string& s1, const std::string& s2 )
      {
        auto it1 = std::find( input_names.begin(), input_names.end(), s1 );
        auto it2 = std::find( input_names.begin(), input_names.end(), s2 );

        if( it1 != input_names.end() && it2 != input_names.end() )
        {
          return std::distance( input_names.begin(), it1 ) < std::distance( input_names.begin(), it2 );
        }

        assert( false );
        return s1 < s2; //lexicographical order
      }

      std::vector<std::string> swap_vars( const std::vector<std::string>& inputs,
                                                    const std::string& var1,
                                                    const std::string& var2 )
      {
        auto result = inputs;
        bool find_var1 = false;
        bool find_var2 = false;
        int idx_var1;
        int idx_var2;

        for( auto i = 0; i < inputs.size(); i++ )
        {
          if( inputs[i] == var1 )
          {
            idx_var1 = i;
            find_var1 = true;
          }
          else if( inputs[i] == var2 )
          {
            idx_var2 = i;
            find_var2 = true;
          }

          if( find_var1 && find_var2 )
          {
            break;
          }
        }

        std::swap( result[idx_var2], result[idx_var1] );
        result.insert( result.begin() + ( idx_var1 > idx_var2 ? idx_var2 : idx_var1 ), "W2" );
        return result;
      }

      std::vector<std::string> sort_vars()
      {
        auto x_tokens = parse_tokens( expr, "x_" );
        auto result = x_tokens;

        for( int i = 1; i < x_tokens.size(); ++i ) 
        {
          std::string key = x_tokens[i];
          int j = i - 1;

          while (j >= 0 && variable_compare( key, x_tokens[j] ) ) 
          {
            auto current_result = result;
            auto vars_right_result = move_vars_to_rightside( current_result );
            result = swap_vars( vars_right_result, x_tokens[j+1], x_tokens[j] );

            if( verbose )
            {
              std::cout << "\nSwap: " << x_tokens[j + 1] << " and " << x_tokens[j] << std::endl;
              std::cout << "Current expr: ";
              print_strings( result );
            }

            std::swap( x_tokens[j], x_tokens[j + 1] );
            --j;
          }
        }
        
        return move_vars_to_rightside( result );
      }

      void run()
      {
        if( is_valid_string() )
        {
          initialization();
          //expr_to_chain();
          //print_chain();
          auto expr_ops  = move_vars_to_rightside( all_tokens );
          auto normal = sort_vars();
          normal.insert( normal.begin(), expr_ops.begin(), expr_ops.end() - input_names.size() - 1 );
          print_strings( normal );
        }
        else
        {
          std::cerr << "[e] Input is not a valid string.\n";
        }

      }

    private:
      std::string expr;
      std::vector<std::string> all_tokens;
      std::vector<std::string> input_names; //the name of variables
      bool verbose = false;
      matrix_chain chain;
      matrix Mr; //power reducing 
      matrix I2; //identity
      matrix Mc; //conjunctive
      matrix Md; //disjunctive
      matrix Mi; //implication
      matrix Mn; //not
  };

  void expr_normalize( const std::string& expr, const std::vector<std::string>& input_names )
  {
    expr_normalize_impl p( expr, input_names );
    p.run();
  }
  

  void test()
  {
    std::string input = "m_i m_c x_1 x_3 m_c x_2 m_n x_3";
    std::vector<std::string> input_names{ "x_1", "x_2", "x_3" };
    expr_normalize( input, input_names );
  }
}
