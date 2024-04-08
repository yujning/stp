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
          if ( ( word.substr(0, 2) != "m_" && word.substr(0, 2) != "x_") 
            || ( word.length() <= 2 )
            || ( word.substr(0, 2) == "m_" && word.length() > 3 ) ) 
          {
            return false;
          }
        }

        return true;
      }

      void matrix_initialization()
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
        auto allTokens = parse_tokens( expr, "" );

        for( const auto& token : allTokens )
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


      void run()
      {
        if( is_valid_string() )
        {
          matrix_initialization();
          expr_to_chain();
          print_chain();
        }
        else
        {
          std::cerr << "[e] Input is not a valid string.\n";
        }

      }

    private:
      std::string expr;
      std::vector<std::string> input_names; //the name of variables
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

    std::string input = "m_i m_c x_1 x_3 m_c x_2 m_n m_i x_3";

    auto mTokens   = parse_tokens( input, "m_" );
    print_strings( mTokens );
    auto xTokens   = parse_tokens( input, "x_" );
    print_strings( xTokens );
    auto allTokens = parse_tokens( input, "" );
    print_strings( allTokens );
    
    std::vector<std::string> mTokens_detail;
    std::vector<std::string> xTokens_detail;
    std::unordered_map<std::string, int> tokenCounts;
    std::unordered_map<std::string, int> tokenIndices;

    parse_tokens_with_details( input, "m_", mTokens_detail, tokenCounts, tokenIndices );
    parse_tokens_with_details( input, "x_", xTokens_detail, tokenCounts, tokenIndices );
    print_strings( mTokens_detail );
    print_strings( xTokens_detail );
    std::vector<std::string> input_names{ "x_1", "x_2", "x_3" };
    expr_normalize( input, input_names );
  }
}
