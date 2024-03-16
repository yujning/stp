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


#include "stp/stp_utils.hpp"


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

  void test()
  {

    std::string input = "m_i m_c x_1 x_3 m_e x_2 m_n x_3";

    auto mTokens   = parse_tokens( input, "m_" );
    print_strings( mTokens );
    auto xTokens   = parse_tokens( input, "x_" );
    print_strings( xTokens );
    auto allTokens = parse_tokens( input, "" );
    print_strings( allTokens );
  }
}
