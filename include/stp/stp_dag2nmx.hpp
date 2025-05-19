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
  \brief header file for dag to normalize
  \author Zhufei Chu
  \author Ruibing Zhang
*/

#pragma once

#include <Eigen/Dense>
#include <cmath>

#include "stp/stp_circuit.hpp"
#include "stp/stp_logic_expr.hpp"
#include "stp/stp_utils.hpp"
#ifdef STP_ENABLE_CUDA
#include "stp/stp_cuda.hpp"
#endif

using stp_expr = std::vector<uint32_t>;

namespace stp
{
class circuit_normalize_impl
{
 public:
  circuit_normalize_impl( const stp_circuit& circuit, const bool& verbose)
      : circuit( circuit ), verbose( verbose )
  {
  }

  matrix run()
  {
    initialization();

    assert( circuit.get_outputs().size() == 1 );

    const id po = circuit.get_outputs()[ 0 ];

    get_chain_code( po );

    if ( verbose )
      {
        print_expr( mc );
      }

    get_chain();

    result = normalize_matrix( chain );
    return result;
  }

  matrix new_run()
  {
    initialization();

    assert( circuit.get_outputs().size() == 1 );

    const id po = circuit.get_outputs()[ 0 ];

    get_chain_code( po );

    if ( verbose )
      {
        print_expr( mc );
      }

    int num_vars_in_expr = 0;
    for ( const uint32_t t : mc )
      {
        if ( circuit.is_pi( t ) )
          {
            num_vars_in_expr++;
          }
      }

    stp_expr expr_ops = move_vars_to_rightside( mc );

    if ( verbose )
      {
        print_expr( expr_ops );
      }

    stp_expr normal = sort_vars();

    normal.insert( normal.begin(), expr_ops.begin(),
                   expr_ops.end() - num_vars_in_expr );

    if ( verbose )
      {
        print_expr( normal );
      }

    chain = expr_to_chain( normal );

    result = matrix_chain_multiply( chain, verbose );

    return result;
  }

#ifdef STP_ENABLE_CUDA
  std::vector<int32_t> run_cuda1(stp_expr& all_pi)
  {
    initialization();

    assert( circuit.get_outputs().size() == 1 );

    const id po = circuit.get_outputs()[ 0 ];

    get_chain_code( po );

    if ( verbose )
	  {
		  print_expr( mc );
	  }

	  //Merge identical variables
    stp_expr expr_ops = move_vars( mc );

	  //Variable right shift
    stp_expr normal = move_vars_to_rightside_1( expr_ops );
       
    all_pi = collect_pi(normal);
	
	  //matrix encoding
    std::vector<int32_t> result_vec;

	  bool I_flag = false; //The previous one is an identity matrix
	  int32_t I_dim = 0;  //dim of Identity matrix

	//读取字符串信息 进行半张量积
	for ( int32_t i = 0; i < normal.size() - circuit.get_inputs().size(); i++ )
	{
		uint32_t n = normal[i];

		//Special matrix
		if( n >= circuit.get_nodes().size() )
		{
			assert( other_matrix.find( n ) != other_matrix.end() );
			std::string str = other_matrix[ n ];

			//Power matrix
			if ( str == "Mr" )
			{
				std::vector<int32_t> Mr = {4, 0, 3};

				//The chain has only one matrix, return
				if(i == 0)
				{
					result_vec = Mr;
				}
				else
				{
					if(I_flag)
					{
						std::vector<int32_t> temp = In_KR_Matrix(I_dim,Mr);
						result_vec = my_semi_tensor_product(result_vec, temp);

						I_dim = 0;
						I_flag = false;
					}
					else
					{
						result_vec = my_semi_tensor_product(result_vec, Mr);
					}
				}
			}
			//Identity matrix
			if ( str[ 0 ] == 'I' )
			{
				//str.erase( str.begin() );
				I_flag = true;
				//I_dim = std::stoi( str );  //string to decimal
				I_dim = std::stoi( str.substr(1) );  //string to decimal
			}
			//Swap matrix
			if ( str[ 0 ] == 'W' )
			{
				//Delete the first character of the string
				str.erase( str.begin() );
				assert( std::stoi( str ) == 2 );

				std::vector<int32_t> M = {4, 0, 2, 1, 3};

				//The chain has only one matrix, return
				if(i == 0)
				{
					result_vec = M;
				}
				else
				{
					if(I_flag)
					{
						std::vector<int32_t> temp = In_KR_Matrix(I_dim,M);
						result_vec = my_semi_tensor_product(result_vec, temp);

						I_dim = 0;
						I_flag = false;
					}
					else
					{
						result_vec = my_semi_tensor_product(result_vec, M);
					}
				}
			}
		}
		else
		{
			matrix matrix_temp = circuit.get_node( n ).get_mtx();
			std::vector<int32_t> Vec = Matrix_to_Vec(matrix_temp);

			if(i == 0)
			{
				result_vec = Vec;
			}
			else
			{
				if(I_flag)
				{
					std::vector<int32_t> temp = In_KR_Matrix(I_dim,Vec);
					result_vec = my_semi_tensor_product(result_vec, temp);

					I_dim = 0;
					I_flag = false; 
				}
				else
				{
					result_vec = my_semi_tensor_product(result_vec, Vec);
				}
			}
		}
	}

	assert( I_flag == false ); 

    return result_vec;
  }


  std::vector<int32_t> run_cuda()
  {
    initialization();
    assert(circuit.get_outputs().size() == 1);
    const id po = circuit.get_outputs()[0];
    
    get_chain_code(po);

    if(verbose)
    {
      print_expr(mc);
    }

    get_chain();

    matrix_chain temp = cuda_normalize_matrix( chain );

    std::vector<std::vector<int32_t>> result_cuda_vec = Matrix_to_Vec_chain(temp);

    std::vector<int32_t> result_cuda = my_chain_multiply_by_multi_thread(result_cuda_vec);

    return result_cuda;
  }
#endif

  std::string run_str( bool new_ = false, bool _using_CUDA = false )
  {
    std::string str = "";
    //using cuda
#ifdef STP_ENABLE_CUDA
    if(_using_CUDA)
    {
        stp_expr all_pi;
        std::vector<int32_t> m = run_cuda1(all_pi);

        std::vector<int32_t> out_vec(m.size() - 1);
        for(int i = 1; i < m.size(); i++)
        {
          int32_t temp = m.size() - 1 -i;
          int32_t input_num = all_pi.size();
          int32_t out_temp=0;

          while(input_num!=0)
          {
            out_temp += temp % 2 * pow(2 , all_pi[all_pi.size() - input_num]);
            temp = temp / 2;
            input_num--;
          }

          out_vec[out_temp] = m[i];

        }
        for(int i = out_vec.size()- 1; i >=0; i--)
        {
          if(out_vec[i] == 1) 
          {
            str += '0';
          }
          else if(out_vec[i] == 0) 
          {
            str += '1';
          }
          else
          {
            std::cout<<"str err"<<std::endl;
          }
        }
        return str;
    }
#endif

      matrix m;

      if(new_)   
      {
        m = new_run();
      }
      else
      {
        m = run();
      }

      for(int i = 0; i < m.cols(); i++)
      {
        if(m(0, i) == 0) 
        {
          str += '0';
        }
        else
        {
          str += '1';
        }
      }
    
    return str;
  }

 private:
  void initialization()
  {
    mc.clear();
    vars_order.clear();
    chain.clear();
    other_matrix.clear();

    Mr = matrix::Zero( 4, 2 );
    Mr << 1, 0, 0, 0, 0, 0, 0, 1;

    num_vars = 0u;
    for ( const id& pi : circuit.get_inputs() )
      {
        num_vars++;
        vars_order[ pi ] = circuit.get_inputs().size() - num_vars + 1;
      }

    other = circuit.get_nodes().size();
  }

  void get_chain_code( const id n, bool use_new = false )
  {
    mc.push_back( n );

    if ( circuit.is_pi( n ) )
      {
        return;
      }

    for ( const auto input : circuit.get_node( n ).get_inputs() )
      {
        get_chain_code( input.index );
      }
  }

  matrix get_matritx( const uint32_t n )
  {
    if ( circuit.is_pi( n ) )
      {
        uint32_t order = vars_order[ n ];
        matrix result( 2, 1 );
        result << order, order;
        return result;
      }
    if ( n >= circuit.get_nodes().size() )
      {
        assert( other_matrix.find( n ) != other_matrix.end() );
        std::string str = other_matrix[ n ];
        if ( str == "Mr" )
          {
            return Mr;
          }
        if ( str[ 0 ] == 'I' )
          {
            str.erase( str.begin() );
            int c = std::stoi( str );
            matrix identity_matrix;
            identity_matrix.setIdentity( c, c );
            return identity_matrix;
          }
        if ( str[ 0 ] == 'W' )
          {
            str.erase( str.begin() );
            assert( std::stoi( str ) == 2 );
            return generate_swap_matrix( 2, 2 );
          }
      }
    return circuit.get_node( n ).get_mtx();
  }

  void print_chain()
  {
    for ( const matrix& m : chain )
      {
        std::cout << m << "\n";
      }
  }

  void get_chain()
  {
    for ( const id& n : mc )
      {
        chain.push_back( get_matritx( n ) );
      }
  }

  /******new method******/
stp_expr move_vars_to_rightside_1(const stp_expr& inputs) {
    stp_expr new_expr;
    stp_expr temp_vars;
    int variable_count = 0;  

    for (size_t i = 0; i < inputs.size(); i++) {
      uint32_t t = inputs[i]; 
        if (!circuit.is_pi(t)) {
            if (t >= circuit.get_nodes().size()) {
                if (variable_count == 0) {
                    new_expr.push_back(inputs[i]);
                } 
                else{
                  std::string str = other_matrix[t];
                  if (str[0] == 'I') {
                      uint32_t dim;

                      size_t digit_start = str.find_first_of("0123456789");

                      if (digit_start != std::string::npos) 
                      {
                        std::string number_part = str.substr( digit_start );

                       dim = std::stoi( number_part );
                      } 
                      else 
                      {
                        assert( false && "No number found." );
                      }
                                        //uint32_t dim = (str[1] - '0');  
                    dim *= std::pow(2, variable_count);
                    std::string new_str = "I" + std::to_string(dim);
                    other_matrix[other] = new_str;
                    new_expr.push_back(other);
                    new_expr.push_back(inputs[i + 1]);
                    ++i; 
                    other++;
                  }
                  else {
                    int dim = std::pow(2, variable_count);
                    std::string new_str = "I" + std::to_string(dim);
                    other_matrix[other] = new_str;
                    new_expr.push_back(other);
                    new_expr.push_back(inputs[i]);
                    other++;
                  }
                }
            }else {
                if (variable_count == 0) {
                    new_expr.push_back(inputs[i]);
                } else {
                    int dim = std::pow(2, variable_count);
                    std::string new_str = "I" + std::to_string(dim);
                    other_matrix[other] = new_str;
                    new_expr.push_back(other);
                    new_expr.push_back(inputs[i]);
                    other++;
                }
            }
        } else {
            variable_count++;
            temp_vars.push_back(inputs[i]);
        }
    }

    new_expr.insert(new_expr.end(), temp_vars.begin(), temp_vars.end());
    return new_expr;
}

  stp_expr move_vars(stp_expr inputs) {
    stp_expr result_expr;

    while (!inputs.empty()) {
        bool pair_found = false; 
        for (size_t i = 0; i < inputs.size(); ++i) {
            uint32_t t = inputs[i];  

            if (circuit.is_pi(t)) { 
                size_t first_var = i;
                size_t second_var = std::distance(inputs.begin(), std::find(inputs.begin() + i + 1, inputs.end(), t));

                if (second_var != inputs.size()) {
                    result_expr.insert(result_expr.end(), inputs.begin(), inputs.begin() + first_var);

                    for (size_t k = first_var + 1; k < second_var; ++k) {
                        uint32_t temp_t = inputs[k];
                        if (circuit.is_pi(temp_t)) {
                            other_matrix[other] = "W2";
                            result_expr.push_back(other);
                            result_expr.push_back(inputs[k]);
                            other++;
                        }
                        else {
                          if (temp_t >= circuit.get_nodes().size()) {
                            std::string str = other_matrix[temp_t];
                            if (str[0] == 'I') {
                              uint32_t dim;

                                size_t digit_start = str.find_first_of("0123456789");

                                if (digit_start != std::string::npos) 
                                {
                                  std::string number_part = str.substr( digit_start );

                                   dim = std::stoi( number_part ) * 2;
                                } 
                                else 
                                {
                                  assert( false && "No number found." );
                                }
                                std::string new_str = "I" + std::to_string(dim);
                                other_matrix[other] = new_str;
                                result_expr.push_back(other);
                                result_expr.push_back(inputs[k + 1]);
                                ++k; 
                                other++;
                            }
                            else {
                                other_matrix[other] = "I2";
                                result_expr.push_back(other);
                                result_expr.push_back(inputs[k]);
                                other++;
                            }
                          }
                          else {
                                other_matrix[other] = "I2";
                                result_expr.push_back(other);
                                result_expr.push_back(inputs[k]);
                                other++;
                          }   
                        }
                    }

                    other_matrix[other] = "Mr";
                    result_expr.push_back(other);
                    result_expr.push_back(t);
                    other++;

                    result_expr.insert(result_expr.end(), inputs.begin() + second_var + 1, inputs.end());

                    inputs = result_expr;
                    result_expr.clear();
                    pair_found = true;
                    break;
                }
            }
        }

        if (!pair_found) {
            result_expr.insert(result_expr.end(), inputs.begin(), inputs.end());
            break;
        }
    }

    return result_expr;
  }

  stp_expr collect_pi(const stp_expr& inputs)
{
    stp_expr collected_pi;

    for (int i = inputs.size() - 1; i > 0; i--)
    {
      uint32_t t = inputs[i];
      if (circuit.is_pi(t))
      {
        collected_pi.push_back(inputs[i]);
      }
      else
      {
          break;
      }
    }

    //std::reverse(collected_pi.begin(), collected_pi.end());
    return collected_pi;
}

  stp_expr move_vars_to_rightside( const stp_expr& inputs )
  {
    stp_expr new_expr;
    stp_expr temp_vars;

    for ( int i = 0; i < inputs.size(); i++ )
      {
        uint32_t t = inputs[ i ];
        if ( circuit.is_pi( t ) )
          {
            temp_vars.push_back( t );
          }
        else
          {
            int count = get_the_number_vars_before_operation( inputs, i );
            if ( count == 0 )
              {
                new_expr.push_back( t );
              }
            else
              {
                int dim = 1 << count;
                std::string str = "I" + std::to_string( dim );
                other_matrix[ other ] = str;
                new_expr.push_back( other );
                new_expr.push_back( t );
                other++;
              }
          }
      }
    new_expr.insert( new_expr.end(), temp_vars.begin(), temp_vars.end() );
    return new_expr;
  }

  int get_the_number_vars_before_operation( const stp_expr& inputs, int op_idx )
  {
    int count = 0;
    for ( int i = 0; i < op_idx; i++ )
      {
        if ( circuit.is_pi( inputs[ i ] ) )
          {
            count++;
          }
      }
    return count;
  }

  stp_expr sort_vars()
  {
    stp_expr var_tokens = parse_vars( mc );
    stp_expr result = var_tokens;

    for ( int i = 1; i < var_tokens.size(); i++ )
      {
        uint32_t key = var_tokens[ i ];
        int j = i - 1;
        while ( j >= 0 && variable_compare( key, var_tokens[ j ] ) )
          {
            stp_expr cur_result = result;
            stp_expr vars_right_result = move_vars_to_rightside( cur_result );

            result = swap_vars( vars_right_result, j );

            std::swap( var_tokens[ j ], var_tokens[ j + 1 ] );
            --j;
          }
      }

    stp_expr temp_result =
        vars_power_reducing( move_vars_to_rightside( result ) );
    result = move_vars_to_rightside( temp_result );
    return result;
  }

  stp_expr parse_vars( const stp_expr& inputs )
  {
    stp_expr vars;
    for ( const uint32_t t : inputs )
      {
        if ( circuit.is_pi( t ) )
          {
            vars.push_back( t );
          }
      }
    return vars;
  }

  bool variable_compare( uint32_t v1, uint32_t v2 )
  {
    return vars_order[ v1 ] < vars_order[ v2 ];
  }

  stp_expr swap_vars( const stp_expr& inputs, int idx )
  {
    stp_expr result = inputs;
    bool find_var = false;
    int idx_var;
    for ( uint32_t i = 0; i < inputs.size(); i++ )
      {
        if ( circuit.is_pi( inputs[ i ] ) )
          {
            idx_var = i;
            find_var = true;
            break;
          }
      }
    assert( find_var );
    std::swap( result[ idx_var + idx ], result[ idx_var + idx + 1 ] );
    other_matrix[ other ] = "W2";
    result.insert( result.begin() + idx_var + idx, other );
    other++;
    return result;
  }

  stp_expr vars_power_reducing( const stp_expr& inputs )
  {
    stp_expr new_expr;
    stp_expr temp_vars;

    for ( int i = 0; i < inputs.size(); i++ )
      {
        if ( !circuit.is_pi( inputs[ i ] ) )
          {
            new_expr.push_back( inputs[ i ] );
          }
        else
          {
            if ( temp_vars.size() > 0 && inputs[ i ] == temp_vars.back() )
              {
                other_matrix[ other ] = "Mr";
                temp_vars.insert( temp_vars.end() - 1, other );
                other++;
              }
            else
              {
                temp_vars.push_back( inputs[ i ] );
              }
          }
      }

    new_expr.insert( new_expr.end(), temp_vars.begin(), temp_vars.end() );
    return new_expr;
  }

  matrix_chain expr_to_chain( const stp_expr& normal )
  {
    for ( const uint32_t token : normal )
      {
        chain.push_back( get_matritx( token ) );
      }

    // for the identity matrix, we should first calculate the kronecker product
    matrix_chain new_chain;

    for ( int i = 0; i < normal.size() - circuit.get_inputs().size(); i++ )
      {
        if ( other_matrix[ normal[ i ] ].substr( 0, 1 ) == "I" )
          {
            new_chain.push_back(
                kronecker_product( chain[ i ], chain[ i + 1 ] ) );
            i++;
          }
        else
          {
            new_chain.push_back( chain[ i ] );
          }
      }

    return new_chain;
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

  void print_expr( const stp_expr& m_c )
  {
    std::cout << "[stp_expr]: ";
    for ( const uint32_t t : m_c )
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
            std::cout << other_matrix[ t ] << " ";
          }
      }
    std::cout << "\n";
  }

 private:
  matrix result;
  const stp_circuit& circuit;
  bool verbose;
  matrix_chain chain;
  matrix Mr;  // power reducing
  stp_expr mc;

  uint32_t other = 0u;
  std::unordered_map<uint32_t, std::string>
      other_matrix;  // record the intermediate results W(m,n) In Mr

  uint32_t num_vars = 0u;
  std::unordered_map<uint32_t, uint32_t> vars_order;
};

}  // namespace stp
