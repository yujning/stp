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

#include "stp/stp_utils.hpp"
#include "stp/stp_circuit.hpp"
#include "stp/stp_logic_expr.hpp"

namespace stp
{
class stp_simulator
{
public:
  stp_simulator(const stp_circuit& circuit) : circuit(circuit) {} 

  //full emulation, obtaining a truth table represented by a string
  std::string run()
  {
    init();
    sim(circuit.get_outputs()[0]);
    // print_info();
    return get_tt();
  }

private: 
  void init()
  {
    num_vars = circuit.get_inputs().size();
    patterns_num = 1 << num_vars;
    sim_info.resize(circuit.get_nodes().size(), std::vector<uint8_t>(patterns_num, 0));

    had_sim.resize(circuit.get_nodes().size(), false);

    for(int i = 0; i < num_vars; i++)
    {
      std::vector<uint8_t>& patterns = sim_info[i];
      uint32_t turn = patterns_num >> (num_vars - i);
      uint8_t val = 1;
      uint32_t count = 0;
      for(int j = 0; j < patterns_num; j++)
      {
        patterns[j] = val;
        count++;
        if(count == turn)
        {
          count = 0;
          val ^= 1;
        }
      }
    }

    for(id input : circuit.get_inputs())
    {
      had_sim[input] = true;
    }
  }

  void sim(const id n)
  {
    if(had_sim[n]) return;

    for(const auto input : circuit.get_node(n).get_inputs())
    {
      sim(input.index);
    }

    compute(n);
    had_sim[n] = true;
  }

  void compute(const id n)
  {
    const matrix& mtx = circuit.get_node(n).get_mtx();
    int bits = 1 << circuit.get_node(n).get_inputs().size();
    for(int count = 0; count < patterns_num; count++)
    {
      uint32_t idx = 0;
      for(const auto& input : circuit.get_node(n).get_inputs())
      {
        idx = (idx << 1) + sim_info[input.index][count];
      }
      sim_info[n][count] = mtx(0, bits - idx - 1);
    }
  }

  std::string get_tt()
  {
    const id po = circuit.get_outputs()[0];
    std::string str = "";
    for(const uint8_t val : sim_info[po])
    {
      if(val) str += "1";
      else    str += "0";
    }
    return str;
  } 

  void print_info()
  {
    for(const auto& input_id : circuit.get_inputs())
    {
      std::cout << circuit.get_node(input_id).get_name() << " "; 
    }
    std::cout << ": ";
    for(const auto& output_id : circuit.get_outputs())
    {
      std::cout << circuit.get_node(output_id).get_name() << " "; 
    }
    std::cout << std::endl;
    for(int i = 0; i < patterns_num; i++)
    {
      for(const auto& input_id : circuit.get_inputs())
      {
        std::cout << int(sim_info[input_id][i]) << "  ";
      }
      std::cout << ":  ";
      for(const auto& output_id : circuit.get_outputs())
      {
        std::cout << int(sim_info[output_id][i]) << " "; 
      }
      std::cout << std::endl;
    }
  }
  
private:
  const stp_circuit& circuit;
  std::vector<std::vector<uint8_t>> sim_info;
  std::vector<bool> had_sim;
  uint32_t num_vars = 0u;
  uint32_t patterns_num = 0u;
};

}
