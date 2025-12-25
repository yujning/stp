#ifndef DSD_HPP
#define DSD_HPP

#include <iostream>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <alice/alice.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>

// 引入算法
#include "../include/algorithms/stp_dsd.hpp"
#include "../include/algorithms/reorder.hpp"   // all_reorders_char(raw)
#include "../include/algorithms/strong_dsd.hpp"
#include "../include/algorithms/mix_dsd.hpp"
namespace alice
{
    class dsd_command : public command
    {
    public:
        explicit dsd_command(const environment::ptr &env)
            : command(env, "Run STP decomposition / raw reorder")
        {
            add_option("-f, --factor", hex_input,
                       "hexadecimal number (must map to 2^n bits)");

            add_option("-x, --raw", raw_input,
                       "raw truth-table with don't care (x), length must be 2^n");
            
            
            add_flag("-s, --strong",
                     "use strong DSD: find first L=2^k with exactly two block types，ACD ");


            add_flag("-m, --mix",
             "prefer DSD (-f) per layer, fallback to strong DSD when needed");

             
            add_flag("-e, --else",
                     "enable exact synthesis + Shannon fallback for -f");
        }

    protected:
  void execute() override
{
    using clk = std::chrono::high_resolution_clock;

    bool use_raw = is_set("raw");
    bool use_hex = is_set("factor");
    bool use_strong = is_set("strong");
    bool use_mix = is_set("mix");
    bool use_else_dec = is_set("else");

    if (use_raw && use_hex)
    {
        std::cout << "❌ Options -f and -x cannot be used together.\n";
        return;
    }

    if (use_else_dec && use_raw)
    {
        std::cout << "❌ Option -e only supports hexadecimal mode (-f).\n";
        return;
    }

        if (use_strong && use_mix)
    {
        std::cout << "❌ Options -s and -m cannot be used together.\n";
        return;
    }


    // ------------ RAW MODE ------------
    if (use_raw)
    {
        std::string raw = raw_input;

        if (!is_power_of_two(raw.size()))
        {
            std::cout << "❌ Error: length (" << raw.size() << ") is not power of 2\n";
            return;
        }
         
        if (raw.size() == 4) {
           std::cout << "⚠ 输入函数已是 2-LUT，无需分解。\n";
        return;
        }
        std::cout << "➡ Raw truth-table mode (-x)\n";
        std::cout << "Input = " << raw << "\n";

        if (use_strong) {
            if (raw.find_first_not_of("01") != std::string::npos) {
                std::cout << "❌ Strong DSD requires raw input of only 0/1.\n";
                return;
            }
            auto t1 = clk::now();
    if (use_else_dec) {
                run_dsd_recursive(raw, true);
            } else {
                RESET_NODE_GLOBAL();
                ORIGINAL_VAR_COUNT = static_cast<int>(std::log2(raw.size()));
                std::vector<int> order;
                for (int i = ORIGINAL_VAR_COUNT; i >= 1; --i) {
                    order.push_back(i);
                }
                for (int v = 1; v <= ORIGINAL_VAR_COUNT; ++v) {
                    new_in_node(v);
                }
                build_strong_dsd_nodes(raw, order, 0);
            }
 
            auto t2 = clk::now();
            auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
            std::cout << "⏱ Strong DSD time = " << us << " us\n";
            return;
        }


        if (use_mix) {
            if (raw.find_first_not_of("01") != std::string::npos) {
                std::cout << "❌ Mixed DSD requires raw input of only 0/1.\n";
                return;
            }
            auto t1 = clk::now();
                        if (use_else_dec) {
                run_dsd_recursive(raw, true);
            } else {
                run_dsd_recursive_mix(raw);
            }
            auto t2 = clk::now();
            auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
            std::cout << "⏱ Mixed DSD time = " << us << " us\n";
            return;
        }


        auto t1 = clk::now();
        all_reorders(raw);
        auto t2 = clk::now();

        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        std::cout << "⏱ RAW decomposition time = " << us << " us\n";

        return;
    }

    // ------------ HEX MODE ------------
    if (!use_hex)
    {
        std::cout << "❌ Please use -f <hex> or -x <raw>\n";
        return;
    }

    std::string hex = hex_input;

    if (hex.rfind("0x",0)==0 || hex.rfind("0X",0)==0)
        hex = hex.substr(2);

    unsigned bit_count = hex.size() * 4;
    unsigned num_vars = 0;
    while ((1u << num_vars) < bit_count) num_vars++;

    if ((1u << num_vars) != bit_count)
    {
        std::cout << "❌ Hex length is not 2^n bits\n";
        return;
    }

        if (bit_count == 4) {
        std::cout << "⚠ 输入函数已是 2-LUT，无需分解。\n";
        return;
    }

    kitty::dynamic_truth_table tt(num_vars);
    kitty::create_from_hex_string(tt, hex);

    std::ostringstream oss;
    kitty::print_binary(tt, oss);
    std::string binary = oss.str();

    std::cout << "📘 Hex " << hex << " => binary " << binary
              << " (len = " << binary.size() << " vars = " << num_vars << ")\n";

    
    if (use_strong) {
        auto t1 = clk::now();
  if (use_else_dec) {
            run_dsd_recursive(binary, true);
        } else {
            RESET_NODE_GLOBAL();
            ORIGINAL_VAR_COUNT = static_cast<int>(std::log2(binary.size()));
            std::vector<int> order;
            for (int i = ORIGINAL_VAR_COUNT; i >= 1; --i) {
                order.push_back(i);
            }
            for (int v = 1; v <= ORIGINAL_VAR_COUNT; ++v) {
                new_in_node(v);
            }
            build_strong_dsd_nodes(binary, order, 0);
        }

        auto t2 = clk::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        std::cout << "⏱ Strong DSD time = " << us << " us\n";
        return;
    }

    if (use_mix) {
        auto t1 = clk::now();
        if (use_else_dec) {
            run_dsd_recursive(binary, true);
        } else {
            run_dsd_recursive_mix(binary);
        }
        auto t2 = clk::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        std::cout << "⏱ Mixed DSD time = " << us << " us\n";
        return;
    }


    // 计时
    auto t1 = clk::now();
    //all_reorders(binary);
     run_dsd_recursive(binary, use_else_dec);
    auto t2 = clk::now();

    auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "⏱ DSD execution time = " << us << " us\n";
}

    private:
        std::string hex_input{};
        std::string raw_input{};
        bool use_else_dec = false;
    };

    ALICE_ADD_COMMAND(dsd, "STP")
}

#endif
