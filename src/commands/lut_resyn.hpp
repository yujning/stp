#ifndef LUT_RESYN_HPP
#define LUT_RESYN_HPP

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <alice/alice.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>

#include "../include/algorithms/bench_lut.hpp"
#include "../include/algorithms/bi_decomposition.hpp"
#include "../include/algorithms/node_global.hpp"
#include "../include/algorithms/truth_table.hpp"

namespace alice
{


static int run_bi_decomp_for_resyn(const std::string& binary01, bool enable_else_dec)
{
    bool prev_minimal_output = BD_MINIMAL_OUTPUT;
    bool prev_enable_else = ENABLE_ELSE_DEC;

    RESET_NODE_GLOBAL();
    ENABLE_ELSE_DEC = enable_else_dec;
    BD_MINIMAL_OUTPUT = true;

    if (!is_power_of_two(binary01.size()))
        throw std::runtime_error("input length must be power of two");

    int n = static_cast<int>(std::log2(binary01.size()));
    ORIGINAL_VAR_COUNT = n;

    TT root;
    root.f01 = binary01;
    root.order.resize(n);
    for (int i = 0; i < n; ++i)
        root.order[i] = n - i;

    for (int v = 1; v <= n; ++v)
        new_in_node(v);

    TT root_shrunk = shrink_to_support(root);
    int root_id = bi_decomp_recursive(root_shrunk, 0);

    BD_MINIMAL_OUTPUT = prev_minimal_output;
    ENABLE_ELSE_DEC = prev_enable_else;

    return root_id;
}

class lut_resyn_command : public command
{
public:
    explicit lut_resyn_command(const environment::ptr& env)
        : command(env, "Resynthesize LUT netlist to 2-LUT")
    {
        add_option("file", input_file, "BENCH file");
        add_option("-o,--output", output_file, "output BENCH file")->required();
        add_flag("-e,--else_dec", use_else_dec,
                 "enable else_dec fallback when BD fails");
    }

protected:
    void execute() override
    {
        if (output_file.empty())
        {
             std::cout << "❌ Output file missing\n";
            return;
        }

               BenchNetlist net;
        if (input_file.empty())
        {
            if (!BENCH_LOADED)
            {
                std::cout << "❌ No BENCH loaded. Run read_bench or provide file.\n";
                return;
            }
            net = BENCH_NETLIST;
        }
        else
        {
            try
            {
                net = read_bench_lut(input_file);
            }
            catch (const std::exception& e)
            {
                std::cout << "❌ " << e.what() << "\n";
                return;
            }
            BENCH_NETLIST = net;
            BENCH_LOADED = true;
            BENCH_SOURCE = input_file;
        }

        std::ofstream fout(output_file);
        if (!fout)
        {
            std::cout << "❌ Cannot open " << output_file << "\n";
            return;
        }

        for (const auto& in : net.inputs)
            fout << "INPUT(" << in << ")\n";
        for (const auto& out : net.outputs)
            fout << "OUTPUT(" << out << ")\n";
        fout << "\n";

        std::vector<std::string> lut_names;
        lut_names.reserve(net.luts.size());
        for (const auto& kv : net.luts)
            lut_names.push_back(kv.first);
        std::sort(lut_names.begin(), lut_names.end());

        int unique_id = 0;

        for (const auto& name : lut_names)
        {
            const auto& lut = net.luts.at(name);

            if (lut.fanins.size() <= 2)
            {
                fout << name << " = LUT " << lut.hex << " (";
                for (size_t i = 0; i < lut.fanins.size(); ++i)
                {
                    fout << lut.fanins[i];
                    if (i + 1 < lut.fanins.size())
                        fout << ", ";
                }
                fout << ")\n";
                continue;
            }

            std::string binary01;
            try
            {
                binary01 = hex_to_binary(lut.hex);
            }
            catch (const std::exception& e)
            {
                std::cout << "❌ " << name << ": " << e.what() << "\n";
                return;
            }

            int root_id = 0;
            try
            {
                root_id = run_bi_decomp_for_resyn(binary01, use_else_dec);
            }
            catch (const std::exception& e)
            {
                std::cout << "❌ " << name << ": " << e.what() << "\n";
                return;
            }

            std::map<int, std::string> name_of;
            for (const auto& node : NODE_LIST)
            {
                if (node.func == "in")
                {
                    int index = node.var_id - 1;
                    if (index < 0 || index >= static_cast<int>(lut.fanins.size()))
                    {
                        std::cout << "❌ " << name << ": invalid fanin index\n";
                        return;
                    }
                    name_of[node.id] = lut.fanins[index];
                }
            }

            for (const auto& node : NODE_LIST)
            {
                if (node.func == "in")
                    continue;

                if (node.id == root_id)
                {
                    name_of[node.id] = name;
                }
                else
                {
                    ++unique_id;
                    name_of[node.id] = name + "_d" + std::to_string(unique_id);
                }
            }

            for (const auto& node : NODE_LIST)
            {
                if (node.func == "in")
                    continue;

                fout << name_of[node.id] << " = LUT 0x"
                     << bin_to_hex(node.func) << " (";

                std::vector<int> rev = node.child;
                std::reverse(rev.begin(), rev.end());

                for (size_t i = 0; i < rev.size(); ++i)
                {
                    fout << name_of[rev[i]];
                    if (i + 1 < rev.size())
                        fout << ", ";
                }
                fout << ")\n";
            }

            fout << "\n";
        }

        std::cout << "✅ LUT resynthesis written to " << output_file << "\n";
    }

private:
    std::string input_file{};
    std::string output_file{};
    bool use_else_dec = false;
};

ALICE_ADD_COMMAND(lut_resyn, "STP")

} // namespace alice

#endif // LUT_RESYN_HPP