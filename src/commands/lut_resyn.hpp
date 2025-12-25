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


#include "../include/algorithms/mix_dsd.hpp"
#include "../include/algorithms/strong_dsd.hpp"
#include "../include/algorithms/stp_dsd.hpp"
namespace alice
{
enum class resyn_strategy
{
    bi_dec,
    dsd,
    strong_dsd,
    mix_dsd
};


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

static int run_dsd_for_resyn(const std::string& binary01, bool enable_else_dec)
{
    if (!is_power_of_two(binary01.size()))
        throw std::runtime_error("input length must be power of two");

    int n = static_cast<int>(std::log2(binary01.size()));

    bool prev_enable_else = ENABLE_ELSE_DEC;
    ENABLE_ELSE_DEC = enable_else_dec;

    RESET_NODE_GLOBAL();
    ORIGINAL_VAR_COUNT = n;

    TT root;
    root.f01 = binary01;
    root.order.resize(n);

    for (int i = 0; i < n; ++i)
        root.order[i] = n - i;

    for (int v = 1; v <= n; ++v)
        new_in_node(v);

    TT root_shrunk = shrink_to_support(root);
    int root_id = dsd_factor(root_shrunk);

    ENABLE_ELSE_DEC = prev_enable_else;

    return root_id;
}

static int run_strong_dsd_for_resyn(const std::string& binary01)
{
    if (!is_power_of_two(binary01.size()))
        throw std::runtime_error("input length must be power of two");

    int n = static_cast<int>(std::log2(binary01.size()));

    RESET_NODE_GLOBAL();
    ORIGINAL_VAR_COUNT = n;

    TT root;
    root.f01 = binary01;
    root.order.resize(n);

    for (int i = 0; i < n; ++i)
        root.order[i] = n - i;

    for (int v = 1; v <= n; ++v)
        new_in_node(v);

    TT root_shrunk = shrink_to_support(root);
    return build_strong_dsd_nodes(root_shrunk.f01, root_shrunk.order, 0);
}

static int run_mix_dsd_for_resyn(const std::string& binary01, bool enable_else_dec)
{
    if (!is_power_of_two(binary01.size()))
        throw std::runtime_error("input length must be power of two");

    int n = static_cast<int>(std::log2(binary01.size()));

        if (enable_else_dec)
    {
        return run_dsd_for_resyn(binary01, true);
    }
    RESET_NODE_GLOBAL();
    ORIGINAL_VAR_COUNT = n;

    TT root;
    root.f01 = binary01;
    root.order.resize(n);

    for (int i = 0; i < n; ++i)
        root.order[i] = n - i;

    for (int v = 1; v <= n; ++v)
        new_in_node(v);

    TT root_shrunk = shrink_to_support(root);
    return dsd_factor_mix_impl(root_shrunk, 0, nullptr, nullptr);
}

class lut_resyn_command : public command
{
public:
    explicit lut_resyn_command(const environment::ptr& env)
        : command(env, "Resynthesize LUT netlist to 2-LUT")
    {
        add_option("file", input_file, "BENCH file");
        add_option("-o,--output", output_file, "output BENCH file")->required();

        add_flag("-b,--bi_dec", use_bi_dec,
                 "use bi-decomposition (default)");
        add_flag("-d,--dsd", use_dsd, "use DSD-based decomposition");
        add_flag("-s,--strong", use_strong_dsd,
                 "use strong DSD (requires -d)");
        add_flag("-m,--mix", use_mix_dsd,
                 "use mixed DSD (requires -d)");
        add_flag("-e,--else_dec", use_else_dec,
                "enable else_dec fallback (BD) or DSD exact+Shannon (-d)");
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

               bool algorithm_selected = use_bi_dec || use_dsd;
        if (use_strong_dsd || use_mix_dsd)
            algorithm_selected = true;

        if (!algorithm_selected)
            use_bi_dec = true;

        if (use_bi_dec && use_dsd)
        {
            std::cout << "❌ Options -b and -d cannot be used together.\n";
            return;
        }

         if (!use_dsd && (use_strong_dsd || use_mix_dsd))
        {
             std::cout << "❌ Options -s or -m require -d.\n";
            return;
        }

         if (use_else_dec && use_strong_dsd)
        {
              std::cout << "❌ Option -e cannot be combined with -s.\n";
            return;
        }

        if (use_strong_dsd && use_mix_dsd)
        {
            std::cout << "❌ Options -s and -m cannot be used together.\n";
            return;
        }

        resyn_strategy strategy = resyn_strategy::bi_dec;
        if (use_dsd)
        {
            strategy = resyn_strategy::dsd;
            if (use_strong_dsd)
                strategy = resyn_strategy::strong_dsd;
            else if (use_mix_dsd)
                strategy = resyn_strategy::mix_dsd;
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

        // 3️⃣ 非 LUT 行（vdd / gnd / assign 等）
        for (const auto& line : net.passthrough_lines)
        {
            fout << line << "\n";
        }

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
                                switch (strategy)
                {
                case resyn_strategy::bi_dec:
                    root_id = run_bi_decomp_for_resyn(binary01, use_else_dec);
                    break;
                case resyn_strategy::dsd:
                    root_id = run_dsd_for_resyn(binary01, use_else_dec);
                    break;
                case resyn_strategy::strong_dsd:
                    root_id = run_strong_dsd_for_resyn(binary01);
                    break;
                case resyn_strategy::mix_dsd:
                    root_id = run_mix_dsd_for_resyn(binary01, use_else_dec);
                    break;
                }
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
    bool use_bi_dec = false;
    bool use_dsd = false;
    bool use_strong_dsd = false;
    bool use_mix_dsd = false;
    bool use_else_dec = false;
};

ALICE_ADD_COMMAND(lut_resyn, "STP")

} // namespace alice

#endif // LUT_RESYN_HPP