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

#include "../include/algorithms/stp_dsd.hpp"
#include "../include/algorithms/strong_dsd.hpp"
#include "../include/algorithms/mix_dsd.hpp"

#include "../include/algorithms/66lut_bidec.hpp"
#include "../include/algorithms/66lut_dsd.hpp"


namespace alice
{

enum class resyn_strategy
{
    bi_dec,
    dsd,
    strong_dsd,
    mix_dsd
};

/*============================================================*
 * BI-DECOMPOSITION（保持原样）
 *============================================================*/
static int run_bi_decomp_for_resyn(const std::string& binary01, bool enable_else_dec,
                                   bool enable_dsd_mix_fallback)
{
    bool prev_minimal_output = BD_MINIMAL_OUTPUT;
    bool prev_enable_else    = ENABLE_ELSE_DEC;
    bool prev_dsd_mix        = BD_ENABLE_DSD_MIX_FALLBACK;

    RESET_NODE_GLOBAL();
    ENABLE_ELSE_DEC   = enable_else_dec;
    BD_ENABLE_DSD_MIX_FALLBACK = enable_dsd_mix_fallback;
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
    ENABLE_ELSE_DEC   = prev_enable_else;
    BD_ENABLE_DSD_MIX_FALLBACK = prev_dsd_mix;

    return root_id;
}

/*============================================================*
 * STRONG DSD（保持原样）
 *============================================================*/
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

static bool run_lut66_for_resyn(const std::string& binary01, int nvars, int& root_id)
{
    TT root;
    root.f01 = binary01;
    root.order.resize(nvars);
    for (int i = 0; i < nvars; ++i)
        root.order[i] = nvars - i;

    TT root_shrunk = shrink_to_support(root);
    const unsigned shrunk_vars = static_cast<unsigned>(root_shrunk.order.size());

    if (shrunk_vars <= 6)
    {
        RESET_NODE_GLOBAL();
        ORIGINAL_VAR_COUNT = static_cast<int>(shrunk_vars);

        std::vector<int> sorted_vars = root_shrunk.order;
        std::sort(sorted_vars.begin(), sorted_vars.end());
        sorted_vars.erase(std::unique(sorted_vars.begin(), sorted_vars.end()), sorted_vars.end());
        for (int var_id : sorted_vars)
            new_in_node(var_id);

        for (int var_id : root_shrunk.order)
        {
            if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), var_id) == FINAL_VAR_ORDER.end())
                FINAL_VAR_ORDER.push_back(var_id);
        }

        std::vector<int> children;
        children.reserve(root_shrunk.order.size());

        for (int var_id : root_shrunk.order)
            children.push_back(new_in_node(var_id));

        new_node(root_shrunk.f01, children);

        root_id = NODE_LIST.empty() ? 0 : NODE_LIST.back().id;
        return true;
    }

    bool success = run_66lut_dsd_and_build_dag(root_shrunk);
    if (!success)
        success = run_strong_bi_dec_and_build_dag(root_shrunk);

    if (success)
        root_id = NODE_LIST.empty() ? 0 : NODE_LIST.back().id;

    return success;
}

/*============================================================*
 * LUT RESYN COMMAND
 *============================================================*/
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
        add_flag("-d,--dsd", use_dsd,
                 "use DSD-based decomposition");
        add_flag("-s,--strong", use_strong_dsd,
                 "use strong DSD (requires -d)");
        add_flag("-m,--mix", use_mix_dsd,
                 "use mixed DSD (requires -d)");
        add_flag("-e,--else_dec", use_else_dec,
                 "exact synthesis + Shannon fallback (DSD)");
        
        add_flag("--dm,--dsd_mix", use_dsd_mix_fallback,
                 "mixed DSD (-m) fallback when BD cannot proceed (-b)");

        
        add_flag("--lut66", use_lut66,
                 "per-LUT 66-LUT decomposition (same as `lut66 -f`)");
    }

protected:
    void execute() override
    {
        if (output_file.empty())
        {
            std::cout << "❌ Output file missing\n";
            return;
        }

        /*---------------- BENCH load ----------------*/
        BenchNetlist net;
        if (input_file.empty())
        {
            if (!BENCH_LOADED)
            {
                std::cout << "❌ No BENCH loaded\n";
                return;
            }
            net = BENCH_NETLIST;
        }
        else
        {
            net = read_bench_lut(input_file);
            BENCH_NETLIST = net;
            BENCH_LOADED  = true;
            BENCH_SOURCE  = input_file;
        }

        /*---------------- option check ----------------*/
        if (use_lut66 && (use_bi_dec || use_dsd || use_strong_dsd || use_mix_dsd|| use_else_dec || use_dsd_mix_fallback))
        {
             std::cout << "❌ --lut66 cannot be combined with other options\n";
            return;
        }

        resyn_strategy strategy = resyn_strategy::bi_dec;
        if (!use_lut66)
        {
            bool algorithm_selected = use_bi_dec || use_dsd;
            if (use_strong_dsd || use_mix_dsd)
                algorithm_selected = true;
            if (!algorithm_selected)
                use_bi_dec = true;

            if (use_dsd_mix_fallback && !use_bi_dec)
            {
                std::cout << "❌ --dm requires -b (bi-decomposition)\n";
                return;
            }

 if (use_bi_dec && use_dsd)
            {
                std::cout << "❌ -b and -d cannot be used together\n";
                return;
            }
            if (!use_dsd && (use_strong_dsd || use_mix_dsd))
            {
                std::cout << "❌ -s / -m require -d\n";
                return;
            }

            if (use_strong_dsd && use_mix_dsd)
            {
                std::cout << "❌ -s and -m cannot be used together\n";
                return;
            }

            if (use_dsd)
            {
                strategy = resyn_strategy::dsd;
                if (use_strong_dsd)
                    strategy = resyn_strategy::strong_dsd;
                else if (use_mix_dsd)
                    strategy = resyn_strategy::mix_dsd;
            }
        }

        /*---------------- output file ----------------*/
        std::ofstream fout(output_file);
        if (!fout)
        {
            std::cout << "❌ Cannot open " << output_file << "\n";
            return;
        }

        /*---------------- write IO ----------------*/
        for (const auto& in : net.inputs)
            fout << "INPUT(" << in << ")\n";
        for (const auto& out : net.outputs)
            fout << "OUTPUT(" << out << ")\n";
        fout << "\n";

        for (const auto& line : net.passthrough_lines)
            fout << line << "\n";
        fout << "\n";

        /*---------------- LUT order ----------------*/
        std::vector<std::string> lut_names;
        for (const auto& kv : net.luts)
            lut_names.push_back(kv.first);
        std::sort(lut_names.begin(), lut_names.end());

        int unique_id = 0;

        /*============================================================*
         * MAIN LOOP
         *============================================================*/
        for (const auto& name : lut_names)
        {
            const auto& lut = net.luts.at(name);

/*------------------------------------------------------------*
 * RAW passthrough for --lut66 (same as other commands)
 *------------------------------------------------------------*/
if (use_lut66 && lut.fanins.size() <= 6)

{
    fout << name << " = LUT " << lut.hex << " (";
    for (size_t i = 0; i < lut.fanins.size(); ++i)
    {
        if (i) fout << ", ";
        fout << lut.fanins[i];
    }
    fout << ")\n";
    continue;
}

/* already 2-LUT (non-lut66 paths) */
if (lut.fanins.size() <= 2)
{
    fout << name << " = LUT " << lut.hex << " (";
    for (size_t i = 0; i < lut.fanins.size(); ++i)
    {
        if (i) fout << ", ";
        fout << lut.fanins[i];
    }
    fout << ")\n";
    continue;
}


            std::string binary01 = hex_to_binary(lut.hex);
            int root_id = 0;

            /*---------------- choose algorithm ----------------*/
      if (use_lut66)
            {
                         bool success = run_lut66_for_resyn(binary01,
                                                static_cast<int>(lut.fanins.size()),
                                                root_id);
                if (!success)
                {
                    std::cout << "❌ 66-LUT decomposition failed for " << name << "\n";
                    return;
                }
            }
            else
            {
                switch (strategy)
                {
                case resyn_strategy::bi_dec:
                    root_id = run_bi_decomp_for_resyn(binary01, use_else_dec,use_dsd_mix_fallback);
                    break;

                case resyn_strategy::dsd:
                    root_id = run_dsd_recursive(binary01, use_else_dec);
                    break;

                case resyn_strategy::strong_dsd:
                     if (use_else_dec)
                        root_id = run_dsd_recursive(binary01, true);
                    else
                        root_id = run_strong_dsd_for_resyn(binary01);
                    break;

                case resyn_strategy::mix_dsd:
                    if (use_else_dec)
                        root_id = run_dsd_recursive(binary01, true);   // -m -e
                    else
                        root_id = run_dsd_recursive_mix(binary01);    // -m
                    break;
                }
            }

            /*---------------- name binding ----------------*/
            std::map<int, std::string> name_of;

            for (const auto& node : NODE_LIST)
            {
                if (node.func == "in")
                {
                    int index = node.var_id - 1;
                    name_of[node.id] = lut.fanins[index];
                }
            }

            for (const auto& node : NODE_LIST)
            {
                if (node.func == "in") continue;

                if (node.id == root_id)
                    name_of[node.id] = name;
                else
                    name_of[node.id] = name + "_d" + std::to_string(++unique_id);
            }

            /*---------------- write LUTs (格式不变) ----------------*/
            for (const auto& node : NODE_LIST)
            {
                if (node.func == "in") continue;

                fout << name_of[node.id] << " = LUT 0x"
                     << bin_to_hex(node.func) << " (";

                std::vector<int> rev = node.child;
                std::reverse(rev.begin(), rev.end());

                for (size_t i = 0; i < rev.size(); ++i)
                {
                    if (i) fout << ", ";
                    fout << name_of[rev[i]];
                }
                fout << ")\n";
            }

            fout << "\n";
        }

        std::cout << "✅ LUT resynthesis written to " << output_file << "\n";
    }

private:
    std::string input_file;
    std::string output_file;
    bool use_bi_dec = false;
    bool use_dsd = false;
    bool use_strong_dsd = false;
    bool use_mix_dsd = false;
    bool use_else_dec = false;
    bool use_dsd_mix_fallback = false;
    bool use_lut66 = false;
};

ALICE_ADD_COMMAND(lut_resyn, "STP")

} // namespace alice

#endif // LUT_RESYN_HPP
