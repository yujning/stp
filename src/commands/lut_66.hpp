#ifndef LUT_66_HPP
#define LUT_66_HPP

#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>

#include <alice/alice.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>

#include "../include/algorithms/node_global.hpp"
#include "../include/algorithms/66lut_bidec.hpp"
#include "../include/algorithms/66lut_dsd.hpp"   // ★ 新增
#include "../include/algorithms/stp_dsd.hpp"
namespace alice
{

class lut_66_command : public command
{
public:
    explicit lut_66_command(const environment::ptr& env)
        : command(env, "66-LUT decomposition (disjoint-first)")
    {
        add_option("-f, --factor", hex_input,
                   "truth table as hex string")->required();

        add_flag("-d, --dsd",
                        "only run 66-LUT Strong DSD (disjoint detection)");

        add_flag("-b, --bidec",
                 "force 66-LUT bi-decomposition (legacy -f behavior)");
    }

protected:
    void execute() override
    {
        using clk = std::chrono::high_resolution_clock;

        std::string hex = hex_input;
        if (hex.rfind("0x", 0) == 0 || hex.rfind("0X", 0) == 0)
            hex = hex.substr(2);

        unsigned bits = hex.size() * 4;
        if (bits == 0)
        {
            std::cout << "❌ Empty truth table.\n";
            return;
        }

        unsigned nvars = 0;
        while ((1u << nvars) < bits) nvars++;

        if ((1u << nvars) != bits)
        {
            std::cout << "❌ TT size is not 2^n.\n";
            return;
        }

        kitty::dynamic_truth_table tt(nvars);
        kitty::create_from_hex_string(tt, hex);

        std::ostringstream oss;
        kitty::print_binary(tt, oss);
        std::string binF = oss.str();

        std::cout << "📘 TT = " << binF
                  << "  (vars=" << nvars << ")\n";

                      TT root;
        root.f01 = binF;
        root.order.resize(nvars);
        for (unsigned i = 0; i < nvars; ++i)
            root.order[i] = static_cast<int>(nvars - i);

        TT root_shrunk = shrink_to_support(root);
        const unsigned shrunk_vars = static_cast<unsigned>(root_shrunk.order.size());

        auto t1 = clk::now();

 if (shrunk_vars <= 6)
{
    std::cout << "🔀 Mode: single 6-LUT (no decomposition needed)\n";

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

    // Preserve the shrunk variable order when emitting the single 6-LUT.
    // write_bench() will reverse the child list again, so we keep the
    // original order here (MSB->LSB) to ensure the final BENCH matches
    // the printed variable mapping.
    for (int var_id : root_shrunk.order)
        children.push_back(new_in_node(var_id));

    new_node(root_shrunk.f01, children);

    auto t2 = clk::now();
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "⏱ time = " << us << " us\n";
            return;
        }


        bool success = false;

        const bool only_dsd   = is_set("dsd");
        const bool only_bidec = is_set("bidec");

        if (only_dsd && only_bidec)
        {
             std::cout << "❌ Options -d and -b cannot be used together.\n";
            return;
        }
        if (only_dsd || !only_bidec)
        {
            std::cout << "🔀 Mode: 66-LUT Strong DSD (disjoint detection)\n";
            success = run_66lut_dsd_and_build_dag(root_shrunk);

            if (success || only_dsd)
            {
                auto t2 = clk::now();
                auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

                if (!success)
                {
                    std::cout << "❌ Decomposition failed\n";
                    return;
                }

                std::cout << "⏱ time = " << us << " us\n";
                return;
            }

            std::cout << "⚠️ DSD failed, falling back to 66-LUT bi-decomposition...\n";
        }

        std::cout << "🔀 Mode: 66-LUT Bi-Decomposition (-b legacy)\n";
        success = run_strong_bi_dec_and_build_dag(root_shrunk);

        auto t2 = clk::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

        if (!success)
        {
            std::cout << "❌ Decomposition failed\n";
            return;
        }

        std::cout << "⏱ time = " << us << " us\n";
    }

private:
    std::string hex_input{};
};

struct lut_66_command_init
{
    lut_66_command_init()
    {
        alice_globals::get().command_names.emplace_back("66l", "STP");
    }
};
static lut_66_command_init _lut_66_command_init;

_ALICE_ADD_TO_LIST(alice_commands, lut_66_command)

} // namespace alice

#endif
