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
#include "../include/algorithms/66lut_else_dec.hpp"

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
                  "force 66-LUT bi-decomposition (no else-dec fallback)");

        add_flag("-e, --else-dec",
                 "use legacy bi-decomposition flow with else-dec fallback");
        
        add_flag("--only",
                 "run 66-LUT decomposition without strong DSD fallback");

    }

protected:
        static void build_truth_table_as_single_lut(const TT& tt, unsigned nvars)
        {
            RESET_NODE_GLOBAL();
            ORIGINAL_VAR_COUNT = static_cast<int>(nvars);

            std::vector<int> sorted_vars = tt.order;
            std::sort(sorted_vars.begin(), sorted_vars.end());
            sorted_vars.erase(std::unique(sorted_vars.begin(), sorted_vars.end()), sorted_vars.end());
            for (int var_id : sorted_vars)
                new_in_node(var_id);

            for (int var_id : tt.order)
            {
                if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), var_id) == FINAL_VAR_ORDER.end())
                    FINAL_VAR_ORDER.push_back(var_id);
            }

            std::vector<int> children;
            children.reserve(tt.order.size());

            // Preserve the variable order; write_bench will reverse the child list.
            for (int var_id : tt.order)
                children.push_back(new_in_node(var_id));

            new_node(tt.f01, children);
        }
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

            build_truth_table_as_single_lut(root_shrunk, nvars);

            auto t2 = clk::now();
            auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
            std::cout << "⏱ time = " << us << " us\n";
            return;
        }


        bool success = false;

        const bool only_dsd   = is_set("dsd");
        const bool only_bidec = is_set("bidec");
        const bool use_else_dec = is_set("else-dec");
        const bool only_lut66 = is_set("only");
        if (only_dsd && (only_bidec || use_else_dec))
        {
             std::cout << "❌ Option -d cannot be used together with -b or -e.\n";
            return;
        }
        if (only_bidec && use_else_dec)
        {
              std::cout << "❌ Options -b and -e cannot be used together.\n";
            return;
        }
         if (only_dsd || (!only_bidec && !use_else_dec))
        {
            std::cout << "🔀 Mode: 66-LUT Strong DSD (disjoint detection)\n";
            success = run_66lut_dsd_and_build_dag(root_shrunk);

            if (success || only_dsd)
            {
                auto t2 = clk::now();
                auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

                if (!success)
                {
                    std::cout << "⚠️ Decomposition failed, outputting original truth table\n";
                    build_truth_table_as_single_lut(root, nvars);
                    return;

                }

                std::cout << "⏱ time = " << us << " us\n";
                return;
            }

            std::cout << "⚠️ DSD failed, falling back to 66-LUT bi-decomposition...\n";
        }

        std::cout << "🔀 Mode: 66-LUT Bi-Decomposition\n";
        success = run_strong_bi_dec_and_build_dag(root_shrunk);

        if (!success && use_else_dec && !only_lut66)
        {
            std::cout << "⚠️ Bi-Decomposition failed, trying DSD (-f -s) fallback...\n";
            success = run_66lut_else_dec_and_build_dag(root_shrunk);
        }

           if (!success && (only_lut66 || only_bidec))
        {
            std::cout << "⚠️ Decomposition failed, outputting original truth table\n";
            build_truth_table_as_single_lut(root, nvars);
            return;
        }

        if (!success)
        {
            std::cout << "❌ Decomposition failed\n";
            return;
        }

                auto end = clk::now();
        auto final_us = std::chrono::duration_cast<std::chrono::microseconds>(end - t1).count();
        std::cout << "⏱ time = " << final_us << " us\n";
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
