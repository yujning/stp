#ifndef BI_D_HPP
#define BI_D_HPP

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <sstream>

#include <alice/alice.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>

#include "../include/algorithms/node_global.hpp"
#include "../include/algorithms/bi_decomposition.hpp"
#include "../include/algorithms/bi_dec_else_dec.hpp"

// 注意：这里不再 include run_dsd_recursive_mix
// DSD 只通过 bi_decomposition 内部的 one-layer split 使用

namespace alice
{

class bd_command : public command
{
public:
    explicit bd_command(const environment::ptr& env)
        : command(env, "Bi-decomposition (recursive, BD主体)")
    {
        add_option("-f, --factor", hex_input,
                   "truth table as hex string")->required();

        add_flag("-d, --k2_zero", only_k2_zero,
                 "only try bi-decomposition with k2 = 0");

        add_flag("-e, --else_dec", use_else_dec,
                 "enable else_dec fallback (Shannon / exact)");

        add_flag("--dm, --dsd_mix", use_dsd_mix,
                 "enable DSD one-layer fallback inside BD");
    }

protected:
    void execute() override
    {
        using clk = std::chrono::high_resolution_clock;

        use_else_dec = is_set("else_dec");
        use_dsd_mix  = is_set("dsd_mix") || is_set("dm");
        only_k2_zero = is_set("k2_zero");

        // ------------------------------------------------------
        // Parse hex → binary truth table
        // ------------------------------------------------------
        std::string hex = hex_input;
        if (hex.rfind("0x", 0) == 0 || hex.rfind("0X", 0) == 0)
            hex = hex.substr(2);

        unsigned bits = static_cast<unsigned>(hex.size() * 4);
        if (bits == 0)
        {
            std::cout << "❌ Empty truth table\n";
            return;
        }

        unsigned nvars = 0;
        while ((1u << nvars) < bits) ++nvars;

        if ((1u << nvars) != bits)
        {
            std::cout << "❌ TT size is not 2^n\n";
            return;
        }

        kitty::dynamic_truth_table tt(nvars);
        kitty::create_from_hex_string(tt, hex);

        std::ostringstream oss;
        kitty::print_binary(tt, oss);
        std::string binF = oss.str();

        std::cout << "📘 TT = " << binF
                  << "  (vars=" << nvars << ")\n";

        // ------------------------------------------------------
        // Set global control flags (used inside BD recursion)
        // ------------------------------------------------------
        ENABLE_ELSE_DEC            = use_else_dec;
        BD_ENABLE_DSD_MIX_FALLBACK = use_dsd_mix;
        BD_ONLY_K2_EQ_0            = only_k2_zero;

        // ------------------------------------------------------
        // Run BD (single entry point)
        // ------------------------------------------------------
        auto t1 = clk::now();
        bool success = run_bi_decomp_recursive(binF);
        auto t2 = clk::now();

        auto elapsed =
            std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

        if (success)
            std::cout << "⏱ time = " << elapsed << " us\n";
        else
            std::cout << "❌ Decomposition failed\n";
    }

private:
    std::string hex_input{};
    bool use_else_dec  = false;
    bool only_k2_zero  = false;
    bool use_dsd_mix   = false;
};

ALICE_ADD_COMMAND(bd, "STP")

} // namespace alice

#endif // BI_D_HPP
