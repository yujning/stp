#ifndef LUT_66_HPP
#define LUT_66_HPP

#include <chrono>
#include <iostream>
#include <sstream>
#include <string>

#include <alice/alice.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>

#include "../include/algorithms/node_global.hpp"
#include "../include/algorithms/66lut_bidec.hpp"
#include "../include/algorithms/66lut_dsd.hpp"   // ★ 新增

namespace alice
{

class lut_66_command : public command
{
public:
    explicit lut_66_command(const environment::ptr& env)
        : command(env, "66-LUT decomposition (bi-dec or strong DSD)")
    {
        add_option("-f, --factor", hex_input,
                   "truth table as hex string")->required();

        add_flag("-d, --dsd",
                 "use 66-LUT Strong DSD instead of bi-decomposition");
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

        auto t1 = clk::now();

        bool success = false;

        if (is_set("dsd"))
        {
            std::cout << "🔀 Mode: 66-LUT Strong DSD\n";
            success = run_66lut_dsd_and_build_dag(binF);
        }
        else
        {
            std::cout << "🔀 Mode: 66-LUT Bi-Decomposition\n";
            success = run_strong_bi_dec_and_build_dag(binF);
        }

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
