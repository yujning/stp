#ifndef BI_D_HPP
#define BI_D_HPP

#include <iostream>
#include <iomanip>
#include <cmath>
#include <alice/alice.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>

#include "../include/algorithms/node_global.hpp"
#include "../include/algorithms/stp_dsd.hpp"
#include "../include/algorithms/bi_decomposition.hpp"
#include "../include/algorithms/bi_dec_else_dec.hpp"   
#include "../include/algorithms/mix_dsd.hpp"
namespace alice
{

class bd_command : public command
{
public:
    explicit bd_command(const environment::ptr &env)
        : command(env, "Bi-decomposition (recursive)")
    {
        add_option("-f, --factor", hex_input,
                   "truth table as hex string")->required();

        // Only decompose k2=0 cases
        add_flag("-d, --k2_zero", only_k2_zero,
                 "only try bi-decomposition with k2=0,dsd,顶层为2输入");

        // ⭐ 加 else dec fallback 开关
        add_flag("-e, --else_dec", use_else_dec,
                 "enable else_dec fallback when BD fails");

        
        add_flag("--dm, --dsd_mix", use_dsd_mix,
                 "enable mixed DSD (-m) fallback when BD cannot proceed");
    }

protected:
    void execute() override
    {
        using clk = std::chrono::high_resolution_clock;
        
        use_else_dec = is_set("else_dec");
        use_dsd_mix = is_set("dsd_mix") || is_set("dm");
        only_k2_zero = is_set("k2_zero");
        std::string hex = hex_input;

        if (!is_set("factor"))
        {
            std::cout << "❌ Usage: bd -f <hex> [-s]\n";
            return;
        }

        // ------------------------------------------------------
        // 解析 Hex 变成真值表（binary）
        // ------------------------------------------------------
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

        std::cout << "📘 TT = " << binF << "  (vars=" << nvars << ")\n";

        // ------------------------------------------------------
        // 设置 其他分解模式（全局变量）
        // ------------------------------------------------------
        ENABLE_ELSE_DEC = use_else_dec;  // ⭐ 关键一步
        BD_ENABLE_DSD_MIX_FALLBACK = use_dsd_mix;
        BD_ONLY_K2_EQ_0 = only_k2_zero;

        auto t1 = clk::now();

        
        bool success = run_bi_decomp_recursive(binF);

        auto t2 = clk::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

        if (!success)
          {
            if (use_dsd_mix)
            {
                std::cout << "⚠️ BD 失败，启动 DSD -m 回退…\n";
                auto t_mix_start = clk::now();
                run_dsd_recursive_mix(binF);
                auto t_mix_end = clk::now();
                auto mix_us = std::chrono::duration_cast<std::chrono::microseconds>(t_mix_end - t_mix_start).count();
                std::cout << "⏱ DSD -m fallback time = " << mix_us << " us\n";
                return;
            }
            std::cout << "❌ Decomposition failed\n";
                 return;
        }

        std::cout << "⏱ time = " << us << " us\n";
    }

private:
    std::string hex_input{};
    bool use_else_dec = false;  // ⭐ 是否启用 其他分解
    bool only_k2_zero = false;
    bool use_dsd_mix = false;  
};

ALICE_ADD_COMMAND(bd, "STP")

} // namespace alice

#endif
