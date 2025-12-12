#pragma once
// ⭐ 防止 std::is_trivial 污染 percy
#ifdef is_trivial
#undef is_trivial
#endif

#include <iostream>
#include <vector>
#include <cmath>
#include <utility>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>

#include "stp_dsd.hpp"
#include "node_global.hpp"   // new_node / new_in_node


int new_node(const std::string&, const std::vector<int>&);


inline int else_decompose(
    const TT& f,
    const std::vector<int>& orig_children,
    int depth
)
{
   const int n = static_cast<int>(std::log2(f.f01.size()));

    if (n > 4)
    {
        std::cout << "⚠️ depth " << depth
                  << ": exact 2-LUT only supports n ≤ 4\n";
        return build_small_tree(f);
    }

    std::cout << "⚠️ depth " << depth
              << ": EXACT 2-LUT refine (n=" << n << ")\n";

    // =====================================================
    // 1) TT → kitty TT（严格保持变量顺序）
    // =====================================================
    kitty::dynamic_truth_table tt(n);
    const int N = 1 << n;

    for (int i = 0; i < N; ++i)
    {
        // STP → kitty 索引反序（你系统里是固定对的）
        int ki = N - 1 - i;
        if (f.f01[i] == '1')
            kitty::set_bit(tt, ki);
    }

    // =====================================================
    // 2) exact 2-LUT synthesis
    // =====================================================
    mockturtle::klut_network klut;
    std::vector<mockturtle::klut_network::signal> pis;
    pis.reserve(n);

    for (int i = 0; i < n; ++i)
        pis.push_back(klut.create_pi());

    mockturtle::exact_resynthesis<mockturtle::klut_network> resyn(2);
    resyn(klut, tt, pis.begin(), pis.end(),
          [&](auto const& g) { klut.create_po(g); });

    // =====================================================
    // 3) klut → NODE_LIST（关键：PI 不 new）
    // =====================================================
    std::unordered_map<mockturtle::klut_network::node, int> id_map;

    // ⭐ PI 映射：优先复用传入的 child，如果编号不匹配则按当前变量→原始变量映射新建
    {
        int idx = 0;
        klut.foreach_pi([&](auto pi_node) {
            const int var = f.order[idx];
            const int mapped_var = (ORIGINAL_VAR_COUNT > 0)
                                        ? ORIGINAL_VAR_COUNT - var + 1
                                        : var;

            int node_id = (idx < static_cast<int>(orig_children.size()))
                               ? orig_children[idx]
                               : -1;


            if (node_id <= 0)
                node_id = new_in_node(mapped_var);

            id_map[pi_node] = node_id;

            std::cout << "      PI" << idx
                      << " : local var " << var
                      << " (mapped=" << mapped_var << ")"
                      << " → node " << node_id << "\n";
            ++idx;
        });
    }

    // =====================================================
    // 4) gate 展开
    // =====================================================
    klut.foreach_gate([&](auto n_node) {

        auto ltt = klut.node_function(n_node);

        std::ostringstream oss;
        kitty::print_binary(ltt, oss);
        std::string func01 = oss.str();

        std::vector<int> children(ltt.num_vars());

        klut.foreach_fanin(n_node, [&](auto fin_sig, auto i) {
            auto fin_node = klut.get_node(fin_sig);
            children[i] = id_map.at(fin_node);
        });

                // exact 2-LUT 的左右输入与 STP 约定相反，统一在此交换
        if (children.size() == 2)
            std::swap(children[0], children[1]);


        int my_id = new_node(func01, children);
        id_map[n_node] = my_id;
    });

    // =====================================================
    // 5) 输出
    // =====================================================
    auto po_node = klut.get_node(klut.po_at(0));
    return id_map.at(po_node);
}
