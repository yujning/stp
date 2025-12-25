#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cstdint>

#include "stp_dsd.hpp"
#include "node_global.hpp"

// =====================================================
// Debug switch
// =====================================================
inline bool STRONG_DSD_DEBUG_PRINT = true;

// =====================================================
// Pretty print: TT + order (must be paired)
// =====================================================
inline void print_tt_with_order(
    const std::string& title,
    const std::string& tt,
    const std::vector<int>& order,
    int depth = 0)
{
    if (!STRONG_DSD_DEBUG_PRINT) return;

    std::string indent((size_t)depth * 2, ' ');
    std::cout << indent << "📌 " << title << "\n";
    std::cout << indent << "   TT    = " << tt << "\n";
    std::cout << indent << "   order = { ";
    for (int v : order) std::cout << v << " ";
    std::cout << "}\n";
}


// =====================================================
// StrongDsdResult
// =====================================================
struct StrongDsdResult {
    bool found = false;
    size_t L = 0;          // = 2^{|Mx|}
    std::string Mx;        // = block0 + block1, length 2*L
    std::string My;        // length 2^{|My|}
};

// =====================================================
// Split meta
// =====================================================
struct StrongDsdSplit {
    bool found = false;
    StrongDsdResult dsd;

    // 0-based positions in current order
    std::vector<int> mx_pos;
    std::vector<int> my_pos;

    // var IDs in MSB->LSB order (by position)
    std::vector<int> mx_vars_msb2lsb;
    std::vector<int> my_vars_msb2lsb;

    // optional debug info
    std::string block0;
    std::string block1;

    std::string reordered_tt;
};

// =====================================================
// next_combination: comb is strictly increasing indices in [0..m-1]
// =====================================================
inline bool next_combination(std::vector<int>& comb, int m)
{
    int k = (int)comb.size();
    for (int i = k - 1; i >= 0; --i) {
        if (comb[i] < m - k + i) {
            ++comb[i];
            for (int j = i + 1; j < k; ++j) {
                comb[j] = comb[j - 1] + 1;
            }
            return true;
        }
    }
    return false;
}

inline bool is_all_same(const std::string& s, char c)
{
    for (char x : s) if (x != c) return false;
    return true;
}

// =====================================================
// Extract block (Mx subfunction) for a fixed My assignment
// - order positions: 0..n-1 (0 is MSB position)
// - mx_pos_sorted / my_pos_sorted must be sorted ascending by position (MSB->LSB)
// - my_assignment encoded in My's MSB->LSB order
// - mx_index encoded in Mx's MSB->LSB order
// =====================================================
inline std::string extract_block_for_mx(
    const std::string& mf,
    int n,
    const std::vector<int>& mx_pos_sorted,
    const std::vector<int>& my_pos_sorted,
    uint64_t my_assignment)
{
    int k = (int)mx_pos_sorted.size();
    int m = (int)my_pos_sorted.size();

    const size_t sub_size = 1ull << k;
    std::string sub(sub_size, '0');

    for (size_t mx_index = 0; mx_index < sub_size; ++mx_index)
    {
        uint64_t full_index = 0;

        // 1) place My bits (MSB->LSB)
        for (int t = 0; t < m; ++t)
        {
            int pos = my_pos_sorted[t];
            int bit = (int)((my_assignment >> (m - 1 - t)) & 1ull);
            full_index |= (uint64_t(bit) << (n - 1 - pos));
        }

        // 2) place Mx bits (MSB->LSB)
        for (int t = 0; t < k; ++t)
        {
            int pos = mx_pos_sorted[t];
            int bit = (int)((mx_index >> (k - 1 - t)) & 1ull);
            full_index |= (uint64_t(bit) << (n - 1 - pos));
        }

        sub[mx_index] = mf[(size_t)full_index];
    }

    return sub;
}

// =====================================================
// Helper: print current candidate split info
// =====================================================
inline void print_candidate_info(
    int depth,
    int k,
    const std::vector<int>& mx_vars_msb2lsb,
    const std::vector<int>& my_vars_msb2lsb)
{
    if (!STRONG_DSD_DEBUG_PRINT) return;

    std::string indent((size_t)depth * 2, ' ');
    std::cout << indent << "🔎 尝试 |Mx|=" << k << " : Mx={ ";
    for (int v : mx_vars_msb2lsb) std::cout << v << " ";
    std::cout << "}  My={ ";
    for (int v : my_vars_msb2lsb) std::cout << v << " ";
    std::cout << "}\n";
}

// =====================================================
// ★ Subset-enumeration Strong DSD
// Search order: |Mx| = ceil(n/2) down to 1
// Combination order: by var_id ascending (a,b,c...)
// =====================================================
inline StrongDsdSplit run_strong_dsd_by_mx_subset(
    const std::string& mf,
    const std::vector<int>& order,
    int depth_for_print = 0)
{
    StrongDsdSplit out;

    int n = (int)order.size();
    if (n <= 1) return out;

    const size_t total = mf.size();
    if ((size_t(1) << n) != total) return out;

    // (var_id, pos) sorted by var_id (ascending)
    struct VarPos { int var; int pos; };
    std::vector<VarPos> vp;
    vp.reserve(n);
    for (int pos = 0; pos < n; ++pos) vp.push_back({order[pos], pos});
    std::sort(vp.begin(), vp.end(), [](const VarPos& x, const VarPos& y){
        return x.var < y.var;
    });

    //int max_k = (n + 1) / 2;


    int min_k = 1;          // |Mx| >= 1
    int max_k = n - 2;      // |My| >= 2

    for (int k = max_k; k >= min_k; --k)   

    {
        std::vector<int> comb(k);
        for (int i = 0; i < k; ++i) comb[i] = i;

        while (true)
        {
            // mx positions from comb
            std::vector<int> mx_pos;
            mx_pos.reserve(k);
            for (int idx : comb) mx_pos.push_back(vp[idx].pos);

            // my positions = rest
            std::vector<int> my_pos;
            my_pos.reserve(n - k);
            {
                std::vector<char> is_mx(n, 0);
                for (int p : mx_pos) is_mx[p] = 1;
                for (int p = 0; p < n; ++p) if (!is_mx[p]) my_pos.push_back(p);
            }

            // sort by position (MSB->LSB)
            std::sort(mx_pos.begin(), mx_pos.end());
            std::sort(my_pos.begin(), my_pos.end());

            // build vars list in MSB->LSB by position
            std::vector<int> mx_vars_msb2lsb;
            std::vector<int> my_vars_msb2lsb;
            mx_vars_msb2lsb.reserve(mx_pos.size());
            my_vars_msb2lsb.reserve(my_pos.size());
            for (int p : mx_pos) mx_vars_msb2lsb.push_back(order[p]);
            for (int p : my_pos) my_vars_msb2lsb.push_back(order[p]);

            // print candidate split (variables only)
            print_candidate_info(depth_for_print, k, mx_vars_msb2lsb, my_vars_msb2lsb);

            

            int m = n - k;
            uint64_t my_count = 1ull << m;
            uint64_t L = 1ull << k;

            if (STRONG_DSD_DEBUG_PRINT)
{
    std::string reordered;
    reordered.reserve((size_t)my_count * (size_t)L);

    for (uint64_t y = 0; y < my_count; ++y)
    {
        std::string block = extract_block_for_mx(
            mf, n, mx_pos, my_pos, y
        );
        reordered += block;
    }

    std::vector<int> reordered_order;
    reordered_order.reserve(n);
    reordered_order.insert(
        reordered_order.end(),
        my_vars_msb2lsb.begin(),
        my_vars_msb2lsb.end()
    );
    reordered_order.insert(
        reordered_order.end(),
        mx_vars_msb2lsb.begin(),
        mx_vars_msb2lsb.end()
    );

    print_tt_with_order(
        "候选 split 的重排 TT (My|Mx)",
        reordered,
        reordered_order,
        depth_for_print
    );
}

            std::unordered_map<std::string, int> block_index;
            std::vector<std::string> blocks;
            blocks.reserve(2);

            std::string My;
            My.reserve((size_t)my_count);

            bool too_many = false;

            for (uint64_t y = 0; y < my_count; ++y)
            {
                std::string block = extract_block_for_mx(mf, n, mx_pos, my_pos, y);

                auto it = block_index.find(block);
                if (it == block_index.end())
                {
                    if (blocks.size() >= 2) { too_many = true; break; }
                    int id = (int)blocks.size();
                    block_index.emplace(block, id);
                    blocks.push_back(block);
                    My.push_back(id == 0 ? '1' : '0');
                }
                else
                {
                    My.push_back(it->second == 0 ? '1' : '0');
                }
            }

            // non-degenerate: exactly 2 blocks
            if (!too_many && blocks.size() == 2)
            {

                //  // ★ 新增：拒绝 |My| == 1 的 split
                // if (my_vars_msb2lsb.size() == 1) {
                //     // 不接受这个 split，继续枚举
                //     goto NEXT_COMBINATION;
                // }
                out.found = true;
                out.dsd.found = true;
                out.dsd.L = (size_t)L;
                out.dsd.Mx = blocks[0] + blocks[1];
                out.dsd.My = My;

                out.mx_pos = mx_pos;
                out.my_pos = my_pos;
                out.mx_vars_msb2lsb = mx_vars_msb2lsb;
                out.my_vars_msb2lsb = my_vars_msb2lsb;
                

                std::string reordered;
                reordered.reserve((size_t)my_count * (size_t)L);

                for (uint64_t y = 0; y < my_count; ++y)
                {
                    std::string block = extract_block_for_mx(
                        mf, n, mx_pos, my_pos, y
                    );
                    reordered += block;
                }

                out.reordered_tt = reordered;

                out.block0 = blocks[0];
                out.block1 = blocks[1];

                if (STRONG_DSD_DEBUG_PRINT)
                {
                    std::string indent((size_t)depth_for_print * 2, ' ');
                    std::cout << indent << "✅ 命中 Strong DSD split\n";
                    std::cout << indent << "   block0 = " << blocks[0] << "\n";
                    std::cout << indent << "   block1 = " << blocks[1] << "\n";
                    // 打印 My 真值表 + order
                    print_tt_with_order("当前 split 的 My", My, my_vars_msb2lsb, depth_for_print);
                }

                return out;
            }

            NEXT_COMBINATION:
            if (!next_combination(comb, (int)vp.size())) break;
        }
    }

    return out;
}

// =====================================================
// children
// =====================================================
inline std::vector<int>
make_children_from_order_with_placeholder(
    const std::vector<int>& order,
    const std::unordered_map<int, int>* placeholder_nodes,
    const std::vector<int>* local_to_global)
{
    std::vector<int> children;
    children.reserve(order.size());

    // 关键：按 MSB -> LSB 直接放
    for (int var_id : order) {

        if (placeholder_nodes) {
            auto ph = placeholder_nodes->find(var_id);
            if (ph != placeholder_nodes->end()) {
                children.push_back(ph->second);
                continue;
            }
        }

        int global_var_id = var_id;
        if (local_to_global &&
            var_id >= 0 &&
            var_id < (int)local_to_global->size())
        {
            int mapped = (*local_to_global)[var_id];
            if (mapped != 0) global_var_id = mapped;
        }

        children.push_back(new_in_node(global_var_id));
    }

    return children;
}

inline int resolve_global_var_id(
    int local_id,
    const std::vector<int>* local_to_global)
{
    if (local_to_global && local_id >= 0 && local_id < (int)local_to_global->size()) {
        int mapped = (*local_to_global)[local_id];
        if (mapped != 0) return mapped;
    }
    return local_id;
}

// =====================================================
// ★ Recursive Strong DSD (subset enumeration)
// =====================================================
inline int build_strong_dsd_nodes_impl(
    const std::string& mf,
    const std::vector<int>& order,
    int depth,
    const std::vector<int>* local_to_global,
    const std::unordered_map<int, int>* placeholder_nodes)
{
    // 入口：打印当前 TT + order
    print_tt_with_order("进入 Strong DSD", mf, order, depth);

    if (mf.size() <= 4) {
        print_tt_with_order("⏹ Stop (size <= 4)", mf, order, depth);
               auto children = make_children_from_order_with_placeholder(
            order, placeholder_nodes, local_to_global);
        return new_node(mf, children);
    }

    // ① subset-enum split, pass depth for aligned prints
    StrongDsdSplit split = run_strong_dsd_by_mx_subset(mf, order, depth);

    // ===== 关键：拒绝 |My| == 1 的 split（必须在使用 split 之前）=====
// if (split.found && split.my_vars_msb2lsb.size() == 1) {
//     std::string indent((size_t)depth * 2, ' ');
//     std::cout << indent
//               << "⚠️ Skip split: |My| == 1 (not accepted)\n";

//     auto children = make_children_from_order_with_placeholder(
//         order, placeholder_nodes, local_to_global
//     );
//     return new_node(mf, children);
// }


    if (!split.found) {
        std::string indent((size_t)depth * 2, ' ');
        std::cout << indent << "❌ Strong DSD: no valid split\n";
        auto children = make_children_from_order_with_placeholder(
        order, placeholder_nodes, local_to_global);
        return new_node(mf, children);
    }

    const auto& result = split.dsd;

    {
        std::string indent((size_t)depth * 2, ' ');
        std::cout << indent << "✅ L = " << result.L << "\n";
        std::cout << indent << "Mx = " << result.Mx << "\n";
        std::cout << indent << "My = " << result.My << "\n";

        std::cout << indent << "My 使用变量（MSB->LSB）：{ ";
        for (int v : split.my_vars_msb2lsb) std::cout << v << " ";
        std::cout << "}\n";
        std::cout << indent << "Mx 使用变量（MSB->LSB）：{ ";
        for (int v : split.mx_vars_msb2lsb) std::cout << v << " ";
        std::cout << "}\n";
    }

    // ② recurse on My
    // IMPORTANT: your existing children construction reverses order,
    // so keep the same convention: reverse before passing down.
// ===== ① recurse on My (MSB->LSB, no reverse) =====
const std::vector<int>& order_my = split.my_vars_msb2lsb;

print_tt_with_order("递归进入 My", result.My, order_my, depth);

int my_id = build_strong_dsd_nodes_impl(
    result.My,
    order_my,
    depth + 1,
    local_to_global,
    placeholder_nodes
);

// ===== ② recurse on Mx =====
// |Mx| = k
int k = (int)split.mx_vars_msb2lsb.size();
int my_local_id = k + 1;

// order: (my_local, k, k-1, ..., 1)
std::vector<int> order_mx;
order_mx.reserve(k + 1);
order_mx.push_back(my_local_id);
for (int i = k; i >= 1; --i) order_mx.push_back(i);

// placeholder + local_to_global
std::unordered_map<int, int> placeholder_nodes_mx;
std::vector<int> local_to_global_mx(my_local_id + 1, 0);

// My 占据 local_id = k+1
placeholder_nodes_mx[my_local_id] = my_id;

// 继承父 placeholder
const std::unordered_map<int, int> empty_ph;
const auto& parent_ph = placeholder_nodes ? *placeholder_nodes : empty_ph;

// Mx 变量映射：MSB->LSB
for (int i = 0; i < k; ++i) {
    int old_id = split.mx_vars_msb2lsb[i];
    int new_id = k - i;

    auto it = parent_ph.find(old_id);
    if (it != parent_ph.end()) {
        placeholder_nodes_mx[new_id] = it->second;
    } else {
        local_to_global_mx[new_id] =
            resolve_global_var_id(old_id, local_to_global);
    }
}

print_tt_with_order("递归进入 Mx", result.Mx, order_mx, depth);

return build_strong_dsd_nodes_impl(
    result.Mx,
    order_mx,
    depth + 1,
    &local_to_global_mx,
    &placeholder_nodes_mx
);

}

inline int build_strong_dsd_nodes(
    const std::string& mf,
    const std::vector<int>& order,
    int depth = 0)
{
    return build_strong_dsd_nodes_impl(mf, order, depth, nullptr, nullptr);
}
