#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "stp_dsd.hpp"
#include "strong_dsd.hpp"

// =====================================================
// Mixed DSD: prefer normal DSD, fallback to strong DSD per layer
// =====================================================
inline void add_final_var_order(int var_id)
{
    if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), var_id) ==
        FINAL_VAR_ORDER.end())
    {
        FINAL_VAR_ORDER.push_back(var_id);
    }
}

inline int resolve_leaf_node(
    int var_id,
    const std::unordered_map<int, int>* placeholder_nodes,
    const std::vector<int>* local_to_global)
{
    if (placeholder_nodes) {
        auto ph = placeholder_nodes->find(var_id);
        if (ph != placeholder_nodes->end()) {
            return ph->second;
        }
    }

    int global_var_id = resolve_global_var_id(var_id, local_to_global);
    return new_in_node(global_var_id);
}

inline void record_final_var_order(
    int var_id,
    const std::unordered_map<int, int>* placeholder_nodes,
    const std::vector<int>* local_to_global)
{
    if (placeholder_nodes) {
        auto ph = placeholder_nodes->find(var_id);
        if (ph != placeholder_nodes->end()) {
            return;
        }
    }

    int global_var_id = resolve_global_var_id(var_id, local_to_global);
    add_final_var_order(global_var_id);
}

inline int build_small_tree_mix(
    const TT& t,
    const std::vector<int>* local_to_global,
    const std::unordered_map<int, int>* placeholder_nodes)
{
    int nv = t.order.size();

    if (nv == 1)
    {
        int var_id = t.order[0];
        int leaf = resolve_leaf_node(var_id, placeholder_nodes, local_to_global);
        record_final_var_order(var_id, placeholder_nodes, local_to_global);

        if (t.f01 == "10") return leaf;                  // identity
        if (t.f01 == "01") return new_node("01", {leaf}); // NOT
        if (t.f01 == "00") return new_node("0", {});      // const 0
        if (t.f01 == "11") return new_node("1", {});      // const 1
        return leaf;
    }

    if (nv == 2)
    {
        for (int var_id : t.order)
            record_final_var_order(var_id, placeholder_nodes, local_to_global);

        int a = resolve_leaf_node(t.order[0], placeholder_nodes, local_to_global);
        int b = resolve_leaf_node(t.order[1], placeholder_nodes, local_to_global);
        return new_node(t.f01, {a, b});
    }

    std::vector<int> child_ids;
    child_ids.reserve(nv);
    for (int var_id : t.order)
    {
        record_final_var_order(var_id, placeholder_nodes, local_to_global);
        child_ids.push_back(resolve_leaf_node(var_id, placeholder_nodes, local_to_global));
    }

    return new_node(t.f01, child_ids);
}

inline int dsd_factor_mix_impl(
    const TT& f,
    int depth,
    const std::vector<int>* local_to_global,
    const std::unordered_map<int, int>* placeholder_nodes)
{
    int len = f.f01.size();
    if (len <= 4)
        return build_small_tree_mix(f, local_to_global, placeholder_nodes);

    std::string MF12;
    TT phi_tt, psi_tt;

    if (factor_once_with_reorder_01(f, depth, MF12, phi_tt, psi_tt)) {
        std::vector<int> phi_original_vars = phi_tt.order;
        std::vector<int> psi_original_vars = psi_tt.order;

        int n_phi = phi_tt.order.size();
        int n_psi = psi_tt.order.size();

        std::cout << "📌 递归分解 Φ：原始变量 { ";
        for (int v : phi_original_vars) std::cout << v << " ";
        std::cout << "} → 局部编号 { ";
        for (int i = 1; i <= n_phi; i++) std::cout << i << " ";
        std::cout << "}\n";
        std::cout << "   映射关系：";
        for (int i = 0; i < n_phi; i++)
            std::cout << "位置" << (i + 1) << "→变量" << phi_original_vars[i] << " ";
        std::cout << "\n";

        std::cout << "📌 递归分解 Ψ：原始变量 { ";
        for (int v : psi_original_vars) std::cout << v << " ";
        std::cout << "} → 局部编号 { ";
        for (int i = 1; i <= n_psi; i++) std::cout << i << " ";
        std::cout << "}\n";
        std::cout << "   映射关系：";
        for (int i = 0; i < n_psi; i++)
            std::cout << "位置" << (i + 1) << "→变量" << psi_original_vars[i] << " ";
        std::cout << "\n\n";

        int L = dsd_factor_mix_impl(phi_tt, depth + 1, local_to_global, placeholder_nodes);
        int R = dsd_factor_mix_impl(psi_tt, depth + 1, local_to_global, placeholder_nodes);

        return new_node(MF12, {L, R});
    }

    std::cout << "⚠️ DSD -f failed at depth " << depth
              << ", fallback to strong DSD (one layer).\n";

    StrongDsdSplit split = run_strong_dsd_by_mx_subset(f.f01, f.order, depth);
    if (!split.found) {
        return build_small_tree_mix(f, local_to_global, placeholder_nodes);
    }

    const auto& result = split.dsd;

    TT my_tt;
    my_tt.f01 = result.My;
    my_tt.order = split.my_vars_msb2lsb;

    int my_id = dsd_factor_mix_impl(my_tt, depth + 1, local_to_global, placeholder_nodes);

    int k = (int)split.mx_vars_msb2lsb.size();
    int my_local_id = k + 1;

    std::vector<int> order_mx;
    order_mx.reserve(k + 1);
    order_mx.push_back(my_local_id);
    for (int i = k; i >= 1; --i) order_mx.push_back(i);

    std::unordered_map<int, int> placeholder_nodes_mx;
    std::vector<int> local_to_global_mx(my_local_id + 1, 0);

    placeholder_nodes_mx[my_local_id] = my_id;

    const std::unordered_map<int, int> empty_ph;
    const auto& parent_ph = placeholder_nodes ? *placeholder_nodes : empty_ph;

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

    TT mx_tt;
    mx_tt.f01 = result.Mx;
    mx_tt.order = order_mx;

    return dsd_factor_mix_impl(
        mx_tt,
        depth + 1,
        &local_to_global_mx,
        &placeholder_nodes_mx);
}

inline bool run_dsd_recursive_mix(const std::string& binary01)
{
    RESET_NODE_GLOBAL();
    if (!is_power_of_two(binary01.size())) {
        std::cout << "输入长度必须为 2^n\n";
        return false;
    }

    int n = static_cast<int>(std::log2(binary01.size()));
    ORIGINAL_VAR_COUNT = n;

    TT root;
    root.f01 = binary01;
    root.order.resize(n);

    for (int i = 0; i < n; ++i)
        root.order[i] = n - i;

    std::cout << "输入 = " << binary01 << " (n=" << n << ")\n";
    std::cout << "初始映射：";
    for (int i = 0; i < n; i++)
        std::cout << "位置" << (i + 1) << "→变量" << root.order[i] << " ";
    std::cout << "\n\n";

    NODE_LIST.clear();
    NODE_ID = 1;
    STEP_ID = 1;
    FINAL_VAR_ORDER.clear();

    for (int v = 1; v <= n; ++v)
        new_in_node(v);

    TT root_shrunk = shrink_to_support(root);
    int root_id = dsd_factor_mix_impl(root_shrunk, 0, nullptr, nullptr);

    std::cout << "===== 最终 DSD 节点列表 =====\n";
    for (auto& nd : NODE_LIST)
    {
        std::cout << nd.id << " = " << nd.func;

        if (nd.func == "in")
        {
            std::cout << "(var=" << nd.var_id << ")";
        }
        else if (!nd.child.empty())
        {
            std::cout << "(";
            for (size_t i = 0; i < nd.child.size(); ++i)
            {
                std::cout << nd.child[i];
                if (i + 1 < nd.child.size())
                    std::cout << ",";
            }
            std::cout << ")";
        }

        std::cout << "\n";
    }

    std::cout << "Root = " << root_id << "\n";

    std::cout << "FINAL_VAR_ORDER = { ";
    for (int v : FINAL_VAR_ORDER) std::cout << v << " ";
    std::cout << "}\n";

    return true;
}