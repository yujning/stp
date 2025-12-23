#pragma once

#include "stp_dsd.hpp"
#include "strong_dsd.hpp"

// =====================================================
// Mixed DSD: prefer normal DSD, fallback to strong DSD
// =====================================================
inline int dsd_factor_mix(const TT& f, int depth = 0)
{
    int len = f.f01.size();
    if (len <= 4)
        return build_small_tree(f);

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

        int L = dsd_factor_mix(phi_tt, depth + 1);
        int R = dsd_factor_mix(psi_tt, depth + 1);

        return new_node(MF12, {L, R});
    }

    std::cout << "⚠️ DSD -f failed at depth " << depth
              << ", fallback to strong DSD.\n";
    return build_strong_dsd_nodes(f.f01, f.order, depth);
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
    int root_id = dsd_factor_mix(root_shrunk);

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