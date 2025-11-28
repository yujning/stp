#pragma once
#include <bits/stdc++.h>
using namespace std;

#include "excute.hpp"
#include "reorder.hpp"

// ------------------------------
// 额外引入 kitty（仅新增这两行）
// ------------------------------
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

// =====================================================
// DSD 节点结构
// =====================================================
struct DSDNode {
    int id;
    string func;           // 真值表 / MF 串（仅 0/1）
    vector<int> child;
};

static vector<DSDNode> NODE_LIST;
static int NODE_ID = 1;
static int STEP_ID = 1;

// =====================================================
// 01 串 <-> kitty truth table
// =====================================================
static kitty::dynamic_truth_table make_tt_from01(const string& f01)
{
    const size_t len = f01.size();
    assert((len & (len - 1)) == 0);
    unsigned n = static_cast<unsigned>(log2(len));

    kitty::dynamic_truth_table tt(n);
    for (uint64_t i = 0; i < len; i++)
        if (f01[i] == '1')
            kitty::set_bit(tt, i);

    return tt;
}

// 检查变量依赖（support）
static bool has_var_kitty(const string& f01, uint8_t var)
{
    auto tt = make_tt_from01(f01);
    return kitty::has_var(tt, var);
}

// 获取支持集变量下标
static vector<int> get_support_indices(const string& f01)
{
    auto tt = make_tt_from01(f01);
    vector<int> support;

    for (uint8_t i = 0; i < tt.num_vars(); i++)
        if (kitty::has_var(tt, i))
            support.push_back(i);

    return support;
}

// =====================================================
// ⭐ 最终正确版：按支持集压缩真值表（mockturtle 同款）
// =====================================================
static string shrink_to_support(const string& f01)
{
    auto tt = make_tt_from01(f01);

    // 1. 找支持集
    vector<int> support;
    for (uint8_t i = 0; i < tt.num_vars(); i++)
        if (kitty::has_var(tt, i))
            support.push_back(i);

    unsigned new_vars = support.size();
    if (new_vars == tt.num_vars())
        return f01; // 无需缩减

    // 2. 重新构造 truth table，仅保留 support 变量
    kitty::dynamic_truth_table new_tt(new_vars);

    for (uint64_t x = 0; x < (1ull << new_vars); x++)
    {
        uint64_t old_index = 0;

        for (unsigned b = 0; b < new_vars; b++)
        {
            uint64_t bit = (x >> b) & 1;
            old_index |= (bit << support[b]);  // 映射回旧表
        }

        if (kitty::get_bit(tt, old_index))
            kitty::set_bit(new_tt, x);
    }

    // 3. 转回 01 串
    string result;
    result.resize(1ull << new_vars);

    for (uint64_t i = 0; i < result.size(); i++)
        result[i] = kitty::get_bit(new_tt, i) ? '1' : '0';

    return result;
}


static inline string mul_ui(const string& ui, const string& w)
{
    if (ui == "10") return w;
    if (ui == "01") {
        string r; r.reserve(w.size());
        for (size_t i = 0; i + 1 < w.size(); i += 2) {
            r.push_back(w[i + 1]);
            r.push_back(w[i]);
        }
        return r;
    }
    if (ui == "11") return string(w.size(), '1');
    if (ui == "00") return string(w.size(), '0');
    return w;
}

// =====================================================
// 模板计算
// =====================================================
struct TemplateResult {
    string MF;    // 4 bits
    string Mphi;  // block 标签 01 串
    string Mpsi;  // 子函数 01 串
};

static TemplateResult run_case_once(
    const vector<string>& blocks01,
    int s,
    const string& S0,
    const string& S1)
{
    int m = blocks01.size();
    vector<string> W = blocks01;

    string MF = S0 + S1;
    string Mpsi;

    // 尝试反推
    auto try_u = [&](const string& u)->bool {
        if (u == "10" || u == "01") {
            vector<string> cand;
            for (int i = 0; i < m; i++) {
                if (is_constant_block(W[i])) continue;
                string c = mul_ui(u, W[i]);
                if (mul_ui(u, c) == W[i])
                    cand.push_back(c);
            }
            if (cand.empty()) return false;
            for (size_t k = 1; k < cand.size(); k++)
                if (cand[k] != cand[0]) return false;
            Mpsi = cand[0];
            return true;
        }
        return false;
    };

    if (!try_u(S0) && !try_u(S1))
    {
        int pick = -1;
        for (int i = 0; i < m; i++)
            if (!is_constant_block(W[i])) { pick = i; break; }
        if (pick < 0) pick = 0;
        Mpsi = W[pick];
    }

    // Φ
    string exp0 = mul_ui(S0, Mpsi);
    string exp1 = mul_ui(S1, Mpsi);

    string Mphi; Mphi.reserve(m);
    for (int i = 0; i < m; i++)
        Mphi.push_back(W[i] == exp0 ? '1' : '0');

    return { MF, Mphi, Mpsi };
}

// =====================================================
// 构造节点
// =====================================================
static int new_node(const string& func, const vector<int>& child)
{
    int id = NODE_ID++;
    NODE_LIST.push_back({ id, func, child });
    return id;
}

static int build_small_tree_from01(const string& f01)
{
    if (f01.size() == 2) {
        int a = new_node("in", {});
        return new_node(f01, { a });
    }
    if (f01.size() == 4) {
        int a = new_node("in", {});
        int b = new_node("in", {});
        return new_node(f01, { a, b });
    }
    return new_node(f01, {});
}

// =====================================================
// 单步分解（框架不动）
// =====================================================
static bool factor_once_with_reorder_01(
    const string& binary01,
    int depth,
    string& MF12,
    string& phi01,
    string& psi01)
{
    int len = binary01.size();
    if (!is_power_of_two(len) || len <= 4) return false;

    int n = log2(len);
    int r = n / 2;

    vector<stp_data> Mf = binary_to_vec(binary01);

    vector<int> s_order;
    if (r >= 2) {
        s_order.push_back(2);
        for (int s = 1; s <= r; s++)
            if (s != 2) s_order.push_back(s);
    }
    else {
        for (int s = 1; s <= r; s++)
            s_order.push_back(s);
    }

    for (int s : s_order)
    {
        vector<bool> v(n);
        fill(v.begin(), v.begin() + s, true);

        do {
            vector<int> Lambda;
            for (int i = 0; i < n; i++)
                if (v[i]) Lambda.push_back(i + 1);

            // swap-chain 计算
            vector<vector<stp_data>> chain;

            for (int k = s; k >= 1; k--) {
                int j_k = Lambda[k - 1];
                int exp = j_k + (s - 1) - k;
                chain.push_back(generate_swap_vec(2, 1 << exp));
            }
            chain.push_back(generate_swap_vec(1 << (n - s), 1 << s));

            auto Mperm = Vec_chain_multiply(chain, false);
            auto result = Vec_semi_tensor_product(Mf, Mperm);

            string reordered;
            reordered.reserve(len);
            for (size_t i = 1; i < result.size(); i++)
                reordered.push_back(result[i] ? '1' : '0');

            cout << "Λ = { ";
            for (int j : Lambda) cout << j << " ";
            cout << "} -> reordered = " << reordered << "\n";

            int cid = theorem33_case_id(reordered, s);
            if (cid == 0) continue;

            // 分块
            int bl = 1 << s;
            int nb = len / bl;

            vector<string> blocks(nb);
            for (int i = 0; i < nb; i++)
                blocks[i] = reordered.substr(i * bl, bl);

            // S0,S1
            bool has1 = false, has0 = false;
            for (auto& b : blocks) {
                if (is_constant_block(b)) {
                    if (b[0] == '1') has1 = true;
                    if (b[0] == '0') has0 = true;
                }
            }

            vector<pair<string,string>> S_list;
            switch (cid) {
                case 1: S_list = { {"11","00"}, {"00","11"} }; break;
                case 2:
                    if (has1)
                        S_list = { {"11","10"}, {"11","01"},
                                   {"10","11"}, {"01","11"} };
                    else
                        S_list = { {"00","10"}, {"00","01"},
                                   {"10","00"}, {"01","00"} };
                    break;
                case 3: S_list = { {"10","10"}, {"01","01"} }; break;
                case 4: S_list = { {"10","01"}, {"01","10"} }; break;
                case 5: return false;
            }

            auto R = run_case_once(blocks, s, S_list[0].first, S_list[0].second);

            cout << STEP_ID++ << ". MF  = [" << R.MF << "]\n";
            cout << "   MΦ  = [" << R.Mphi << "]\n";
            cout << "   Mψ  = [" << R.Mpsi << "]\n\n";

            MF12 = R.MF;
            phi01 = R.Mphi;
            psi01 = R.Mpsi;
            return true;

        } while (prev_permutation(v.begin(), v.end()));
    }

    return false;
}

// =====================================================
// 主递归 DSD
// =====================================================
static int dsd_factor(const string& f01_raw, int depth=0)
{
    // ⭐ 关键：首先缩减支持集（mockturtle 同款逻辑）
    string f01 = shrink_to_support(f01_raw);

    int len = f01.size();
    if (!is_power_of_two(len) || len <= 4)
        return build_small_tree_from01(f01);

    string MF12, phi01, psi01;
    if (!factor_once_with_reorder_01(f01, depth, MF12, phi01, psi01))
        return build_small_tree_from01(f01);

    int L = dsd_factor(phi01, depth + 1);
    int R = dsd_factor(psi01, depth + 1);
    return new_node(MF12, { L, R });
}

// =====================================================
// 顶层入口
// =====================================================
inline bool run_dsd_recursive(const string& binary01)
{
    if (!is_power_of_two(binary01.size())) {
        cout << "输入长度必须为 2^n\n";
        return false;
    }

    // 打印支持集（未缩减）
    auto support = get_support_indices(binary01);
    cout << "support vars = { ";
    for (auto v : support) cout << v << " ";
    cout << "}\n";

    NODE_LIST.clear();
    NODE_ID = 1;
    STEP_ID = 1;

    cout << "输入 = " << binary01
         << " (len = " << binary01.size() << ")\n";

    int root = dsd_factor(binary01);

    cout << "===== 最终 DSD 节点列表 =====\n";
    for (auto& nd : NODE_LIST) {
        cout << nd.id << " = " << nd.func;
        if (nd.child.size() == 1)
            cout << "(" << nd.child[0] << ")";
        else if (nd.child.size() == 2)
            cout << "(" << nd.child[0] << "," << nd.child[1] << ")";
        cout << "\n";
    }

    cout << "Root = " << root << "\n";
    return true;
}
