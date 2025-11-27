#pragma once
#include <bits/stdc++.h>
using namespace std;

#include "excute.hpp"
#include "reorder.hpp"

struct DSDNode {
    int id;                 // 节点编号
    string func;            // 真值表 / MF 串（只含 '0' / '1'）
    vector<int> child;      // 子节点编号
};

static vector<DSDNode> NODE_LIST;
static int NODE_ID = 1;
static int STEP_ID = 1;     // 分解步骤编号：1,2,3,...


static inline string mul_ui(const string& ui, const string& w)
{
    if (ui == "10") return w;              // 恒等 (identity)

    if (ui == "01") {                      // swap：每 2bit 交换
        string r;
        r.reserve(w.size());
        for (size_t i = 0; i + 1 < w.size(); i += 2) {
            r.push_back(w[i + 1]);
            r.push_back(w[i]);
        }
        return r;
    }

    if (ui == "11") {                      // 常量 1
        return string(w.size(), '1');
    }

    if (ui == "00") {                      // 常量 0（原来常量 2）
        return string(w.size(), '0');
    }

    // 理论上不会走到这里，保险起见当恒等
    return w;
}

static void indent(int d)
{
    while (d--) cout << "   ";
}

static string case_note(int cid)
{
    switch (cid) {
        case 1: return "两个常量块（如全0与全1）";
        case 2: return "一种常量块 + 一种非常量块";
        case 3: return "只有一种非常量块";
        case 4: return "两种非常量块且互补";
        case 5: return "只有一种常量块（函数恒定）";
        default: return "未知";
    }
}

// =============================
// 通用模板分解：给定 blocks01 + S0,S1
// 说明：
//   - blocks01：每个块是 0/1 串
//   - S0, S1 ∈ {"00","01","10","11"}
//   - MF = S0 + S1（4 位 0/1 串）
// =============================
struct TemplateResult {
    string MF;    // 长度 4，仅含 '0' / '1'
    string Mphi;  // 标签串：每位 '0' 或 '1'
    string Mpsi;  // 子函数真值表：0/1 串
};

static TemplateResult run_case_once(
    const vector<string>& blocks01,
    int s,
    const string& S0,
    const string& S1)
{
    // 不再做 0→2 映射，直接在 01 域工作
    vector<string> W = blocks01;

    int m = (int)W.size();
    string MF = S0 + S1;   // 4-bit 0/1 串

    // ----------------- 求 Mψ -----------------
    string Mpsi;
    bool found = false;

    auto try_u = [&](const string& u)->bool {
        // 只对“恒等 / swap”尝试反推
        if (u == "10" || u == "01") {
            vector<string> cand;

            for (int i = 0; i < m; i++) {
                // 跳过常量块（全 0 或全 1）
                if (W[i] == string(W[i].size(), W[i][0]))
                    continue;

                string c = mul_ui(u, W[i]);

                // 检查是否 c 反推回去能得到 W[i]
                if (mul_ui(u, c) == W[i])
                    cand.push_back(c);
            }

            if (cand.empty()) return false;

            // 多个非常量块反推结果必须一致
            for (size_t k = 1; k < cand.size(); k++)
                if (cand[k] != cand[0])
                    return false;

            Mpsi = cand[0];    // 0/1 串
            return true;
        }
        return false;
    };

    if (!found && try_u(S0)) found = true;
    if (!found && try_u(S1)) found = true;

    if (!found) {
        int pick = -1;
        for (int i = 0; i < m; i++) {
            if (!is_constant_block(blocks01[i])) {
                pick = i; break;
            }
        }
        if (pick < 0) pick = 0;
        Mpsi = W[pick];
    }

    // ----------------- 求 MΦ 标签（0/1 串） -----------------
    string exp0 = mul_ui(S0, Mpsi);
    string exp1 = mul_ui(S1, Mpsi);
    string Mphi; Mphi.reserve(m);

    for (int i = 0; i < m; i++) {
        // 这里沿用原逻辑：匹配 S0 的记为 '1'，其余记为 '0'
        if (W[i] == exp0) Mphi.push_back('1');
        else              Mphi.push_back('0');
    }

    return { MF, Mphi, Mpsi };
}

// =============================
// 创建节点 + 小 LUT 叶子结构
// =============================
static int new_node(const string& func, const vector<int>& child)
{
    int id = NODE_ID++;
    NODE_LIST.push_back({ id, func, child });
    return id;
}

// 长度 ≤4：构造 2-LUT 叶子结构（直接使用 0/1 真值表）
static int build_small_tree_from01(const string& f01)
{
    int len = (int)f01.size();

    if (len == 2) {
        int in = new_node("in", {});
        return new_node(f01, { in });
    }
    else if (len == 4) {
        int inL = new_node("in", {});
        int inR = new_node("in", {});
        return new_node(f01, { inL, inR });
    }
    // 理论上不会进来
    return new_node(f01, {});
}

// =============================
// 单步：对一个 01 串做 “重排 + case + 模板”
// 输入：binary01（当前要分解的函数，01 域）
// 输出：MF12（4位 0/1 串），phi01, psi01（下一轮 01 域真值表）
// 返回：是否成功分解（false = 无法按定理3.3分解）
// =============================
static bool factor_once_with_reorder_01(
    const string& binary01,
    int depth,
    string& MF12,
    string& phi01,
    string& psi01)
{
    int len = (int)binary01.size();
    if (!is_power_of_two(len) || len <= 4) return false;

    int n = (int)log2(len);
    int r = n / 2;

    vector<stp_data> Mf = binary_to_vec(binary01);

    // --------- s 的优先级：n>=4 时优先考虑 s=2（2-LUT 风格） ---------
    vector<int> s_order;
    if (r >= 2) {
        s_order.push_back(2);           // 先 s = 2
        for (int s = 1; s <= r; ++s)
            if (s != 2) s_order.push_back(s); // 再 1,3,4,5,...
    }
    else {
        for (int s = 1; s <= r; ++s)
            s_order.push_back(s);       // 例如 n=3 时，只会有 s=1
    }

    // --------- 按 s_order 遍历 ---------
    for (int idx_s = 0; idx_s < (int)s_order.size(); ++idx_s)
    {
        int s = s_order[idx_s];

        vector<bool> v(n);
        fill(v.begin(), v.begin() + s, true);

        do {
            // Λ
            vector<int> Lambda;
            for (int i = 0; i < n; i++)
                if (v[i]) Lambda.push_back(i + 1);

            // swap_chain
            vector<vector<stp_data>> swap_chain;
            for (int k = s; k >= 1; k--) {
                int j_k = Lambda[k - 1];
                int exp = j_k + (s - 1) - k;
                swap_chain.push_back(generate_swap_vec(2, pow(2, exp)));
            }
            swap_chain.push_back(generate_swap_vec(pow(2, n - s), pow(2, s)));

            vector<stp_data> Mperm =
                Vec_chain_multiply(swap_chain, false);

            vector<stp_data> result =
                Vec_semi_tensor_product(Mf, Mperm);

            string reordered;
            reordered.reserve(len);
            for (size_t i = 1; i < result.size(); ++i)
                reordered.push_back(result[i] ? '1' : '0');

            // ---- 打印所有重排 ----
            cout << "Λ = { ";
            for (int j : Lambda) cout << j << " ";
            cout << "}  => reordered = " << reordered << "\n";

            int cid = theorem33_case_id(reordered, s);
            if (cid != 0) {
                // ---------- 分块 ----------
                int bl = 1 << s;
                int nb = len / bl;
                vector<string> blocks01(nb);
                for (int i = 0; i < nb; i++)
                    blocks01[i] = reordered.substr(i * bl, bl);

                // ---------- 按 case 选 S0,S1 ----------
                bool has_const1 = false, has_const0 = false;
                for (auto& b : blocks01) {
                    if (is_constant_block(b)) {
                        if (b[0] == '1') has_const1 = true;
                        if (b[0] == '0') has_const0 = true;
                    }
                }

                vector<pair<string, string>> S_list;
                switch (cid) {
                    case 1:
                        // 两个常量块（全0 & 全1）
                        // 原来 {("11","22"),("22","11")}
                        // 现在 const1="11", const0="00"
                        S_list = { {"11","00"}, {"00","11"} };
                        break;

                    case 2:
                        if (has_const1) {
                            S_list = {
                                {"11","10"}, {"11","01"},
                                {"10","11"}, {"01","11"}
                            };
                        } else {
                            S_list = {
                                {"00","10"}, {"00","01"},
                                {"10","00"}, {"01","00"}
                            };
                        }
                        break;

                    case 3:
                        S_list = { {"10","10"}, {"01","01"} };
                        break;

                    case 4:
                        S_list = { {"10","01"}, {"01","10"} };
                        break;

                    case 5:
                        return false;   // 恒定函数，不再分解
                }

                string S0 = S_list[0].first;
                string S1 = S_list[0].second;

                auto R = run_case_once(blocks01, s, S0, S1);

                // ---------- 打印这一步分解 ----------
                cout << STEP_ID++ << ". MF  = [" << R.MF  << "]\n";
                cout << "   MΦ  = [" << R.Mphi << "]\n";
                cout << "   Mψ  = [" << R.Mpsi << "]\n\n";

                MF12 = R.MF;        // 4-bit 0/1 串
                phi01 = R.Mphi;     // 0/1 标签串
                psi01 = R.Mpsi;     // 0/1 真值表
                return true;
            }

        } while (prev_permutation(v.begin(), v.end()));
    }

    return false;
}

// =============================
// 主递归：对 01 域真值表做 DSD
// =============================
static int dsd_factor(const string& f01, int depth=0)
{
    int len = (int)f01.size();

    // 小函数：直接构造 2-LUT 结构
    if (len <= 4 || !is_power_of_two(len)) {
        return build_small_tree_from01(f01);
    }

    string MF12, phi01, psi01;
    if (!factor_once_with_reorder_01(f01, depth, MF12, phi01, psi01)) {
        // 不能按定理 3.3 分解，当成 2-LUT
        return build_small_tree_from01(f01);
    }

    int left_id  = dsd_factor(phi01, depth + 1);
    int right_id = dsd_factor(psi01, depth + 1);

    return new_node(MF12, { left_id, right_id });
}

// =============================
// 顶层入口：递归 DSD
// =============================
inline bool run_dsd_recursive(const string& binary01)
{
    if (!is_power_of_two(binary01.size())) {
        cout << "输入长度必须是 2 的整数次幂\n";
        return false;
    }

    NODE_LIST.clear();
    NODE_ID = 1;
    STEP_ID = 1;

    cout << "输入 = " << binary01
         << " (len = " << binary01.size() << ")\n";

    int root = dsd_factor(binary01, 0);

    cout << "===== 最终 DSD 节点列表 =====\n";
    for (auto& nd : NODE_LIST) {
        cout << nd.id << " = " << nd.func;
        if (nd.child.size() == 1) {
            cout << "(" << nd.child[0] << ")";
        }
        else if (nd.child.size() == 2) {
            cout << "(" << nd.child[0] << "," << nd.child[1] << ")";
        }
        cout << "\n";
    }
    cout << "Root = " << root << "\n";
    return true;
}
