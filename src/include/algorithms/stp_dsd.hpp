#pragma once
#include <bits/stdc++.h>
using namespace std;

// bench 输出用：最终输入变量顺序（0=a,1=b,...）
inline std::vector<int> FINAL_VAR_ORDER;

#include "excute.hpp"
#include "reorder.hpp"

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

// =====================================================
// DSD 节点结构：输入节点带 var_id (0=a,1=b,...)
// =====================================================
struct DSDNode {
    int id;                 // 节点编号（1,2,3,...）
    std::string func;       // "in" / "0010" / "10" / "0001" / ...
    std::vector<int> child; // 子节点 id
    int var_id = -1;        // 仅输入节点有效：0=a,1=b,...
};

static std::vector<DSDNode> NODE_LIST;
static int NODE_ID = 1;
static int STEP_ID = 1;

// =====================================================
// Truth table + 变量顺序
//   f01        : 长度 2^n 的 01 串
//   order[i]   : bit i 对应的“变量编号”（0=a,1=b,...）
// =====================================================
struct TT {
    std::string f01;
    std::vector<int> order;
};

// =====================================================
// 工具：01 串 → kitty truth table
// =====================================================
static kitty::dynamic_truth_table make_tt_from01(const std::string& f01)
{
    const size_t len = f01.size();
    assert((len & (len - 1)) == 0);
    unsigned n = static_cast<unsigned>(std::log2(len));

    kitty::dynamic_truth_table tt(n);
    for (uint64_t i = 0; i < len; i++)
        if (f01[i] == '1')
            kitty::set_bit(tt, i);
    return tt;
}

// 支持集 bit 下标（0-based）
static std::vector<int> get_support_bits(const std::string& f01)
{
    auto tt = make_tt_from01(f01);
    std::vector<int> support;
    for (uint8_t i = 0; i < tt.num_vars(); i++)
        if (kitty::has_var(tt, i))
            support.push_back(i);
    return support;
}

// =====================================================
// 按支持集缩减 TT：f01 + order 同步缩减
// =====================================================
static TT shrink_to_support(const TT& in)
{
    auto tt = make_tt_from01(in.f01);

    std::vector<int> supp_bits;
    for (uint8_t i = 0; i < tt.num_vars(); i++)
        if (kitty::has_var(tt, i))
            supp_bits.push_back(i);

    unsigned new_vars = supp_bits.size();
    if (new_vars == tt.num_vars())
        return in;  // 不需要缩减

    kitty::dynamic_truth_table new_tt(new_vars);

    for (uint64_t x = 0; x < (1ull << new_vars); x++)
    {
        uint64_t old_index = 0;
        for (unsigned b = 0; b < new_vars; b++)
        {
            uint64_t bit = (x >> b) & 1;
            old_index |= (bit << supp_bits[b]);   // 映射回原 bit 位置
        }
        if (kitty::get_bit(tt, old_index))
            kitty::set_bit(new_tt, x);
    }

    TT out;
    out.f01.resize(1ull << new_vars);
    for (uint64_t i = 0; i < out.f01.size(); i++)
        out.f01[i] = kitty::get_bit(new_tt, i) ? '1' : '0';

    // 同步缩减变量顺序：新第 b 个变量 = 原来 bit = supp_bits[b] 对应的变量编号
    out.order.reserve(new_vars);
    for (unsigned b = 0; b < new_vars; b++)
    {
        int bit_pos = supp_bits[b];       // 0-based bit index
        out.order.push_back(in.order[bit_pos]);
    }

    return out;
}

// =====================================================
// mul_ui：你原来的 UI 乘法，不动
// =====================================================
static inline std::string mul_ui(const std::string& ui, const std::string& w)
{
    if (ui == "10") return w;
    if (ui == "01") {
        std::string r; r.reserve(w.size());
        for (size_t i = 0; i + 1 < w.size(); i += 2) {
            r.push_back(w[i + 1]);
            r.push_back(w[i]);
        }
        return r;
    }
    if (ui == "11") return std::string(w.size(), '1');
    if (ui == "00") return std::string(w.size(), '0');
    return w;
}

// =====================================================
// 模板计算：保持你原先的逻辑
// =====================================================
struct TemplateResult {
    std::string MF;    // 4 bits
    std::string Mphi;  // block 标签 01 串
    std::string Mpsi;  // 子函数 01 串
};

static TemplateResult run_case_once(
    const std::vector<std::string>& blocks,
    int /*s*/,
    const std::string& S0,
    const std::string& S1)
{
    int m = blocks.size();
    auto W = blocks;

    std::string MF = S0 + S1;
    std::string Mpsi;

    auto try_u = [&](const std::string& u)->bool {
        if (u == "10" || u == "01") {
            std::vector<std::string> cand;
            for (int i = 0; i < m; i++) {
                if (is_constant_block(W[i])) continue;
                auto c = mul_ui(u, W[i]);
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

    std::string exp0 = mul_ui(S0, Mpsi);
    std::string exp1 = mul_ui(S1, Mpsi);

    std::string Mphi; Mphi.reserve(m);
    for (int i = 0; i < m; i++)
        Mphi.push_back(W[i] == exp0 ? '1' : '0');

    return { MF, Mphi, Mpsi };
}

// =====================================================
// 创建节点
// =====================================================
static int new_node(const std::string& func, const std::vector<int>& child)
{
    int id = NODE_ID++;
    NODE_LIST.push_back({ id, func, child, -1 });
    return id;
}

static int new_in_node(int var_id)   // var_id: 0=a,1=b,...
{
    int id = NODE_ID++;
    NODE_LIST.push_back({ id, "in", {}, var_id });
    return id;
}

// =====================================================
// 小规模构树（n<=2）
// =====================================================
static int build_small_tree(const TT& t)
{
    int nv = t.order.size();

    // ======================================================
    //  case nv = 1：单变量函数（长度 2 的真值表）
    // ======================================================
    if (nv == 1)
    {
        int a = new_in_node(t.order[0]);  // 变量节点

        // 01 → NOT x
        if (t.f01 == "01")
            return new_node("01", {a});

        // 10 → identity → 直接返回输入
        if (t.f01 == "10")
            return a;

        // 00 → 常 0
        if (t.f01 == "00")
            return new_node("0", {});

        // 11 → 常 1
        if (t.f01 == "11")
            return new_node("1", {});

        // 不可能走到这里
        return a;
    }

    // ======================================================
    // case nv = 2：普通两变量节点，照常建结构
    // ======================================================
    if (nv == 2)
    {
        int a = new_in_node(t.order[0]);
        int b = new_in_node(t.order[1]);
        return new_node(t.f01, {a, b});
    }

    // ======================================================
    // 兜底：>2 不会在这里触发
    // ======================================================
    return new_node(t.f01, {});
}

// =====================================================
// factor_once_with_reorder_01
// 现在：
//   - Λ 打印为“变量编号集合”（0=a,1=b,...）
//   - 内部 swap-chain 仍使用 bit 位置（单独的 Lambda_bits）
// =====================================================
// =====================================================
//      ⭐ 完整修正版：支持 LSB=a 的变量语义
//      Λ 显示用 0-based 变量编号（0=a）
//      STP 公式内部用 j=1..n 映射：j = n - i
// =====================================================
static bool factor_once_with_reorder_01(
    const TT& in,
    int depth,
    std::string& MF12,
    TT& phi_tt,
    TT& psi_tt)
{
    const string& bin = in.f01;
    int len = bin.size();
    if (!is_power_of_two(len) || len <= 4)
        return false;

    int n = log2(len);
    int r = n / 2;

    auto Mf = binary_to_vec(bin);

    // s 优先级与原版保持一致
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
            // ============================================================
            // 1) Lambda_view 按你习惯的语义：0=a,1=b,...(LSB→MSB)
            // ============================================================
            vector<int> Lambda_view;      // 用于显示
            for (int i = 0; i < n; i++)
                if (v[i]) Lambda_view.push_back(i);

            // 打印 Λ：对应你的语义
            cout << "Λ = { ";
            for (int x : Lambda_view) cout << x << " ";
            cout << "}";

            // ============================================================
            // 2) 映射到 STP 公式的坐标：j = n - i
            //    i = 0..n-1 (LSB→MSB)
            //    j = 1..n   (MSB→LSB)
            // ============================================================
            vector<int> Lambda_j;
            for (int i : Lambda_view)
                Lambda_j.push_back(n - i);   // n-i: 1..n

            // 排序：STP 公式要求 j 升序
            sort(Lambda_j.begin(), Lambda_j.end());

            // ============================================================
            // 3) 生成 swap-chain（使用 STP 论文公式）
            // ============================================================
            vector<vector<stp_data>> chain;

            for (int k = s; k >= 1; k--) {
                int j_k = Lambda_j[k - 1];        // j_k ∈ [1..n]
                int exp = j_k + (s - 1) - k;
                chain.push_back(generate_swap_vec(2, 1 << exp));
            }
            chain.push_back(generate_swap_vec(1 << (n - s), 1 << s));

            auto Mperm  = Vec_chain_multiply(chain, false);
            auto result = Vec_semi_tensor_product(Mf, Mperm);

            // ============================================================
            // 4) 取得重排后的真值表
            // ============================================================
            string reordered;
            reordered.reserve(len);
            for (size_t i = 1; i < result.size(); i++)
                reordered.push_back(result[i] ? '1' : '0');

            cout << " -> reordered = " << reordered << "\n";

            // ============================================================
            // 5) 判断分解类型
            // ============================================================
            int cid = theorem33_case_id(reordered, s);
            if (cid == 0) continue;

            // ============================================================
            // 6) 基于 block 分割（按原逻辑）
            // ============================================================
            int bl = 1 << s;
            int nb = len / bl;
            vector<string> blocks(nb);

            for (int i = 0; i < nb; i++)
                blocks[i] = reordered.substr(i * bl, bl);

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

            // ============================================================
            // 7) 计算新变量顺序（关键！映射保持 0→a）
            // 
            // STP 的 newPos 是按 “j 的排序”
            // 我们要把 j 转回 i：i = n - j
            // ============================================================
            int n_phi = n - s;
            int n_psi = s;

            // STP 认为 newPos = [Ω_j , Λ_j]（都是 j 的排序）
            vector<int> j_Omega;
            vector<bool> inLambda_j(n + 1,false);
            for (int j : Lambda_j) inLambda_j[j] = true;

            for (int j = 1; j <= n; j++)
                if (!inLambda_j[j]) j_Omega.push_back(j);

            // 拼接：Ω + Λ
            vector<int> newPos_j = j_Omega;
            newPos_j.insert(newPos_j.end(), Lambda_j.begin(), Lambda_j.end());

            // 映射到 i （0-based LSB→MSB）
            vector<int> newPos_i;
            for (int j : newPos_j)
                newPos_i.push_back(n - j);    // i = n - j

            // 映射到真正的变量编号
            vector<int> newOrder;
            for (int i : newPos_i)
                newOrder.push_back(in.order[i]);

            vector<int> phi_order(newOrder.begin(), newOrder.begin()+n_phi);
            vector<int> psi_order(newOrder.begin()+n_phi, newOrder.end());

            // 打印
            cout << STEP_ID++ << ". MF = [" << R.MF << "]\n";
            cout << "   MΦ = [" << R.Mphi << "]\n";
            cout << "   Mψ = [" << R.Mpsi << "]\n";
            cout << "   var_order = { ";
            for (int x: newOrder) cout << x << " ";
            cout << "}\n";
            cout << "   Φ vars = { ";
            for (int x: phi_order) cout << x << " ";
            cout << "}  Ψ vars = { ";
            for (int x: psi_order) cout << x << " ";
            cout << "}\n\n";

            // 返回
            MF12          = R.MF;
            phi_tt.f01    = R.Mphi;
            phi_tt.order  = phi_order;
            psi_tt.f01    = R.Mpsi;
            psi_tt.order  = psi_order;

            return true;

        } while (prev_permutation(v.begin(), v.end()));
    }

    return false;
}

// =====================================================
// dsd_factor：递归 DSD
// =====================================================
static int dsd_factor(const TT& f_raw, int depth = 0)
{
    TT f = shrink_to_support(f_raw);

    int len = (int)f.f01.size();
    if (!is_power_of_two(len) || len <= 4)
        return build_small_tree(f);

    std::string MF12;
    TT phi_tt, psi_tt;

    if (!factor_once_with_reorder_01(f, depth, MF12, phi_tt, psi_tt))
        return build_small_tree(f);

    int L = dsd_factor(phi_tt, depth + 1);
    int R = dsd_factor(psi_tt, depth + 1);

    return new_node(MF12, { L, R });
}

// =====================================================
// 打印 DSD（变量数字版）
// =====================================================
static void print_dsd_pretty(int id)
{
    const auto& nd = NODE_LIST[id - 1];

    if (nd.func == "in") {
        std::cout << nd.var_id;
        return;
    }

    std::cout << nd.func << "(";
    for (size_t i = 0; i < nd.child.size(); i++)
    {
        print_dsd_pretty(nd.child[i]);
        if (i + 1 < nd.child.size()) std::cout << ",";
    }
    std::cout << ")";
}

// =====================================================
// 顶层入口：run_dsd_recursive
//  - binary01: 真值表 01 串，按 (a,b,c,...) = (var0,var1,...) 排序
//  - 初始 order = {0,1,2,...,n-1}
//  - 最终把 cur.order 复制到 FINAL_VAR_ORDER（目前就是 {0..n-1}）
// =====================================================
inline bool run_dsd_recursive(const std::string& binary01)
{
    if (!is_power_of_two(binary01.size())) {
        std::cout << "输入长度必须为 2^n\n";
        return false;
    }

    int n = static_cast<int>(std::log2(binary01.size()));

    auto supp = get_support_bits(binary01);
    std::cout << "support bit indices = { ";
    for (auto v : supp) std::cout << v << " ";
    std::cout << "}\n";

    TT root;
    root.f01 = binary01;
    root.order.resize(n);
    for (int i = 0; i < n; i++)
        root.order[i] = i;        // 0=a,1=b,...

    std::cout << "初始变量顺序 (0=a,1=b,...) = { ";
    for (int x : root.order) std::cout << x << " ";
    std::cout << "}\n";

    NODE_LIST.clear();
    NODE_ID = 1;
    STEP_ID = 1;

    std::cout << "输入 = " << binary01
              << " (len = " << binary01.size() << ")\n";

    TT cur = root;      // dsd_factor 不会改 cur（参数按 const&）
    int root_id = dsd_factor(cur);

    std::cout << "===== 最终 DSD 节点列表 =====\n";
    for (auto& nd : NODE_LIST)
    {
        std::cout << nd.id << " = " << nd.func;
        if (nd.func == "in")
            std::cout << "(var=" << nd.var_id << ")";
        else if (nd.child.size() == 1)
            std::cout << "(" << nd.child[0] << ")";
        else if (nd.child.size() == 2)
            std::cout << "(" << nd.child[0] << "," << nd.child[1] << ")";
        std::cout << "\n";
    }

    // 现在 FINAL_VAR_ORDER = {0,1,...,n-1}
    FINAL_VAR_ORDER = cur.order;

    std::cout << "Root = " << root_id << "\n";
    return true;
}
