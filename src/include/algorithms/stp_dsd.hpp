#pragma once
#include <bits/stdc++.h>
using namespace std;
#include <set>
#include <algorithm>
inline std::vector<int> FINAL_VAR_ORDER;
inline int ORIGINAL_VAR_COUNT = 0;
#include "excute.hpp"
#include "reorder.hpp"

// ================================================
// kitty truth table
// ================================================
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

// ================================================
// DSD Node
// ================================================
struct DSDNode {
    int id;
    string func;
    vector<int> child;
    int var_id = -1;
};

static vector<DSDNode> NODE_LIST;
static int NODE_ID = 1;
static int STEP_ID = 1;

// ================================================
// TT = truth table + variable order
// order[i] = x(i+1)
// ================================================
struct TT {
    string f01;
    vector<int> order;
};

// ================================================
// make_tt_from01
// ================================================
static kitty::dynamic_truth_table make_tt_from01(const string& f01)
{
    size_t len = f01.size();
    unsigned n = log2(len);
    kitty::dynamic_truth_table tt(n);

    for (uint64_t i = 0; i < len; i++)
        if (f01[i] == '1')
            kitty::set_bit(tt, i);

    return tt;
}

// ================================================
// get support bits
// ================================================
static vector<int> get_support_bits(const string& f01)
{
    auto tt = make_tt_from01(f01);
    vector<int> supp;

    for (int i = 0; i < tt.num_vars(); i++){
        if (kitty::has_var(tt, i))
            supp.push_back(i);
    std::cout <<"i="<<i<<endl;
    }
    return supp;
}

// ================================================
// shrink_to_support（按 support 缩减 TT）
// ================================================
static TT shrink_to_support(const TT& in)
{
    auto tt = make_tt_from01(in.f01);

    vector<int> supp;
    for (int i = 0; i < tt.num_vars(); i++)
        if (kitty::has_var(tt, i))
            supp.push_back(i);

    if (supp.size() == tt.num_vars())
        return in;

    unsigned nv = supp.size();
    kitty::dynamic_truth_table new_tt(nv);

    for (uint64_t x = 0; x < (1ull << nv); x++)
    {
        uint64_t old = 0;
        for (int b = 0; b < nv; b++)
        {
            uint64_t bit = (x >> b) & 1;
            old |= (bit << supp[b]);
        }
        if (kitty::get_bit(tt, old))
            kitty::set_bit(new_tt, x);
    }

    TT out;
    out.f01.resize(1ull << nv);

    for (uint64_t i = 0; i < out.f01.size(); i++)
        out.f01[i] = kitty::get_bit(new_tt, i) ? '1' : '0';

    out.order.reserve(nv);
    for (int b : supp)
        out.order.push_back(in.order[b]);

    return out;
}

// ================================================
// mul_ui
// ================================================
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

// ================================================
// TemplateResult
// ================================================
struct TemplateResult {
    string MF;
    string Mphi;
    string Mpsi;
};

// ================================================
// run_case_once
// ================================================
static TemplateResult run_case_once(
    const vector<string>& blocks,
    int s,
    const string& S0,
    const string& S1)
{
    int m = blocks.size();
    vector<string> W = blocks;

    string MF = S0 + S1;
    string Mpsi;

    auto try_u = [&](const string& u)->bool {
        if (u=="10" || u=="01") {
            vector<string> cand;
            for (int i = 0; i < m; i++) {
                if (is_constant_block(W[i])) continue;
                string c = mul_ui(u, W[i]);
                if (mul_ui(u, c) == W[i])
                    cand.push_back(c);
            }
            if (cand.empty()) return false;
            for (int k=1;k<cand.size();k++)
                if (cand[k] != cand[0]) return false;
            Mpsi = cand[0];
            return true;
        }
        return false;
    };

    if (!try_u(S0) && !try_u(S1))
    {
        int p = -1;
        for (int i = 0; i < m; i++)
            if (!is_constant_block(W[i])) { p=i; break; }
        if (p < 0) p = 0;
        Mpsi = W[p];
    }

    string exp0 = mul_ui(S0, Mpsi);
    string exp1 = mul_ui(S1, Mpsi);

    string Mphi;
    Mphi.reserve(m);
    for (int i = 0; i < m; i++)
        Mphi.push_back(W[i] == exp0 ? '1' : '0');

    return { MF, Mphi, Mpsi };
}
// =====================================================
// new_node / new_in_node
// =====================================================
static int new_node(const string& func, const vector<int>& child)
{
    int id = NODE_ID++;
    NODE_LIST.push_back({ id, func, child, -1 });
    return id;
}

static int new_in_node(int var_id)  // var_id = 1..n
{
    int id = NODE_ID++;
    NODE_LIST.push_back({ id, "in", {}, var_id });
    return id;
}

// =====================================================
// build_small_tree
// =====================================================
// =====================================================
// build_small_tree - 修正版：记录变量到 FINAL_VAR_ORDER
// =====================================================
static int build_small_tree(const TT& t)
{
    int nv = t.order.size();

    if (nv == 1)
    {
        int var_id = t.order[0];  // 原始变量编号
        int a = new_in_node(var_id);
        
        // 🔥 记录到 FINAL_VAR_ORDER（如果还没记录）
        if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), var_id) 
            == FINAL_VAR_ORDER.end())
        {
            FINAL_VAR_ORDER.push_back(var_id);
        }
        
        if (t.f01 == "10") return a;                  // identity
        if (t.f01 == "01") return new_node("01",{a}); // NOT
        if (t.f01 == "00") return new_node("0",{});   // const 0
        if (t.f01 == "11") return new_node("1",{});   // const 1
        return a;
    }

    if (nv == 2)
    {
        // 🔥 记录两个变量
        for (int var_id : t.order)
        {
            if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), var_id) 
                == FINAL_VAR_ORDER.end())
            {
                FINAL_VAR_ORDER.push_back(var_id);
            }
        }
        
        int a = new_in_node(t.order[0]);
        int b = new_in_node(t.order[1]);
        return new_node(t.f01, { a, b });
    }

    return new_node(t.f01, {});
}
// =====================================================
// factor_once_with_reorder_01
// 完整正确 STP 重排：
//   - TT.order[i] = 变量编号（1-based）
//   - Lambda 显示用变量编号
//   - Lambda_j = STP 论文中 j = n-i
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
            vector<int> Lambda_bits;
            for (int i = 0; i < n; i++)
                if (v[i]) Lambda_bits.push_back(i);

            // 🔥 关键：STP 中 j 的定义
            // bit i 对应局部编号 j = n - i
            // 但在 in.order 中：
            //   in.order[0] 是局部编号 1 的原始变量
            //   in.order[1] 是局部编号 2 的原始变量
            //   ...
            // 所以：局部编号 j 对应 in.order[j-1]
            
            vector<int> Lambda_j;
            for (int bit : Lambda_bits)
                Lambda_j.push_back(n - bit);
            
            sort(Lambda_j.begin(), Lambda_j.end());

            cout << "Λ = { ";
            for (int j : Lambda_j) cout << j << " ";
            cout << "}";

            // 生成 swap-chain
            vector<vector<stp_data>> chain;

            for (int k = s; k >= 1; k--) {
                int j_k = Lambda_j[k - 1];
                int exp = j_k + (s - 1) - k;
                chain.push_back(generate_swap_vec(2, 1 << exp));
            }
            chain.push_back(generate_swap_vec(1 << (n - s), 1 << s));

            auto Mperm  = Vec_chain_multiply(chain, false);
            auto result = Vec_semi_tensor_product(Mf, Mperm);

            string reordered;
            reordered.reserve(len);
            for (size_t i = 1; i < result.size(); i++)
                reordered.push_back(result[i] ? '1' : '0');

            cout << " -> reordered = " << reordered << "\n";

            int cid = theorem33_case_id(reordered, s);
            if (cid == 0) continue;

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
                        S_list = {
                            {"11","10"}, {"11","01"},
                            {"10","11"}, {"01","11"}
                        };
                    else
                        S_list = {
                            {"00","10"}, {"00","01"},
                            {"10","00"}, {"01","00"}
                        };
                    break;
                case 3: S_list = { {"10","10"}, {"01","01"} }; break;
                case 4: S_list = { {"10","01"}, {"01","10"} }; break;
                case 5: return false;
            }

            auto R = run_case_once(blocks, s, S_list[0].first, S_list[0].second);

            // 🔥 计算新的局部编号顺序
            int n_phi = n - s;

            vector<bool> inLam_j(n + 1, false);
            for (int j : Lambda_j) inLam_j[j] = true;

            vector<int> Omega_j;
            for (int j = 1; j <= n; j++)
                if (!inLam_j[j]) Omega_j.push_back(j);

            vector<int> newPos_j = Omega_j;
            newPos_j.insert(newPos_j.end(), Lambda_j.begin(), Lambda_j.end());

            // 🔥🔥🔥 关键修正：
            // in.order[i] 存储的是位置 (i+1) 的原始变量编号
            // 即：局部编号 j 对应 in.order[j-1]
            
            vector<int> newOrder_original;
            for (int j : newPos_j) {
                newOrder_original.push_back(in.order[j - 1]);  // j从1开始，数组从0开始
            }

            vector<int> phi_order_original(newOrder_original.begin(), 
                                          newOrder_original.begin() + n_phi);
            vector<int> psi_order_original(newOrder_original.begin() + n_phi, 
                                          newOrder_original.end());

            cout << STEP_ID++ << ". MF = [" << R.MF << "]\n";
            cout << "   MΦ = [" << R.Mphi << "]\n";
            cout << "   MΨ = [" << R.Mpsi << "]\n";
            
            cout << "   重排详情：\n";
            for (int i = 0; i < newPos_j.size(); i++) {
                int j = newPos_j[i];
                int orig = in.order[j - 1];
                cout << "     新位置" << (i+1) << " = 局部编号" << j 
                     << " → 原始变量" << orig << "\n";
            }
            
            cout << "   新局部顺序 = { ";
            for (int j : newPos_j) cout << j << " ";
            cout << "}\n";
            
            cout << "   新原始变量顺序 = { ";
            for (int v : newOrder_original) cout << v << " ";
            cout << "}\n";
            
            cout << "   Φ 原始变量 = { ";
            for (int v : phi_order_original) cout << v << " ";
            cout << "}  Ψ 原始变量 = { ";
            for (int v : psi_order_original) cout << v << " ";
            cout << "}\n\n";

            MF12 = R.MF;
            phi_tt.f01 = R.Mphi;
            psi_tt.f01 = R.Mpsi;
            
            // 🔥 order[i] 存储位置 (i+1) 的原始变量编号
            phi_tt.order = phi_order_original;
            psi_tt.order = psi_order_original;

            return true;

        } while (prev_permutation(v.begin(), v.end()));
    }

    return false;
}
// dsd_factor - 修正版
// =====================================================
static int dsd_factor(const TT& f_raw, int depth=0)
{
    TT f = shrink_to_support(f_raw);

    int len = f.f01.size();
    if(len <= 4)  
        return build_small_tree(f);

    string MF12;
    TT phi_tt, psi_tt;

    if(!factor_once_with_reorder_01(f, depth, MF12, phi_tt, psi_tt))
        return build_small_tree(f);

    vector<int> phi_original_vars = phi_tt.order;
    vector<int> psi_original_vars = psi_tt.order;
    
    int n_phi = phi_tt.order.size();
    int n_psi = psi_tt.order.size();
    
    cout << "📌 递归分解 Φ：原始变量 { ";
    for (int v : phi_original_vars) cout << v << " ";
    cout << "} → 局部编号 { ";
    for (int i = 1; i <= n_phi; i++) cout << i << " ";
    cout << "}\n";
    cout << "   映射关系：";
    for (int i = 0; i < n_phi; i++)
        cout << "位置" << (i+1) << "→变量" << phi_original_vars[i] << " ";
    cout << "\n";
    
    cout << "📌 递归分解 Ψ：原始变量 { ";
    for (int v : psi_original_vars) cout << v << " ";
    cout << "} → 局部编号 { ";
    for (int i = 1; i <= n_psi; i++) cout << i << " ";
    cout << "}\n";
    cout << "   映射关系：";
    for (int i = 0; i < n_psi; i++)
        cout << "位置" << (i+1) << "→变量" << psi_original_vars[i] << " ";
    cout << "\n\n";

    int L = dsd_factor(phi_tt, depth+1);
    int R = dsd_factor(psi_tt, depth+1);

    return new_node(MF12,{L,R});
}
inline bool run_dsd_recursive(const std::string& binary01)
{
    if (!is_power_of_two(binary01.size())) {
        std::cout << "输入长度必须为 2^n\n";
        return false;
    }

    int n = static_cast<int>(std::log2(binary01.size()));
    ORIGINAL_VAR_COUNT = n;
    
    TT root;
    root.f01 = binary01;
    root.order.resize(n);

    // 🔥 修正：order[i] 存储位置 (i+1) 的原始变量编号
    // 初始时：位置1→变量1, 位置2→变量2, ...
    for (int i = 0; i < n; ++i)
        root.order[i] = i + 1;  // 位置 (i+1) 对应变量 (i+1)

    std::cout << "输入 = " << binary01 << " (n=" << n << ")\n";
    std::cout << "初始映射：";
    for (int i = 0; i < n; i++)
        std::cout << "位置" << (i+1) << "→变量" << root.order[i] << " ";
    std::cout << "\n\n";

    NODE_LIST.clear();
    NODE_ID = 1;
    STEP_ID = 1;

    int root_id = dsd_factor(root);

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

    std::cout << "Root = " << root_id << "\n";
    return true;
}