#pragma once
#include <bits/stdc++.h>
using namespace std;

// bench 输出会用到最终变量顺序
inline std::vector<int> FINAL_VAR_ORDER;

#include "excute.hpp"
#include "reorder.hpp"

// kitty 真值表
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

// =====================================================
// DSD 节点结构：输入节点带 var_id (1=a,2=b,...)
// =====================================================
struct DSDNode {
    int id;                 // 节点编号（1,2,3,...）
    std::string func;       // 真值表 / MF 串（"in" / "0010" / "10" ...）
    std::vector<int> child; // 子节点编号（指向 NODE_LIST 中的 id）
    int var_id = -1;        // 对输入节点："in" 时存变量编号；其它节点 = -1
};

static std::vector<DSDNode> NODE_LIST;
static int NODE_ID = 1;
static int STEP_ID = 1;

// =====================================================
// Truth table + 变量顺序 封装
//   f01  : 01 串，长度 2^n
//   order[i] : 第 i 个 bit 变量（低位）对应的原始变量编号(1=a,2=b,...)
// =====================================================
struct TT {
    std::string f01;
    std::vector<int> order;
};

// =====================================================
// 工具：01 串 <-> kitty truth table
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
// 按支持集缩减 TT（f01 + order 同步缩减）
// =====================================================
static TT shrink_to_support(const TT& in)
{
    auto tt = make_tt_from01(in.f01);

    std::vector<int> support_bits;
    for (uint8_t i = 0; i < tt.num_vars(); i++)
        if (kitty::has_var(tt, i))
            support_bits.push_back(i);

    unsigned new_vars = support_bits.size();
    if (new_vars == tt.num_vars())
        return in;  // 支持集已经覆盖全部变量，不用缩减

    kitty::dynamic_truth_table new_tt(new_vars);

    for (uint64_t x = 0; x < (1ull << new_vars); x++)
    {
        uint64_t old_index = 0;
        for (unsigned b = 0; b < new_vars; b++)
        {
            uint64_t bit = (x >> b) & 1;
            old_index |= (bit << support_bits[b]);  // 映射回原来的 bit 位置
        }
        if (kitty::get_bit(tt, old_index))
            kitty::set_bit(new_tt, x);
    }

    TT out;
    out.f01.resize(1ull << new_vars);
    for (uint64_t i = 0; i < out.f01.size(); i++)
        out.f01[i] = kitty::get_bit(new_tt, i) ? '1' : '0';

    // 同步缩减变量顺序：新第 b 个变量对应旧的 support_bits[b]
    out.order.reserve(new_vars);
    for (unsigned b = 0; b < new_vars; b++)
    {
        int bit_pos = support_bits[b];      // 0-based bit index
        out.order.push_back(in.order[bit_pos]); // 映射到原始变量编号
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
// 模板计算
// =====================================================
struct TemplateResult {
    std::string MF;    // 4 bit (u1,u2)
    std::string Mphi;  // block 标签 01 串，长度 = #blocks = 2^(n-s)
    std::string Mpsi;  // 子函数 01 串，长度 = block 大小 = 2^s
};

static TemplateResult run_case_once(
    const std::vector<std::string>& blocks01,
    int /*s*/,
    const std::string& S0,
    const std::string& S1)
{
    int m = blocks01.size();
    std::vector<std::string> W = blocks01;

    std::string MF = S0 + S1;
    std::string Mpsi;

    auto try_u = [&](const std::string& u)->bool {
        if (u == "10" || u == "01") {
            std::vector<std::string> cand;
            for (int i = 0; i < m; i++) {
                if (is_constant_block(W[i])) continue;
                std::string c = mul_ui(u, W[i]);
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
// 构造节点
// =====================================================
static int new_node(const std::string& func, const std::vector<int>& child)
{
    int id = NODE_ID++;
    NODE_LIST.push_back({ id, func, child, -1 });
    return id;
}

// 构造“输入”节点：绑定变量编号 var_id
static int new_in_node(int var_id)
{
    int id = NODE_ID++;
    DSDNode nd;
    nd.id     = id;
    nd.func   = "in";
    nd.child  = {};
    nd.var_id = var_id;
    NODE_LIST.push_back(nd);
    return id;
}

// 小规模真值表构树（n<=2），根据 TT.order 绑定变量编号
static int build_small_tree(const TT& t)
{
    int nv = (int)t.order.size();

    if (nv == 1)
    {
        int a = new_in_node(t.order[0]);        // 唯一变量
        return new_node(t.f01, { a });
    }

    if (nv == 2)
    {
        int a = new_in_node(t.order[0]);        // 第一个变量
        int b = new_in_node(t.order[1]);        // 第二个变量
        return new_node(t.f01, { a, b });
    }

    // 安全兜底（几乎不会到这里）
    return new_node(t.f01, {});
}

// =====================================================
// 单步分解（输入 / 输出都是 TT）
// 负责：
//   - 计算重排后的真值表
//   - 根据 Λ 得到新的变量顺序 newOrder
//   - 划分 Φ / Ψ 对应变量集合，填到 phi_tt.order / psi_tt.order
// =====================================================
static bool factor_once_with_reorder_01(
    const TT& in,
    int /*depth*/,
    std::string& MF12,
    TT& phi_tt,
    TT& psi_tt)
{
    const std::string& binary01 = in.f01;

    int len = binary01.size();
    if (!is_power_of_two(len) || len <= 4) return false;

    int n = (int)std::log2(len);  // 当前变量个数
    int r = n / 2;

    std::vector<stp_data> Mf = binary_to_vec(binary01);

    // s 的尝试顺序（优先 2）
    std::vector<int> s_order;
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
        std::vector<bool> v(n);
        std::fill(v.begin(), v.begin() + s, true);

        do {
            // 选中的 bit 位置（1-based）
            std::vector<int> Lambda;
            for (int i = 0; i < n; i++)
                if (v[i]) Lambda.push_back(i + 1);

            // 构造 swap-chain
            std::vector<std::vector<stp_data>> chain;

            for (int k = s; k >= 1; k--) {
                int j_k = Lambda[k - 1];
                int exp = j_k + (s - 1) - k;
                chain.push_back(generate_swap_vec(2, 1 << exp));
            }
            chain.push_back(generate_swap_vec(1 << (n - s), 1 << s));

            auto Mperm  = Vec_chain_multiply(chain, false);
            auto result = Vec_semi_tensor_product(Mf, Mperm);

            std::string reordered;
            reordered.reserve(len);
            for (size_t i = 1; i < result.size(); i++)
                reordered.push_back(result[i] ? '1' : '0');

            std::cout << "Λ = { ";
            for (int j : Lambda) std::cout << j << " ";
            std::cout << "} -> reordered = " << reordered << "\n";

            int cid = theorem33_case_id(reordered, s);
            if (cid == 0) continue;

            // === 计算新的变量顺序 newOrder ===
            // Lambda 是选中的变量位置（1-based）
            std::vector<bool> inLambda(n + 1, false);
            for (int p : Lambda) inLambda[p] = true;

            std::vector<int> OmegaPos;
            for (int p = 1; p <= n; ++p)
                if (!inLambda[p]) OmegaPos.push_back(p);

            // 新 bit 位置顺序：先 Ω，再 Λ
            std::vector<int> newPos;
            newPos.reserve(n);
            for (int p : OmegaPos) newPos.push_back(p);
            for (int p : Lambda)   newPos.push_back(p);

            std::vector<int> newOrder;
            newOrder.reserve(n);
            for (int p : newPos)
                newOrder.push_back(in.order[p - 1]); // 映射到变量编号

            // === 分块 ===
            int bl = 1 << s;
            int nb = len / bl;
            std::vector<std::string> blocks(nb);
            for (int i = 0; i < nb; i++)
                blocks[i] = reordered.substr(i * bl, bl);

            bool has1 = false, has0 = false;
            for (auto& b : blocks) {
                if (is_constant_block(b)) {
                    if (b[0] == '1') has1 = true;
                    if (b[0] == '0') has0 = true;
                }
            }

            std::vector<std::pair<std::string,std::string>> S_list;
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

            // Mphi 长度 = 2^(n-s)，Mpsi 长度 = 2^s
            int n_phi = n - s;
            int n_psi = s;

            std::vector<int> phi_order(newOrder.begin(), newOrder.begin() + n_phi);
            std::vector<int> psi_order(newOrder.begin() + n_phi, newOrder.end());

            std::cout << STEP_ID++ << ". MF  = [" << R.MF  << "]\n";
            std::cout << "   MΦ  = [" << R.Mphi << "]\n";
            std::cout << "   Mψ  = [" << R.Mpsi << "]\n";
            std::cout << "   var_order = { ";
            for (auto x : newOrder) std::cout << x << " ";
            std::cout << "}\n";
            std::cout << "   Φ vars = { ";
            for (auto x : phi_order) std::cout << x << " ";
            std::cout << "}   Ψ vars = { ";
            for (auto x : psi_order) std::cout << x << " ";
            std::cout << "}\n\n";

            MF12         = R.MF;
            phi_tt.f01   = R.Mphi;
            phi_tt.order = std::move(phi_order);
            psi_tt.f01   = R.Mpsi;
            psi_tt.order = std::move(psi_order);
            return true;

        } while (std::prev_permutation(v.begin(), v.end()));
    }

    return false;
}

// =====================================================
// 主递归 DSD （用 TT）
// =====================================================
static int dsd_factor(const TT& f_raw, int depth = 0)
{
    // 先按支持集缩减（同步变量顺序）
    TT f = shrink_to_support(f_raw);

    int len = f.f01.size();
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
// 可选：用数字递归打印 DSD 结构（变量用 1,2,3,... 表示）
// =====================================================
static void print_dsd_pretty(int id)
{
    const DSDNode& nd = NODE_LIST[id - 1];

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
//   - 输入 binary01 (如从 hex 转出的 01 串，按 a,b,c,d,... 顺序)
//   - 初始化变量顺序为 {1,2,...,n} 对应 {a,b,...}
//   - 递归分解，构造 NODE_LIST
//   - 把最终变量顺序保存到 FINAL_VAR_ORDER
// =====================================================
inline bool run_dsd_recursive(const std::string& binary01)
{
    if (!is_power_of_two(binary01.size())) {
        std::cout << "输入长度必须为 2^n\n";
        return false;
    }

    int n = static_cast<int>(std::log2(binary01.size()));

    // 打印 bit-level 支持集（原始）
    auto support_bits = get_support_bits(binary01);
    std::cout << "support bit indices = { ";
    for (auto v : support_bits) std::cout << v << " ";
    std::cout << "}\n";

    // 初始化 TT：变量顺序为 1..n （a,b,c,...）
    TT root0;
    root0.f01 = binary01;
    root0.order.resize(n);
    for (int i = 0; i < n; ++i)
        root0.order[i] = i + 1;

    std::cout << "初始变量顺序 (1=a,2=b,...) = { ";
    for (auto x : root0.order) std::cout << x << " ";
    std::cout << "}\n";

    NODE_LIST.clear();
    NODE_ID = 1;
    STEP_ID = 1;

    std::cout << "输入 = " << binary01
              << " (len = " << binary01.size() << ")\n";

    // 用一个可修改副本 top 做递归（支持集缩减会改它的 order）
    TT top = root0;
    int root_id = dsd_factor(top);

    std::cout << "===== 最终 DSD 节点列表 =====\n";
    for (auto& nd : NODE_LIST) {
        std::cout << nd.id << " = " << nd.func;
        if (nd.func == "in") {
            std::cout << "(var=" << nd.var_id << ")";
        } else if (nd.child.size() == 1) {
            std::cout << "(" << nd.child[0] << ")";
        } else if (nd.child.size() == 2) {
            std::cout << "(" << nd.child[0] << "," << nd.child[1] << ")";
        }
        std::cout << "\n";
    }

    // 保存最终变量顺序（可能已被 shrink_to_support 缩减）
    FINAL_VAR_ORDER = top.order;

    std::cout << "Root = " << root_id << "\n";
    return true;
}
