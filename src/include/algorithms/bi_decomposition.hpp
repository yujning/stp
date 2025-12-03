#pragma once
#include <vector>
#include <string>
#include <iostream>

// 直接包含 stp_dsd.hpp，里面已经有 struct TT 定义
#include "stp_dsd.hpp"

using std::string;
using std::vector;

// =====================================================
// Bi-decomposition 结果结构
// =====================================================
struct BiDecompResult {
    int k1, k2, k3;             // |Gamma|, |Theta|, |Lambda|
    vector<int> Gamma;          // Γ 变量编号
    vector<int> Theta;          // Θ 变量编号
    vector<int> Lambda;         // Λ 变量编号

    string F01;                 // F(u,v) 的 4bit 真值表，顺序 00,01,10,11

    TT phi_tt;                  // φ 的真值表 + 变量顺序（Gamma ∪ Theta）
    TT psi_tt;                  // ψ 的真值表 + 变量顺序（Theta ∪ Lambda）
};

// =====================================================
// 小工具：本文件用自己的函数名，避免和 reorder.hpp 冲突
// =====================================================
static inline bool bi_is_power_of_two(size_t x)
{
    return x && ((x & (x - 1)) == 0);
}

// (u,v) in {0,1}^2, idx = (u<<1)|v
static inline int bi_F_eval(int Fmask, int u, int v)
{
    int idx = (u << 1) | v;
    return (Fmask >> idx) & 1;
}

// 是否常数块（01 域）
static inline bool bi_is_constant_block(const string& b)
{
    if (b.empty()) return false;
    char c0 = b[0];
    for (char c : b)
        if (c != c0) return false;
    return true;
}

// =====================================================
// 一条约束：φ(row) 与 ψ(col) 在 F 下要等于 f
// =====================================================
struct BiConstraint {
    int r;   // row index  → φ 的输入索引 (Gamma,Theta)
    int c;   // col index  → ψ 的输入索引 (Theta,Lambda)
    int f;   // f 的输出 0/1
};

static inline bool bi_check_constraint_partial(
    const BiConstraint& con,
    const vector<int>& U,
    const vector<int>& V,
    int Fmask)
{
    int ur = U[con.r];
    int vc = V[con.c];
    // 这里 U 始终是已定的（0/1），V 可能是 -1
    if (vc == -1) return true;
    return bi_F_eval(Fmask, ur, vc) == con.f;
}

// =====================================================
// 只回溯 V（ψ），U 已经由块模式固定好
// =====================================================
static bool bi_backtrack_V(
    int c_idx,
    vector<int>& V,
    const vector<BiConstraint>& cons,
    const vector<vector<int>>& consByCol,
    const vector<int>& U,
    int Nc,
    int Fmask)
{
    if (c_idx == Nc)
        return true;  // 所有列都赋值完成

    for (int val = 0; val <= 1; ++val)
    {
        V[c_idx] = val;
        bool ok = true;
        for (int cid : consByCol[c_idx])
        {
            if (!bi_check_constraint_partial(cons[cid], U, V, Fmask)) {
                ok = false;
                break;
            }
        }
        if (ok && bi_backtrack_V(c_idx + 1, V, cons, consByCol, U, Nc, Fmask))
            return true;

        V[c_idx] = -1;
    }
    return false;
}

// =====================================================
// 调试：按 (k1,k2,k3) 分块打印真值表
// =====================================================
static inline void bi_print_blocks_debug(const string& f01,
                                         int k1, int k2, int k3)
{
    int R = 1 << k1;
    int C = 1 << k2;
    int B = 1 << k3;

    std::cout << "[Block view] k1=" << k1
              << " k2=" << k2
              << " k3=" << k3
              << "  (R=" << R
              << ", C=" << C
              << ", B=" << B << ")\n";

    for (int r = 0; r < R; ++r) {
        std::cout << "    row " << r << " : ";
        for (int c = 0; c < C; ++c) {
            string blk;
            blk.reserve(B);
            for (int p = 0; p < B; ++p) {
                int idx = (r << (k2 + k3)) | (c << k3) | p;
                blk.push_back(f01[idx]);
            }
            std::cout << blk << "  ";
        }
        std::cout << "\n";
    }
}

// =====================================================
// 主函数：固定顺序、k2=1、φ 按“常数=1、非常数=0”固定，
//        然后搜索所有 (F, ψ) 使得 f = F(φ, ψ)
// =====================================================
static std::vector<BiDecompResult>
enumerate_bi_decomposition_k2_eq1_no_reorder(const TT& in)
{
    std::vector<BiDecompResult> results;

    const string& f01 = in.f01;
    size_t len = f01.size();
    if (!bi_is_power_of_two(len))
        return results;

    int n = (int)std::log2((double)len);
    if ((int)in.order.size() != n)
        return results;
    if (n < 3)
        return results;

    int k2 = 1;
    int max_k3 = (n - k2) / 2;  // k3 = 1..floor((n-k2)/2)

    for (int k3 = 1; k3 <= max_k3; ++k3)
    {
        int k1 = n - k2 - k3;
        if (k1 <= 0) continue;

        // ===== 打印块视图（方便你对照论文） =====
        bi_print_blocks_debug(f01, k1, k2, k3);

        // ===== Γ, Θ, Λ 变量集合 =====
        vector<int> Gamma, Theta, Lambda;
        Gamma.reserve(k1);
        for (int i = 0; i < k1; ++i)
            Gamma.push_back(in.order[i]);

        Theta.push_back(in.order[k1]);   // k2 = 1

        Lambda.reserve(k3);
        for (int i = 0; i < k3; ++i)
            Lambda.push_back(in.order[k1 + 1 + i]);

        // φ 输入维度: k1 + k2
        // ψ 输入维度: k2 + k3
        int Nr = 1 << (k1 + k2);   // #rows = 2^(k1+1)
        int Nc = 1 << (k2 + k3);   // #cols = 2^(1+k3)

        // ===== 用连续分块来固定 φ（常数=1，非常数=0） =====
        int B = 1 << k3;
        vector<int> U(Nr, 0);      // U[r] = φ(Γ,Θ)
        for (int r = 0; r < Nr; ++r)
        {
            int start = r * B;
            string blk = f01.substr(start, B);
            U[r] = bi_is_constant_block(blk) ? 1 : 0;
        }
        // 注意：对应 11100110，k1=k2=k3=1：
        // blocks = [11,10,01,10] → U = [1,0,0,0] → φ = 1000

        // ===== 构造所有约束 =====
        vector<BiConstraint> cons;
        cons.reserve(len);

        vector<vector<int>> consByCol(Nc);

        for (int a = 0; a < (1 << k1); ++a)
        for (int b = 0; b < (1 << k2); ++b)       // b = 0/1
        for (int c = 0; c < (1 << k3); ++c)
        {
            int idx = (a << (k2 + k3)) | (b << k3) | c;
            int fv  = (f01[idx] == '1');

            int r   = (a << k2) | b;   // φ 输入索引 (Γ,Θ)
            int col = (b << k3) | c;   // ψ 输入索引 (Θ,Λ)

            int cid = (int)cons.size();
            cons.push_back({ r, col, fv });
            consByCol[col].push_back(cid);
        }

        // ===== 遍历所有 16 种 F，找 ψ 的可行解 =====
        for (int Fmask = 0; Fmask < 16; ++Fmask)
        {
            vector<int> V(Nc, -1); // ψ 的值

            if (!bi_backtrack_V(0, V, cons, consByCol, U, Nc, Fmask))
                continue;   // 这个 F 不行

            // --- 找到一组 (F, ψ) 与固定 φ 相容，记录结果 ---
            BiDecompResult R;
            R.k1 = k1;
            R.k2 = k2;
            R.k3 = k3;
            R.Gamma = Gamma;
            R.Theta = Theta;
            R.Lambda = Lambda;

            // F01
            R.F01.assign(4, '0');
            for (int u = 0; u <= 1; ++u)
            for (int v = 0; v <= 1; ++v) {
                int idx = (u << 1) | v;
                R.F01[idx] = bi_F_eval(Fmask, u, v) ? '1' : '0';
            }

            // φ truth table：直接用 U
            R.phi_tt.f01.assign(Nr, '0');
            for (int r = 0; r < Nr; ++r)
                R.phi_tt.f01[r] = U[r] ? '1' : '0';

            R.phi_tt.order.clear();
            R.phi_tt.order.insert(R.phi_tt.order.end(),
                                  Gamma.begin(), Gamma.end());
            R.phi_tt.order.insert(R.phi_tt.order.end(),
                                  Theta.begin(), Theta.end());

            // ψ truth table：用 V
            R.psi_tt.f01.assign(Nc, '0');
            for (int c = 0; c < Nc; ++c)
                R.psi_tt.f01[c] = V[c] ? '1' : '0';

            R.psi_tt.order.clear();
            R.psi_tt.order.insert(R.psi_tt.order.end(),
                                  Theta.begin(), Theta.end());
            R.psi_tt.order.insert(R.psi_tt.order.end(),
                                  Lambda.begin(), Lambda.end());

            results.push_back(std::move(R));
        }
    }

    return results;
}

