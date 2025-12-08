#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include "reorder.hpp"
#include "stp_dsd.hpp"

using std::string;
using std::vector;

// =====================================================
// BiDecompResult
// =====================================================
struct BiDecompResult {
    int k1, k2, k3;
    vector<int> Gamma;
    vector<int> Theta;
    vector<int> Lambda;
    string F01;
    TT phi_tt;
    TT psi_tt;
};

// 判断块是不是常数（全0或全1）
static bool is_const_block(const string& s)
{
    if (s.empty()) return false;
    char b0 = s[0];
    for (char b : s)
        if (b != b0) return false;
    return true;
}

// 判断两块是否按位互补
static bool is_complement_block(const string& a, const string& b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
    {
        if (a[i] == b[i]) return false;       // 互补必须处处相反
    }
    return true;
}

// 打印结构矩阵（调试用，保持简单）
static void print_structure_matrix(
        int k1, int k2, int k3,
        const vector<vector<string>> &blk,
        const vector<int> &phi_bits,
        const vector<int> &psi_bits)
{
    int R = 1 << k1;
    int C = 1 << k2;

    std::cout << "Mf blocks:\n";
    for (int r = 0; r < R; r++)
    {
        std::cout << "  Row " << r << " : ";
        for (int c = 0; c < C; c++)
            std::cout << blk[r][c] << " ";
        std::cout << "\n";
    }

    std::cout << "\nφ table bits:\n";
    int idx = 0;
    for (int r = 0; r < R; r++)
    {
        std::cout << "  row " << r << " : ";
        for (int c = 0; c < C; c++)
            std::cout << phi_bits[idx++] << " ";
        std::cout << "\n";
    }

    std::cout << "\nψ bits:\n  ";
    for (int b : psi_bits) std::cout << b;
    std::cout << "\n\n";
}



// 根据公式(34)进行真值表重排
static string apply_variable_reordering(
    const string& f01,
    int n,
    const vector<int>& Gamma_indices,  // 1-based
    const vector<int>& Theta_indices,  // 1-based
    const vector<int>& Lambda_indices, // 1-based
    int k1, int k2, int k3)
{
    // 将 f01 转换为向量形式
    vector<stp_data> Mf = binary_to_vec(f01);

    // 构造交换矩阵链（从右到左按公式34的顺序）
    vector<vector<stp_data>> swap_chain;

    // 第一部分：W[2^k1, 2^k2]
    swap_chain.push_back(generate_swap_vec(std::pow(2, k1), std::pow(2, k2)));

    // 第二部分：⊗_{i=k2}^1 W[2, 2^{j_i+(k2-1)-i}]
    for (int i = k2; i >= 1; --i)
    {
        int j_i = Theta_indices[i - 1];  // Θ 的第 i 个变量（1-based）
        int exp = j_i + (k2 - 1) - i;
        swap_chain.push_back(generate_swap_vec(2, std::pow(2, exp)));
    }

    // 第三部分：W[2^{k1+k2}, 2^k3]
    swap_chain.push_back(generate_swap_vec(std::pow(2, k1 + k2), std::pow(2, k3)));

    // 第四部分：⊗_{i=k3}^1 W[2, 2^{j_i+(k3-1)-i}]
    for (int i = k3; i >= 1; --i)
    {
        int j_i = Lambda_indices[i - 1];  // Λ 的第 i 个变量（1-based）
        int exp = j_i + (k3 - 1) - i;
        swap_chain.push_back(generate_swap_vec(2, std::pow(2, exp)));
    }

    // 矩阵链乘法
    vector<stp_data> Mperm = Vec_chain_multiply(swap_chain, false);
    vector<stp_data> result = Vec_semi_tensor_product(Mf, Mperm);

    // 转换回字符串（跳过第一个维度元素）
    string reordered;
    for (size_t i = 1; i < result.size(); ++i)
        reordered.push_back(result[i] ? '1' : '0');

    return reordered;
}

// u 作用在一个长度为 B 的块 P 上：
// "10" -> 恒等  (返回 P)
// "01" -> 取反  (返回 ~P)
// "00" -> 全 0
// "11" -> 全 1
static string apply_u_to_block(const string& P, const string& u)
{
    if (u == "10")
    {
        // 单位矩阵：uM = M
        return P;
    }
    else if (u == "01")
    {
        // 反对角：uM = ~M
        string Q = P;
        for (char &ch : Q)
            ch = (ch == '0' ? '1' : '0');
        return Q;
    }
    else if (u == "00")
    {
        // 选第二行：全 0
        return string(P.size(), '0');
    }
    else if (u == "11")
    {
        // 两行都是 1：全 1
        return string(P.size(), '1');
    }
    // 理论上不会到这里
    return P;
}


static vector<BiDecompResult>
enumerate_one_case(const TT& in, int k1, int k2, int k3)
{
    vector<BiDecompResult> results;

    const string &f01 = in.f01;
    if (f01.empty()) return results;

    int n = (int)std::log2((double)f01.size());
    if ((int)in.order.size() != n) return results;
    if (k1 + k2 + k3 != n) return results;

    if (k1 == k3)
    {
        int gamma_first = in.order[0];           
        int lambda_first = in.order[k1 + k2];    
        
        if (gamma_first > lambda_first)
        {
            return results;
        }
    }

    int R = 1 << k1;
    int C = 1 << k2;
    int B = 1 << k3;

    auto is_const_block = [&](const string& b){
        if (b.empty()) return false;
        for (char c : b)
            if (c != b[0]) return false;
        return true;
    };

    auto is_complement = [&](const string& a, const string& b){
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i)
            if (a[i] == b[i]) return false;
        return true;
    };

    // 1) 构造块矩阵 Mf: blk[r][c]
    vector<vector<string>> blk(R, vector<string>(C));
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < C; ++c)
        {
            string s;
            for (int l = 0; l < B; ++l)
            {
                int idx = (r << (k2+k3)) | (c << k3) | l;
                s.push_back(f01[idx]);
            }
            blk[r][c] = s;
        }

    // 2) 每列做"模式分类"
    struct ColType
    {
        bool ok;
        vector<string> const_blocks;
        vector<string> nonconst_blocks;
    };

    vector<ColType> cols(C);

    for (int c = 0; c < C; ++c)
    {
        ColType ct{true, {}, {}};
        
        auto push_unique = [&](vector<string>& v, const string& s){
            for (auto &x : v) if (x == s) return;
            v.push_back(s);
        };

        for (int r = 0; r < R; ++r)
        {
            const string &b = blk[r][c];
            if (is_const_block(b)) 
                push_unique(ct.const_blocks, b);
            else                   
                push_unique(ct.nonconst_blocks, b);
        }

        int nc = (int)ct.nonconst_blocks.size();
        int cc = (int)ct.const_blocks.size();

        // 合法情况判断
        if (nc == 0 && cc == 1)
        {
            // 纯常数列（单一常数）✓
        }
        else if (nc == 0 && cc == 2)
        {
            // 纯常数列（两种常数）✓
        }
        else if (nc == 1 && cc == 0)
        {
            // 纯非常数列 ✓
        }
        else if (nc == 1 && cc == 1)
        {
            // 一种常数 + 一种非常数 ✓
        }
        else if (nc == 2 && cc == 0 && is_complement(ct.nonconst_blocks[0], ct.nonconst_blocks[1]))
        {
            // 互补非常数对 ✓
        }
        else
        {
            // 其他情况不合法
            ct.ok = false;
        }

        if (!ct.ok) return results;
        cols[c] = ct;
    }

    // 3) 收集必需的 u 类型（优先处理强约束列）
    set<string> required_u;

    // 第一遍：处理强约束列（互补对、混合列、纯非常数列）
    for (int c = 0; c < C; ++c)
    {
        const auto &ct = cols[c];
        int nc = (int)ct.nonconst_blocks.size();
        int cc = (int)ct.const_blocks.size();
        
        if (nc == 2 && cc == 0)
        {
            // 互补对：强制需要 10 和 01
            required_u.insert("10");
            required_u.insert("01");
        }
        else if (nc == 1 && cc == 1)
        {
            // 混合列：需要 非常数u + 对应常数u
            required_u.insert("10");  // 默认用10处理非常数块
            
            string const_val = ct.const_blocks[0];
            if (const_val == string(B, '1'))
                required_u.insert("11");
            else if (const_val == string(B, '0'))
                required_u.insert("00");
        }
        else if (nc == 1 && cc == 0)
        {
            // 纯非常数列
            required_u.insert("10");
        }
    }

    // 第二遍：处理纯常数列
    for (int c = 0; c < C; ++c)
    {
        const auto &ct = cols[c];
        int nc = (int)ct.nonconst_blocks.size();
        int cc = (int)ct.const_blocks.size();
        
        if (nc == 0 && cc == 1)
        {
            // 单一常数块的纯常数列
            // 任意2种u都可以覆盖，如果已经有2种u就不需要添加
            if (required_u.size() < 2)
            {
                string const_val = ct.const_blocks[0];
                if (const_val == string(B, '1'))
                    required_u.insert("11");
                else
                    required_u.insert("00");
            }
        }
        else if (nc == 0 && cc == 2)
        {
            // 两种常数块的纯常数列
            // 如果已经有2种u，可以尝试用现有u覆盖
            // 否则需要添加常数u
            if (required_u.size() < 2)
            {
                // 需要补充u
                for (const string& b : ct.const_blocks)
                {
                    if (b == string(B, '0'))
                        required_u.insert("00");
                    else if (b == string(B, '1'))
                        required_u.insert("11");
                    
                    if (required_u.size() >= 2) break;
                }
            }
            // 如果已经有2种u，纯常数列可以被任意2种u覆盖
        }
    }

    vector<string> u_types(required_u.begin(), required_u.end());
    
    // 4) 检查是否可分解
    if (u_types.size() > 2)
    {
        std::cout << "  ⚠️  需要 " << u_types.size() << " 种 u，不可分解（";
        for (const string &u : u_types) std::cout << u << " ";
        std::cout << "）\n";
        return results;
    }
    
    if (u_types.empty())
    {
        // 不应该发生
        return results;
    }

    // 5) 确定全局 u1 和 u2
    string global_u1, global_u2;
    if (u_types.size() == 2)
    {
        global_u1 = u_types[0];
        global_u2 = u_types[1];
    }
    else if (u_types.size() == 1)
    {
        global_u1 = u_types[0];
        global_u2 = global_u1;  // 只有1种u，重复使用
    }
    else
    {
        return results;
    }

    // 6) 构造 F
    string F01;
    if (u_types.size() == 2)
    {
        F01 = global_u1 + global_u2;
    }
    else
    {
        F01 = global_u1 + global_u1;
    }

     // 7+8) 利用全局 u1, u2，同步构造 φ-table 和 ψ-table
    vector<int> phi_bits(R * C, 0);
    vector<int> psi_bits(C * B, 0);

    auto blocks_equal = [&](const string& a, const string& b) {
        return a == b;
    };

    for (int c = 0; c < C; ++c)
    {
        // 收集这一列所有不同的块 W_r
        vector<string> distinct_blocks;
        auto add_unique = [&](const string& s) {
            if (std::find(distinct_blocks.begin(), distinct_blocks.end(), s) 
                == distinct_blocks.end())
            {
                distinct_blocks.push_back(s);
            }
        };

        for (int r = 0; r < R; ++r)
            add_unique(blk[r][c]);

        // 为这一列构造候选的 ψ 列向量 P
        vector<string> candidate_P;
        auto add_candidate = [&](const string& s) {
            if (std::find(candidate_P.begin(), candidate_P.end(), s) 
                == candidate_P.end())
            {
                candidate_P.push_back(s);
            }
        };

        // 把出现过的块加进候选
        for (const string& w : distinct_blocks)
            add_candidate(w);

        // 加上全 0 / 全 1
        string all0(B, '0'), all1(B, '1');
        add_candidate(all0);
        add_candidate(all1);

        // 非常数块再加上它的补
        for (const string& w : distinct_blocks)
        {
            if (!is_const_block(w))
            {
                string comp = w;
                for (char &ch : comp)
                    ch = (ch == '0' ? '1' : '0');
                add_candidate(comp);
            }
        }

        bool foundP = false;
        string chosenP;

        // 尝试每一个候选 P，看是否能用 (u1, u2) 解释这一列的所有块
        for (const string& P : candidate_P)
        {
            bool ok = true;

            // 预先算好 u1P, u2P
            string u1P = apply_u_to_block(P, global_u1);
            string u2P = apply_u_to_block(P, global_u2);

            for (const string& W : distinct_blocks)
            {
                if (!blocks_equal(W, u1P) && !blocks_equal(W, u2P))
                {
                    ok = false;
                    break;
                }
            }

            if (ok)
            {
                foundP = true;
                chosenP = P;
                break;
            }
        }

        if (!foundP)
        {
            std::cout << "  ✘ 列 " << c << " 在给定全局 u 下无法找到合适的 ψ 和 φ\n";
            return results;   // 这个 (k1,k2,k3) / 变量划分失败
        }

        // 填 ψ：这一列对应的 ψ 列向量就是 chosenP
        for (int l = 0; l < B; ++l)
            psi_bits[c*B + l] = (chosenP[l] == '1') ? 1 : 0;

        // 填 φ：对于每一行，看它的块等于 u1P 还是 u2P 来决定 φ 位
        string u1P = apply_u_to_block(chosenP, global_u1);
        string u2P = apply_u_to_block(chosenP, global_u2);

        for (int r = 0; r < R; ++r)
        {
            const string& W = blk[r][c];

            if (blocks_equal(W, u1P))
            {
                // 这一行这一列用 u1
                phi_bits[r*C + c] = 1;   // 约定 1 -> 选 global_u1
            }
            else if (blocks_equal(W, u2P))
            {
                // 用 u2
                phi_bits[r*C + c] = 0;   // 0 -> 选 global_u2
            }
            else
            {
                // 理论上不该发生，因为上面已检查所有块 ∈ {u1P, u2P}
                std::cout << "  ✘ 行 " << r << " 列 " << c 
                          << " 无法匹配到 u1 或 u2\n";
                return results;
            }
        }
    }


    // 9) 变量集合
    vector<int> Gamma, Theta, Lambda;
    for (int i = 0; i < k1; ++i) Gamma.push_back(in.order[i]);
    for (int i = 0; i < k2; ++i) Theta.push_back(in.order[k1 + i]);
    for (int i = 0; i < k3; ++i) Lambda.push_back(in.order[k1 + k2 + i]);

    BiDecompResult Rst;
    Rst.k1 = k1;
    Rst.k2 = k2;
    Rst.k3 = k3;
    Rst.Gamma = Gamma;
    Rst.Theta = Theta;
    Rst.Lambda = Lambda;
    Rst.F01 = F01;

    Rst.phi_tt.f01.resize(R*C);
    for (int i = 0; i < R*C; ++i)
        Rst.phi_tt.f01[i] = phi_bits[i] ? '1' : '0';
    Rst.phi_tt.order = Gamma;
    Rst.phi_tt.order.insert(Rst.phi_tt.order.end(), Theta.begin(), Theta.end());

    Rst.psi_tt.f01.resize(C*B);
    for (int i = 0; i < C*B; ++i)
        Rst.psi_tt.f01[i] = psi_bits[i] ? '1' : '0';
    Rst.psi_tt.order = Theta;
    Rst.psi_tt.order.insert(Rst.psi_tt.order.end(), Lambda.begin(), Lambda.end());

    std::cout << "\n===== Matrix Form (Theorem 4.2) =====\n";
    std::cout << "k1=" << k1 << "  k2=" << k2 << "  k3=" << k3 << "\n";
    std::cout << "F = " << Rst.F01 << "\n";
    std::cout << "u_types: ";
    for (const string &u : u_types) std::cout << u << " ";
    std::cout << "\n";
    std::cout << "global_u1 = " << global_u1 << ", global_u2 = " << global_u2 << "\n\n";
    print_structure_matrix(k1, k2, k3, blk, phi_bits, psi_bits);

    results.push_back(Rst);
    return results;
}

//重排变量，所有解版本
static vector<BiDecompResult>
enumerate_bi_decomposition_all_permutations(const TT& in)
{
    vector<BiDecompResult> results;

    const string &f01 = in.f01;
    if (f01.empty()) return results;

    int n = (int)std::log2((double)f01.size());
    if ((int)in.order.size() != n) return results;

    // 枚举 k2 和 k3 的大小
    for (int k2 = 1; k2 <= n - 2; ++k2)
    {
        int max_k3 = (n - k2) / 2;

        for (int k3 = 1; k3 <= max_k3; ++k3)
        {
            int k1 = n - k2 - k3;
            if (k1 <= 0) continue;

            std::cout << "\n========== 枚举 k1=" << k1 << ", k2=" << k2 << ", k3=" << k3 << " ==========\n";

            // 枚举 Θ 的所有 C(n, k2) 种组合
            vector<bool> theta_mask(n, false);
            std::fill(theta_mask.begin(), theta_mask.begin() + k2, true);

            do {
                // 获取 Θ 的变量索引（1-based用于公式计算）
                vector<int> Theta_indices;
                for (int i = 0; i < n; ++i)
                    if (theta_mask[i])
                        Theta_indices.push_back(i + 1);

                // 剩余变量用于 Γ 和 Λ
                vector<int> remaining;
                for (int i = 0; i < n; ++i)
                    if (!theta_mask[i])
                        remaining.push_back(i + 1);

                // 枚举 Λ 的所有 C(n-k2, k3) 种组合
                vector<bool> lambda_mask(remaining.size(), false);
                std::fill(lambda_mask.begin(), lambda_mask.begin() + k3, true);

                do {
                    vector<int> Lambda_indices;
                    vector<int> Gamma_indices;

                    for (size_t i = 0; i < remaining.size(); ++i)
                    {
                        if (lambda_mask[i])
                            Lambda_indices.push_back(remaining[i]);
                        else
                            Gamma_indices.push_back(remaining[i]);
                    }

                    // 避免对称重复：当 k1 == k3 时，要求 Γ[0] < Λ[0]
                    if (k1 == k3 && Gamma_indices[0] > Lambda_indices[0])
                        continue;

                    // 打印当前尝试的分组
                    std::cout << "  尝试 Γ={";
                    for (int g : Gamma_indices) std::cout << g << " ";
                    std::cout << "}, Θ={";
                    for (int t : Theta_indices) std::cout << t << " ";
                    std::cout << "}, Λ={";
                    for (int l : Lambda_indices) std::cout << l << " ";
                    std::cout << "}\n";

                    // 构造重排后的真值表
                    string reordered_f01 = apply_variable_reordering(
                        f01, n, Gamma_indices, Theta_indices, Lambda_indices, k1, k2, k3);

                    // 构造重排后的 TT 对象
                    TT reordered_tt;
                    reordered_tt.f01 = reordered_f01;
                    
                    // 新的变量顺序：Γ, Θ, Λ（转回0-based）
                    reordered_tt.order.clear();
                    for (int idx : Gamma_indices) reordered_tt.order.push_back(idx );
                    for (int idx : Theta_indices) reordered_tt.order.push_back(idx );
                    for (int idx : Lambda_indices) reordered_tt.order.push_back(idx );

                    // 在重排后的真值表上尝试分解
                    auto sub = enumerate_one_case(reordered_tt, k1, k2, k3);
                    
                    if (!sub.empty())
                    {
                        std::cout << "    ✓ 找到解！\n";
                    }
                    
                    results.insert(results.end(), sub.begin(), sub.end());

                } while (std::prev_permutation(lambda_mask.begin(), lambda_mask.end()));

            } while (std::prev_permutation(theta_mask.begin(), theta_mask.end()));
        }
    }

    return results;
}

// =====================================================
// 顶层：枚举所有 k2,k3（不重排变量）不重排变量版本
//   k2 = 1..n-2
//   k3 = 1..(n-k2)/2   （保证 k3 <= k1）
//   一旦某个 k2 有解，就返回该 k2 的所有解
// =====================================================
static vector<BiDecompResult>
enumerate_bi_decomposition_no_reorder(const TT& in)
{
    vector<BiDecompResult> results;

    const string &f01 = in.f01;
    if (f01.empty()) return results;

    int n = (int)std::log2((double)f01.size());
    if ((int)in.order.size() != n) return results;

    // 不再在第一个 k2 有结果就 return，
    // 而是把所有 k2,k3 的分解都收集起来。
    for (int k2 = 1; k2 <= n - 2; ++k2)
    {
        int r = (n - k2) / 2;   // 最大 k3，保证 k3 <= k1

        for (int k3 = 1; k3 <= r; ++k3)
        {
            int k1 = n - k2 - k3;
            if (k1 <= 0) continue;

            auto sub = enumerate_one_case(in, k1, k2, k3);
            results.insert(results.end(), sub.begin(), sub.end());
        }
    }

    return results;
}

// =====================================================
// 边重排边找，找到第一个解就立即返回
// =====================================================
static bool
find_first_bi_decomposition(const TT& in, BiDecompResult& out)
{
    const string &f01 = in.f01;
    if (f01.empty()) return false;

    int n = (int)std::log2((double)f01.size());
    if ((int)in.order.size() != n) return false;

    // 枚举 k2 和 k3 的大小
    for (int k2 = 1; k2 <= n - 2; ++k2)
    {
        int max_k3 = (n - k2) / 2;

        for (int k3 = 1; k3 <= max_k3; ++k3)
        {
            int k1 = n - k2 - k3;
            if (k1 <= 0) continue;

            // 先试试不重排的情况
            auto sub = enumerate_one_case(in, k1, k2, k3);
            if (!sub.empty())
            {
                out = sub[0];
                return true;
            }

            // 枚举 Θ 的所有 C(n, k2) 种组合
            vector<bool> theta_mask(n, false);
            std::fill(theta_mask.begin(), theta_mask.begin() + k2, true);

            do {
                // 获取 Θ 的变量索引（1-based用于公式计算）
                vector<int> Theta_indices;
                for (int i = 0; i < n; ++i)
                    if (theta_mask[i])
                        Theta_indices.push_back(i + 1);

                // 剩余变量用于 Γ 和 Λ
                vector<int> remaining;
                for (int i = 0; i < n; ++i)
                    if (!theta_mask[i])
                        remaining.push_back(i + 1);

                // 枚举 Λ 的所有 C(n-k2, k3) 种组合
                vector<bool> lambda_mask(remaining.size(), false);
                std::fill(lambda_mask.begin(), lambda_mask.begin() + k3, true);

                do {
                    vector<int> Lambda_indices;
                    vector<int> Gamma_indices;

                    for (size_t i = 0; i < remaining.size(); ++i)
                    {
                        if (lambda_mask[i])
                            Lambda_indices.push_back(remaining[i]);
                        else
                            Gamma_indices.push_back(remaining[i]);
                    }

                    // 避免对称重复：当 k1 == k3 时，要求 Γ[0] < Λ[0]
                    if (k1 == k3 && Gamma_indices[0] > Lambda_indices[0])
                        continue;

                    // 构造重排后的真值表
                    string reordered_f01 = apply_variable_reordering(
                        f01, n, Gamma_indices, Theta_indices, Lambda_indices, k1, k2, k3);

                    // 构造重排后的 TT 对象
                    TT reordered_tt;
                    reordered_tt.f01 = reordered_f01;
                    
                    // 新的变量顺序：Γ, Θ, Λ（转回0-based）
                    reordered_tt.order.clear();
                    for (int idx : Gamma_indices) reordered_tt.order.push_back(idx);
                    for (int idx : Theta_indices) reordered_tt.order.push_back(idx);
                    for (int idx : Lambda_indices) reordered_tt.order.push_back(idx);

                    // 在重排后的真值表上尝试分解
                    sub = enumerate_one_case(reordered_tt, k1, k2, k3);
                    
                    // 只要找到一个解，立即返回
                    if (!sub.empty())
                    {
                        out = sub[0];
                        return true;
                    }

                } while (std::prev_permutation(lambda_mask.begin(), lambda_mask.end()));

            } while (std::prev_permutation(theta_mask.begin(), theta_mask.end()));
        }
    }

    return false;
}

// =====================================================
// 递归双分解（参照 DSD 的编号和递归方式）
// =====================================================

// 递归双分解主函数
static int bi_decomp_recursive(const TT& f, int depth = 0)
{
    int len = f.f01.size();
    
    // 🔥 基本情况：2输入或更少，直接建小树
    if (len <= 4)
        return build_small_tree(f);
    
    // 🔥 尝试找到第一个双分解
    BiDecompResult result;
    bool found = find_first_bi_decomposition(f, result);
    
    if (!found)
    {
        // 找不到双分解，回退到建小树
        cout << "⚠️  深度 " << depth << "：无法双分解，回退到直接建树\n";
        return build_small_tree(f);
    }
    
    // 🔥 找到了双分解，打印信息
    cout << "\n" << string(depth*2, ' ') << "📌 深度 " << depth << " 双分解成功：\n";
    cout << string(depth*2, ' ') << "   k1=" << result.k1 
         << "  k2=" << result.k2 << "  k3=" << result.k3 << "\n";
    
    cout << string(depth*2, ' ') << "   Γ = { ";
    for (int v : result.Gamma) cout << v << " ";
    cout << "}\n";
    
    cout << string(depth*2, ' ') << "   Θ = { ";
    for (int v : result.Theta) cout << v << " ";
    cout << "}\n";
    
    cout << string(depth*2, ' ') << "   Λ = { ";
    for (int v : result.Lambda) cout << v << " ";
    cout << "}\n";
    
    cout << string(depth*2, ' ') << "   F(u,v) = " << result.F01 << "\n";
    
    // 🔥 记录变量到 FINAL_VAR_ORDER（像 DSD 那样）
    for (int v : result.Gamma)
    {
        if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), v) 
            == FINAL_VAR_ORDER.end())
        {
            FINAL_VAR_ORDER.push_back(v);
        }
    }
    for (int v : result.Theta)
    {
        if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), v) 
            == FINAL_VAR_ORDER.end())
        {
            FINAL_VAR_ORDER.push_back(v);
        }
    }
    for (int v : result.Lambda)
    {
        if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), v) 
            == FINAL_VAR_ORDER.end())
        {
            FINAL_VAR_ORDER.push_back(v);
        }
    }
    
    // 🔥 准备 φ 和 ψ 的递归
    TT phi_tt = result.phi_tt;
    TT psi_tt = result.psi_tt;
    
    int n_phi = phi_tt.order.size();
    int n_psi = psi_tt.order.size();
    
    cout << string(depth*2, ' ') << "📌 递归分解 φ：原始变量 { ";
    for (int v : phi_tt.order) cout << v << " ";
    cout << "} → 局部编号 { ";
    for (int i = 1; i <= n_phi; i++) cout << i << " ";
    cout << "}\n";
    cout << string(depth*2, ' ') << "   映射关系：";
    for (int i = 0; i < n_phi; i++)
        cout << "位置" << (i+1) << "→变量" << phi_tt.order[i] << " ";
    cout << "\n";
    
    cout << string(depth*2, ' ') << "📌 递归分解 ψ：原始变量 { ";
    for (int v : psi_tt.order) cout << v << " ";
    cout << "} → 局部编号 { ";
    for (int i = 1; i <= n_psi; i++) cout << i << " ";
    cout << "}\n";
    cout << string(depth*2, ' ') << "   映射关系：";
    for (int i = 0; i < n_psi; i++)
        cout << "位置" << (i+1) << "→变量" << psi_tt.order[i] << " ";
    cout << "\n\n";
    
    // 🔥 递归分解 φ 和 ψ
    int L = bi_decomp_recursive(phi_tt, depth + 1);
    int R = bi_decomp_recursive(psi_tt, depth + 1);
    
    // 🔥 创建当前节点（用 F 作为函数）
    return new_node(result.F01, {L, R});
}

// =====================================================
// 顶层调用入口（参照 run_dsd_recursive）
// =====================================================
inline bool run_bi_decomp_recursive(const std::string& binary01)
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

    for (int i = 0; i < n; ++i)
        root.order[i] = i + 1;  // 位置 (i+1) 对应变量 (i+1)

    std::cout << "======= 双分解递归开始 =======\n";
    std::cout << "输入 = " << binary01 << " (n=" << n << ")\n";
    std::cout << "初始映射：";
    for (int i = 0; i < n; i++)
        std::cout << "位置" << (i+1) << "→变量" << root.order[i] << " ";
    std::cout << "\n\n";

    NODE_LIST.clear();
    NODE_ID = 1;
    STEP_ID = 1;
    FINAL_VAR_ORDER.clear();

    // 🔥 可选：先缩减到 support
    TT root_shrunk = shrink_to_support(root);
    int root_id = bi_decomp_recursive(root_shrunk, 0);
    
    // 🔥 打印最终节点列表
    std::cout << "\n===== 最终双分解节点列表 =====\n";
    for (auto& nd : NODE_LIST)
    {
        std::cout << nd.id << " = " << nd.func;

        if (nd.func == "in")
        {
            // 输入节点：显示原始变量编号
            std::cout << "(var=" << nd.var_id << ")";
        }
        else if (!nd.child.empty())
        {
            // 打印所有子节点
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
