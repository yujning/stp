#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include <numeric>

#include "node_global.hpp"

// =====================================================
// 工具
// =====================================================
inline uint64_t pow2(int k) { return 1ull << k; }

// =====================================================
// 变量重排（编码映射）
// new_order: MSB -> LSB（1-based 原变量编号）
// 原始 TT 顺序假定：n,n-1,...,1
// =====================================================
inline std::string reorder_tt_by_var_order(
    const std::string& mf,
    int n,
    const std::vector<int>& new_order_msb2lsb)
{
    const size_t N = mf.size();
    std::string out(N, '0');

    for (size_t new_idx = 0; new_idx < N; ++new_idx)
    {
        uint64_t old_idx = 0;
        for (int i = 0; i < n; ++i)
        {
            int bit = (new_idx >> (n - 1 - i)) & 1;
            int var = new_order_msb2lsb[i];   // 1-based
            int pos = var - 1;                // 原 MSB 位置
            old_idx |= (uint64_t(bit) << (n - 1 - pos));
        }
        out[new_idx] = mf[old_idx];
    }
    return out;
}

// =====================================================
// 枚举 (x,y,z)，按 y 从小到大
// 约束：x+y ≤ 6, y+z+1 ≤ 6, x+y+z = n
// =====================================================
inline std::vector<std::tuple<int,int,int>>
enumerate_xyz(int n)
{
    std::vector<std::tuple<int,int,int>> out;
    for (int y = 1; y <= 4; ++y)
    for (int x = 0; x <= 6; ++x)
    {
        int z = n - x - y;
        if (z < 0) continue;
        if (x + y > 6) continue;
        if (y + z + 1 > 6) continue;
        out.emplace_back(x,y,z);
    }
    return out;
}

// =====================================================
// 枚举变量划分 (A,B,C)，保持 MSB->LSB
// =====================================================
inline std::vector<
    std::tuple<std::vector<int>,std::vector<int>,std::vector<int>>>
enumerate_variable_partitions(int n, int x, int y, int z)
{
    std::vector<
        std::tuple<std::vector<int>,std::vector<int>,std::vector<int>>> parts;

    std::vector<int> vars(n);
    std::iota(vars.begin(), vars.end(), 1);

    std::vector<bool> sel_b(n,false);
    std::fill(sel_b.begin(), sel_b.begin()+y, true);

    do {
        std::vector<int> B, rest;
        for (int i=0;i<n;i++)
            (sel_b[i]?B:rest).push_back(vars[i]);

        if (x == 0) {
            parts.emplace_back(std::vector<int>{}, B, rest);
            continue;
        }

        std::vector<bool> sel_a(rest.size(), false);
        std::fill(sel_a.begin(), sel_a.begin()+x, true);

        do {
            std::vector<int> A,C;
            for (size_t i=0;i<rest.size();i++)
                (sel_a[i]?A:C).push_back(rest[i]);
            parts.emplace_back(A,B,C);
        } while (std::prev_permutation(sel_a.begin(), sel_a.end()));

    } while (std::prev_permutation(sel_b.begin(), sel_b.end()));

    return parts;
}

// =====================================================
// A|B|C → MF index (MSB→LSB)
// =====================================================
inline uint64_t mf_index(uint64_t a, uint64_t b, uint64_t c, int y, int z)
{
    return (a << (y + z)) | (b << z) | c;
}

// =====================================================
// 打印 Mr / MZ
// =====================================================
inline void print_MZ_delta(int x, int y)
{
    uint64_t Ix = pow2(x);
    uint64_t t  = pow2(y);

    std::cout << "🟩 MZ = I_{2^" << x << "} ⊗ Mr_{2^" << y << "}\n";
    std::cout << "🟩 Mr = δ_" << (t*t) << "[ ";
    for (uint64_t i=0;i<t;i++)
        std::cout << (i*t+i+1) << " ";
    std::cout << "]\n";

    std::cout << "🟩 MZ = δ_" << (Ix*t*t) << "[ ";
    for (uint64_t a=0;a<Ix;a++)
        for (uint64_t i=0;i<t;i++)
            std::cout << (a*t*t + i*t + i + 1) << " ";
    std::cout << "]\n";
}

// =====================================================
// Step 1：MF′ + MZ → MXY
// =====================================================
inline std::string compute_MXYX_from_MF(
    const std::string& MF,
    int x,int y,int z)
{
    uint64_t HN = pow2(x);
    uint64_t MN = pow2(y);
    uint64_t LN = pow2(z);

    uint64_t MXY_cols = pow2(x+2*y+z);
    std::string MXY(MXY_cols,'x');

    for (uint64_t a=0;a<HN;a++)
    for (uint64_t b=0;b<MN;b++)
    {
        uint64_t r = b*MN + b;
        uint64_t block = a*MN*MN + r;

        for (uint64_t c=0;c<LN;c++)
        {
            uint64_t mf_col  = mf_index(a,b,c,y,z);
            uint64_t mxy_col = block*LN + c;
            MXY[mxy_col] = MF[mf_col];
        }
    }
    return MXY;
}

// =====================================================
// Step 2：求 MX / MY
// =====================================================
inline bool solve_MX_MY_from_MXY(
    const std::string& MXY,
    int x,int y,int z,
    std::string& MX,
    std::string& MY)
{
    uint64_t MYN = pow2(x+y);
    uint64_t block_size = pow2(y+z);

    std::vector<std::string> blocks;
    MY.resize(MYN);

    for (uint64_t i=0;i<MYN;i++)
    {
        std::string blk = MXY.substr(i*block_size, block_size);

        bool hit=false;
        for (size_t k=0;k<blocks.size();k++)
        {
            bool ok=true;
            for (uint64_t j=0;j<block_size;j++)
                if (blk[j]!='x' && blocks[k][j]!='x'
                    && blk[j]!=blocks[k][j]) { ok=false; break; }

            if (ok)
            {
                for (uint64_t j=0;j<block_size;j++)
                    if (blocks[k][j]=='x') blocks[k][j]=blk[j];
                MY[i]=(k==0?'1':'0');
                hit=true; break;
            }
        }

        if (!hit)
        {
            if (blocks.size()==2) return false;
            blocks.push_back(blk);
            MY[i]=(blocks.size()==1?'1':'0');
        }
    }

    if (blocks.size()!=2) return false;

    MX.clear();
    for (auto& b:blocks)
        for (char c:b)
            MX.push_back(c=='x'?'0':c);

    return true;
}

// =====================================================
// children 构造（含占位符）
// =====================================================
inline std::vector<int> make_children_with_placeholder(
    const std::vector<int>& order_msb2lsb,
    int placeholder_var_id,
    int placeholder_node_id)
{
    std::vector<int> ch;
    ch.reserve(order_msb2lsb.size());

    for (int var_id : order_msb2lsb)
    {
        if (var_id == placeholder_var_id) continue;
        new_in_node(var_id);
        if (std::find(FINAL_VAR_ORDER.begin(),
                      FINAL_VAR_ORDER.end(),
                      var_id) == FINAL_VAR_ORDER.end())
            FINAL_VAR_ORDER.push_back(var_id);
    }

    for (auto it = order_msb2lsb.rbegin();
         it != order_msb2lsb.rend(); ++it)
    {
        int var_id = *it;
        if (var_id == placeholder_var_id)
            ch.push_back(placeholder_node_id);
        else
            ch.push_back(new_in_node(var_id));
    }

    return ch;
}

// =====================================================
// ★ 入口函数（最终）
// =====================================================
inline bool run_strong_bi_dec_and_build_dag(const std::string& MF)
{
    int n=0;
    while ((1u<<n)<MF.size()) n++;
    if ((1u<<n)!=MF.size()) return false;

    RESET_NODE_GLOBAL();

    for (auto [x,y,z] : enumerate_xyz(n))
    {
        auto parts = enumerate_variable_partitions(n,x,y,z);
        for (const auto& [A,B,C] : parts)
        {
            std::vector<int> new_order;
            new_order.insert(new_order.end(), A.begin(), A.end());
            new_order.insert(new_order.end(), B.begin(), B.end());
            new_order.insert(new_order.end(), C.begin(), C.end());

            std::string MFp =
                reorder_tt_by_var_order(MF,n,new_order);

            std::cout<<"📐 Try "<<x<<"+"<<y<<"+"<<z
                     <<"  A={ "; for(int v:A) std::cout<<v<<" ";
            std::cout<<"} B={ "; for(int v:B) std::cout<<v<<" ";
            std::cout<<"} C={ "; for(int v:C) std::cout<<v<<" ";
            std::cout<<"}\n";

            print_MZ_delta(x,y);

            std::string MXY =
                compute_MXYX_from_MF(MFp,x,y,z);
            std::cout<<"🟨 MXY = "<<MXY<<"\n";

            std::string MX,MY;
            if (!solve_MX_MY_from_MXY(MXY,x,y,z,MX,MY))
            {
                std::cout<<"❌ MX/MY unsat\n";
                continue;
            }

            std::cout<<"✅ Strong Bi-Decomposition found\n";
            std::cout<<"🟦 MY = "<<MY<<"\n";
            std::cout<<"🟥 MX = "<<MX<<"\n";

            // ===== 构造 DAG =====
            ORIGINAL_VAR_COUNT = n;

            TT tt_my;
            tt_my.f01 = MY;
            tt_my.order.insert(tt_my.order.end(), A.begin(), A.end());
            tt_my.order.insert(tt_my.order.end(), B.begin(), B.end());

            auto children_my = make_children_from_order(tt_my);
            std::reverse(children_my.begin(), children_my.end());
            int my_node = new_node(MY, children_my);

            int placeholder_var_id = n + 1;
            std::vector<int> order_mx;
            order_mx.push_back(placeholder_var_id);
            order_mx.insert(order_mx.end(), B.begin(), B.end());
            order_mx.insert(order_mx.end(), C.begin(), C.end());

            auto children_mx = make_children_with_placeholder(
                order_mx, placeholder_var_id, my_node);
            std::reverse(children_mx.begin(), children_mx.end());

            new_node(MX, children_mx);
            return true;
        }
    }

    std::cout<<"❌ No valid strong bi-decomposition\n";
    return false;
}
