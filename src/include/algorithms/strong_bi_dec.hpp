#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <cstdint>

// =====================================================
// 工具
// =====================================================
inline uint64_t pow2(int k) { return 1ull << k; }

// A|B|C → MF index (MSB→LSB)
inline uint64_t mf_index(uint64_t a, uint64_t b, uint64_t c, int y, int z)
{
    return (a << (y + z)) | (b << z) | c;
}

// =====================================================
// 分区枚举（与你描述完全一致）
// 前两个 ≤6，后两个 ≤5
// =====================================================
inline std::vector<std::tuple<int,int,int>>
enumerate_partitions(int n)
{
    std::vector<std::tuple<int,int,int>> ps;
    for (int y = 1; y <= 4; ++y)
    for (int x = 0; x <= 6; ++x)
    {
        int z = n - x - y;
        if (z < 0) continue;
        if (x + y > 6) continue;
        if (y + z + 1 > 6) continue;
        ps.emplace_back(x,y,z);
    }

    std::sort(ps.begin(), ps.end(),
        [](auto&a,auto&b){
            int xa,ya,za; std::tie(xa,ya,za)=a;
            int xb,yb,zb; std::tie(xb,yb,zb)=b;
            if (ya!=yb) return ya<yb;
            return xa>xb;
        });
    ps.erase(std::unique(ps.begin(), ps.end()), ps.end());
    return ps;
}

// =====================================================
// MZ = I ⊗ Mr 的 δ 打印
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
// Step 1：由 MF + MZ 反推 MXY（未知填 'x'）
// =====================================================
inline std::string compute_MXYX_from_MF(
    const std::string& MF,
    int x,int y,int z)
{
    uint64_t HN = pow2(x);
    uint64_t MN = pow2(y);
    uint64_t LN = pow2(z);

    uint64_t block_size = pow2(y+z);          // ★ B,C
    uint64_t MXY_cols   = pow2(x+2*y+z);

    std::string MXY(MXY_cols,'x');

    for (uint64_t a=0;a<HN;a++)
    for (uint64_t b=0;b<MN;b++)
    {
        uint64_t r = b*MN + b;                // Mr 行
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
// Step 2：strong_dsd 风格求 MX / MY（允许 x）
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
            {
                if (blk[j]!='x' && blocks[k][j]!='x'
                    && blk[j]!=blocks[k][j]) { ok=false; break; }
            }
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
// bd -s 入口
// =====================================================
inline bool run_strong_bi_dec_and_build_dag(const std::string& MF)
{
    int n=0;
    while ((1u<<n)<MF.size()) n++;
    if ((1u<<n)!=MF.size()) return false;

    std::cout<<"🔀 Try strong bi-decomposition (shared vars)...\n";

    for (auto [x,y,z]: enumerate_partitions(n))
    {
        std::cout<<"📐 Try "<<x<<"+"<<y<<"+"<<z<<"\n";
        print_MZ_delta(x,y);

        std::string MXY = compute_MXYX_from_MF(MF,x,y,z);
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

        
        return true;
    }


    std::cout<<"❌ No valid strong bi-decomposition\n";
    return false;
}
