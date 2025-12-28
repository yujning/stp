#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <iostream>

// =====================================================
// 工具
// =====================================================
inline uint64_t pow2(int k)
{
    return 1ull << k;
}

// MF 列索引：A|B|C（MSB->LSB）
inline uint64_t mf_index(uint64_t a, uint64_t b, uint64_t c, int y, int z)
{
    return (a << (y + z)) | (b << z) | c;
}

// =====================================================
// 只做一件事：
// 根据 MF 和 (x,y,z)，通过
//   MZ = I_{2^x} ⊗ Mr_{2^y}
// 反推出部分确定的 MXY，其余列填 'x'
// =====================================================
inline std::string compute_MXYX_from_MF(
    const std::string& MF, // 二进制真值表，长度 2^(x+y+z)
    int x,
    int y,
    int z)
{
    const uint64_t HN = pow2(x);
    const uint64_t MN = pow2(y);
    const uint64_t LN = pow2(z);

    const uint64_t t = MN;

    const uint64_t MXY_cols = pow2(x + 2*y + z);
    const uint64_t block_size = LN;

    // 初始化为未知
    std::string MXYX(MXY_cols, 'x');

    // 展开 MZ = I ⊗ Mr 的选择规则
    for (uint64_t a = 0; a < HN; ++a)
    {
        for (uint64_t b = 0; b < MN; ++b)
        {
            // Mr,t = δ_{t^2}[ (i-1)t+i ]，0-based
            uint64_t r = b * t + b;

            // 在 MXY 中的块号
            uint64_t block_id = a * t * t + r;

            for (uint64_t c = 0; c < LN; ++c)
            {
                uint64_t mf_col  = mf_index(a, b, c, y, z);
                uint64_t mxy_col = block_id * block_size + c;

                if (mxy_col < MXY_cols)
                    MXYX[mxy_col] = MF[mf_col];
            }
        }
    }

    return MXYX;
}

inline void print_MZ_delta(int x, int y)
{
    uint64_t Ix = 1ull << x;
    uint64_t t  = 1ull << y;

    uint64_t rows = 1ull << (x + 2*y);
    uint64_t cols = 1ull << (x + y);

    std::cout << "🟩 MZ = I_{2^" << x << "} ⊗ Mr_{2^" << y << "}\n";
    std::cout << "   size = " << rows << " x " << cols << "\n";
    std::cout << "   δ_" << rows << "[ ";

    for (uint64_t a = 0; a < Ix; ++a)
    {
        for (uint64_t b = 0; b < t; ++b)
        {
            // Mr,t 的 (i−1)t+i，0-based
            uint64_t r = b * t + b;

            // 在 I ⊗ Mr 中的行号（0-based）
            uint64_t row = a * t * t + r;

            // δ 是 1-based
            std::cout << (row + 1) << " ";
        }
    }

    std::cout << "]\n";
}


// =====================================================
// bd -s 调用的入口函数（只打印 MXY）
// =====================================================
inline bool run_strong_bi_dec_and_build_dag(const std::string& MF)
{
    // -------- 自动推变量数 --------
    int n = 0;
    while ((1u << n) < MF.size()) ++n;
    if ((1u << n) != MF.size()) return false;

    if (n < 7 || n > 10)
    {
        std::cout << "⚠️ strong bi-dec only supports 7~10 vars\n";
        return false;
    }

    // --------------------------------------------------
    // 这里用一个固定策略（你之后可以枚举）
    // 例：优先 y=1，其次 x 尽量大
    // --------------------------------------------------
    int x = n - 2;
    int y = 1;
    int z = n - x - y;

    if (x < 0 || z < 0) return false;
    if (x + y > 6 || y + z + 1 > 6) return false;

    std::cout << "🔀 Try strong bi-decomposition (shared vars)...\n";
    std::cout << "📐 Partition: x=" << x
              << " y=" << y
              << " z=" << z << "\n";

    print_MZ_delta(x, y);


    // -------- 计算 MXYX --------
    std::string MXYX = compute_MXYX_from_MF(MF, x, y, z);

    std::cout << "🟨 MXY = " << MXYX << "\n";

    return true;
}
