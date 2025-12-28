#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <cstdint>
#include <algorithm>
#include <iostream>

#include "node_global.hpp"

// =====================================================
// Strong Bi-Dec result
// =====================================================
struct StrongBiDecResult
{
    bool found = false;
    int x = 0; // |A|  (MY-only, high)
    int y = 0; // |B|  (shared)
    int z = 0; // |C|  (MX-only, low)

    // TT encoding (same as your DSD/BD code):
    // MY: length 2^(x+y), output bit used as selector: '1'->choose block0, '0'->choose block1
    // MX: length 2^(1+y+z) = 2 * 2^(y+z), interpreted as block0 || block1
    std::string MY;
    std::string MX;

    // For debug:
    // MXMY_whole: length 2^(x+2y+z) in variable order [A | B_my | B_mx | C] (MSB->LSB)
    std::string MXMY_whole;
};

// =====================================================
// Helpers
// =====================================================
inline int tt_num_vars_from_size(size_t sz)
{
    int n = 0;
    size_t v = 1;
    while (v < sz) { v <<= 1; ++n; }
    return (v == sz) ? n : -1;
}

// Index in MSB->LSB order for a bit-vector of length n
inline uint64_t pack_bits_msb(const std::vector<int>& bits_msb) {
    uint64_t idx = 0;
    for (int b : bits_msb) idx = (idx << 1) | (uint64_t)(b & 1);
    return idx;
}

// MF order (no reorder): [A(x) | B(y) | C(z)] MSB->LSB
inline uint64_t mf_index_ABC(uint64_t a, uint64_t b, uint64_t c, int y, int z)
{
    return (a << (y + z)) | (b << z) | c;
}

// MXMY_whole order: [A(x) | B_my(y) | B_mx(y) | C(z)] MSB->LSB
inline uint64_t g_index_ABBC(uint64_t a, uint64_t bmy, uint64_t bmx, uint64_t c, int y, int z)
{
    // total bits = x + 2y + z, but x is not used in formula directly here
    // idx = a<<(2y+z) | bmy<<(y+z) | bmx<<z | c
    return (a << (2*y + z)) | (bmy << (y + z)) | (bmx << z) | c;
}

// Given MY index over [A|B_my] : idx = a<<y | bmy
inline uint64_t my_index_AB(uint64_t a, uint64_t bmy, int y)
{
    return (a << y) | bmy;
}

// Given MX index over [s | B_mx | C] : idx = s<<(y+z) | bmx<<z | c
inline uint64_t mx_index_sBC(uint64_t s, uint64_t bmx, uint64_t c, int y, int z)
{
    return (s << (y + z)) | (bmx << z) | c;
}

// =====================================================
// Partition enumeration (7~10 vars, 2 LUTs constraints)
//   MY uses A+B : x+y <= 6
//   MX uses sel + B + C : 1+y+z <= 6
// =====================================================
inline std::vector<std::tuple<int,int,int>> enumerate_partitions_7_to_10(int n)
{
    std::vector<std::tuple<int,int,int>> parts;
    if (n < 7 || n > 10) return parts;

    for (int y = 1; y <= 4; ++y)
    for (int x = 5; x >= 0; --x)
    {
        int z = n - x - y;
        if (z < 0) continue;
        if (x + y > 6) continue;
        if (1 + y + z > 6) continue;
        parts.emplace_back(x,y,z);
    }

    std::sort(parts.begin(), parts.end(),
        [](auto const& a, auto const& b){
            int xa,ya,za; std::tie(xa,ya,za)=a;
            int xb,yb,zb; std::tie(xb,yb,zb)=b;
            if (ya != yb) return ya < yb;
            return xa > xb;
        });

    parts.erase(std::unique(parts.begin(), parts.end()), parts.end());
    return parts;
}

// =====================================================
// Build Mr (diagonal selection) in delta form
//   t = 2^y
//   Mr = δ_{t^2}[ (i-1)*t + i ] for i=1..t
// Build MZ = I_{2^x} ⊗ Mr in delta form indices (1-based for printing)
// =====================================================
inline void build_Mr_MZ_delta_indices(int x, int y,
                                     std::vector<int>& mr_idx_1based,
                                     std::vector<int>& mz_idx_1based)
{
    const int t = 1 << y;       // t = 2^y
    const int Ix = 1 << x;      // 2^x

    mr_idx_1based.clear();
    mr_idx_1based.reserve(t);
    for (int i = 1; i <= t; ++i)
        mr_idx_1based.push_back((i - 1) * t + i);  // diagonal positions, 1-based

    // MZ = I_{Ix} ⊗ Mr
    mz_idx_1based.clear();
    mz_idx_1based.reserve((size_t)Ix * (size_t)t);
    const int block = t * t; // Mr rows
    for (int k = 0; k < Ix; ++k) {
        int base = k * block;
        for (int i = 1; i <= t; ++i) {
            mz_idx_1based.push_back(base + (i - 1) * t + i);
        }
    }
}

// =====================================================
// Core solver (YOUR model):
// 1) Treat unknown whole matrix G = MXMY as function over [A|B_my|B_mx|C]
//    Size: 2^(x+2y+z)
// 2) MZ selects only diagonal where B_my == B_mx to produce MF
//    => We can fill G on those diagonal assignments using MF.
// 3) On this expanded G, do strong_dsd-like split:
//    My vars = [A|B_my] (x+y bits)
//    Mx vars = [B_mx|C] (y+z bits)
//    For each My assignment, extract block over Mx-vars; allow only 2 unique blocks.
//    MY is selector of which block, MX = block0||block1
// 4) Rebuild a full G' from (MY,MX) and verify that diagonal extraction equals MF.
// =====================================================
inline bool solve_strong_bi_fixed_stp(
    const std::string& mf, // MF truth table string in MSB->LSB order [A|B|C]
    int x, int y, int z,
    StrongBiDecResult& out)
{
    const int n = x + y + z;
    if ((int)mf.size() != (1 << n)) return false;

    const uint64_t AN = 1ull << x;
    const uint64_t BN = 1ull << y;
    const uint64_t CN = 1ull << z;

    const int nG = x + 2*y + z;
    const uint64_t GN = 1ull << nG;          // size of whole G = MXMY
    const uint64_t MYN = 1ull << (x + y);    // MY length
    const uint64_t BLK = 1ull << (y + z);    // block length over [B_mx|C]
    const uint64_t MXN = 1ull << (1 + y + z);// MX length = 2*BLK

    // -------------------------------------------------
    // Step 1: fill diagonal of G using MF, others = '0'
    // diagonal means: B_my == B_mx
    // -------------------------------------------------
    std::string G(GN, '0');

    for (uint64_t a = 0; a < AN; ++a)
    for (uint64_t b = 0; b < BN; ++b)
    for (uint64_t c = 0; c < CN; ++c)
    {
        uint64_t idxF = mf_index_ABC(a, b, c, y, z);
        char fv = mf[idxF];

        uint64_t idxG = g_index_ABBC(a, b, b, c, y, z);
        G[idxG] = fv;
    }

    // -------------------------------------------------
    // Step 2: strong_dsd-like extraction on expanded G
    // My assignment enumerates [A|B_my] in MSB->LSB:
    //   my = a<<y | bmy
    // For each my, the block is over [B_mx|C] (MSB->LSB):
    //   mx = bmx<<z | c
    // -------------------------------------------------
    std::unordered_map<std::string, int> block_id;
    std::vector<std::string> blocks;
    blocks.reserve(2);

    std::string MY(MYN, '0');

    for (uint64_t my = 0; my < MYN; ++my)
    {
        uint64_t a  = (y == 0) ? my : (my >> y);
        uint64_t bmy = (y == 0) ? 0  : (my & (BN - 1ull));

        std::string blk(BLK, '0');

        // enumerate [B_mx|C] as MSB->LSB
        for (uint64_t mx = 0; mx < BLK; ++mx)
        {
            uint64_t bmx = (z == 0) ? mx : (mx >> z);
            uint64_t c   = (z == 0) ? 0  : (mx & (CN - 1ull));

            uint64_t idxG = g_index_ABBC(a, bmy, bmx, c, y, z);
            blk[mx] = G[idxG];
        }

        auto it = block_id.find(blk);
        if (it == block_id.end())
        {
            if (blocks.size() >= 2) {
                return false; // more than 2 blocks -> not 2-LUT strong split under this partition
            }
            int id = (int)blocks.size();
            block_id.emplace(blk, id);
            blocks.push_back(blk);
            it = block_id.find(blk);
        }

        // selector encoding: 1 -> block0, 0 -> block1
        MY[my] = (it->second == 0) ? '1' : '0';
    }

    // Need non-degenerate: must really have 2 distinct blocks
    if (blocks.size() != 2) return false;

    std::string MX = blocks[0] + blocks[1]; // length 2*BLK == 2^(1+y+z)

    // -------------------------------------------------
    // Step 3: rebuild a full G' from (MY,MX) for debug/verification
    // G'(a,bmy,bmx,c) = MX( sel, bmx, c ) where sel = MY(a,bmy)
    // -------------------------------------------------
    std::string G2(GN, '0');
    for (uint64_t a = 0; a < AN; ++a)
    for (uint64_t bmy = 0; bmy < BN; ++bmy)
    {
        uint64_t my = my_index_AB(a, bmy, y);
        uint64_t sel = (MY[my] == '1') ? 0ull : 1ull;

        for (uint64_t bmx = 0; bmx < BN; ++bmx)
        for (uint64_t c = 0; c < CN; ++c)
        {
            uint64_t idxMX = mx_index_sBC(sel, bmx, c, y, z);
            char v = MX[idxMX];

            uint64_t idxG = g_index_ABBC(a, bmy, bmx, c, y, z);
            G2[idxG] = v;
        }
    }

    // -------------------------------------------------
    // Step 4: apply diagonal selection (equiv to MZ) and verify equals MF
    // Extract only where bmy==bmx:
    // MF'(a,b,c) = G2(a,b,b,c)
    // -------------------------------------------------
    for (uint64_t a = 0; a < AN; ++a)
    for (uint64_t b = 0; b < BN; ++b)
    for (uint64_t c = 0; c < CN; ++c)
    {
        uint64_t idxF = mf_index_ABC(a,b,c,y,z);
        uint64_t idxG = g_index_ABBC(a,b,b,c,y,z);
        if (G2[idxG] != mf[idxF]) return false;
    }

    // success
    out.found = true;
    out.x = x; out.y = y; out.z = z;
    out.MY = std::move(MY);
    out.MX = std::move(MX);
    out.MXMY_whole = std::move(G2);
    return true;
}

// =====================================================
// Main entry: run strong bi-dec and build DAG
// =====================================================
inline bool run_strong_bi_dec_and_build_dag(const std::string& mf)
{
    const int n = tt_num_vars_from_size(mf.size());
    if (n < 7 || n > 10) return false;

    ORIGINAL_VAR_COUNT = n;

    // var IDs in your system: MSB->LSB = n,n-1,...,1
    std::vector<int> vars_msb2lsb;
    vars_msb2lsb.reserve(n);
    for (int v = n; v >= 1; --v) vars_msb2lsb.push_back(v);

    FINAL_VAR_ORDER = vars_msb2lsb;

    StrongBiDecResult res;
    res.found = false;

    auto parts = enumerate_partitions_7_to_10(n);
    for (auto const& p : parts)
    {
        int x,y,z; std::tie(x,y,z) = p;
        if (solve_strong_bi_fixed_stp(mf, x, y, z, res)) break;
    }
    if (!res.found) return false;

    // ---------------- prints ----------------
    std::cout << "✅ Strong Bi-Decomposition found\n";
    std::cout << "📐 Partition: " << res.x << "+" << res.y << "+" << res.z << "\n";

    std::cout << "  H (A, size " << res.x << "): { ";
    for (int i = 0; i < res.x; ++i) std::cout << vars_msb2lsb[i] << " ";
    std::cout << "}\n";

    std::cout << "  M (B shared, size " << res.y << "): { ";
    for (int i = res.x; i < res.x + res.y; ++i) std::cout << vars_msb2lsb[i] << " ";
    std::cout << "}\n";

    std::cout << "  L (C, size " << res.z << "): { ";
    for (int i = res.x + res.y; i < n; ++i) std::cout << vars_msb2lsb[i] << " ";
    std::cout << "}\n";

    std::cout << "🟦 MY (len=2^(" << (res.x + res.y) << ")) = " << res.MY << "\n";
    std::cout << "🟥 MX (len=2^(1+" << (res.y + res.z) << ")) = " << res.MX << "\n";
    std::cout << "🟨 MXMY_whole (vars=[A|B_my|B_mx|C], len=2^(" << (res.x + 2*res.y + res.z) << "))\n";
    std::cout << "🟨 MXMY = " << res.MXMY_whole << "\n";

    // Mr / MZ delta print
    {
        std::vector<int> mr_idx, mz_idx;
        build_Mr_MZ_delta_indices(res.x, res.y, mr_idx, mz_idx);

        const int t = 1 << res.y;
        const int Ix = 1 << res.x;

        std::cout << "🟩 Mr = δ_" << (t*t) << "[ ";
        for (int v : mr_idx) std::cout << v << " ";
        std::cout << "]\n";

        std::cout << "🟩 MZ = I_{2^" << res.x << "} ⊗ Mr = δ_" << (Ix * t * t) << "[ ";
        for (int v : mz_idx) std::cout << v << " ";
        std::cout << "]\n";
    }

    // ---------------- build DAG (MY -> MX) ----------------
    // MY inputs: [A|B] = first (x+y) vars in MSB->LSB
    std::vector<int> order_my;
    order_my.reserve(res.x + res.y);
    for (int i = 0; i < res.x + res.y; ++i) order_my.push_back(vars_msb2lsb[i]);

    int my_node = new_node(
        res.MY,
        make_children_from_order_with_placeholder(order_my, nullptr, nullptr)
    );

    // MX inputs: [sel | B | C]  (B and C are original vars, shared B used here too)
    const int placeholder_id = n + 1;
    std::unordered_map<int,int> ph;
    ph[placeholder_id] = my_node;

    std::vector<int> order_mx;
    order_mx.reserve(1 + res.y + res.z);
    order_mx.push_back(placeholder_id);

    // B shared vars
    for (int i = res.x; i < res.x + res.y; ++i) order_mx.push_back(vars_msb2lsb[i]);
    // C vars
    for (int i = res.x + res.y; i < n; ++i) order_mx.push_back(vars_msb2lsb[i]);

    new_node(
        res.MX,
        make_children_from_order_with_placeholder(order_mx, &ph, nullptr)
    );

    return true;
}
