#pragma once
#include <bits/stdc++.h>
using namespace std;


#include <chrono>
#include <fstream>



#include "excute.hpp"




//-----------------------------------------
// 生成单位矩阵表示（用你的容器格式）
//-----------------------------------------
vector<stp_data> identity_vec(int dim)
{
    vector<stp_data> I(dim * dim + 1);
    I[0] = dim;
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            I[1 + i * dim + j] = (i == j) ? 1 : 0;
    return I;
}

//-----------------------------------------
// 将输入二进制串转为列向量矩阵格式
//-----------------------------------------
vector<stp_data> binary_to_vec(const string &bin)
{
    int rows = bin.size(); // 列向量
    vector<stp_data> M(rows + 1);
    M[0] = rows; // 行数
    for (int i = 0; i < rows; ++i)
        M[i + 1] = bin[i] - '0';
    return M;
}

//-----------------------------------------
// 打印矩阵(调试)
//-----------------------------------------
void print_vec(const vector<stp_data> &M)
{
    cout << "[";
    for (size_t i = 0; i < M.size(); ++i)
    {
        cout << M[i];
        if (i + 1 < M.size()) cout << " ";
    }
    cout << "]\n";
}

//-----------------------------------------
// 主过程：根据定理3.8生成所有重排
//-----------------------------------------
void all_reorders(const string &binary)
{
    int len = binary.size();
    int n = log2(len);  // 变量数
    int r = n / 2;      // floor(n/2)

    cout << "Input binary = " << binary << " (n=" << n << ")\n";

    vector<stp_data> Mf = binary_to_vec(binary);

    // 遍历 s = 1..r
    for (int s = 1; s <= r; ++s)
    {
        vector<int> vars(n);
        iota(vars.begin(), vars.end(), 1);

        vector<int> comb(s);
        vector<bool> v(n);
        fill(v.begin(), v.begin() + s, true);
        cout << "\n===== s=" << s << " =====\n";

        // 遍历所有Λ组合
        do {
            vector<int> Lambda;
            for (int i = 0; i < n; ++i)
                if (v[i]) Lambda.push_back(i + 1);

            vector<vector<stp_data>> swap_chain;

            // 构造 ∏_{k=s→1} W[2, 2^(j_k+(s-1)-k)]
            for (int k = s; k >= 1; --k)
            {
                int j_k = Lambda[k - 1];
                int exp = j_k + (s - 1) - k;
                auto Wmat = generate_swap_vec(2, pow(2, exp));
                swap_chain.push_back(Wmat);
            }
            // 最后乘 W[2^(n-s), 2^s]
            swap_chain.push_back(generate_swap_vec(pow(2, n - s), pow(2, s)));

            // 链乘得到最终重排矩阵
            vector<stp_data> Mperm = Vec_chain_multiply(swap_chain, false);

            // 应用重排矩阵
            vector<stp_data> result = Vec_semi_tensor_product(Mf, Mperm);

            cout << "Λ = { ";
            for (int j : Lambda) cout << j << " ";
            cout << "}  => reordered binary: ";
            for (size_t i = 1; i < result.size(); ++i)
                cout << result[i];
            cout << "\n";

        } while (prev_permutation(v.begin(), v.end()));
    }
}
