#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include "stp_dsd.hpp"

struct StrongDsdResult {
    bool found = false;
    size_t L = 0;
    std::string Mx;
    std::string My;
};

inline std::string reorder_high_bit_first(const std::string& mf)
{
    const size_t n = mf.size();
    if (n == 0) {
        return mf;
    }

    size_t vars = static_cast<size_t>(std::log2(n));
    if ((1ull << vars) != n) {
        return mf;
    }

    std::string out(n, '0');
    for (size_t i = 0; i < n; ++i) {
        size_t rev = 0;
        for (size_t b = 0; b < vars; ++b) {
            rev = (rev << 1) | ((i >> b) & 1ull);
        }
        out[rev] = mf[i];
    }

    return out;
}

inline StrongDsdResult run_strong_dsd(const std::string& mf)
{
    StrongDsdResult result;
    const size_t n = mf.size();
    if (n == 0) {
        return result;
    }

    for (size_t L = 2; L <= n / 2; L <<= 1) {
        if (n % L != 0) {
            continue;
        }

        const size_t block_count = n / L;
        std::unordered_map<std::string, int> block_index;
        std::vector<std::string> ordered_blocks;
        ordered_blocks.reserve(2);
        std::string My;
        My.reserve(block_count);
        bool too_many_blocks = false;

        for (size_t i = 0; i < block_count; ++i) {
            std::string block = mf.substr(i * L, L);
            auto it = block_index.find(block);
            if (it == block_index.end()) {
                if (ordered_blocks.size() >= 2) {
                    too_many_blocks = true;
                    break;
                }
                int idx = static_cast<int>(ordered_blocks.size());
                block_index.emplace(block, idx);
                ordered_blocks.push_back(block);
                My.push_back(idx == 0 ? '1' : '0');
            } else {
                My.push_back(it->second == 0 ? '1' : '0');
            }
        }

        if (!too_many_blocks && ordered_blocks.size() <= 2) {
            result.found = true;
            result.L = L;
            if (ordered_blocks.size() == 1) {
                result.Mx = ordered_blocks[0] + ordered_blocks[0];
            } else {
                result.Mx = ordered_blocks[0] + ordered_blocks[1];
            }
            result.My = My;
            return result;
        }
    }

    return result;
}

inline std::vector<int> make_children_from_order_with_placeholder(
    const std::vector<int>& order,
    int placeholder_id)
{
    std::vector<int> children;
    children.reserve(order.size());
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        int var_id = *it;
        if (var_id == 0) {
            children.push_back(placeholder_id);
        } else {
            children.push_back(new_in_node(var_id));
        }
    }
    return children;
}

inline int build_strong_dsd_nodes(
    const std::string& mf,
    const std::vector<int>& order,
    int depth = 0)
{
    std::string indent(static_cast<size_t>(depth) * 2, ' ');
    std::cout << indent << "变量顺序：{ ";
    for (int v : order) {
        std::cout << v << " ";
    }
    std::cout << "}\n";

    if (mf.size() <= 4) {
        std::cout << indent << "⏹ Stop (size <= 4): " << mf << "\n";
        auto children = make_children_from_order_with_placeholder(order, 0);
        return new_node(mf, children);
    }

    auto result = run_strong_dsd(mf);
    if (!result.found) {
        std::cout << indent << "❌ Strong DSD: no valid L found for " << mf << "\n";
        auto children = make_children_from_order_with_placeholder(order, 0);
        return new_node(mf, children);
    }

    std::cout << indent << "✅ L = " << result.L << "\n";
    std::cout << indent << "Mx = " << result.Mx << "\n";
    std::cout << indent << "My = " << result.My << "\n";

    size_t vars_my = static_cast<size_t>(std::log2(result.My.size()));
    std::vector<int> order_my(order.begin(), order.begin() + vars_my);
    std::vector<int> order_rest(order.begin() + vars_my, order.end());
    std::reverse(order_my.begin(), order_my.end());
    std::reverse(order_rest.begin(), order_rest.end());

    std::cout << indent << "My 使用变量：{ ";
    for (int v : order_my) {
        std::cout << v << " ";
    }
    std::cout << "}\n";
    std::cout << indent << "Mx 使用变量：{ ";
    for (int v : order_rest) {
        std::cout << v << " ";
    }
    std::cout << "My }\n";

    int my_id = build_strong_dsd_nodes(result.My, order_my, depth + 1);

    std::vector<int> order_mx;
    order_mx.reserve(1 + order_rest.size());
    order_mx.insert(order_mx.end(), order_rest.begin(), order_rest.end());
    order_mx.push_back(0);

    auto children = make_children_from_order_with_placeholder(order_mx, my_id);
    return new_node(result.Mx, children);
}