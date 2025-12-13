#pragma once

// ⭐ 防止 std::is_trivial 污染 percy
#ifdef is_trivial
#undef is_trivial
#endif

#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <stdexcept>
#include <string>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/print.hpp>
#include <kitty/constructors.hpp>

#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>

#include "stp_dsd.hpp"
#include "node_global.hpp"   // new_node / new_in_node

// 你已有的接口（外部实现）
int new_node(const std::string&, const std::vector<int>&);

// 如果你系统里有常量节点，建议你接上（可选）
// int const0_node();
// int const1_node();

// 如果你已经有 build_small_tree(f) 就保留，否则你自己实现一个 fallback
int build_small_tree(const TT& f);



inline int else_decompose(
    const TT& f,
    const std::vector<int>& orig_children,
    int depth
)
{
  // f.f01 是二进制字符串，长度应该是 2^n
  const uint32_t bits = static_cast<uint32_t>(f.f01.size());
  const uint32_t n = static_cast<uint32_t>(std::log2(bits));

  if ((1u << n) != bits)
    throw std::runtime_error("else_decompose: f01 length is not power-of-two");

  if (n > 4u)
  {
    std::cout << "⚠️ depth " << depth
              << ": exact 2-LUT only supports n ≤ 4 (got n=" << n << ")\n";
    return build_small_tree(f);
  }

  // 可选：检查字符串只含 0/1
  for (char c : f.f01)
    if (c != '0' && c != '1')
      throw std::runtime_error("else_decompose: f01 contains non-binary char");

  std::cout << "⚠️ depth " << depth
            << ": EXACT 2-LUT refine (n=" << n << ")\n";
  std::cout << "f=" << f.f01 << "\n";

  // 1) binary string -> kitty TT
  kitty::dynamic_truth_table tt(n);
  kitty::create_from_binary_string(tt, f.f01);
  std::cout << "[DEBUG] kitty hex = " << kitty::to_hex(tt) << "\n";

  // 2) build klut and create EXACTLY n PIs
  mockturtle::klut_network klut;
  std::vector<mockturtle::klut_network::signal> pis;
  pis.reserve(n);

  for (uint32_t i = 0; i < n; ++i)
    pis.push_back(klut.create_pi());

  // 防御：必须保证输入数量对
  if (pis.size() != n)
    throw std::runtime_error("else_decompose: PI creation failed");

  // 3) exact 2-LUT
  mockturtle::exact_resynthesis<mockturtle::klut_network> resyn(2);
  resyn(klut, tt, pis.begin(), pis.end(),
        [&](auto const& s) {
          klut.create_po(s);
        });

  std::cout << "Exact 2-LUT count = " << klut.num_gates() << "\n";

  // 这里先返回一个占位（如果你要接回自己的 node 图，再加翻译逻辑）
  // 最简单：把 PO 视为结果
  auto po_sig = klut.po_at(0);
  (void)po_sig;

  return 0;
}

