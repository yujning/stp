#pragma once

// ⭐ 防止 std::is_trivial 污染 percy
#ifdef is_trivial
#undef is_trivial
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>


#include "node_global.hpp"

// 前向声明（定义在 stp_dsd.hpp 中）
struct TT;
static int build_small_tree(const TT& t);
static int dsd_factor(const TT& f, int depth);

inline int dsd_else_decompose(
    const TT& f,
    int depth)
{
  const uint32_t bits = static_cast<uint32_t>(f.f01.size());
  const uint32_t n = static_cast<uint32_t>(std::log2(bits));

  if ((1u << n) != bits)
    throw std::runtime_error("dsd_else_decompose: f01 length is not power-of-two");

  auto orig_children = make_children_from_order(f);

  std::cout << "⚠️ depth " << depth
            << ": fallback Shannon decomposition (n=" << n << ")\n";
  std::cout << "f=" << f.f01 << "\n";
  std::cout << "输入节点 = { ";
  for (size_t i = 0; i < orig_children.size(); ++i)
  {
    std::cout << orig_children[i];
    if (i + 1 < orig_children.size())
      std::cout << ",";
  }
  std::cout << " } 变量 = { ";
  for (size_t i = 0; i < f.order.size(); ++i)
  {
    std::cout << f.order[i];
    if (i + 1 < f.order.size())
      std::cout << ",";
  }
  std::cout << " }\n";

  const auto pivot_node = orig_children.back();
  const auto pivot_var = f.order.empty() ? -1 : f.order.front();

  const auto half_bits = bits / 2u;
  TT f_pos{ f.f01.substr(0, half_bits), {} };
  TT f_neg{ f.f01.substr(half_bits), {} };

  if (!f.order.empty())
  {
    f_pos.order.assign(f.order.begin() + 1, f.order.end());
    f_neg.order.assign(f.order.begin() + 1, f.order.end());
  }

  std::cout << "  split depth " << depth
            << " pivot var=" << pivot_var
            << " pos f=" << f_pos.f01 << "\n";
  std::cout << "  split depth " << depth
            << " pivot var=" << pivot_var
            << " neg f=" << f_neg.f01 << "\n" << std::flush;

  const auto pos_node = dsd_factor(f_pos, depth + 1);
  const auto neg_node = dsd_factor(f_neg, depth + 1);

  const auto pos_term = new_node("1000", { pivot_node, pos_node });
  const auto neg_term = new_node("0010", { pivot_node, neg_node });

  return new_node("1110", { pos_term, neg_term });
}