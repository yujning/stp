#pragma once
#include <vector>
#include <map>
#include <tuple>
#include <string>
#include "bi_decomposition.hpp"
#include "bi_dec_else_dec.hpp"
#include "stp_dsd.hpp"
// 前向声明，避免循环依赖
struct DSDNode;

// C++17 允许 inline 全局变量，只需定义一次即可，全工程自动共享
inline int NODE_ID = 1;
inline int STEP_ID = 1;
inline bool ENABLE_ELSE_DEC = false;

inline int ORIGINAL_VAR_COUNT = 0;

inline std::vector<int> FINAL_VAR_ORDER;

inline std::vector<DSDNode> NODE_LIST;

inline std::map<int, int> INPUT_NODE_CACHE;

inline std::map<std::tuple<std::string, std::vector<int>>, int> NODE_HASH;

// ======================================================
// 🔥 重置所有全局节点状态（每次运行 bd/dsd 前必须调用）
// ======================================================
inline void RESET_NODE_GLOBAL()
{
    NODE_ID = 1;
    STEP_ID = 1;
    ORIGINAL_VAR_COUNT = 0;

    NODE_LIST.clear();
    FINAL_VAR_ORDER.clear();

    INPUT_NODE_CACHE.clear();
    NODE_HASH.clear();
}
static inline std::vector<int> make_children_from_order(const TT& t)
{
  std::vector<int> ch;
  ch.reserve(t.order.size());

  for (int var_id : t.order)
  {
    // 复用缓存输入节点（不会重复 new）
    ch.push_back(new_in_node(var_id));

    // 顺便维护 FINAL_VAR_ORDER（如果你需要）
    if (std::find(FINAL_VAR_ORDER.begin(), FINAL_VAR_ORDER.end(), var_id) == FINAL_VAR_ORDER.end())
      FINAL_VAR_ORDER.push_back(var_id);
  }
  return ch;
}
