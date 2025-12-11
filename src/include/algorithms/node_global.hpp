#pragma once
#include <vector>
#include <map>
#include <tuple>
#include <string>

// 前向声明，避免循环依赖
struct DSDNode;

// C++17 允许 inline 全局变量，只需定义一次即可，全工程自动共享
inline int NODE_ID = 1;
inline int STEP_ID = 1;
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
