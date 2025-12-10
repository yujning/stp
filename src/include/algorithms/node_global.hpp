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

