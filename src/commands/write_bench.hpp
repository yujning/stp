#ifndef WRITE_BENCH_HPP
#define WRITE_BENCH_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <alice/alice.hpp>

extern std::vector<DSDNode> NODE_LIST;
extern inline std::vector<int> FINAL_VAR_ORDER;

namespace alice
{

// ========================================================
// 变量编号反向映射：max→a, 次 max→b, ..., 1→最后一个字母
//   例如变量数=7： 7→a, 6→b, 5→c, 4→d, 3→e, 2→f, 1→g
// ========================================================
static std::string varname_from_id(int v)
{
    int max_id = FINAL_VAR_ORDER.size();
    char name = 'a' + (max_id - v);  // 反向映射
    return std::string(1, name);
}

// ========================================================
// 二进制转 hex
// ========================================================
static std::string bin_to_hex(const std::string& bin)
{
    int v = std::stoi(bin, nullptr, 2);
    char buf[50];
    sprintf(buf, "%x", v);
    return std::string(buf);
}

// ========================================================
// write_bench
// ========================================================
class write_bench_command : public command
{
public:
    explicit write_bench_command(const environment::ptr& env)
        : command(env, "Write DSD result as LUT BENCH file")
    {
        add_option("file", filename,
            "output benchmark filename")->required();
        add_option("--out, -o", filename,
            "output benchmark filename");
    }

protected:
    void execute() override
    {
        if (filename.empty())
        {
            std::cout << "❌ No output file provided\n";
            return;
        }

        std::ofstream fout(filename);
        if (!fout)
        {
            std::cout << "❌ Cannot open " << filename << "\n";
            return;
        }

        if (NODE_LIST.empty())
        {
            std::cout << "❌ NODE_LIST is empty! (DSD not run?)\n";
            return;
        }

        // ====================================================
        // 1) 输入节点：按 DSDNode.var_id 输出
        // ====================================================
        for (auto &n : NODE_LIST)
        {
            if (n.func == "in")
            {
                std::string name = varname_from_id(n.var_id);
                fout << "INPUT(" << name << ")\n";
            }
        }

        // ====================================================
        // 2) 输出节点 F0
        // ====================================================
        fout << "OUTPUT(F0)\n\n";

        // ====================================================
        // 3) 节点命名：输入用字母，非输入用 new_nX
        // ====================================================
        std::map<int, std::string> name_of;

        // 输入节点：变量名
        for (auto &n : NODE_LIST)
        {
            if (n.func == "in")
                name_of[n.id] = varname_from_id(n.var_id);
        }

        // 内部节点：new_nX
        for (auto &n : NODE_LIST)
        {
            if (n.func != "in")
                name_of[n.id] = "new_n" + std::to_string(n.id);
        }

        // Root 重命名为 F0
        int root_id = NODE_LIST.back().id;
        name_of[root_id] = "F0";

        // ====================================================
        // 4) 输出 LUT，其中 child 左右互换（只对 2 输入）
        // ====================================================
        for (auto &n : NODE_LIST)
        {
            if (n.func == "in") continue;

            fout << name_of[n.id] << " = LUT 0x"
                 << bin_to_hex(n.func) << " (";

            int cnt = n.child.size();

            if (cnt == 2)
            {
                // ⭐ 左右互换
                fout << name_of[n.child[1]] << ", "
                     << name_of[n.child[0]];
            }
            else
            {
                // 单输入或多输入：顺序不动
                for (int i = 0; i < cnt; i++)
                {
                    fout << name_of[n.child[i]];
                    if (i + 1 < cnt) fout << ", ";
                }
            }

            fout << ")\n";
        }

        std::cout << "✅ BENCH written to " << filename << "\n";
    }

private:
    std::string filename{};
};

ALICE_ADD_COMMAND(write_bench, "STP")

} // namespace alice

#endif
