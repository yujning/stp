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
// 变量编号 0=a,1=b,... 直接映射为字母
// ========================================================
static std::string varname_from_id(int v)
{
    return std::string(1, char('a' + v));
}

// ========================================================
// 二进制转 hex
// ========================================================
static std::string bin_to_hex(const std::string& bin)
{
    int v = std::stoi(bin, nullptr, 2);
    char buf[20];
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

        if (FINAL_VAR_ORDER.empty())
        {
            std::cout << "❌ FINAL_VAR_ORDER empty! (DSD not run?)\n";
            return;
        }

        int nvars = FINAL_VAR_ORDER.size();

        // ====================================================
        // 1) 输出全部 INPUT() —— 即使某些变量未在 NODE_LIST 内出现
        // ====================================================
        for (int v = 0; v < nvars; v++)
        {
            fout << "INPUT(" << varname_from_id(v) << ")\n";
        }

        // ====================================================
        // 2) 输出 F0
        // ====================================================
        fout << "OUTPUT(F0)\n\n";

        // ====================================================
        // 3) 节点命名
        // ====================================================
        std::map<int, std::string> name_of;

        // 输入节点：根据 var_id 映射字母
        for (auto &n : NODE_LIST)
        {
            if (n.func == "in")
                name_of[n.id] = varname_from_id(n.var_id);
        }

        // 非输入：new_nX
        for (auto &n : NODE_LIST)
        {
            if (n.func != "in")
                name_of[n.id] = "new_n" + std::to_string(n.id);
        }

        // Root 改名 F0
        int root_id = NODE_LIST.back().id;
        name_of[root_id] = "F0";

        // ====================================================
        // 4) 输出 LUT: 二输入交换 child[1], child[0]
        // ====================================================
        for (auto &n : NODE_LIST)
        {
            if (n.func == "in") continue;

            fout << name_of[n.id] << " = LUT 0x"
                 << bin_to_hex(n.func) << " (";

            int k = n.child.size();

            if (k == 2)
            {
                // ⭐ 两输入：交换顺序
                fout << name_of[n.child[1]] << ", "
                     << name_of[n.child[0]];
            }
            else
            {
                // 单输入或其他：保持原顺序
                for (int i = 0; i < k; i++)
                {
                    fout << name_of[n.child[i]];
                    if (i + 1 < k) fout << ", ";
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
