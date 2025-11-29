#ifndef WRITE_BENCH_HPP
#define WRITE_BENCH_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <alice/alice.hpp>

// 来自 DSD 的节点 & 变量顺序
extern std::vector<DSDNode> NODE_LIST;
extern inline std::vector<int> FINAL_VAR_ORDER;
extern inline int ORIGINAL_VAR_COUNT;

namespace alice
{

// ========================================================
// 映射变量编号 → 字母：最大编号 → 'a'
// 变量1是最高位，所以编号越大，字母越靠前
// 例：n=4 时，4→'a', 3→'b', 2→'c', 1→'d'
// ========================================================
static std::string varname_from_id(int v)
{
    int maxv = ORIGINAL_VAR_COUNT;
    int offset = maxv - v;      // v=max → offset=0 → 'a'
    return std::string(1, char('a' + offset));
}


// ========================================================
// 二进制转 hex
// ========================================================
static std::string bin_to_hex(const std::string& bin)
{
    int v = std::stoi(bin, nullptr, 2);
    char buf[32];
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

        if (ORIGINAL_VAR_COUNT == 0)
        {
            std::cout << "❌ ORIGINAL_VAR_COUNT is 0! (DSD not run?)\n";
            return;
        }

        // ====================================================
        // 1) 输入变量：输出所有变量 1 到 n
        //    按从大到小输出（最大编号 → 'a' 先输出）
        // ====================================================
        for (int v = ORIGINAL_VAR_COUNT; v >= 1; v--)
            fout << "INPUT(" << varname_from_id(v) << ")\n";

        fout << "OUTPUT(F0)\n\n";

        // ====================================================
        // 2) 节点命名
        // ====================================================
        std::map<int,std::string> name_of;

        for (auto &n : NODE_LIST)
        {
            if (n.func == "in")
                name_of[n.id] = varname_from_id(n.var_id);
        }

        for (auto &n : NODE_LIST)
        {
            if (n.func != "in")
                name_of[n.id] = "new_n" + std::to_string(n.id);
        }

        int root_id = NODE_LIST.back().id;
        name_of[root_id] = "F0";

        // ====================================================
        // 3) 输出 LUT
        //    双输入：左右互换
        // ====================================================
        for (auto &n : NODE_LIST)
        {
            if (n.func == "in") continue;

            fout << name_of[n.id] << " = LUT 0x"
                 << bin_to_hex(n.func) << " (";

            int cnt = n.child.size();

            if (cnt == 2)
            {
                fout << name_of[n.child[1]] << ", "
                     << name_of[n.child[0]];
            }
            else
            {
                for (int i = 0; i < cnt; i++)
                {
                    fout << name_of[n.child[i]];
                    if (i + 1 < cnt) fout << ", ";
                }
            }

            fout << ")\n";
        }

        std::cout << "✅ BENCH written to " << filename << "\n";
        
        // 🔥 打印变量映射
        std::cout << "📋 变量映射（最高位→'a'）：\n";
        for (int v = ORIGINAL_VAR_COUNT; v >= 1; v--)
            std::cout << "   变量" << v << " → '" << varname_from_id(v) << "'\n";
    }

private:
    std::string filename{};
};

ALICE_ADD_COMMAND(write_bench, "STP")

} // namespace alice

#endif
