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

namespace alice
{
    // a,b,c,...,z,aa,ab...
    static std::string index_to_name(int index)
    {
        std::string name;
        index--;
        while (index >= 0)
        {
            name.push_back('a' + (index % 26));
            index = index / 26 - 1;
        }
        std::reverse(name.begin(), name.end());
        return name;
    }

    static std::string bin_to_hex(const std::string& bin)
    {
        int v = std::stoi(bin, nullptr, 2);
        char buf[10];
        sprintf(buf, "%x", v);
        return std::string(buf);
    }

    class write_bench_command : public command
    {
    public:
        explicit write_bench_command(const environment::ptr& env)
            : command(env, "Write DSD result as LUT BENCH file")
        {
            // ⭐ 支持写法：write_bench 11.bench
            add_option("file", filename,
                       "output benchmark filename")
                ->required();

            // ⭐ 同时支持写法：write_bench -o 11.bench
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

            // --- 节点名映射 ---
            std::map<int,std::string> name_of;
            int in_count = 0;

            for (auto &n : NODE_LIST)
            {
                if (n.func == "in")
                    name_of[n.id] = index_to_name(++in_count);
                else
                    name_of[n.id] = "new_n" + std::to_string(n.id);
            }

            // --- OUTPUT (force the output to be F0) ---
            int root = NODE_LIST.back().id;

            // override root name
            name_of[root] = "F0";

            // print output section
            fout << "OUTPUT(F0)\n\n";


            // --- INPUT ---
            for (auto &n : NODE_LIST)
                if (n.func == "in")
                    fout << "INPUT(" << name_of[n.id] << ")\n";

            fout << "OUTPUT(" << name_of[root] << ")\n\n";

            // --- LUTs ---
            for (auto &n : NODE_LIST)
            {
                if (n.func == "in") continue;

                fout << name_of[n.id] << " = LUT 0x"
                     << bin_to_hex(n.func) << " (";

                // for (size_t i=0; i<n.child.size(); i++)
                // {
                //     fout << name_of[n.child[i]];
                //     if (i+1 < n.child.size()) fout << ", ";
                // }
                // --- print children in reverse order ---
                for (int i = (int)n.child.size() - 1; i >= 0; i--)
                {
                    fout << name_of[n.child[i]];
                    if (i > 0) fout << ", ";
                }

                fout << ")\n";
            }

            std::cout << "✅ BENCH written to " << filename << "\n";
        }

    private:
        std::string filename{};
    };

    ALICE_ADD_COMMAND(write_bench, "STP")
}

#endif
