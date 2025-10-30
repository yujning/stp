#ifndef DSD_HPP
#define DSD_HPP

#include <iostream>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <alice/alice.hpp>

// 引入算法文件
#include "../include/algorithms/stp_dsd.hpp"

namespace alice
{
    class dsd_command : public command
    {
    public:
        explicit dsd_command(const environment::ptr &env)
            : command(env, "Run STP decomposition on a hexadecimal input")
        {
            add_option("-f, --factor", hex_input, "hexadecimal number (must produce 2^n-length binary)");
        }

    protected:
        void execute() override
        {
            if (hex_input.empty())
            {
                std::cout << "Please specify a hexadecimal number after -f (e.g. -f 0x10)" << std::endl;
                return;
            }

            // 1️⃣ 将输入的16进制转为十进制
            unsigned long value = 0;
            try
            {
                value = std::stoul(hex_input, nullptr, 16);
            }
            catch (...)
            {
                std::cout << "Invalid hexadecimal input: " << hex_input << std::endl;
                return;
            }

            if (value == 0)
            {
                std::cout << "0 is not valid input." << std::endl;
                return;
            }

            // 2️⃣ 转为二进制字符串
            std::string binary;
            {
                unsigned long temp = value;
                while (temp > 0)
                {
                    binary.insert(binary.begin(), (temp & 1) ? '1' : '0');
                    temp >>= 1;
                }
            }

            // 3️⃣ 检查长度是否为2的幂
            size_t len = binary.size();
            bool is_pow2 = (len & (len - 1)) == 0;
            if (!is_pow2)
            {
                std::cout << "❌ Error: binary length (" << len << ") is not a power of 2." << std::endl;
                std::cout << "   Please input a hex value whose binary length = 2^n (e.g. 0x8 -> 1000, len=4)." << std::endl;
                return;
            }

            // 若二进制长度小于2的幂，补齐到完整位宽（例如 3→4）
            size_t next_pow2 = 1;
            while (next_pow2 < len) next_pow2 <<= 1;
            if (next_pow2 != len)
            {
                binary = std::string(next_pow2 - len, '0') + binary;
            }

            std::cout << "✅ Hex " << hex_input << " => binary " << binary 
                      << " (length = " << binary.size() << ")" << std::endl;

            // 4️⃣ 调用 STP 不相交分解算法
            all_reorders(binary);
        }

    private:
        std::string hex_input{};
    };

    ALICE_ADD_COMMAND(dsd, "STP")
}

#endif
