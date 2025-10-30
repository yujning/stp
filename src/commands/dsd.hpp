#ifndef DSD_HPP
#define DSD_HPP

#include <iostream>
#include <iomanip>
#include <cmath>
#include <alice/alice.hpp>

namespace alice
{
    // ./your_tool dsd -f 0x10
    class dsd_command : public command
    {
    public:
        explicit dsd_command(const environment::ptr &env)
            : command(env, "Check if a hex number is a power of two and show its info")
        {
            add_option("-f, --factor", hex_input, "hexadecimal number (must be power of 2)");
        }

    protected:
        void execute() override
        {
            if (hex_input.empty())
            {
                std::cout << "Please specify a hexadecimal number after -f (e.g. -f 0x10)" << std::endl;
                return;
            }

            // 支持 “0x” 或不带 “0x” 的16进制输入
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
                std::cout << "0 is not a valid power of two." << std::endl;
                return;
            }

            // 检查是否为2的次方
            if ((value & (value - 1)) != 0)
            {
                std::cout << "Error: " << hex_input << " (" << value << ") is NOT a power of 2." << std::endl;
                return;
            }

            // 计算指数
            int exponent = static_cast<int>(std::log2(value));

            std::cout << "✅ " << hex_input << " is a power of 2." << std::endl;
            std::cout << "Decimal value : " << value << std::endl;
            std::cout << "Log2(value)   : " << exponent << std::endl;
        }

    private:
        std::string hex_input{};
    };

    ALICE_ADD_COMMAND(dsd, "Utilities");
}

#endif
