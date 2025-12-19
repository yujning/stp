#pragma once

#include <string>
#include <sstream>
#include <stdexcept>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/print.hpp>

namespace alice
{

inline std::string bin_to_hex(const std::string& bin)
{
    std::string b = bin;

    while (b.size() % 4 != 0)
        b = "0" + b;

    static const char* hex_map = "0123456789abcdef";
    std::string hex;
    hex.reserve(b.size() / 4);

    for (size_t i = 0; i < b.size(); i += 4)
    {
        int v = (b[i] - '0') * 8 +
                (b[i+1] - '0') * 4 +
                (b[i+2] - '0') * 2 +
                (b[i+3] - '0');
        hex.push_back(hex_map[v]);
    }

    while (hex.size() > 1 && hex[0] == '0')
        hex.erase(hex.begin());

    return hex;
}

inline std::string hex_to_binary(const std::string& hex)
{
    std::string clean = hex;
    if (clean.rfind("0x", 0) == 0 || clean.rfind("0X", 0) == 0)
        clean = clean.substr(2);

    if (clean.empty())
        throw std::runtime_error("empty hex string");

    unsigned bit_count = clean.size() * 4;
    unsigned num_vars = 0;
    while ((1u << num_vars) < bit_count)
        num_vars++;

    if ((1u << num_vars) != bit_count)
        throw std::runtime_error("hex length is not 2^n bits");

    kitty::dynamic_truth_table tt(num_vars);
    kitty::create_from_hex_string(tt, clean);

    std::ostringstream oss;
    kitty::print_binary(tt, oss);
    return oss.str();
}

} // namespace alice
