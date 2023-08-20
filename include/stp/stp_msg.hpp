#include <iostream>
#include <fmt/format.h>

namespace stp
{
  inline void MyStp::getVersion()
  {
    std::cout << fmt::format( "STP engine version is {}\n", 0.02 ) << std::endl;
  }
}
