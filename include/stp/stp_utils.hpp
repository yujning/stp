#include <iostream>
#include <string>
#include <sstream>
#include <vector>

namespace stp
{
  std::vector<std::string> parse_tokens( const std::string &input, const std::string &prefix ) 
  {
    std::istringstream iss(input);
    std::string token;
    std::vector<std::string> result;

    while (iss >> token) 
    {
      if (token.compare(0, prefix.size(), prefix) == 0) 
      {
        result.push_back(token);
      }
    }

    return result;
  }
}

