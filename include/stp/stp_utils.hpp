#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>

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

  void parse_tokens_with_details( const std::string &input, 
                                  const std::string &prefix, 
                                  std::vector<std::string> &result, 
                                  std::unordered_map<std::string, int> &tokenCounts, 
                                  std::unordered_map<std::string, int> &tokenIndices ) 
  {
    std::istringstream iss(input);
    std::string token;
    int currentIndex = 0;

    while (iss >> token) 
    {
      if (token.compare(0, prefix.size(), prefix) == 0) 
      {
        if (tokenCounts.find(token) == tokenCounts.end()) 
        {
          tokenCounts[token] = 0;
        } 
        else 
        {
          tokenCounts[token]++;
          token = token + "_" + std::to_string(tokenCounts[token]);
        }
        result.push_back(token);

        if (tokenIndices.find(token) == tokenIndices.end()) 
        {
          tokenIndices[token] = currentIndex;
        }
      }

      currentIndex++;
    }
  }

}

