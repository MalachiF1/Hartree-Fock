#pragma once
#include <string>

class Output
{
  public:
    Output(const std::string& filename);
    void write(const std::string& content);
    void writeSeperator(char sep = '-', size_t length = 99);

  private:
    std::string filename;
};
