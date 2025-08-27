#include "Output.hpp"

#include <fstream>
#include <sstream>

Output::Output(const std::string& filename) : filename(filename)
{
    std::ofstream outputFile(this->filename, std::ios::trunc);
    if (!outputFile.is_open())
    {
        throw std::runtime_error("Could not open output file: " + filename);
    }
    outputFile.close();
}

void Output::write(const std::string& content)
{
    std::ofstream outputFile(this->filename, std::ios::app);
    if (!outputFile.is_open())
    {
        throw std::runtime_error("Could not open output file: " + filename);
    }
    outputFile << content;
    outputFile.close();
}

void Output::writeSeperator(char sep, size_t length)
{
    std::stringstream ss;
    for (size_t i = 0; i < length; ++i) { ss << sep; }
    ss << "\n";
    this->write(ss.str());
}
