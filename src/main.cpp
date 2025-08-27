#include "Input.hpp"
#include "Molecule.hpp"
#include "Output.hpp"
#include "SCF.hpp"

#include <fstream>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        return 1;
    }

    std::string inputFile  = argv[1];
    std::string outputFile = argv[2];

    Input input(inputFile);
    std::shared_ptr<Output> output = std::make_shared<Output>(outputFile);

    // Write the user input to the output file.
    std::ifstream inputContents(inputFile);
    if (!inputContents.is_open())
    {
        throw std::runtime_error("Could not open input file: " + inputFile);
    }
    output->writeSeperator();
    output->write("User Input:\n");
    output->writeSeperator();
    std::string inputString((std::istreambuf_iterator<char>(inputContents)), std::istreambuf_iterator<char>());
    while (inputString.front() == '\n') inputString.erase(0, 1); // Remove leading newlines
    while (inputString.back() == '\n') inputString.pop_back(); // Remove trailing newlines
    output->write(inputString + "\n");
    inputContents.close();
    output->writeSeperator();
    output->write("\n");

    std::pair<Molecule, SCFOptions> inputData = input.read();
    SCF scf(inputData.first, inputData.second, output);
    scf.run();

    return 0;
}
