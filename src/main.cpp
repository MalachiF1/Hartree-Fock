#include "Input.hpp"
#include "Molecule.hpp"
#include "SCF.hpp"

#include <iostream>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    std::pair<Molecule, SCFOptions> input = Input::read(argv[1]);
    SCF scf(input.first, input.second);
    scf.run();

    return 0;
}
