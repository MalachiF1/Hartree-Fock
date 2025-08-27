#pragma once
#include "Molecule.hpp"
#include "SCF.hpp"

#include <string>


class Input
{
  public:
    Input(const std::string& filename);
    std::pair<Molecule, SCFOptions> read();

  private:
    std::string filename;
    static Molecule readMoleculeBlock(const std::string& moleculeBlock, const std::string& basisBlock);
    // static std::string readBasisBlock(const std::string& basisBlock);
    static SCFOptions readSCFBlock(const std::string& SCFBlock);
};
