#pragma once
#include "Molecule.hpp"
#include "SCF.hpp"

#include <string>


class Input
{
  public:
    static std::pair<Molecule, SCFOptions> read(const std::string& filename);
  private:
    static Molecule readMoleculeBlock(const std::string& moleculeBlock, const std::string& basisBlock);
    // static std::string readBasisBlock(const std::string& basisBlock);
    static SCFOptions readSCFBlock(const std::string& SCFBlock);
};
