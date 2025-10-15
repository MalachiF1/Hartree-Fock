#pragma once
#include "Molecule.hpp"
#include "SCF.hpp"

#include <array>
#include <string>
#include <vector>

class Input
{
  public:
    Input(const std::string& filename);

    struct InputSettings
    {
        Molecule molecule;
        SCFOptions scfOptions;
        BasisSet basis;
    };
    InputSettings read();

    const std::vector<std::string>& getWarnings() const;

  private:
    std::string filename;
    std::vector<std::string> warnings;

    struct MoleculeSettings
    {
        int charge;
        unsigned int multiplicity;
        std::vector<Atom> geometry;
    };

    struct BasisSettings
    {
        std::string basisName;
    };

    enum OptionName : size_t
    {
        MAX_ITER,
        ENERGY_TOL,
        DENSITY_TOL,
        DIIS,
        DIIS_MAX_SIZE,
        DIIS_ERROR_TOL,
        DIRECT,
        SCHWARTZ_THRESH,
        DENSITY_THRESH,
        SYMMETRY,
        SYMMETRY_TOL,
        PRINT_FULL_MOS,
        UNRESTRICTED,
        GUESS_MIX,
        DAMP,
        MAX_DAMP_CYCLES,
        STOP_DAMP_THRESH,
        LEVEL_SHIFT,
        LSHIFT_GAP_TOL,
        MAX_LSHIFT_CYCLES,
        STOP_LSHIFT_THRESH,
        COUNT // number of options
    };

    struct SCFInputSettings
    {
        SCFOptions scfOptions;
        std::array<bool, Input::COUNT> optionsSet {};
    };

    static MoleculeSettings parseMoleculeBlock(const std::string& moleculeBlock);
    static BasisSettings parseBasisBlock(const std::string& basisBlock);
    static SCFInputSettings parseSCFBlock(const std::string& SCFBlock);
    static void validateSettings(
        const MoleculeSettings& moleculeSettings,
        const BasisSettings& basisSettings,
        const SCFInputSettings& scfSettings,
        std::vector<std::string>& warnings
    );
};
