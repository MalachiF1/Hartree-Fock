#include "Input.hpp"

#include "Utils.hpp"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <ranges>
#include <sstream>
#include <string>
#include <vector>

Input::Input(const std::string& filename) : filename(filename) {}

const std::vector<std::string>& Input::getWarnings() const
{
    return this->warnings;
}

Input::InputSettings Input::read()
{
    this->warnings.clear();

    std::ifstream inputFile(this->filename);
    if (!inputFile.is_open())
    {
        throw std::runtime_error("Could not open input file: " + filename);
    }

    // Flags
    bool inMoleculeBlock = false;
    bool inBasisBlock    = false;
    bool inSCFBlock      = false;

    // Input blocks
    std::stringstream moleculeBlock;
    std::stringstream basisBlock;
    std::stringstream SCFBlock;

    // Remove empty lines and comments, seperate into the blocks.
    std::string line;
    size_t lineNumber = 0;

    auto toLowerString = [](std::string_view sv) -> std::string
    {
        auto view = sv | std::ranges::views::transform(::tolower);
        return std::string(view.begin(), view.end());
    };

    while (std::getline(inputFile, line))
    {
        std::istringstream iss(line);
        std::string token;

        ++lineNumber;

        while (iss >> token)
        {
            std::string lwrToken = toLowerString(token);
            if (token[0] == '#') // Comment
                break;

            if (lwrToken == "$molecule")
            {
                inMoleculeBlock = true;
                continue;
            }
            else if (lwrToken == "$basis")
            {
                inBasisBlock = true;
                continue;
            }
            else if (lwrToken == "$scf")
            {
                inSCFBlock = true;
                continue;
            }
            else if (lwrToken == "$end")
            {
                inMoleculeBlock = false;
                inBasisBlock    = false;
                inSCFBlock      = false;
                continue;
            }
            else if (inSCFBlock)
            {
                SCFBlock << token << " ";
            }
            else if (inMoleculeBlock)
            {
                moleculeBlock << token << " ";
            }
            else if (inBasisBlock)
            {
                basisBlock << token << " ";
            }
            else
            {
                throw std::runtime_error("Unexpected token at line " + std::to_string(lineNumber) + ": " + token);
            }
        }
    }
    inputFile.close();

    MoleculeSettings moleculeSettings = parseMoleculeBlock(moleculeBlock.str());
    BasisSettings basisSettings       = parseBasisBlock(basisBlock.str());
    SCFInputSettings scfInputSettings = parseSCFBlock(SCFBlock.str());

    validateSettings(moleculeSettings, basisSettings, scfInputSettings, this->warnings);

    Molecule molecule(
        moleculeSettings.charge,
        moleculeSettings.multiplicity,
        moleculeSettings.geometry,
        scfInputSettings.scfOptions.useSymmetry,
        scfInputSettings.scfOptions.symmetryTolerance
    );

    // Important to use molecule.getGeometry() and not moleculeSettings.geometry
    // because the Molecule constructor may have modified the geometry
    BasisSet basis(basisSettings.basisName, molecule.getGeometry());

    SCFOptions options = scfInputSettings.scfOptions;
    if (!scfInputSettings.optionsSet[UNRESTRICTED])
    {
        options.unrestricted = molecule.getMultiplicity() != 1;
    }

    if (!scfInputSettings.optionsSet[MAX_DAMP_CYCLES])
    {
        options.maxDampIter = options.maxIter;
    }

    if (!scfInputSettings.optionsSet[MAX_LSHIFT_CYCLES])
    {
        options.maxLshiftIter = options.maxIter;
    }

    return {.molecule = molecule, .scfOptions = options, .basis = basis};
}

Input::MoleculeSettings Input::parseMoleculeBlock(const std::string& moleculeBlock)
{
    std::istringstream molStream(moleculeBlock);
    MoleculeSettings settings;

    molStream >> settings.charge;
    if (molStream.fail())
        throw std::runtime_error("Invalid charge in $molecule block.");

    molStream >> settings.multiplicity;
    if (molStream.fail())
        throw std::runtime_error("Invalid multiplicity in $molecule block.");

    std::string atomStr, xStr, yStr, zStr;
    while (molStream >> atomStr)
    {
        if (!(molStream >> xStr >> yStr >> zStr))
        {
            throw std::runtime_error("Incomplete atom specification in $molecule block for atom '" + atomStr + "'.");
        }

        std::istringstream coordStream(xStr + " " + yStr + " " + zStr);
        double x, y, z;
        coordStream >> x >> y >> z;

        if (coordStream.fail())
        {
            throw std::runtime_error(
                "Invalid coordinates in $molecule block for atom " + atomStr + ": '" + xStr + " " + yStr + " " + zStr
                + ""
            );
        }

        if (std::ranges::all_of(atomStr, ::isdigit))
        {
            unsigned atomicNumber = std::stoi(atomStr);
            settings.geometry.emplace_back(atomicNumber, Vec3(x, y, z));
        }
        else
        {
            auto it = Utils::nameToAtomicNumber.find(atomStr);
            if (it == Utils::nameToAtomicNumber.end())
                throw std::runtime_error("Invalid atom symbol \"" + atomStr + "\" in $molecule block.");
            settings.geometry.emplace_back(it->second, Vec3(x, y, z));
        }
    }
    if (settings.geometry.empty())
        throw std::runtime_error("No atoms specified in $molecule block.");

    // check for atoms with distance < 1e-4
    for (size_t i = 0; i < settings.geometry.size(); ++i)
    {
        for (size_t j = i + 1; j < settings.geometry.size(); ++j)
        {
            double dist = (settings.geometry[i].coords - settings.geometry[j].coords).norm();
            if (dist < 1e-4)
            {
                throw std::runtime_error(
                    "Atoms " + std::to_string(i + 1) + " and " + std::to_string(j + 1)
                    + " in $molecule block are too close together (distance < 1e-4)."
                );
            }
        }
    }


    return settings;
}

Input::BasisSettings Input::parseBasisBlock(const std::string& basisBlock)
{
    std::istringstream basisStream(basisBlock);
    BasisSettings settings;
    basisStream >> settings.basisName;
    std::string extraToken;
    if (basisStream.fail() || (basisStream >> extraToken))
        throw std::runtime_error("Error in basis set specification in $basis block.");
    return settings;
}

Input::SCFInputSettings Input::parseSCFBlock(const std::string& SCFBlock)
{

    auto toLowerString = [](std::string_view sv) -> std::string
    {
        auto view = sv | std::ranges::views::transform(::tolower);
        return std::string(view.begin(), view.end());
    };

    Input::SCFInputSettings settings;

    std::istringstream SCFStream(SCFBlock);
    std::string token;
    while (SCFStream >> token)
    {
        std::string tokenLwr = toLowerString(token);
        if (tokenLwr == "max_iter")
        {
            SCFStream >> settings.scfOptions.maxIter;
            settings.optionsSet[MAX_ITER] = true;
        }
        else if (tokenLwr == "energy_tol")
        {
            SCFStream >> settings.scfOptions.energyTol;
            settings.optionsSet[ENERGY_TOL] = true;
        }
        else if (tokenLwr == "density_tol")
        {
            SCFStream >> settings.scfOptions.densityTol;
            settings.optionsSet[DENSITY_TOL] = true;
        }
        else if (tokenLwr == "diis")
        {
            std::string boolStr;
            SCFStream >> boolStr;
            if (toLowerString(boolStr) == "true" || boolStr == "1")
                settings.scfOptions.useDIIS = true;
            else if (toLowerString(boolStr) == "false" || boolStr == "0")
                settings.scfOptions.useDIIS = false;
            else
                throw std::runtime_error("Invalid boolean value for " + token + " in $scf block: " + boolStr);
            settings.optionsSet[DIIS] = true;
        }
        else if (tokenLwr == "diis_size")
        {
            SCFStream >> settings.scfOptions.DIISmaxSize;
            settings.optionsSet[DIIS_MAX_SIZE] = true;
        }
        else if (tokenLwr == "diis_tol")
        {
            SCFStream >> settings.scfOptions.DIISErrorTol;
            settings.optionsSet[DIIS_ERROR_TOL] = true;
        }
        else if (tokenLwr == "direct")
        {
            std::string boolStr;
            SCFStream >> boolStr;
            if (toLowerString(boolStr) == "true" || boolStr == "1")
                settings.scfOptions.direct = true;
            else if (toLowerString(boolStr) == "false" || boolStr == "0")
                settings.scfOptions.direct = false;
            else
                throw std::runtime_error("Invalid boolean value for " + token + " in $scf block: " + boolStr);
            settings.optionsSet[DIRECT] = true;
        }
        else if (tokenLwr == "schwartz_thresh")
        {
            SCFStream >> settings.scfOptions.schwartzThreshold;
            settings.optionsSet[SCHWARTZ_THRESH] = true;
        }
        else if (tokenLwr == "density_thresh")
        {
            SCFStream >> settings.scfOptions.densityThreshold;
            settings.optionsSet[DENSITY_THRESH] = true;
        }
        else if (tokenLwr == "symmetry")
        {
            std::string boolStr;
            SCFStream >> boolStr;
            if (toLowerString(boolStr) == "true" || boolStr == "1")
                settings.scfOptions.useSymmetry = true;
            else if (toLowerString(boolStr) == "false" || boolStr == "0")
                settings.scfOptions.useSymmetry = false;
            else
                throw std::runtime_error("Invalid boolean value for " + token + " in $scf block: " + boolStr);
            settings.optionsSet[SYMMETRY] = true;
        }
        else if (tokenLwr == "print_full_mos")
        {
            std::string boolStr;
            SCFStream >> boolStr;
            if (toLowerString(boolStr) == "true" || boolStr == "1")
                settings.scfOptions.printFullMOs = true;
            else if (toLowerString(boolStr) == "false" || boolStr == "0")
                settings.scfOptions.printFullMOs = false;
            else
                throw std::runtime_error("Invalid boolean value for " + token + " in $scf block: " + boolStr);
            settings.optionsSet[PRINT_FULL_MOS] = true;
        }
        else if (tokenLwr == "symmetry_tol")
        {
            SCFStream >> settings.scfOptions.symmetryTolerance;
            settings.optionsSet[SYMMETRY_TOL] = true;
        }
        else if (tokenLwr == "unrestricted")
        {
            std::string boolStr;
            SCFStream >> boolStr;
            if (toLowerString(boolStr) == "true" || boolStr == "1")
                settings.scfOptions.unrestricted = true;
            else if (toLowerString(boolStr) == "false" || boolStr == "0")
                settings.scfOptions.unrestricted = false;
            else
                throw std::runtime_error("Invalid boolean value for " + token + " in $scf block: " + boolStr);
            settings.optionsSet[UNRESTRICTED] = true;
        }
        else if (tokenLwr == "guess_mix")
        {
            std::string mixStr;
            SCFStream >> mixStr;
            if (toLowerString(mixStr) == "true")
                settings.scfOptions.guessMix = 1;
            else if (toLowerString(mixStr) == "false")
                settings.scfOptions.guessMix = 0;
            else if (std::ranges::all_of(mixStr, ::isdigit))
            {
                int mixVal = std::stoi(mixStr);
                if (mixVal < 0 || mixVal > 10)
                    throw std::runtime_error(
                        "Invalid integer value for " + token + " in $scf block: \"" + mixStr
                        + ". Must be between 0 and 10."
                    );
                settings.scfOptions.guessMix = mixVal;
            }
            else
            {
                throw std::runtime_error(
                    "Invalid value for " + token + " in $scf block: \"" + mixStr
                    + ". Must be an integer between 0 and 10"
                );
            }
            settings.optionsSet[GUESS_MIX] = true;
        }
        else if (tokenLwr == "damp")
        {
            SCFStream >> settings.scfOptions.damp;
            settings.optionsSet[DAMP] = true;
        }
        else if (tokenLwr == "max_damp_cycles")
        {
            SCFStream >> settings.scfOptions.maxDampIter;
            settings.optionsSet[MAX_DAMP_CYCLES] = true;
        }
        else if (tokenLwr == "stop_damp_thresh")
        {
            SCFStream >> settings.scfOptions.stopDampThresh;
            settings.optionsSet[STOP_DAMP_THRESH] = true;
        }
        else if (tokenLwr == "level_shift")
        {
            SCFStream >> settings.scfOptions.levelShift;
            settings.optionsSet[LEVEL_SHIFT] = true;
        }
        else if (tokenLwr == "lshift_gap_tol")
        {
            SCFStream >> settings.scfOptions.lshiftGapTol;
            settings.optionsSet[LSHIFT_GAP_TOL] = true;
        }
        else if (tokenLwr == "max_lshift_cycles")
        {
            SCFStream >> settings.scfOptions.maxLshiftIter;
            settings.optionsSet[MAX_LSHIFT_CYCLES] = true;
        }
        else if (tokenLwr == "stop_lshift_thresh")
        {
            SCFStream >> settings.scfOptions.stopLshiftThresh;
            settings.optionsSet[STOP_LSHIFT_THRESH] = true;
        }
        else
        {
            throw std::runtime_error("Unknown option in $scf block: " + token);
        }

        if (SCFStream.fail())
        {
            throw std::runtime_error("Error reading value for option: " + token);
        }
    }


    return settings;
}

void Input::validateSettings(
    const MoleculeSettings& moleculeSettings,
    const BasisSettings& basisSettings,
    const SCFInputSettings& scfSettings,
    std::vector<std::string>& warnings
)
{
    std::vector<std::string> errors;

    if (scfSettings.optionsSet[MAX_ITER] && scfSettings.scfOptions.maxIter == 0)
        errors.emplace_back("MAX_ITER in $scf block must be greater than 0.");

    if (scfSettings.optionsSet[ENERGY_TOL] && scfSettings.scfOptions.energyTol < 0)
        errors.emplace_back("ENERGY_TOL in $scf block must be non-negative.");

    if (scfSettings.optionsSet[DENSITY_TOL] && scfSettings.scfOptions.densityTol < 0)
        errors.emplace_back("DENSITY_TOL in $scf block must be non-negative.");

    if (scfSettings.optionsSet[DIIS_MAX_SIZE] && scfSettings.scfOptions.DIISmaxSize < 1)
        errors.emplace_back("DIIS_SIZE in $scf block must be greater than 0.");

    if (scfSettings.optionsSet[DIIS_ERROR_TOL] && scfSettings.scfOptions.DIISErrorTol < 0)
        errors.emplace_back("DIIS_TOL in $scf block must be non-negative.");

    if (scfSettings.optionsSet[SCHWARTZ_THRESH] && scfSettings.scfOptions.schwartzThreshold < 0)
        errors.emplace_back("SCHWARTZ_THRESH in $scf block must be non-negative.");

    if (scfSettings.optionsSet[DENSITY_THRESH] && scfSettings.scfOptions.densityThreshold < 0)
        errors.emplace_back("DENSITY_THRESH in $scf block must be non-negative.");

    if (scfSettings.optionsSet[SYMMETRY_TOL] && scfSettings.scfOptions.symmetryTolerance < 0)
        errors.emplace_back("SYMMETRY_TOL in $scf block must be non-negative.");

    if (scfSettings.optionsSet[DAMP] && (scfSettings.scfOptions.damp < 0 || scfSettings.scfOptions.damp > 100))
        errors.emplace_back("DAMP in $scf block must be between 0 and 100.");

    if (scfSettings.optionsSet[MAX_DAMP_CYCLES] && scfSettings.scfOptions.maxDampIter <= 0)
        errors.emplace_back("MAX_DAMP_CYCLES in $scf block must be a positive integer.");

    if (scfSettings.optionsSet[STOP_DAMP_THRESH] && scfSettings.scfOptions.stopDampThresh < 0)
        errors.emplace_back("STOP_DAMP_THRESH in $scf block must be non-negative.");

    if (scfSettings.optionsSet[LEVEL_SHIFT] && scfSettings.scfOptions.levelShift < 0)
        errors.emplace_back("LEVEL_SHIFT in $scf block must be non-negative.");

    if (scfSettings.optionsSet[LSHIFT_GAP_TOL] && scfSettings.scfOptions.lshiftGapTol <= 0)
        errors.emplace_back("LSHIFT_GAP_TOL in $scf block must be positive.");

    if (scfSettings.optionsSet[MAX_LSHIFT_CYCLES] && scfSettings.scfOptions.maxLshiftIter <= 0)
        errors.emplace_back("MAX_LSHIFT_CYCLES in $scf block must be a positive integer.");

    if (scfSettings.optionsSet[STOP_LSHIFT_THRESH] && scfSettings.scfOptions.stopLshiftThresh < 0)
        errors.emplace_back("STOP_LSHIFT_THRESH in $scf block must be non-negative.");

    size_t electronCount = std::accumulate(
        moleculeSettings.geometry.begin(),
        moleculeSettings.geometry.end(),
        -moleculeSettings.charge,
        [](size_t acc, const Atom& atom) { return acc + atom.atomicNumber; }
    );

    bool unrestricted = scfSettings.scfOptions.unrestricted;
    if (!scfSettings.optionsSet[UNRESTRICTED])
    {
        unrestricted = (moleculeSettings.multiplicity != 1);
    }

    if (scfSettings.optionsSet[UNRESTRICTED])
    {
        if (!scfSettings.scfOptions.unrestricted && moleculeSettings.multiplicity != 1)
            errors.emplace_back(
                "Multiplicity set to " + std::to_string(moleculeSettings.multiplicity)
                + " but the SCF method is restricted (RHF). Set multiplicity to 1 or use unrestricted (UHF) method."
            );
    }

    if (scfSettings.scfOptions.guessMix > 0 && !unrestricted)
        errors.emplace_back("The guess_mix option can only be used with unrestricted calculations.");

    if (electronCount + 1 < moleculeSettings.multiplicity || (electronCount - moleculeSettings.multiplicity) % 2 == 0)
        errors.emplace_back("Invalid combination of charge and multiplicity.");

    // --- Warnings ---
    if (scfSettings.optionsSet[DIIS_MAX_SIZE] && !scfSettings.scfOptions.useDIIS)
        warnings.emplace_back("DIIS_SIZE is set, but it will have no effect because diis is disabled.");
    if (scfSettings.optionsSet[DIIS_ERROR_TOL] && !scfSettings.scfOptions.useDIIS)
        warnings.emplace_back("DIIS_TOL is set, but it will have no effect because diis is disabled.");
    if (scfSettings.optionsSet[MAX_DAMP_CYCLES] && scfSettings.scfOptions.damp == 0)
        warnings.emplace_back("MAX_DAMP_CYCLES is set, but it will have no effect because damping is off.");
    if (scfSettings.optionsSet[STOP_DAMP_THRESH] && scfSettings.scfOptions.damp == 0)
        warnings.emplace_back("STOP_DAMP_THRESH is set, but it will have no effect because damping is off.");
    if (scfSettings.optionsSet[DENSITY_THRESH] && !scfSettings.scfOptions.direct)
        warnings.emplace_back("DENSITY_THRESH is set, but it will have no effect because direct SCF is not enabled.");
    if (scfSettings.optionsSet[SYMMETRY_TOL] && !scfSettings.scfOptions.useSymmetry)
        warnings.emplace_back("SYMMETRY_TOL is set, but it will have no effect because symmetry is not enabled.");
    if (scfSettings.optionsSet[LSHIFT_GAP_TOL] && scfSettings.scfOptions.levelShift == 0)
        warnings.emplace_back("LSHIFT_GAP_TOL is set, but it will have no effect because level shifting is off.");
    if (scfSettings.optionsSet[MAX_LSHIFT_CYCLES] && scfSettings.scfOptions.levelShift == 0)
        warnings.emplace_back("MAX_LSHIFT_CYCLES is set, but it will have no effect because level shifting is off.");
    if (scfSettings.optionsSet[STOP_LSHIFT_THRESH] && scfSettings.scfOptions.levelShift == 0)
        warnings.emplace_back("STOP_LSHIFT_THRESH is set, but it will have no effect because level shifting is off.");

    auto toLowerString = [](std::string_view sv) -> std::string
    {
        auto view = sv | std::ranges::views::transform(::tolower);
        return std::string(view.begin(), view.end());
    };
    std::string lowercaseName = toLowerString(basisSettings.basisName);
    size_t startPos           = 0;
    while ((startPos = lowercaseName.find("*", startPos)) != std::string::npos)
    {
        lowercaseName.replace(startPos, 1, "_st_");
        startPos += 4;
    }
    std::string basisPath = "basis_sets/" + lowercaseName + ".json";
    if (!std::filesystem::exists(basisPath))
    {
        throw std::runtime_error("Basis set '" + basisSettings.basisName + "' not found.");
    }

    if (!errors.empty())
    {
        std::string errorMessage = "INPUT ERROR:\n";
        for (const auto& error : errors) { errorMessage += "- " + error + "\n"; }
        throw std::runtime_error(errorMessage);
    }
}
