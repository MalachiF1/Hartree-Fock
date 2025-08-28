#include "Input.hpp"

#include "Utils.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

Input::Input(const std::string& filename) : filename(filename) {}

std::pair<Molecule, SCFOptions> Input::read()
{
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

    auto toLowerString = [](const std::string& str)
    {
        std::string lowerStr = str;
        std::transform(
            lowerStr.begin(), lowerStr.end(), lowerStr.begin(), [](unsigned char c) { return std::tolower(c); }
        );
        return lowerStr;
    };

    while (std::getline(inputFile, line))
    {
        std::istringstream iss(line);
        std::string token;

        ++lineNumber;

        while (iss >> token)
        {
            if (token[0] == '#') // Comment
                break;

            if (toLowerString(token) == "$molecule")
            {
                inMoleculeBlock = true;
                continue;
            }
            else if (toLowerString(token) == "$basis")
            {
                inBasisBlock = true;
                continue;
            }
            else if (toLowerString(token) == "$scf")
            {
                inSCFBlock = true;
                continue;
            }
            else if (toLowerString(token) == "$end")
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

    Molecule molecule  = readMoleculeBlock(moleculeBlock.str(), basisBlock.str());
    SCFOptions options = readSCFBlock(SCFBlock.str());

    if (!options.unrestricted && molecule.getMultiplicity() != 1)
        throw std::runtime_error(
            "Multiplicity set to" + std::to_string(molecule.getMultiplicity())
            + " but the SCF method is restricted (RHF). Set multiplicity to 1 or use unrestricted (UHF) method."
        );

    return {molecule, options};
}

Molecule Input::readMoleculeBlock(const std::string& moleculeBlock, const std::string& basisBlock)
{
    std::istringstream molStream(moleculeBlock);
    int charge;
    unsigned multiplicity;
    std::vector<Atom> geometry;
    molStream >> charge;
    if (molStream.fail())
        throw std::runtime_error("Invalid charge in $molecule block.");

    molStream >> multiplicity;
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
                + "'"
            );
        }

        if (std::all_of(atomStr.begin(), atomStr.end(), ::isdigit))
        {
            unsigned atomicNumber = std::stoi(atomStr);
            geometry.emplace_back(atomicNumber, Vec3(x, y, z));
        }
        else
        {
            auto it = Utils::nameToAtomicNumber.find(atomStr);
            if (it == Utils::nameToAtomicNumber.end())
                throw std::runtime_error("Invalid atom symbol \"" + atomStr + "\" in $molecule block.");
            geometry.emplace_back(it->second, Vec3(x, y, z));
        }
    }
    if (geometry.empty())
        throw std::runtime_error("No atoms specified in $molecule block.");

    std::istringstream basisStream(basisBlock);
    std::string basisSetName;
    basisStream >> basisSetName;
    if (basisStream.fail())
        throw std::runtime_error("Error in basis set specification in $basis block.");

    // Some basic validation
    Molecule molecule(charge, multiplicity, basisSetName, geometry);

    if (molecule.getElectronCount() + 1 < multiplicity || (molecule.getElectronCount() - multiplicity) % 2 == 0)
        throw std::runtime_error("Invalid combination of charge and multiplicity.");

    return molecule;
}


SCFOptions Input::readSCFBlock(const std::string& SCFBlock)
{

    auto toLowerString = [](const std::string& str)
    {
        std::string lowerStr = str;
        std::transform(
            lowerStr.begin(), lowerStr.end(), lowerStr.begin(), [](unsigned char c) { return std::tolower(c); }
        );
        return lowerStr;
    };

    std::istringstream SCFStream(SCFBlock);
    SCFOptions options;
    std::string token;
    while (SCFStream >> token)
    {
        std::string tokenLwr = toLowerString(token);
        if (tokenLwr == "max_iter")
        {
            SCFStream >> options.maxIter;
            if (options.maxIter == 0)
                throw std::runtime_error(token + " in $scf block must be greater than 0.");
        }
        else if (tokenLwr == "energy_tol")
        {
            SCFStream >> options.energyTol;
            if (options.energyTol < 0)
                throw std::runtime_error(
                    "Invalid value for " + token + " in $scf block: " + std::to_string(options.energyTol)
                    + ". Must be non-negative."
                );
        }
        else if (tokenLwr == "density_tol")
        {
            SCFStream >> options.densityTol;
            if (options.densityTol < 0)
                throw std::runtime_error(
                    "Invalid value for " + token + " in $scf block: " + std::to_string(options.densityTol)
                    + ". Must be non-negative."
                );
        }
        else if (tokenLwr == "diis")
        {
            std::string boolStr;
            SCFStream >> boolStr;
            if (toLowerString(boolStr) == "true" || boolStr == "1")
                options.useDIIS = true;
            else if (toLowerString(boolStr) == "false" || boolStr == "0")
                options.useDIIS = false;
            else
                throw std::runtime_error("Invalid boolean value for " + token + " in $scf block: " + boolStr);
        }
        else if (tokenLwr == "diis_size")
        {
            SCFStream >> options.DIISmaxSize;
            if (options.DIISmaxSize < 1)
                throw std::runtime_error(token + " in $scf block must be greater than 0.");
        }
        else if (tokenLwr == "diis_start")
        {
            SCFStream >> options.DIISstart;
            if (options.DIISstart < 1)
                throw std::runtime_error(token + " in $scf block must be greater than 0.");
        }
        else if (tokenLwr == "diis_tol")
        {
            SCFStream >> options.DIISErrorTol;
            if (options.DIISErrorTol < 0)
                throw std::runtime_error(
                    "Invalid value for " + token + " in $scf block: " + std::to_string(options.DIISErrorTol)
                    + ". Must be non-negative."
                );
        }
        else if (tokenLwr == "direct")
        {
            std::string boolStr;
            SCFStream >> boolStr;
            if (toLowerString(boolStr) == "true" || boolStr == "1")
                options.direct = true;
            else if (toLowerString(boolStr) == "false" || boolStr == "0")
                options.direct = false;
            else
                throw std::runtime_error("Invalid boolean value for " + token + " in $scf block: " + boolStr);
        }
        else if (tokenLwr == "schwartz_thresh")
        {
            SCFStream >> options.schwartzThreshold;
            if (options.schwartzThreshold < 0)
                throw std::runtime_error(
                    "Invalid value for " + token + " in $scf block: " + std::to_string(options.schwartzThreshold)
                    + ". Must be non-negative."
                );
        }
        else if (tokenLwr == "density_thresh")
        {
            SCFStream >> options.densityThreshold;
            if (options.densityThreshold < 0)
                throw std::runtime_error(
                    "Invalid value for " + token + " in $scf block: " + std::to_string(options.densityThreshold)
                    + ". Must be non-negative."
                );
        }
        else if (tokenLwr == "symmetry")
        {
            std::string boolStr;
            SCFStream >> boolStr;
            if (toLowerString(boolStr) == "true" || boolStr == "1")
                options.useSymmetry = true;
            else if (toLowerString(boolStr) == "false" || boolStr == "0")
                options.useSymmetry = false;
            else
                throw std::runtime_error("Invalid boolean value for " + token + " in $scf block: " + boolStr);
        }
        else if (tokenLwr == "print_full_mos")
        {
            std::string boolStr;
            SCFStream >> boolStr;
            if (toLowerString(boolStr) == "true" || boolStr == "1")
                options.printFullMOs = true;
            else if (toLowerString(boolStr) == "false" || boolStr == "0")
                options.printFullMOs = false;
            else
                throw std::runtime_error("Invalid boolean value for " + token + " in $scf block: " + boolStr);
        }
        else if (tokenLwr == "symmetry_tol")
        {
            SCFStream >> options.symmetryTolerance;
            if (options.symmetryTolerance < 0)
                throw std::runtime_error(
                    "Invalid value for " + token + " in $scf block: " + std::to_string(options.symmetryTolerance)
                    + ". Must be non-negative."
                );
        }
        else if (tokenLwr == "unrestricted")
        {
            std::string boolStr;
            SCFStream >> boolStr;
            if (toLowerString(boolStr) == "true" || boolStr == "1")
                options.unrestricted = true;
            else if (toLowerString(boolStr) == "false" || boolStr == "0")
                options.unrestricted = false;
            else
                throw std::runtime_error("Invalid boolean value for " + token + " in $scf block: " + boolStr);
        }
        else if (tokenLwr == "guess_mix")
        {
            std::string mixStr;
            SCFStream >> mixStr;
            if (toLowerString(mixStr) == "true")
                options.guessMix = 1;
            else if (toLowerString(mixStr) == "false")
                options.guessMix = 0;
            else if (std::all_of(mixStr.begin(), mixStr.end(), ::isdigit))
            {
                int mixVal = std::stoi(mixStr);
                if (mixVal < 0 || mixVal > 10)
                    throw std::runtime_error(
                        "Invalid integer value for " + token + " in $scf block: \"" + mixStr
                        + "\". Must be between 0 and 10."
                    );
                options.guessMix = mixVal;
            }
            else
                throw std::runtime_error(
                    "Invalid value for " + token + " in $scf block: \"" + mixStr
                    + "\". Must be an integer between 0 and 10"
                );
        }
        else if (tokenLwr == "damp")
        {
            SCFStream >> options.damp;
            if (options.damp < 0 || options.damp > 100)
                throw std::runtime_error(
                    "Invalid integer value for " + token + " in $scf block: " + std::to_string(options.damp)
                    + ". Must be between 0 and 100."
                );
        }
        else if (tokenLwr == "max_damp_cycles")
        {
            SCFStream >> options.maxDampIter;
            if (options.maxDampIter <= 0)
                throw std::runtime_error(
                    "Invalid value for " + token + " in $scf block: " + std::to_string(options.maxDampIter)
                    + ". Must be positive integer."
                );
        }
        else if (tokenLwr == "stop_damp_thresh")
        {
            SCFStream >> options.stopDampThresh;
            if (options.stopDampThresh < 0)
                throw std::runtime_error(
                    "Invalid value for " + token + " in $scf block: " + std::to_string(options.stopDampThresh)
                    + ". Must be non-negative."
                );
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

    // Make some basic checks on the options.
    if (options.guessMix > 0 && !options.unrestricted)
        throw std::runtime_error("The guess_mix option can only be used with unrestricted calculations.");


    return options;
}
