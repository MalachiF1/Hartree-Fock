#pragma once

#include "Symmetry/PointGroup.hpp"

#include <string>


namespace Symmetry
{

struct CharacterTable
{
    std::vector<std::string> irrepNames;
    std::vector<std::string> classNames;
    std ::vector<std::vector<double>> characters;

    std::string toString() const;
};

const CharacterTable& getCharacterTable(const PointGroup& pg);

} // namespace Symmetry
