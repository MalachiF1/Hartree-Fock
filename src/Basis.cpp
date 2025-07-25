#include "Basis.hpp"

#include <cctype>
#include <stdexcept>

const std::map<std::string, std::map<int, std::vector<Shell>>> Basis::basisSets = {
    {"sto-3g",
     {
         {1, // Hydrogen
          {{0, {3.425250914, 0.6239137298, 0.1688554040}, {0.1543289673, 0.5353281423, 0.4446345422}}}},
         {8, // Oxygen
          {
              {0, {130.7093214, 23.80886605, 6.443608313}, {0.1543289673, 0.5353281423, 0.4446345422}},   // 1s
              {0, {5.033151319, 1.169596125, 0.380389600}, {-0.09996722919, 0.3995128261, 0.7001154689}}, // 2s
              {1, {5.033151319, 1.169596125, 0.380389600}, {0.1559162750, 0.6076837186, 0.3919573931}}    // 2p
          }},
         {9, // Fluorine
          {
              {0, {166.6791340, 30.36081233, 8.216820672}, {0.15432897, 0.53532814, 0.44463454}},          // 1s
              {0, {6.464803249, 1.502281245, 0.4885884864}, {-0.09996722919, 0.3995128261, 0.7001154689}}, // 2s
              {1, {6.464803249, 1.502281245, 0.4885884864}, {0.155916275, 0.6076837186, 0.3919573931}}     // 2p
          }},
     }}
};

std::map<int, std::vector<Shell>> Basis::getBasis(const std::string& name, const std::vector<int>& elements)
{
    // convert name to lowercase for case-insensitive comparison
    std::string lowercaseName = name;
    std::transform(
        lowercaseName.begin(),
        lowercaseName.end(),
        lowercaseName.begin(),
        [](unsigned char c) { return std::tolower(c); }
    );

    if (basisSets.find(lowercaseName) == basisSets.end())
    {
        throw std::runtime_error("Basis set '" + name + "' not found.");
    }

    const auto& fullBasis = basisSets.at(lowercaseName);
    std::map<int, std::vector<Shell>> resultBasis;

    for (int element : elements)
    {
        if (fullBasis.count(element))
        {
            resultBasis[element] = fullBasis.at(element);
        }
        else
        {
            throw std::runtime_error(
                "Element with atomic number " + std::to_string(element) + " not found in basis set '" + name + "'."
            );
        }
    }
    return resultBasis;
}
