#pragma once

#include "Utils.hpp"

#include <map>
#include <string>
#include <vector>

class Basis
{
  public:
    // Retrieves a basis set for a list of elements.
    static std::map<int, std::vector<Shell>> getBasis(const std::string& name, const std::vector<int>& elements);

  private:
    /*
     *  Hardcoded basis sets for simplicity - change later to read from JSON or other format.
     *  The structure is: map<basis_name, map<atomic_number, vector<shell>>>
     */
    static const std::map<std::string, std::map<int, std::vector<Shell>>> basisSets;
};
