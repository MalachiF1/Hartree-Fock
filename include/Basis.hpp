#pragma once

#include "Utils.hpp"

#include <map>
#include <string>
#include <vector>

/**
 * @brief Helper class to retrieve basis sets for molecular calculations.
 *
 * This class provides a method to retrieve a basis set for a list of elements.
 * The basis set is defined in the `basis_set` directory, using the Basis Set Exchange JSON format.
 * The file should be named as `<name>.json` (in lower-case).
 */
class Basis
{
  public:
    /*
     * Retrieves a basis set for a list of elements.
     *
     * @param name The name of the basis set (e.g., "STO-3G", "6-31G").
     * @param elements A vector of atomic numbers for which the basis set is requested.
     * @note The basis set must be defined in the `basis_set` directory, with the Basis Set Exchange JSON format.
     *       The file should be named as `<name>.json` and located in the `basis_set` directory.
     * @throw std::runtime_error If the basis set file cannot be found or parsed (e.g. if the JSON file does not
     *             exist, or does not contain the required elements).
     * @return A map where the key is the atomic number and the value is a vector of Shell objects.
     */
    static std::map<int, std::vector<Shell>> getBasis(const std::string& name, const std::vector<int>& elements);
};
