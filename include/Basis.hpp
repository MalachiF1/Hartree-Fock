#pragma once

#include "Utils.hpp"

#include <map>
#include <string>
#include <vector>

/**
 * Struct to hold shell information from a basis set
 */
struct RawShell
{
    RawShell() = default;
    RawShell(int angularMomentum, const std::vector<double>& exponents, const std::vector<double>& coefficients) :
        angularMomentum(angularMomentum), exponents(exponents), coefficients(coefficients)
    {
    }
    unsigned angularMomentum;
    std::vector<double> exponents;
    std::vector<double> coefficients;
};


/**
 * Struct to hold shell information with respect to the basis set data arrays.
 */
struct Shell
{
    unsigned l;
    size_t nao;
    size_t nprim;
    size_t primOffset;
    size_t aoOffset;
    size_t coeffOffset;

    std::string toString() const;
};

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
    Basis(const std::string& name, const std::vector<Atom>& geometry);
    Basis() = default;
    size_t getShellCount() const { return shells.size(); }
    size_t getPrimitiveCount() const { return exps.size(); }
    size_t getAOCount() const { return lx.size(); }
    const std::vector<Shell>& getShells() const { return shells; }
    const std::vector<double>& getExponents() const { return exps; }
    const std::vector<double>& getCoefficients() const { return coefficients; }
    const std::vector<double>& getCx() const { return cx; }
    const std::vector<double>& getCy() const { return cy; }
    const std::vector<double>& getCz() const { return cz; }
    const std::vector<unsigned>& getLx() const { return lx; }
    const std::vector<unsigned>& getLy() const { return ly; }
    const std::vector<unsigned>& getLz() const { return lz; }


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
    static std::map<int, std::vector<RawShell>> readBasis(const std::string& name, const std::vector<int>& elements);

  private:
    // Name of the basis set (e.g., "STO-3G", "6-31G").
    std::string name;

    // Values per primitives per shell (not per AO).
    // Indexing: primOffset + p  -  primitive p of shell
    std::vector<double> exps;
    std::vector<double> cx;
    std::vector<double> cy;
    std::vector<double> cz;

    // Values per AO.
    // Indexing: aoOffset + a  -  atomic orbital a of shell
    std::vector<unsigned> lx;
    std::vector<unsigned> ly;
    std::vector<unsigned> lz;

    // Values per primitives per AO.
    // Indexing: (coeffOffset + p * shell.nao) + a  -  coefficient for primitive p of atomic orbital a of shell
    std::vector<double> coefficients; // includes the normalization factor

    // Shell information
    std::vector<Shell> shells;
};
