#pragma once

#include "Utils.hpp"

#include <algorithm>
#include <string>
#include <vector>

/**
 * @brief Represents the electron repulsion tensor for a molecule.
 *
 * Flattened 1D vector representing the electron repulsion tensor.
 * We only store the unique elements (8-fold symmetry). The tensor can be indexed normaly as ElectronRepulsionTensor(i, j, k, l).
 * This class stores the electron repulsion integrals in canonical indexing order: (i >= j, k >= l). However, it allows calling
 * the elements with any indices, and it will automatically swap the indices to maintain the canonical indexing.
 */
class ElectronRepulsionTensor : public std::vector<double>
{
  public:
    using std::vector<double>::vector; // inheritance from std::vector<double> allows us to use the vector interface directly

    ElectronRepulsionTensor() = default;

    /**
     * Constructs an electron repulsion tensor for a given basis size.
     * The tensor is initialized with zeros, and is of the same size as the number of unique elements.
     *
     * @param basisSize The number of basis functions.
     */
    ElectronRepulsionTensor(size_t basisSize) : std::vector<double>(tensorSize(basisSize), 0.0) {}

    /**
     * Overload the () operator to access the tensor elements (by reference, so you can set them as well).
     * Indices are swaped automatically to maintain canonical indexing.
     *
     * @param i The first index
     * @param j The second index
     * @param k The third index
     * @param l The fourth index
     * @return A reference to the element at (i, j, k, l) in the tensor.
     */
    double& operator()(size_t i, size_t j, size_t k, size_t l)
    {
        // Enforce i >= j and k >= l
        const auto [j_, i_] = std::minmax(j, i);
        const auto [l_, k_] = std::minmax(l, k);

        const size_t big_I = i_ * (i_ + 1) / 2 + j_;
        const size_t big_K = k_ * (k_ + 1) / 2 + l_;

        // Enforce big_I >= big_K
        const auto [big_K_, big_I_] = std::minmax(big_K, big_I);

        const size_t index = big_I_ * (big_I_ + 1) / 2 + big_K_;
        return std::vector<double>::operator[](index);
    };

    /**
     * Overload the () operator to access the tensor elements (by reference, so you can set them as well).
     * Indices are swaped automatically to maintain canonical indexing.
     *
     * @param i The first index
     * @param j The second index
     * @param k The third index
     * @param l The fourth index
     * @return A reference to the element at (i, j, k, l) in the tensor.
     */
    double const& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        // Enforce i >= j and k >= l
        const auto [j_, i_] = std::minmax(j, i);
        const auto [l_, k_] = std::minmax(l, k);

        const size_t big_I = i_ * (i_ + 1) / 2 + j_;
        const size_t big_K = k_ * (k_ + 1) / 2 + l_;

        // Enforce big_I >= big_K
        const auto [big_K_, big_I_] = std::minmax(big_K, big_I);

        const size_t index = big_I_ * (big_I_ + 1) / 2 + big_K_;
        return std::vector<double>::operator[](index);
    };

  private:
    /**
     * Returns the size of the tensor (number of unique elements) for a given basis size.
     *
     * @param basisSize The number of basis functions.
     * @return The size of the tensor
     */
    size_t tensorSize(size_t basisSize) const
    {
        // The number of unique pairs (i, j) where i >= j
        size_t M = basisSize * (basisSize + 1) / 2;

        // The number of unique quartets (i, j, k, l) where i >= j, k >= l and (i, j) >= (k, l).
        return M * (M + 1) / 2;
    }
};


/**
 * @brief Represents a molecule built from a set of atoms.
 *
 * This class encapsulates the properties of a molecule, including its charge, multiplicity, geometry, and point-group
 * symmetry. It provides methods to obtain these properties.
 *
 */
class Molecule
{
  public:
    Molecule(
        int charge, int multiplicity, const std::vector<Atom>& geometry, bool detectSymmetry = true, double symmetryTolerance = 1e-5
    );

    /**
     * @return A string containing the charge, multiplicity, number of electrons, number of basis functions, and geometry.
     */
    std::string toString() const;

    /**
     * @return The total nuclear repulsion energy.
     */
    double nuclearRepulsion() const;

    // Getters for molecule properties
    int getCharge() const { return charge; }
    size_t getMultiplicity() const { return multiplicity; }
    double getSymmetryTol() const { return symmetryTolerance; }
    size_t getElectronCount() const { return electronCount; }
    const std::vector<Atom>& getGeometry() const { return geometry; }
    const std::string& getPointGroup() const { return pointGroupName; }
    // const PointGroup& getAbelianSubgroup() const { return abelianSubgroup; }
    std::vector<std::string> getAOLabels() const;

  private:
    size_t basisFunctionCount; // number of atomic orbitals in the molecule

    /**
     * Counts the number of electrons in the molecule based on the geometry and charge
     * @return The total number of electrons in the molecule.
     */
    size_t countElectrons() const;

    const int charge;
    const int multiplicity;
    const double symmetryTolerance;
    std::vector<Atom> geometry;
    size_t electronCount;
    std::string pointGroupName;
};
