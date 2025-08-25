#pragma once

#include "AtomicOrbital.hpp"
#include "Symmetry/PointGroup.hpp"
#include "Utils.hpp"

#include <string>
#include <system_error>
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
        // if i < j or k < l, swap i and j or k and l to maintain canonical indexing
        if (i < j || k < l)
        {
            if (i < j)
                std::swap(i, j);
            if (k < l)
                std::swap(k, l);
        }

        // Calculate the index in the flattened vector
        size_t big_I = (i * (i + 1) / 2) + j;
        size_t big_K = (k * (k + 1) / 2) + l;

        // if I < K, swap I and K to l to maintian canonical indexing
        if (big_I < big_K)
            std::swap(big_I, big_K);

        // Calculate the index in the flattened vector
        size_t index = (big_I * (big_I + 1) / 2) + big_K;
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
        // if i < j or k < l, swap i and j or k and l to maintain canonical indexing
        if (i < j || k < l)
        {
            if (i < j)
                std::swap(i, j);
            if (k < l)
                std::swap(k, l);
        }

        // Calculate the index in the flattened vector
        size_t big_I = (i * (i + 1) / 2) + j;
        size_t big_K = (k * (k + 1) / 2) + l;

        // if I < K, swap I and K to l to maintian canonical indexing
        if (big_I < big_K)
            std::swap(big_I, big_K);

        // Calculate the index in the flattened vector
        size_t index = (big_I * (big_I + 1) / 2) + big_K;
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

        // The number of unique quartets (i, j, k, l) where i >= j and k >= l
        return M * (M + 1) / 2; // Total number of unique elements in the tensor
    }
};


/**
 * @brief Represents a molecule built from atomic orbital basis functions (which are themselves built by cartesian primitive gaussians).
 *
 * This class encapsulates the properties of a molecule, including its charge, multiplicity, basis set, and geometry.
 * It provides methods to compute molecular properties such as the overlap, kinetic energy, nuclear attraction matrices,
 * the electron repulsion tensor, and the nuclear-nuclear repulsion energy.
 *
 */
class Molecule
{
  public:
    Molecule(
        int charge,
        int multiplicity,
        const std::string& basisName,
        const std::vector<Atom>& geometry,
        bool detectSymmetry      = true,
        double symmetryTolerance = 1e-5
    );


    /**
     * @return A string containing the charge, multiplicity, number of electrons, number of basis functions, and geometry.
     */
    std::string toString() const;

    /**
     * @return The total nuclear repulsion energy.
     */
    double nuclearRepulsion() const;

    /**
     * @return The overlap matrix S as an Eigen::MatrixXd.
     */
    Eigen::MatrixXd overlapMatrix() const;

    /**
     * @return The kinetic energy matrix T as an Eigen::MatrixXd.
     */

    /**
     * @return The kinetic energy matrix T as an Eigen::MatrixXd.
     */
    Eigen::MatrixXd kineticMatrix() const;

    /**
     * @return The nuclear attraction matrix V as an Eigen::MatrixXd.
     */
    Eigen::MatrixXd nuclearAttractionMatrix() const;

    /**
     * Computes the electron repulsion integral tensor for the molecule.
     *
     * @param threshold The Schwartz screening threshold below which integrals are considered negligible and set to zero.
     * @return A ElectronRepulsionTensor object containing the electron repulsion integrals.
     */
    ElectronRepulsionTensor electronRepulsionTensor(double threshold = 1e-10) const;

    /**
     * Computes the Schwartz screening matrix for the molecule.
     * This matrix is used to screen out negligible electron repulsion integrals.
     * This method is used by the SCF class when building the Fock matrix in the direct method.
     * The electronRepulsionTensor method also calcultes and uses the Schwartz screening matrix,
     * but doesn't store it (so it can set the elements in the tensor directly).
     *
     * @return An Eigen::MatrixXd representing the Schwartz screening matrix.
     */
    Eigen::MatrixXd schwartzScreeningMatrix() const;

    // Getters for molecule properties
    const std::vector<AtomicOrbital>& getAtomicOrbitals() const { return atomicOrbitals; }
    size_t getBasisFunctionCount() const { return basisFunctionCount; }
    int getCharge() const { return charge; }
    int getMultiplicity() const { return multiplicity; }
    double getSymmetryTol() const { return symmetryTolerance; }
    size_t getElectronCount() const { return electronCount; }
    const std::vector<Atom>& getGeometry() const { return geometry; }
    const PointGroup& getPointGroup() const { return pointGroup; }
    std::vector<std::string> getAOLabels() const;

  private:
    size_t basisFunctionCount; // number of atomic orbitals in the molecule

    /**
     * Counts the number of electrons in the molecule based on the geometry and charge
     * @return The total number of electrons in the molecule.
     */
    size_t countElectrons() const;

    /**
     * Builds the atomic orbitals for the molecule based on the provided basis set.
     *
     * @param basisName The name of the basis set to use (case insensitive).
     * @Throws std::runtime_error if the basis set is not found.
     */
    void buildBasis(const std::string& basisName);

    const int charge;
    const int multiplicity;
    const double symmetryTolerance;
    std::vector<Atom> geometry;
    size_t electronCount;
    std::vector<AtomicOrbital> atomicOrbitals;
    PointGroup pointGroup;
};
