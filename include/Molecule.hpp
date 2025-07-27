#pragma once

#include "AtomicOrbital.hpp"
#include "Utils.hpp"

#include <string>
#include <vector>


/*
 * Flattened 1D vector representing the electron repulsion tensor.
 * We only store the unique elements (8-fold symmetry).
 * Canonical indexing is paramnount (i >= j, k >= l).
 */
class ElectronRepulsionTensor : public std::vector<double>
{
  public:
    // inheritance from std::vector<double> allows us to use the vector interface directly
    using std::vector<double>::vector;

    ElectronRepulsionTensor() = default;

    /*
     * Constructs an electron repulsion tensor for a given basis size.
     * The tensor is initialized with zeros.
     * @param basisSize The number of basis functions.
     */
    ElectronRepulsionTensor(size_t basisSize) : std::vector<double>(tensorSize(basisSize), 0.0) {}

    /*
     * overload the () operator to access the tensor elements (by reference, so you can set them as well).
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

    /*
     * overload the () with const reference operator
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
    /*
     * Returns the size of the tensor for a given basis size.
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


class Molecule
{
  public:
    Molecule(int charge, int multiplicity, const std::string& basisName, const std::vector<Atom>& geometry);

    /*
     * Returns a string representation of the molecule.
     * @return A string containing the charge, multiplicity, number of electrons, number of basis functions, and geometry.
     */
    std::string toString() const;

    /*
     * Calculates the total nuclear-nuclear repulsion energy.
     * @return The total nuclear repulsion energy.
     */
    double nuclearRepulsion() const;

    /*
     * Computes the overlap matrix of the molecule.
     * @return The overlap matrix S as an Eigen::MatrixXd.
     */
    Eigen::MatrixXd overlapMatrix() const;

    /*
     * Computes the kinetic energy matrix of the molecule.
     * @return The kinetic energy matrix T as an Eigen::MatrixXd.
     */

    /*
     * Computes the kinetic energy matrix of the molecule.
     * @return The kinetic energy matrix T as an Eigen::MatrixXd.
     */
    Eigen::MatrixXd kineticMatrix() const;

    /*
     * Computes the nuclear attraction matrix of the molecule.
     * @return The nuclear attraction matrix V as an Eigen::MatrixXd.
     */
    Eigen::MatrixXd nuclearAttractionMatrix() const;

    /*
     * Computes the electron repulsion integral tensor for the molecule.
     * @param threshold The Schwartz screening threshold below which integrals are considered negligible and set to zero.
     * @return A flattened 1D vector representing the electron repulsion tensor.
     */
    ElectronRepulsionTensor electronRepulsionTensor(double threshold = 1e-10) const;

    // Getters for molecule properties
    const std::vector<AtomicOrbital>& getAtomicOrbitals() const { return atomicOrbitals; }
    size_t getBasisFunctionCount() const { return basisFunctionCount; }
    int getCharge() const { return charge; }
    int getMultiplicity() const { return multiplicity; }
    size_t getElectronCount() const { return electronCount; }
    const std::vector<Atom>& getGeometry() const { return geometry; }

  private:
    size_t basisFunctionCount; // number of atomic orbitals in the molecule

    /* Counts the number of electrons in the molecule based on the geometry and charge
     * @return The total number of electrons in the molecule.
     */
    size_t countElectrons() const;

    /*
     * Builds the atomic orbitals for the molecule based on the provided basis set.
     * @param basisName The name of the basis set to use (case insensitive).
     * Throws std::runtime_error if the basis set is not found.
     */
    void buildBasis(const std::string& basisName);

    const int charge;                          // charge of the molecule
    const int multiplicity;                    // multiplicity of the molecule (2S + 1)
    const std::vector<Atom> geometry;          // geometry of the molecule, containing atomic numbers and coordinates
    const size_t electronCount;                // total number of electrons in the molecule
    std::vector<AtomicOrbital> atomicOrbitals; // vector of atomic orbitals built from the geometry and basis set
};
