#pragma once

#include "AtomicOrbital.hpp"
#include "Utils.hpp"

#include <string>
#include <vector>

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
    std::vector<double> electronRepulsionTensor(double threshold = 1e-12) const;

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
