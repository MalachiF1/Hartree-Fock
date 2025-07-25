#pragma once
#include "Molecule.hpp"

#include <Eigen/Dense>

class SCF
{
  public:
    SCF(const Molecule& molecule);

    /*
     * Runs the Self-Consistent Field (SCF) calculation.
     * @param max_iter Maximum number of iterations to perform.
     * @param energy_tol Convergence tolerance for electronic energy.
     * @param density_tol Convergence tolerance for density matrix.
     * @param schwartz_threshold Threshold for Schwartz screening of two-electron integrals.
     */
    void run(size_t maxIter = 50, double energyTol = 1e-8, double densityTol = 1e-8, double schwartzThreshold = 1e-10);

  private:
    // A constant reference to the molecule object.
    const Molecule& molecule;

    // SCF state variables that are constant throughout the calculation.
    size_t basisCount;       // number of basis functions
    size_t occupiedCount;    // number of occupied orbitals
    double nuclearEnergy;    // Nuclear repulsion energy
    Eigen::MatrixXd T;       // Kinetic energy matrix
    Eigen::MatrixXd V;       // Nuclear attraction matrix
    Eigen::MatrixXd h;       // Core Hamiltonian (T+V)
    Eigen::MatrixXd S;       // Overlap matrix
    Eigen::MatrixXd X;       // Orthogonalization matrix S^(-1/2)
    std::vector<double> Vee; // Two-electron repulsion integrals

    // SCF state variables that change during iterations.
    double electronicEnergy; // Electronic energy
    Eigen::MatrixXd D;       // Density matrix
    Eigen::MatrixXd F;       // Fock matrix
    Eigen::MatrixXd C;       // MO coefficient matrix

    /*
     * Initializes the SCF calculation by computing necessary matrices and parameters (one-time calculations).
     * This includes the overlap matrix, kinetic energy matrix, nuclear attraction matrix,
     * and two-electron repulsion integrals.
     * This method should be called before running the SCF iterations.
     * @note This method is called once when starting the SCF calculation.
     * @param schwartz_threshold The threshold for Schwartz screening of two-electron integrals.
     */
    void initialize(double schwartzThreshold = 1e-10);

    /*
     * Computes the initial guess for the density matrix and Fock matrix.
     * This is currently done using the core Hamiltonian (T + V).
     * @note This method is called once at the beginning of the SCF iterations.
     */
    void computeInitialGuessDensity();

    /*
     * Builds the Fock matrix using the current density matrix and two-electron integrals.
     */
    void buildFockMatrix();

    /*
     * Diagonalizes the Fock matrix to obtain the molecular orbitals and updates the density matrix and coefficient matrix.
     */
    void diagonalizeAndUpdate();

    /*
     * Prints the status of a single SCF iteration.
     * @param iter The current iteration number.
     * @param dE The change in electronic energy from the previous iteration.
     * @param dD The change in density matrix from the previous iteration.
     */
    void printIteration(int iter, double dE, double dD) const;

    /*
     * Prints the final results of the SCF calculation, including whether it converged and the final electronic energy.
     * @param converged True if the SCF calculation converged, false otherwise.
     */
    void printFinalResults(bool converged) const;
};
