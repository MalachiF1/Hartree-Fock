#pragma once
#include "DIIS.hpp"
#include "Molecule.hpp"

#include <Eigen/Dense>

/**
 * @brief Class for performing Self-Consistent Field (SCF) calculations on a molecule.
 *
 * This class performs the SCF calculation to find the ground state electronic structure of a molecule.
 * Currently, it only implements the restricted Hartree-Fock (RHF) method.
 */
class SCF
{
  public:
    SCF(const Molecule& molecule);

    /**
     * Runs the Self-Consistent Field (SCF) calculation.
     *
     * @param max_iter Maximum number of iterations to perform.
     * @param energy_tol Convergence tolerance for electronic energy.
     * @param density_tol Convergence tolerance for density matrix.
     * @param schwartz_threshold Threshold for Schwartz screening of two-electron integrals.
     * @param use_DIIS Whether to use Direct Inversion in the Iterative Subspace (DIIS) method for convergence acceleration.
     * @param DIISmaxSize Maximum size of the DIIS error vector (and number of previous Fock matrices to extrapolate from).
     * @param DIISstart The iteration number at which to start using DIIS (error vectors/Fock matrices are stored prior to DIISstart
     *        but,the extrapolation only starts at the DIISstart itteration.
     * @param DIISErrorTol Tolerance for the DIIS error norm to determine convergence.
     */
    void run(
        size_t maxIter           = 50,
        double energyTol         = 1e-8,
        double densityTol        = 1e-8,
        double schwartzThreshold = 1e-10,
        bool useDIIS             = true,
        size_t DIISmaxSize       = 8,
        unsigned DIISstart       = 1,
        double DIISErrorTol      = 1e-8
    );

  private:
    // A constant reference to the molecule object.
    const Molecule& molecule;

    // SCF state variables that are constant throughout the calculation.
    size_t basisCount;           // number of basis functions
    size_t occupiedCount;        // number of occupied orbitals
    double nuclearEnergy;        // Nuclear repulsion energy
    Eigen::MatrixXd T;           // Kinetic energy matrix
    Eigen::MatrixXd V;           // Nuclear attraction matrix
    Eigen::MatrixXd h;           // Core Hamiltonian (T+V)
    Eigen::MatrixXd S;           // Overlap matrix
    Eigen::MatrixXd X;           // Orthogonalization matrix S^(-1/2)
    ElectronRepulsionTensor Vee; // Two-electron repulsion integrals

    // SCF state variables that change during iterations.
    double electronicEnergy; // Electronic energy
    Eigen::MatrixXd D;       // Density matrix
    Eigen::MatrixXd F;       // Fock matrix
    Eigen::MatrixXd C;       // MO coefficient matrix

    // Pointer to the DIIS object (if used).
    std::unique_ptr<DIIS> diis_handler;

    /**
     * Initializes the SCF calculation by computing necessary matrices and parameters (one-time calculations).
     * This includes the overlap matrix, kinetic energy matrix, nuclear attraction matrix, and two-electron
     * repulsion integrals. This method should be called before running the SCF iterations.
     *
     * @note This method is called once when starting the SCF calculation.
     * @param schwartz_threshold The threshold for Schwartz screening of two-electron integrals.
     */
    void initialize(double schwartzThreshold = 1e-10);

    /**
     * Computes the initial guess for the density matrix and Fock matrix. This is currently done using the core Hamiltonian (T + V).
     *
     * @note This method is called once at the beginning of the SCF iterations.
     */
    void computeInitialGuessDensity();

    /**
     * Builds the Fock matrix using the current density matrix and two-electron integrals.
     */
    void buildFockMatrix();

    /**
     * Diagonalizes the Fock matrix to obtain the molecular orbitals and updates the density matrix and coefficient matrix.
     */
    void diagonalizeAndUpdate();
    void diagonalizeAndUpdate(const Eigen::MatrixXd& F_prime); // overload for DIIS

    /**
     * Prints the status of a single SCF iteration.
     *
     * @param iter The current iteration number.
     * @param dE The change in electronic energy from the previous iteration.
     * @param dD The change in density matrix from the previous iteration.
     */
    void printIteration(int iter, double dE, double dD) const;

    /**
     * Prints the status of a single SCF iteration.
     *
     * @param iter The current iteration number.
     * @param dE The change in electronic energy from the previous iteration.
     * @param dD The change in density matrix from the previous iteration.
     * @param DIISError The DIIS error norm for the current iteration.
     */
    void printIteration(int iter, double dE, double dD, double DIISError) const; // overload for DIIS

    /**
     * Prints the final results of the SCF calculation, including whether it converged and the final electronic energy.
     *
     * @param converged True if the SCF calculation converged, false otherwise.
     */
    void printFinalResults(bool converged) const;
};
