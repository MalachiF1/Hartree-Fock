#pragma once
#include <Eigen/Dense>

class DIIS
{
  public:
    DIIS(size_t maxSize = 8);

    /*
     * Updates the DIIS storage with the latest error and Fock matrix (in orthogonal basis).
     *
     * @param F The current Fock matrix.
     * @param D The current density matrix.
     * @param S The overlap matrix.
     * @param X The orthogonalization matrix (S^(-1/2)).
     */
    void update(const Eigen::MatrixXd& F, const Eigen::MatrixXd& D, const Eigen::MatrixXd& S, const Eigen::MatrixXd& X);

    /*
     * Extrapolates the Fock matrix using the stored error vectors and Fock matrices.
     *
     * @param orthoF The current Fock matrix in the orthogonalized basis.
     * @return The extrapolated Fock matrix.
     */
    Eigen::MatrixXd extrapolate(const Eigen::MatrixXd& orthoF);

    /*
     * Returns the current error norm, which is the norm of the last error vector.
     *
     * @return The norm of the last error vector.
     */
    double getErrorNorm() const;

  private:
    size_t maxSize;                                 // Maximum number of error vectors to store
    std::vector<Eigen::MatrixXd> errorVectors;      // History of error vectors
    std::vector<Eigen::MatrixXd> fockMatricesOrtho; // History of Fock matrices in orthogonalized basis

    /*
     * Computes the coefficients for the DIIS extrapolation.
     *
     * @return The coefficients for the DIIS extrapolation.
     */
    Eigen::VectorXd computeCoefficients();

    /*
     * Builds the B (Pulay) matrix required for the DIIS linear system.
     *
     * @return The B matrix for the DIIS extrapolation.
     */
    Eigen::MatrixXd buildBMatrix() const;
};
