#pragma once
#include <Eigen/Dense>

class DIIS
{
  public:
    DIIS(size_t maxSize = 8);

    /*
     * Updates the DIIS storage with the latest error and Fock matrix (in orthogonal basis).
     *
     * @param F_alpha The current alpha Fock matrix.
     * @param F_beta The current beta Fock matrix.
     * @param errorVector The current error vector.k
     */
    void update(const Eigen::MatrixXd& F_alpha, const Eigen::MatrixXd& F_beta, const Eigen::MatrixXd& errorVector);

    /*
     * Extrapolates the Fock matrix using the stored error vectors and Fock matrices.
     *
     * @param F_alpha The current alpha Fock matrix.
     * @param F_beta The current beta Fock matrix.
     * @return A pair containing the new extrapolated alpha and beta Fock matrices.
     */
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> extrapolate(const Eigen::MatrixXd& F_alpha, const Eigen::MatrixXd& F_beta);

    /*
     * Returns the current error norm, which is the norm of the last error vector.
     *
     * @return The norm of the last error vector.
     */
    double getErrorNorm() const;

  private:
    size_t maxSize;                                 // Maximum number of error vectors to store
    std::vector<Eigen::MatrixXd> errorVectors;      // History of error vectors
    std::vector<Eigen::MatrixXd> alphaFockMatrices; // History of Fock matrices
    std::vector<Eigen::MatrixXd> betaFockMatrices;

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
