#include "DIIS.hpp"

DIIS::DIIS(size_t maxSize) : maxSize(maxSize) {}


void DIIS::update(const Eigen::MatrixXd& F, const Eigen::MatrixXd& D, const Eigen::MatrixXd& S, const Eigen::MatrixXd& X)
{
    // Add the new error vector and Fock matrix to the history.
    errorVectors.emplace_back(X.transpose() * (F * D * S - S * D * F) * X);
    fockMatricesOrtho.emplace_back(X.transpose() * F * X);

    // If the storage exceeds the maximum size, remove the oldest entry.
    if (errorVectors.size() > maxSize)
    {
        errorVectors.erase(errorVectors.begin());
        fockMatricesOrtho.erase(fockMatricesOrtho.begin());
    }
}


Eigen::MatrixXd DIIS::extrapolate(const Eigen::MatrixXd& orthoF)
{
    if (errorVectors.empty())
    {
        return orthoF; // Cannot extrapolate without history.
    }

    Eigen::MatrixXd B = buildBMatrix();
    Eigen::VectorXd b = Eigen::VectorXd::Zero(B.rows());
    b(B.rows() - 1)   = -1.0;

    // Solve the linear system Bc = b for the DIIS coefficients.
    Eigen::VectorXd coeffs = B.colPivHouseholderQr().solve(b);

    // Use the coefficients to form the new extrapolated Fock matrix.
    Eigen::MatrixXd newF = Eigen::MatrixXd::Zero(orthoF.rows(), orthoF.cols());
    for (size_t i = 0; i < fockMatricesOrtho.size(); ++i) { newF += coeffs(i) * fockMatricesOrtho[i]; }

    return newF;
}


Eigen::MatrixXd DIIS::buildBMatrix() const
{
    size_t size       = errorVectors.size();
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(size + 1, size + 1);

    for (size_t i = 0; i < size; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            // B_ij = <e_i | e_j>
            B(i, j) = B(j, i) = (errorVectors[i].array() * errorVectors[j].array()).sum();
        }
    }

    // The last row and column are -1.
    B.row(size).setConstant(-1.0);
    B.col(size).setConstant(-1.0);
    B(size, size) = 0.0; // The bottom-right element is 0.

    return B;
}


double DIIS::getErrorNorm() const
{
    if (errorVectors.empty())
    {
        return 0.0; // No error to report.
    }
    return errorVectors.back().norm();
}
