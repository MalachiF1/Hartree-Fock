#include "DIIS.hpp"

DIIS::DIIS(size_t maxSize) : maxSize(maxSize)
{
    errorVectors.reserve(maxSize);
    fockMatrices.reserve(maxSize);
}


void DIIS::update(const Eigen::MatrixXd& F, const Eigen::MatrixXd& errorVector)
{
    // If the storage exceeds the maximum size, remove the oldest entry.
    if (errorVectors.size() >= maxSize)
    {
        errorVectors.erase(errorVectors.begin());
        fockMatrices.erase(fockMatrices.begin());
    }

    // Add the new error vector and Fock matrix to the history.
    errorVectors.emplace_back(errorVector);
    fockMatrices.emplace_back(F);
}


Eigen::MatrixXd DIIS::extrapolate(const Eigen::MatrixXd& F)
{
    if (errorVectors.empty())
        return F; // Cannot extrapolate without history.

    Eigen::MatrixXd B = buildBMatrix();
    Eigen::VectorXd b = Eigen::VectorXd::Zero(B.rows());
    b(B.rows() - 1)   = -1.0;

    // Solve the linear system Bc = b for the DIIS coefficients.
    Eigen::VectorXd coeffs = B.colPivHouseholderQr().solve(b);

    // Use the coefficients to form the new extrapolated Fock matrix.
    Eigen::MatrixXd newF = Eigen::MatrixXd::Zero(F.rows(), F.cols());
    for (size_t i = 0; i < fockMatrices.size(); ++i) { newF += coeffs(i) * fockMatrices[i]; }

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
        return 0.0;

    return errorVectors.back().norm();
}
