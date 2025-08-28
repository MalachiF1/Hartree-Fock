#include "DIIS.hpp"

DIIS::DIIS(size_t maxSize) : maxSize(maxSize)
{
    errorVectors.reserve(maxSize);
    alphaFockMatrices.reserve(maxSize);
    betaFockMatrices.reserve(maxSize);
}


void DIIS::update(const Eigen::MatrixXd& F_alpha, const Eigen::MatrixXd& F_beta, const Eigen::MatrixXd& errorVector)
{
    // If the storage exceeds the maximum size, remove the oldest entry.
    if (errorVectors.size() >= maxSize)
    {
        errorVectors.erase(errorVectors.begin());
        alphaFockMatrices.erase(alphaFockMatrices.begin());
        betaFockMatrices.erase(betaFockMatrices.begin());
    }

    // Add the new error vector and Fock matrix to the history.
    errorVectors.emplace_back(errorVector);
    alphaFockMatrices.emplace_back(F_alpha);
    betaFockMatrices.emplace_back(F_beta);
}


std::pair<Eigen::MatrixXd, Eigen::MatrixXd> DIIS::extrapolate(const Eigen::MatrixXd& F_alpha, const Eigen::MatrixXd& F_beta)
{
    if (errorVectors.empty())
        return {F_alpha, F_beta}; // Cannot extrapolate without history.

    Eigen::MatrixXd B = buildBMatrix();
    Eigen::VectorXd b = Eigen::VectorXd::Zero(B.rows());
    b(B.rows() - 1)   = -1.0;

    // Solve the linear system Bc = b for the DIIS coefficients.
    Eigen::VectorXd coeffs = B.colPivHouseholderQr().solve(b);

    // Use the coefficients to form the new extrapolated Fock matrix.
    Eigen::MatrixXd newF_alpha = Eigen::MatrixXd::Zero(F_alpha.rows(), F_alpha.cols());
    Eigen::MatrixXd newF_beta  = Eigen::MatrixXd::Zero(F_beta.rows(), F_beta.cols());
    for (size_t i = 0; i < alphaFockMatrices.size(); ++i)
    {
        newF_alpha += coeffs(i) * alphaFockMatrices[i];
        newF_beta += coeffs(i) * betaFockMatrices[i];
    }
    return {newF_alpha, newF_beta};
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
