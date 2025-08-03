#include "Utils.hpp"

#include <cmath>
#include <unordered_map>

int dfact(int n)
{
    if (n == -1)
    {
        return 1;
    }
    else if (n < 0)
    {
        throw std::invalid_argument("Input must be a non-negative odd integer or -1.");
    }
    int result = 1;
    for (int i = n; i > 0; i -= 2) { result *= i; }
    return result;
}

double boys(int m, double T)
{
    double result;

    if (T < 1e-9)
    {
        // For small T, use the series expansion for stability.
        // F_m(T) = 1/(2m+1) - T/(2m+3) + ...
        result = 1.0 / (2.0 * m + 1.0);
    }
    else
    {
        // Calculate the base case, F_0(T), using the error function.
        // F_0(T) = (sqrt(pi) / 2*sqrt(T)) * erf(sqrt(T))
        double sqrt_T = std::sqrt(T);
        double F0     = (std::sqrt(M_PI) / (2.0 * sqrt_T)) * std::erf(sqrt_T);

        if (m == 0)
        {
            result = F0;
        }
        else
        {
            // Use the upward recurrence relation to get F_m(T) from F_0(T).
            // F_{i+1}(T) = ( (2i+1)F_i(T) - exp(-T) ) / (2T)
            double F_i = F0; // Start with F_0
            for (int i = 0; i < m; ++i)
            {
                double F_i_plus_1 = ((2.0 * i + 1.0) * F_i - std::exp(-T)) / (2.0 * T);
                F_i               = F_i_plus_1;
            }
            result = F_i;
        }
    }

    return result;
}


Eigen::MatrixXd inverseSqrtMatrix(const Eigen::MatrixXd& S)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S);
    const Eigen::VectorXd& eigenvalues  = solver.eigenvalues();
    const Eigen::MatrixXd& eigenvectors = solver.eigenvectors();

    Eigen::VectorXd sqrtInvEigenvalues(eigenvalues.size());
    for (int i = 0; i < eigenvalues.size(); ++i)
    {
        if (eigenvalues(i) > 0)
        {
            sqrtInvEigenvalues(i) = 1.0 / std::sqrt(eigenvalues(i));
        }
        else
        {
            // Handle very small or zero eigenvalues if they occur
            sqrtInvEigenvalues(i) = 0;
        }
    }

    return eigenvectors * sqrtInvEigenvalues.asDiagonal() * eigenvectors.transpose();
}
