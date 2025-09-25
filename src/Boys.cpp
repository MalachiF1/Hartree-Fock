#include "Boys.hpp"

#include "Utils.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

/**
 * This implementation is based on:
 * Weiss, A. K., & Ochsenfeld, C. (2015). A rigorous and optimized strategy for the evaluation of the B oys function
 * kernel in molecular electronic structure theory. Journal of Computational Chemistry, 36(18), 1390-1398.
 */
void Boys::calculateBoys(unsigned m_max, double T, std::span<double> F)
{
    if (T > 30.0)
    {
        // For large T, use the asymptotic expansion.
        F[0] = SQRT_PI_BY_2 / std::sqrt(T);
        for (unsigned m = 0; m < m_max; ++m) { F[m + 1] = (m + 0.5) * F[m] / T; }
        return;
    }

    const auto& gridManager    = GridManager::getInstance();
    bool useDownwardRecursion  = (T <= TRANSITION_VALUES[m_max]);
    auto gridIndex             = static_cast<size_t>(T / GRID_STEP_DELTA);
    double T_grid              = gridIndex * GRID_STEP_DELTA;
    double T_diff              = T - T_grid;
    double half_T_diff_squared = 0.5 * T_diff * T_diff;

    if (m_max == 0 || m_max == 1)
    {
        F[0] = gridManager.getValue(0, gridIndex) - (T_diff * gridManager.getValue(1, gridIndex))
             + (half_T_diff_squared * gridManager.getValue(2, gridIndex));

        if (m_max == 0)
            return;

        F[1] = gridManager.getValue(1, gridIndex) - (T_diff * gridManager.getValue(2, gridIndex))
             + (half_T_diff_squared * gridManager.getValue(3, gridIndex));
        return;
    }

    unsigned startOrder = useDownwardRecursion ? m_max : 0;
    double F_start      = gridManager.getValue(startOrder, gridIndex)
                   - (T_diff * gridManager.getValue(startOrder + 1, gridIndex))
                   + (half_T_diff_squared * gridManager.getValue(startOrder + 2, gridIndex));
    F[startOrder] = F_start;

    // double expMinusT;
    // if (T > EXP_CUTOFF)
    //     expMinusT = schraudolphExp(-T);
    // else
    // expMinusT = std::exp(-T);

    double expMinusT = std::exp(-T);

    if (useDownwardRecursion)
    {
        for (unsigned m = m_max; m >= 1; --m) { F[m - 1] = (2 * T * F[m] + expMinusT) / (2 * (m - 1) + 1); }
    }
    else
    {
        for (unsigned m = 0; m < m_max; ++m) { F[m + 1] = ((2.0 * m + 1.0) * F[m] - expMinusT) / (2.0 * T); }
    }
}


double Boys::analyticalBoys(unsigned m, double T)
{
    if (T < 1e-9)
    {
        return 1.0 / (2.0 * m + 1.0);
    }

    double sum  = 0.0;
    double term = 1.0 / (m + 0.5); // Initial term for the iterative series.
    sum += term;

    for (int i = 1; i < 200; ++i)
    {
        term *= T / (m + i + 0.5);
        sum += term;
        if (std::abs(term) < 1e-15) // Convergence check.
            break;
    }

    return sum * std::exp(-T) / 2.0;
}


GridManager::GridManager()
{
    const size_t gridMaxOrder = MAX_L + 2;

    const size_t totalSize = static_cast<size_t>(gridMaxOrder + 1) * numGridPoints;
    grid.resize(totalSize);

    std::ifstream in("boys_grid.data", std::ios_base::binary);
    if (!in.read((char*)grid.data(), totalSize * sizeof(double)))
    {
        for (size_t i = 0; i <= gridMaxOrder; ++i)
        {
            for (size_t j = 0; j < numGridPoints; ++j)
            {
                double T_j                  = j * GRID_STEP_DELTA;
                grid[i * numGridPoints + j] = Boys::analyticalBoys(i, T_j);
            }
        }

        std::ofstream out("boys_grid.data", std::ios_base::binary);
        out.write((char*)grid.data(), totalSize * sizeof(double));
    }
}
