#include "Utils.hpp"

#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <cmath>

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
    if (T < 1e-10)
    {
        return 1.0 / (2.0 * m + 1.0);
    }
    return boost::math::hypergeometric_1F1(m + 0.5, m + 1.5, -T) / (2.0 * m + 1.0);
}
