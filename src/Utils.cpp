#include "Utils.hpp"

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
