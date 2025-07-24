#include "PrimitiveGaussian.hpp"
#include "Utils.hpp"

#include <iostream>

int main()
{
    PrimitiveGaussian pg(0.5, 1.0, Vec3(1.0, 2.0, 3.0), 1, 0, 2);
    std::cout << pg.toString() << std::endl;
    return 0;
}
