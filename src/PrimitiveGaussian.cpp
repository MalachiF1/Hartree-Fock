#include "PrimitiveGaussian.hpp"

#include "Utils.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

PrimitiveGaussian::PrimitiveGaussian(double exponent, double coeff, const Vec3& coords, const Eigen::Vector3i& angularMomentum) :
    exponent(exponent),
    coeff(coeff),
    coords(coords),
    angularMomentum(angularMomentum),
    normConst(normalizationConst(exponent, angularMomentum))
{
}

double PrimitiveGaussian::normalizationConst(double exponent, const Eigen::Vector3i& angularMomentum)
{
    int l = angularMomentum.x();
    int m = angularMomentum.y();
    int n = angularMomentum.z();

    // Calculate the normalization constant using the formula:
    // N = (2 * alpha / pi)^(3/4) * (4 * alpha)^(l + m + n) / sqrt((2l - 1)!! * (2m - 1)!! * (2n - 1)!!)
    double prefactor1  = std::pow(2.0 * exponent / M_PI, 0.75);
    double prefactor2  = std::pow(4.0 * exponent, (l + m + n) / 2.0);
    double denominator = std::sqrt(dfact((2 * l) - 1) * dfact((2 * m) - 1) * dfact((2 * n) - 1));
    return prefactor1 * prefactor2 / denominator;
}

std::string PrimitiveGaussian::toString() const
{
    std::stringstream ss;
    ss << std::fixed << std::setprecision(8);
    ss << "PrimitiveGaussian {\n"
       << "\talpha:    " << this->exponent << "\n"
       << "\tcoeff:    " << this->coeff << "\n"
       << "\tcoords:   (" << this->coords.x() << ", " << this->coords.y() << ", " << this->coords.z() << ")\n"
       << "\t(l,m,n):  (" << this->angularMomentum.x() << ", " << this->angularMomentum.y() << ", "
       << this->angularMomentum.z() << ")\n"
       << "\tnorm:     " << this->normConst << "\n"
       << "}";
    return ss.str();
}
