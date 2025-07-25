#include "PrimitiveGaussian.hpp"

#include "Utils.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

PrimitiveGaussian::PrimitiveGaussian(double exponent, double coeff, const Vec3& coords, int l, int m, int n) :
    exponent(exponent), coeff(coeff), coords(coords), l(l), m(m), n(n), normConst(normalizationConst(exponent, l, m, n))
{
}

double PrimitiveGaussian::normalizationConst(double exponent, int l, int m, int n)
{
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
       << "\t(l,m,n):  (" << this->l << ", " << this->m << ", " << this->n << ")\n"
       << "\tnorm:     " << this->normConst << "\n"
       << "}";
    return ss.str();
}
