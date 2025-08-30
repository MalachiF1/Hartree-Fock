#include "PrimitiveGaussian.hpp"

#include "Utils.hpp"

#include <cmath>
#include <fmt/core.h>
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
    return fmt::format(
        "{{\n"
        "\talpha:    {:.8f}\n"
        "\tcoeff:    {:.8f}\n"
        "\tcoords:   ({:.8f}, {:.8f}, {:.8f})\n"
        "\t(l,m,n):  ({}, {}, {})\n"
        "\tnorm:     {:.8f}\n"
        "}}",
        this->exponent,
        this->coeff,
        this->coords.x(),
        this->coords.y(),
        this->coords.z(),
        this->angularMomentum.x(),
        this->angularMomentum.y(),
        this->angularMomentum.z(),
        this->normConst
    );
}

std::partial_ordering PrimitiveGaussian::operator<=>(const PrimitiveGaussian& other) const
{
    if (auto cmp = this->exponent <=> other.exponent; cmp != 0)
        return cmp;
    if (auto cmp = this->coeff <=> other.coeff; cmp != 0)
        return cmp;
    if (auto cmp = this->angularMomentum.x() <=> other.angularMomentum.x(); cmp != 0)
        return cmp;
    if (auto cmp = this->angularMomentum.y() <=> other.angularMomentum.y(); cmp != 0)
        return cmp;
    if (auto cmp = this->angularMomentum.z() <=> other.angularMomentum.z(); cmp != 0)
        return cmp;
    if (auto cmp = this->normConst <=> other.normConst; cmp != 0)
        return cmp;

    return std::partial_ordering::equivalent;
}
