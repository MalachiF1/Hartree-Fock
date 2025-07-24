#pragma once

#include "Utils.hpp"

#include <string>

// Represents a single primitive Gaussian function.
class PrimitiveGaussian
{
  public:
    const double exponent;  // Exponent
    const double coeff;     // Contraction coefficient
    const Vec3 coords;      // Center coordinates
    const int l, m, n;      // Angular momentum numbers (x, y, z)
    const double normConst; // Normalization constant

    PrimitiveGaussian(double exponent, double coeff, const Vec3& coords, int l, int m, int n);

    std::string toString() const;

  private:
    /**
     * Computes the normalization constant of a primitive gaussian.
     *
     * @param exponent The exponent of the gaussian
     * @param l Angular momentum number in the x direction
     * @param m Angular momentum number in the y direction
     * @param n Angular momentum number in the z direction
     * @return The nomralization constant of a primitive gaussian with the specified exponent and angular momentum
     */
    static double normalizationConst(double exponent, int l, int m, int n);
};
