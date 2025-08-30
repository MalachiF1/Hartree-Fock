#pragma once

#include "Utils.hpp"

#include <string>

/**
 * @brief Represents a cartesian primitive Gaussian function
 *
 * This class represents a cartesian primitive Gaussian: $ N C x^l y^m z^n e^{-\alpha (r - A)^2} $
 * where: N is the normalization constant, C is the contraction coefficient, l, m, n are the angular
 * momentum numbers in the x, y, z directions respectively, $\alpha$ is the exponent, and $A$ is the
 * center of the Gaussian. This class is used to store the parameters of the primitive Gaussians of
 * an atomic orbital in the AtomicOrbital class.
 */
class PrimitiveGaussian
{
  public:
    double exponent;                 // Exponent
    double coeff;                    // Contraction coefficient
    Vec3 coords;                     // Center coordinates
    Eigen::Vector3i angularMomentum; // Angular momentum numbers (x, y, z)
    double normConst;                // Normalization constant

    PrimitiveGaussian(double exponent, double coeff, const Vec3& coords, const Eigen::Vector3i& angularMomentum);

    /**
     * Returns a string representation of the primitive Gaussian.
     * @return A string containing the exponent, coefficient, coordinates, angular momentum numbers, and normalization constant.
     */
    std::string toString() const;

    std::partial_ordering operator<=>(const PrimitiveGaussian& other) const;
    bool operator==(const PrimitiveGaussian& other) const = default;

  private:
    /**
     * Computes the normalization constant of a primitive gaussian according to the formula:
     * $ N = (2 * alpha / pi)^(3/4) * (4 * alpha)^(l + m + n) / sqrt((2l - 1)!! * (2m - 1)!! * (2n - 1)!!) $
     *
     * @param exponent The exponent of the gaussian
     * @param l Angular momentum number in the x direction
     * @param m Angular momentum number in the y direction
     * @param n Angular momentum number in the z direction
     * @return The nomralization constant of a primitive gaussian with the specified exponent and angular momentum
     */
    static double normalizationConst(double exponent, const Eigen::Vector3i& angularMomentum);
};
