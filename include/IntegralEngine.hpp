#pragma once
#include "Basis.hpp"

#include <span>

class IntegralEngine
{
  public:
    static void overlap(const Basis& basis, const Shell& shellA, const Shell& shellB, Eigen::MatrixXd& S);

  private:
    /**
     * @brief Holds precomputed values for a primitive Gaussian pair (a, b).
     */
    struct PrimitivePairData
    {
        double p;           // Exponent sum: alpha_a + alpha_b
        double K;           // Pre-exponential factor: exp(- (alpha_a*alpha_b/p) * |A-B|^2 )
        Eigen::Vector3d P;  // New center: (alpha_a*A + alpha_b*B) / p
        Eigen::Vector3d PA; // Vector from new center to old center A: P - A
        Eigen::Vector3d PB; // Vector from new center to old center B: P - B
    };

    /**
     * @brief Computes pre-calculated quantities for a pair of primitive Gaussians.
     */
    static PrimitivePairData computePrimitivePairData(
        double alpha1, const Eigen::Vector3d& A, double alpha2, const Eigen::Vector3d& B
    );

    /**
     * @brief Computes 1D Hermite expansion coefficients E_k(i, j, a, b).
     */
    static void computeHermiteCoeffs(unsigned l1, unsigned l2, double PA, double PB, std::span<double> E);

    /**
     * @brief Computes the 1D overlap integral for a given angular momentum pair.
     */
    static double compute1dOverlap(
        unsigned l1, unsigned l2, double PA, double PB, const std::vector<double>& I, std::span<double> E_buffer
    );

    /**
     * @brief Computes the 3D overlap integral between two primitive Gaussians.
     */
    static double computePrimitiveOverlap(
        unsigned l1,
        unsigned m1,
        unsigned n1,
        unsigned l2,
        unsigned m2,
        unsigned n2,
        const PrimitivePairData& primPair,
        const std::vector<double>& I,
        std::span<double> E_buffer
    );
};
