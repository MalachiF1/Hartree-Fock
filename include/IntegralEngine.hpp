#pragma once
#include "Basis.hpp"

#include <span>

class IntegralEngine
{
  public:
    static void overlap(const Basis& basis, const Shell& shellA, const Shell& shellB, Eigen::MatrixXd& S);
    static void kinetic(const Basis& basis, const Shell& shellA, const Shell& shellB, Eigen::MatrixXd& T);

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
     * @brief Computes the Boys function F_m(T).
     */
    static double boys(unsigned m, double T);

    struct EBuffer
    {
        size_t l1;
        size_t l2;
        size_t max_t;
        std::vector<double> data;
        EBuffer() = default;
        EBuffer(unsigned l1, unsigned l2) : l1(l1), l2(l2), max_t(l1 + l2), data((l1 + 1) * (l2 + 1) * (max_t + 2) / 2)
        {
        }

        double& operator()(unsigned i, unsigned j, unsigned t)
        {
            size_t offset_i = (i * (l2 + 1) * (i + l2 + 1)) / 2;
            size_t offset_j = (j * (i + 1)) + ((j * (j - 1)) / 2);
            size_t index    = offset_i + offset_j + t;
            assert(index < data.size());
            return data[index];
        }
    };


    /**
     * @brief Computes 1D Hermite expansion coefficients E_k(i, j, a, b).
     */
    static void computeHermiteCoeffs(unsigned l1, unsigned l2, double p, double PA, double PB, EBuffer& E_buffer);


    struct RBuffer
    {
        size_t max_t;
        size_t max_u;
        size_t max_v;
        size_t max_n;
        std::vector<double> data;

        RBuffer() = default;
        RBuffer(size_t max_t, size_t max_u, size_t max_v) :
            max_t(max_t),
            max_u(max_u),
            max_v(max_v),
            max_n(max_t + max_u + max_v + 1),
            data((max_t + 1) * (max_u + 1) * (max_v + 1) * max_n)
        {
        }

        double& operator()(size_t t, size_t u, size_t v, size_t n)
        {
            size_t index = n + (v * max_n) + (u * max_n * (max_v + 1)) + (t * (max_u + 1) * (max_v + 1) * max_n);
            assert(index < data.size());
            return data[index];
        }

        void clear() { std::ranges::fill(data, -1.0); }
    };

    static void computeAuxiliaryIntegrals(
        unsigned t_max, unsigned u_max, unsigned v_max, double p, const Vec3& PC, std::span<double> F, RBuffer& R_buffer
    );
};
