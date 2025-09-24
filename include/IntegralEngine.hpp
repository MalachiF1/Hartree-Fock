#pragma once
#include "Basis.hpp"
#include "Molecule.hpp"

#include <span>

class IntegralEngine
{
  public:
    static void overlap(const Basis& basis, const Shell& shellA, const Shell& shellB, Eigen::MatrixXd& S);
    static void kinetic(const Basis& basis, const Shell& shellA, const Shell& shellB, Eigen::MatrixXd& T);
    static void nuclearAttraction(
        const std::vector<Atom>& geometry, const Basis& basis, const Shell& shellA, const Shell& shellB, Eigen::MatrixXd& V
    );

    static void electronRepulsion(
        const Basis& basis,
        const Shell& shellA,
        const Shell& shellB,
        const Shell& shellC,
        const Shell& shellD,
        ElectronRepulsionTensor& G
    );

  private:
    /**
     * @brief Holds precomputed values for a primitive Gaussian pair (a, b).
     */
    struct PrimitivePairData
    {
        double p; // Exponent sum: alpha_a + alpha_b
        double K; // Pre-exponential factor: exp(- (alpha_a*alpha_b/p) * |A-B|^2 )

        // New center: (alpha_a*A + alpha_b*B) / p
        double Px;
        double Py;
        double Pz;

        // Vector from new center to old center A: P - A
        double PAx;
        double PAy;
        double PAz;

        // Vector from new center to old center B: P - B
        double PBx;
        double PBy;
        double PBz;
    };

    /**
     * @brief Computes pre-calculated quantities for a pair of primitive Gaussians.
     */
    static PrimitivePairData computePrimitivePairData(
        double alpha1, double xa, double ya, double za, double alpha2, double xb, double yb, double zb
    );

    // /**
    //  * @brief Computes the Boys function F_m(T).
    //  */
    // static void calculateBoys(unsigned m_max, double T, std::span<double> F);

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
            const size_t offset_i = i * (l2 + 1) * (i + l2 + 1) / 2;
            const size_t offset_j = (j * (i + 1)) + j * (j - 1) / 2;
            const size_t index    = offset_i + offset_j + t;
            assert(index < data.size());
            return data[index];
        }

        const double& operator()(unsigned i, unsigned j, unsigned t) const
        {
            const size_t offset_i = i * (l2 + 1) * (i + l2 + 1) / 2;
            const size_t offset_j = (j * (i + 1)) + j * (j - 1) / 2;
            const size_t index    = offset_i + offset_j + t;
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
            // const size_t index = n + (v * max_n) + (u * max_n * (max_v + 1)) + (t * (max_u + 1) * (max_v + 1) * max_n);
            const size_t index = t + u * (max_t + 1) + v * (max_t + 1) * (max_u + 1)
                               + n * (max_t + 1) * (max_u + 1) * (max_v + 1);
            assert(index < data.size());
            return data[index];
        }

        const double& operator()(size_t t, size_t u, size_t v, size_t n) const
        {
            const size_t index = t + u * (max_t + 1) + v * (max_t + 1) * (max_u + 1)
                               + n * (max_t + 1) * (max_u + 1) * (max_v + 1);
            assert(index < data.size());
            return data[index];
        }


        void clear() { std::ranges::fill(data, -1.0); }
    };

    static void computeAuxiliaryIntegrals(
        unsigned t_max, unsigned u_max, unsigned v_max, double p, const Vec3& PC, std::span<double> F, RBuffer& R_buffer
    );
};
