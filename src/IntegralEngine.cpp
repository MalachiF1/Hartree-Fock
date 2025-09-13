#include "IntegralEngine.hpp"

#include <cmath>
#include <numeric>
#include <span>
#include <vector>

void IntegralEngine::overlap(const Basis& basis, const Shell& shellA, const Shell& shellB, Eigen::MatrixXd& S)
{
    const auto& exps   = basis.getExponents();
    const auto& coeffs = basis.getCoefficients();
    const auto& lx     = basis.getLx();
    const auto& ly     = basis.getLy();
    const auto& lz     = basis.getLz();

    const size_t nprimA = shellA.nprim;
    const size_t nprimB = shellB.nprim;
    const size_t naoA   = shellA.nao;
    const size_t naoB   = shellB.nao;

    const Eigen::Vector3d& centerA = shellA.center;
    const Eigen::Vector3d& centerB = shellB.center;

    // Allocate scratch space for Hermite coefficients once per shell pair.
    const unsigned max_L = shellA.l + shellB.l;
    std::vector<double> E_buffer(max_L + 1);

    // Loop over primitive pairs
    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double alphaA = exps[shellA.primOffset + pA];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double alphaB = exps[shellB.primOffset + pB];

            // Compute data for this primitive pair
            const PrimitivePairData primPair = computePrimitivePairData(alphaA, centerA, alphaB, centerB);

            // Pre-calculate fundamental integrals I_t for this primitive pair.
            std::vector<double> I(max_L + 1, 0.0);
            const double I0 = std::sqrt(M_PI / primPair.p);
            I[0]            = I0;
            if (max_L > 0)
            {
                double I_prev = I0;
                for (size_t t = 2; t <= max_L; t += 2)
                {
                    I[t]   = (t - 1.0) / (2.0 * primPair.p) * I_prev;
                    I_prev = I[t];
                }
            }

            // This primitive pair contributes to all AO pairs in the shell pair
            for (size_t a = 0; a < naoA; ++a)
            {
                const double coeffA = coeffs[(shellA.coeffOffset + pA * naoA) + a];

                for (size_t b = 0; b < naoB; ++b)
                {
                    const double coeffB = coeffs[(shellB.coeffOffset + pB * naoB) + b];

                    const double S_prim = computePrimitiveOverlap(
                        lx[shellA.aoOffset + a],
                        ly[shellA.aoOffset + a],
                        lz[shellA.aoOffset + a],
                        lx[shellB.aoOffset + b],
                        ly[shellB.aoOffset + b],
                        lz[shellB.aoOffset + b],
                        primPair,
                        I,
                        E_buffer
                    );

                    S(shellA.aoOffset + a, shellB.aoOffset + b) += coeffA * coeffB * S_prim;
                }
            }
        }
    }
}

IntegralEngine::PrimitivePairData IntegralEngine::computePrimitivePairData(
    double alpha1, const Eigen::Vector3d& A, double alpha2, const Eigen::Vector3d& B
)
{
    double p    = alpha1 + alpha2;
    double oo_p = 1.0 / p; // one over p

    Eigen::Vector3d P  = (alpha1 * A + alpha2 * B) * oo_p;
    Eigen::Vector3d PA = P - A;
    Eigen::Vector3d PB = P - B;

    double AmB2 = (A - B).squaredNorm();
    double K    = std::exp(-alpha1 * alpha2 * oo_p * AmB2);

    return {.p = p, .K = K, .P = P, .PA = PA, .PB = PB};
}

void IntegralEngine::computeHermiteCoeffs(unsigned l1, unsigned l2, double PA, double PB, std::span<double> E)
{
    const unsigned L = l1 + l2;
    assert(E.size() >= L + 1);
    std::fill(E.begin(), E.begin() + L + 1, 0.0);
    E[0] = 1.0;

    // Build up E coefficients for (i, 0) from i=1 to l1
    for (size_t i = 1; i <= l1; ++i)
    {
        E[i] = E[i - 1];
        for (int k = i - 1; k >= 1; --k) { E[k] = E[k - 1] + PA * E[k]; }
        E[0] = PA * E[0];
    }

    // Build up E coefficients for (l1, j) from j=1 to l2
    for (size_t j = 1; j <= l2; ++j)
    {
        E[l1 + j] = E[l1 + j - 1];
        for (int k = l1 + j - 1; k >= 1; --k) { E[k] = E[k - 1] + PB * E[k]; }
        E[0] = PB * E[0];
    }
}

double IntegralEngine::compute1dOverlap(
    unsigned l1, unsigned l2, double PA, double PB, const std::vector<double>& I, std::span<double> E_buffer
)
{
    const int L = l1 + l2;
    computeHermiteCoeffs(l1, l2, PA, PB, E_buffer);

    double sum = 0.0;
    for (int i = 0; i <= L; i += 2) { sum += E_buffer[i] * I[i]; }
    return sum;
}

double IntegralEngine::computePrimitiveOverlap(
    unsigned l1,
    unsigned m1,
    unsigned n1,
    unsigned l2,
    unsigned m2,
    unsigned n2,
    const PrimitivePairData& primPair,
    const std::vector<double>& I,
    std::span<double> E_buffer
)
{
    double Sx = compute1dOverlap(l1, l2, primPair.PA.x(), primPair.PB.x(), I, E_buffer);
    double Sy = compute1dOverlap(m1, m2, primPair.PA.y(), primPair.PB.y(), I, E_buffer);
    double Sz = compute1dOverlap(n1, n2, primPair.PA.z(), primPair.PB.z(), I, E_buffer);

    return primPair.K * Sx * Sy * Sz;
}
