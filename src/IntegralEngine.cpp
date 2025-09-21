#include "IntegralEngine.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fmt/core.h>
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
    EBuffer Ex(shellA.l, shellB.l);
    EBuffer Ey(shellA.l, shellB.l);
    EBuffer Ez(shellA.l, shellB.l);

    // Loop over primitive pairs
    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double alphaA = exps[shellA.primOffset + pA];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double alphaB = exps[shellB.primOffset + pB];

            // Compute data for this primitive pair
            const PrimitivePairData primPair = computePrimitivePairData(alphaA, centerA, alphaB, centerB);
            double prefactor                 = std::pow(M_PI / primPair.p, 1.5) * primPair.K;

            computeHermiteCoeffs(shellA.l, shellB.l, primPair.p, primPair.PA.x(), primPair.PB.x(), Ex);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair.p, primPair.PA.y(), primPair.PB.y(), Ey);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair.p, primPair.PA.z(), primPair.PB.z(), Ez);

            // This primitive pair contributes to all AO pairs in the shell pair
            for (size_t a = 0; a < naoA; ++a)
            {
                const double coeffA = coeffs[(shellA.coeffOffset + pA * naoA) + a];

                for (size_t b = 0; b < naoB; ++b)
                {
                    const double coeffB = coeffs[(shellB.coeffOffset + pB * naoB) + b];

                    unsigned l1 = lx[shellA.aoOffset + a], m1 = ly[shellA.aoOffset + a], n1 = lz[shellA.aoOffset + a];
                    unsigned l2 = lx[shellB.aoOffset + b], m2 = ly[shellB.aoOffset + b], n2 = lz[shellB.aoOffset + b];

                    S(shellA.aoOffset + a,
                      shellB.aoOffset + b) += coeffA * coeffB * prefactor * Ex(l1, l2, 0) * Ey(m1, m2, 0) * Ez(n1, n2, 0);
                }
            }
        }
    }
}

void IntegralEngine::kinetic(const Basis& basis, const Shell& shellA, const Shell& shellB, Eigen::MatrixXd& T)
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

    // Allocate scratch space. Max angular momentum needed is L+2 for kinetic integrals.
    EBuffer Ex(shellA.l, shellB.l + 2);
    EBuffer Ey(shellA.l, shellB.l + 2);
    EBuffer Ez(shellA.l, shellB.l + 2);

    // Loop over primitive pairs
    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double alphaA = exps[shellA.primOffset + pA];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double alphaB = exps[shellB.primOffset + pB];

            const PrimitivePairData primPair = computePrimitivePairData(alphaA, centerA, alphaB, centerB);
            double prefactor                 = std::pow(M_PI / primPair.p, 1.5) * primPair.K;

            computeHermiteCoeffs(shellA.l, shellB.l + 2, primPair.p, primPair.PA.x(), primPair.PB.x(), Ex);
            computeHermiteCoeffs(shellA.l, shellB.l + 2, primPair.p, primPair.PA.y(), primPair.PB.y(), Ey);
            computeHermiteCoeffs(shellA.l, shellB.l + 2, primPair.p, primPair.PA.z(), primPair.PB.z(), Ez);

            for (size_t a = 0; a < naoA; ++a)
            {
                const unsigned l1 = lx[shellA.aoOffset + a], m1 = ly[shellA.aoOffset + a], n1 = lz[shellA.aoOffset + a];
                const double coeffA = coeffs[(shellA.coeffOffset + pA * naoA) + a];

                for (size_t b = 0; b < naoB; ++b)
                {
                    const unsigned l2 = lx[shellB.aoOffset + b], m2 = ly[shellB.aoOffset + b], n2 = lz[shellB.aoOffset + b];
                    const double coeffB = coeffs[(shellB.coeffOffset + pB * naoB) + b];

                    double term1 = alphaB * (2 * (l2 + m2 + n2) + 3) * Ex(l1, l2, 0) * Ey(m1, m2, 0) * Ez(n1, n2, 0);

                    double term2 = -2.0 * alphaB * alphaB
                                 * (Ex(l1, l2 + 2, 0) * Ey(m1, m2, 0) * Ez(n1, n2, 0)
                                    + Ex(l1, l2, 0) * Ey(m1, m2 + 2, 0) * Ez(n1, n2, 0)
                                    + Ex(l1, l2, 0) * Ey(m1, m2, 0) * Ez(n1, n2 + 2, 0));

                    double term3 = 0.0;
                    if (l2 >= 2)
                        term3 -= l2 * (l2 - 1) * Ex(l1, l2 - 2, 0) * Ey(m1, m2, 0) * Ez(n1, n2, 0);
                    if (m2 >= 2)
                        term3 -= m2 * (m2 - 1) * Ex(l1, l2, 0) * Ey(m1, m2 - 2, 0) * Ez(n1, n2, 0);
                    if (n2 >= 2)
                        term3 -= n2 * (n2 - 1) * Ex(l1, l2, 0) * Ey(m1, m2, 0) * Ez(n1, n2 - 2, 0);
                    term3 *= 0.5;

                    T(shellA.aoOffset + a, shellB.aoOffset + b) += coeffA * coeffB * prefactor * (term1 + term2 + term3);
                }
            }
        }
    }
}

void IntegralEngine::nuclearAttraction(
    const std::vector<Atom>& geometry, const Basis& basis, const Shell& shellA, const Shell& shellB, Eigen::MatrixXd& V
)
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

    const unsigned max_L_total = shellA.l + shellB.l;

    // Allocate scratch space
    RBuffer R_buffer(max_L_total, max_L_total, max_L_total);
    std::vector<double> F(max_L_total + 1);
    EBuffer Ex(shellA.l, shellB.l);
    EBuffer Ey(shellA.l, shellB.l);
    EBuffer Ez(shellA.l, shellB.l);

    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double alphaA = exps[shellA.primOffset + pA];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double alphaB = exps[shellB.primOffset + pB];

            const PrimitivePairData primPair = computePrimitivePairData(alphaA, centerA, alphaB, centerB);
            double prefactor                 = (2.0 * M_PI / primPair.p) * primPair.K;

            // Hermite coefficients are independent of nuclear centers, compute once per AO pair
            computeHermiteCoeffs(shellA.l, shellB.l, primPair.p, primPair.PA.x(), primPair.PB.x(), Ex);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair.p, primPair.PA.y(), primPair.PB.y(), Ey);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair.p, primPair.PA.z(), primPair.PB.z(), Ez);

            for (const auto& atom : geometry)
            {
                const Eigen::Vector3d& C = atom.coords;
                const auto Z             = static_cast<double>(atom.atomicNumber);

                const Eigen::Vector3d PC = primPair.P - C;
                const double T           = primPair.p * PC.squaredNorm();

                for (unsigned i = 0; i <= max_L_total; ++i) { F[i] = boys(i, T); }
                computeAuxiliaryIntegrals(max_L_total, max_L_total, max_L_total, primPair.p, PC, F, R_buffer);

                for (size_t a = 0; a < naoA; ++a)
                {
                    const unsigned l1   = lx[shellA.aoOffset + a];
                    const unsigned m1   = ly[shellA.aoOffset + a];
                    const unsigned n1   = lz[shellA.aoOffset + a];
                    const double coeffA = coeffs[(shellA.coeffOffset + pA * naoA) + a];

                    for (size_t b = 0; b < naoB; ++b)
                    {
                        const unsigned l2   = lx[shellB.aoOffset + b];
                        const unsigned m2   = ly[shellB.aoOffset + b];
                        const unsigned n2   = lz[shellB.aoOffset + b];
                        const double coeffB = coeffs[(shellB.coeffOffset + pB * naoB) + b];

                        double sum = 0.0;
                        for (unsigned t = 0; t <= l1 + l2; ++t)
                        {
                            double sum_u = 0.0;
                            for (unsigned u = 0; u <= m1 + m2; ++u)
                            {
                                double sum_v = 0.0;
                                for (unsigned v = 0; v <= n1 + n2; ++v)
                                {
                                    sum_v += Ez(n1, n2, v) * R_buffer(t, u, v, 0);
                                }
                                sum_u += Ey(m1, m2, u) * sum_v;
                            }
                            sum += Ex(l1, l2, t) * sum_u;
                        }

                        V(shellA.aoOffset + a, shellB.aoOffset + b) -= coeffA * coeffB * prefactor * Z * sum;
                    }
                }
            }
        }
    }
}





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

namespace
{
double SQRT_PI_BY_2 = std::sqrt(M_PI) / 2.0;
}

double IntegralEngine::boys(unsigned m, double T)
{
    if (T < 1e-9)
    {
        // For small T, use the series expansion for stability.
        // F_m(T) = 1/(2m+1) - T/(2m+3) + ...
        return 1.0 / (2.0 * m + 1.0);
    }

    if (T > 30.0)
    {
        // For large T, use the asymptotic expansion.
        double F_i = SQRT_PI_BY_2 / std::sqrt(T);

        for (unsigned i = 0; i < m; ++i) { F_i = (i + 0.5) * F_i / T; }
        return F_i;
    }

    // Use the upward recursion relation.
    double sqrt_T = std::sqrt(T);
    double F_i    = SQRT_PI_BY_2 * std::erf(sqrt_T) / sqrt_T;
    for (unsigned i = 0; i < m; ++i) { F_i = ((2.0 * i + 1.0) * F_i - std::exp(-T)) / (2.0 * T); }
    return F_i;
}

void IntegralEngine::computeHermiteCoeffs(unsigned l1, unsigned l2, double p, double PA, double PB, EBuffer& E_buffer)
{
    // fmt::println("Computing Hermite coefficients for l1 = {}, l2 = {}", l1, l2);
    // Base case: E(0,0,0) = 1.0 (K factor is applied by the caller)
    E_buffer(0, 0, 0) = 1.0;

    const double oo_2p = 0.5 / p;

    for (unsigned i = 1; i <= l1; ++i)
    {
        size_t t = i + 1;
        while (t-- > 0)
        {
            double val = 0.0;
            if (t <= i - 1)
                val = PA * E_buffer(i - 1, 0, t);
            if (t >= 1)
                val += oo_2p * E_buffer(i - 1, 0, t - 1);
            if (t + 1 <= i - 1)
                val += (t + 1) * E_buffer(i - 1, 0, t + 1);
            E_buffer(i, 0, t) = val;
        }
    }

    for (unsigned j = 1; j <= l2; ++j)
    {
        for (unsigned i = 0; i <= l1; ++i)
        {
            size_t t = i + j + 1;
            while (t-- > 0)
            {
                double val = 0.0;
                if (t <= i + j - 1)
                    val = PB * E_buffer(i, j - 1, t);
                if (t >= 1)
                    val += oo_2p * E_buffer(i, j - 1, t - 1);
                if (t + 1 <= i + j - 1)
                    val += (t + 1) * E_buffer(i, j - 1, t + 1);
                E_buffer(i, j, t) = val;
            }
        }
    }
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
    // This function implements the recurrence relations for R integrals:
    // R_{t,u,v}^{n} = t-1 * R_{t-2,u,v}^{n+1} + (P-C)_x * R_{t-1,u,v}^{n+1} (if t>0)
    // R_{t,u,v}^{n} = u-1 * R_{t,u-2,v}^{n+1} + (P-C)_y * R_{t,u-1,v}^{n+1} (if u>0)
    // R_{t,u,v}^{n} = v-1 * R_{t,u,v-2}^{n+1} + (P-C)_z * R_{t,u,v-1}^{n+1} (if v>0)
    // Base case: R_{0,0,0}^{n} = (-2p)^n * F[n]

    const unsigned n_max = t_max + u_max + v_max;

    for (unsigned n = 0; n <= n_max; ++n) { R_buffer(0, 0, 0, n) = std::pow(-2.0 * p, n) * F[n]; }

    for (int n = n_max; n >= 0; --n)
    {
        // Build up t dimension from t = 1 up to max_t
        for (unsigned t = 1; t <= t_max; ++t)
        {
            // Handle the edge case for t=1 to prevent out-of-bounds access
            double term1         = (t > 1) ? (t - 1) * R_buffer(t - 2, 0, 0, n + 1) : 0.0;
            R_buffer(t, 0, 0, n) = term1 + (PC.x() * R_buffer(t - 1, 0, 0, n + 1));
        }

        // Build up u dimension
        for (unsigned t = 0; t <= t_max; ++t)
        {
            for (unsigned u = 1; u <= u_max; ++u)
            {
                double term1         = (u > 1) ? (u - 1) * R_buffer(t, u - 2, 0, n + 1) : 0.0;
                R_buffer(t, u, 0, n) = term1 + (PC.y() * R_buffer(t, u - 1, 0, n + 1));
            }
        }

        // Build up v dimension
        for (unsigned t = 0; t <= t_max; ++t)
        {
            for (unsigned u = 0; u <= u_max; ++u)
            {
                for (unsigned v = 1; v <= v_max; ++v)
                {
                    double term1         = (v > 1) ? (v - 1) * R_buffer(t, u, v - 2, n + 1) : 0.0;
                    R_buffer(t, u, v, n) = term1 + (PC.z() * R_buffer(t, u, v - 1, n + 1));
                }
            }
        }
    }
}
