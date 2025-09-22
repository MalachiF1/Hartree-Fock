#include "IntegralEngine.hpp"

#include "Boys.hpp"

#include <cassert>
#include <cmath>
#include <fmt/core.h>
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

                Boys::calculateBoys(max_L_total, T, F);
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

void IntegralEngine::electronRepulsion(
    const Basis& basis, const Shell& shellA, const Shell& shellB, const Shell& shellC, const Shell& shellD, ElectronRepulsionTensor& G
)
{
    const auto& exps   = basis.getExponents();
    const auto& coeffs = basis.getCoefficients();
    const auto& lx     = basis.getLx();
    const auto& ly     = basis.getLy();
    const auto& lz     = basis.getLz();

    const size_t nprimA = shellA.nprim, naoA = shellA.nao;
    const size_t nprimB = shellB.nprim, naoB = shellB.nao;
    const size_t nprimC = shellC.nprim, naoC = shellC.nao;
    const size_t nprimD = shellD.nprim, naoD = shellD.nao;

    const Eigen::Vector3d& centerA = shellA.center;
    const Eigen::Vector3d& centerB = shellB.center;
    const Eigen::Vector3d& centerC = shellC.center;
    const Eigen::Vector3d& centerD = shellD.center;

    // Allocate scratch space for Hermite coefficients
    EBuffer Ex_ab(shellA.l, shellB.l), Ey_ab(shellA.l, shellB.l), Ez_ab(shellA.l, shellB.l);
    EBuffer Ex_cd(shellC.l, shellD.l), Ey_cd(shellC.l, shellD.l), Ez_cd(shellC.l, shellD.l);

    const unsigned max_L_total = shellA.l + shellB.l + shellC.l + shellD.l;
    RBuffer R_buffer(max_L_total, max_L_total, max_L_total);
    std::vector<double> F(max_L_total + 1);

    // Loop over primitive quartets
    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double alphaA = exps[shellA.primOffset + pA];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double alphaB = exps[shellB.primOffset + pB];

            const PrimitivePairData primPair_ab = computePrimitivePairData(alphaA, centerA, alphaB, centerB);

            computeHermiteCoeffs(shellA.l, shellB.l, primPair_ab.p, primPair_ab.PA.x(), primPair_ab.PB.x(), Ex_ab);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair_ab.p, primPair_ab.PA.y(), primPair_ab.PB.y(), Ey_ab);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair_ab.p, primPair_ab.PA.z(), primPair_ab.PB.z(), Ez_ab);

            for (size_t pC = 0; pC < nprimC; ++pC)
            {
                const double alphaC = exps[shellC.primOffset + pC];
                for (size_t pD = 0; pD < nprimD; ++pD)
                {
                    const double alphaD = exps[shellD.primOffset + pD];

                    const PrimitivePairData primPair_cd = computePrimitivePairData(alphaC, centerC, alphaD, centerD);

                    const double prefactor = (2.0 * std::pow(M_PI, 2.5)
                                              / (primPair_ab.p * primPair_cd.p * std::sqrt(primPair_ab.p + primPair_cd.p)))
                                           * primPair_ab.K * primPair_cd.K;

                    computeHermiteCoeffs(shellC.l, shellD.l, primPair_cd.p, primPair_cd.PA.x(), primPair_cd.PB.x(), Ex_cd);
                    computeHermiteCoeffs(shellC.l, shellD.l, primPair_cd.p, primPair_cd.PA.y(), primPair_cd.PB.y(), Ey_cd);
                    computeHermiteCoeffs(shellC.l, shellD.l, primPair_cd.p, primPair_cd.PA.z(), primPair_cd.PB.z(), Ez_cd);

                    const double delta         = (primPair_ab.p * primPair_cd.p) / (primPair_ab.p + primPair_cd.p);
                    const Eigen::Vector3d P_ab = primPair_ab.P;
                    const Eigen::Vector3d P_cd = primPair_cd.P;
                    const double T             = delta * (P_ab - P_cd).squaredNorm();

                    Boys::calculateBoys(max_L_total, T, F);
                    computeAuxiliaryIntegrals(max_L_total, max_L_total, max_L_total, delta, P_ab - P_cd, F, R_buffer);

                    for (size_t a = 0; a < naoA; ++a)
                    {
                        const unsigned i  = shellA.aoOffset + a;
                        const unsigned l1 = lx[shellA.aoOffset + a];
                        const unsigned m1 = ly[shellA.aoOffset + a];
                        const unsigned n1 = lz[shellA.aoOffset + a];
                        const double cA   = coeffs[(shellA.coeffOffset + pA * naoA) + a];

                        for (size_t b = 0; b < naoB; ++b)
                        {
                            const unsigned j = shellB.aoOffset + b;
                            if (&shellA == &shellB && j > i) // Exploit permutational symmetry
                                continue;

                            const unsigned l2 = lx[shellB.aoOffset + b];
                            const unsigned m2 = ly[shellB.aoOffset + b];
                            const unsigned n2 = lz[shellB.aoOffset + b];
                            const double cB   = coeffs[(shellB.coeffOffset + pB * naoB) + b];

                            for (size_t c = 0; c < naoC; ++c)
                            {
                                const unsigned k  = shellC.aoOffset + c;
                                const unsigned l3 = lx[shellC.aoOffset + c];
                                const unsigned m3 = ly[shellC.aoOffset + c];
                                const unsigned n3 = lz[shellC.aoOffset + c];
                                const double cC   = coeffs[(shellC.coeffOffset + pC * naoC) + c];

                                for (size_t d = 0; d < naoD; ++d)
                                {
                                    const unsigned l = shellD.aoOffset + d;

                                    // Exploit permutational symmetry
                                    if (&shellC == &shellD && l > k)
                                        continue;
                                    if (((&shellA == &shellC && &shellB == &shellD)
                                         || (&shellA == &shellD && &shellB == &shellC))
                                        && ((i * (i + 1) / 2) + j) < ((k * (k + 1) / 2) + l))
                                        continue;

                                    const unsigned l4 = lx[shellD.aoOffset + d];
                                    const unsigned m4 = ly[shellD.aoOffset + d];
                                    const unsigned n4 = lz[shellD.aoOffset + d];
                                    const double cD   = coeffs[(shellD.coeffOffset + pD * naoD) + d];


                                    double prim_int = 0.0;
                                    for (unsigned t1 = 0; t1 <= l1 + l2; ++t1)
                                    {
                                        double sum_u1 = 0.0;
                                        for (unsigned u1 = 0; u1 <= m1 + m2; ++u1)
                                        {
                                            double sum_v1 = 0.0;
                                            for (unsigned v1 = 0; v1 <= n1 + n2; ++v1)
                                            {
                                                double tmp = 0.0;
                                                for (unsigned t2 = 0; t2 <= l3 + l4; ++t2)
                                                {
                                                    double sum_u2 = 0.0;
                                                    for (unsigned u2 = 0; u2 <= m3 + m4; ++u2)
                                                    {
                                                        double sum_v2 = 0.0;
                                                        for (unsigned v2 = 0; v2 <= n3 + n4; ++v2)
                                                        {
                                                            int signFactor = (t2 + u2 + v2) % 2 == 0 ? 1 : -1;
                                                            sum_v2 += signFactor * Ez_cd(n3, n4, v2)
                                                                    * R_buffer(t1 + t2, u1 + u2, v1 + v2, 0);
                                                        }
                                                        sum_u2 += Ey_cd(m3, m4, u2) * sum_v2;
                                                    }
                                                    tmp += Ex_cd(l3, l4, t2) * sum_u2;
                                                }
                                                sum_v1 += Ez_ab(n1, n2, v1) * tmp;
                                            }
                                            sum_u1 += Ey_ab(m1, m2, u1) * sum_v1;
                                        }
                                        prim_int += Ex_ab(l1, l2, t1) * sum_u1;
                                    }

                                    G(i, j, k, l) += cA * cB * cC * cD * prefactor * prim_int;
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

// namespace
// {
// double SQRT_PI_BY_2 = std::sqrt(M_PI) / 2.0;
//
// std::array<double, 82> TRANSITION_VALUES = {
//     2.35,  2.35,  2.35,  2.35,  2.38,  2.82,  3.05,  3.72,  3.76,  4.01,  4.61,  4.92,  5.15,  5.63,
//     6.02,  6.33,  6.83,  7.15,  7.43,  7.94,  8.27,  8.59,  8.88,  9.21,  9.84,  9.99,  10.40, 10.71,
//     11.15, 11.54, 11.88, 12.37, 12.69, 13.21, 13.42, 13.81, 14.25, 14.65, 14.97, 15.31, 15.84, 15.96,
//     16.22, 16.69, 17.15, 17.49, 17.92, 18.28, 18.61, 18.85, 19.18, 19.68, 19.94, 20.53, 20.73, 21.25,
//     21.54, 21.89, 22.29, 22.72, 23.07, 23.44, 23.90, 24.17, 24.56, 24.85, 25.23, 25.63, 25.96, 26.43,
//     26.78, 27.09, 27.56, 27.89, 28.27, 28.71, 29.00, 29.38, 29.84, 30.22, 30.56, 30.85
// };
//
// } // namespace
//
// /**
//  * This implementation is based on:
//  * Weiss, A. K., & Ochsenfeld, C. (2015). A rigorous and optimized strategy for the evaluation of the B oys function
//  * kernel in molecular electronic structure theory. Journal of Computational Chemistry, 36(18), 1390-1398.
//  */
// void IntegralEngine::calculateBoys(unsigned m_max, double T, std::span<double> F)
// {
//     if (T > 30.0)
//     {
//         // For large T, use the asymptotic expansion.
//         F[0] = SQRT_PI_BY_2 / std::sqrt(T);
//         for (unsigned m = 1; m <= m_max; ++m) { F[m] = (m + 0.5) * F[m - 1] / T; }
//         return;
//     }
//
//     if (m_max == 0)
//     {
//         double sqrt_T = std::sqrt(T);
//         F[0]          = SQRT_PI_BY_2 * std::erf(sqrt_T) / sqrt_T;
//         return;
//     }
//
//     if (T <= TRANSITION_VALUES[m_max - 1])
//     {
//         // use upwards recursion
//         double sqrt_T = std::sqrt(T);
//         F[0]          = SQRT_PI_BY_2 * std::erf(sqrt_T) / sqrt_T;
//         for (unsigned m = 0; m < m_max; ++m) { F[m + 1] = ((2.0 * m + 1.0) * F[m] - std::exp(-T)) / (2.0 * T); }
//         return;
//     }
//
//     // use downward recursion
//     F[m_max] = dfact((2 * m_max) - 1) * std::sqrt(M_PI) / (std::pow(2, m_max + 1) * std::pow(T, m_max + 0.5));
//     for (unsigned m = m_max; m >= 1; --m) { F[m - 1] = ((2 * T * F[m]) + std::exp(-T)) / ((2 * (m - 1)) + 1); }
// }

void IntegralEngine::computeHermiteCoeffs(unsigned l1, unsigned l2, double p, double PA, double PB, EBuffer& E_buffer)
{
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

// namespace
// {
//
// double E(int i, int l1, int l2, double Q, double exponentA, double exponentB)
// {
//     if (i < 0 || i > (l1 + l2) || l1 < 0 || l2 < 0)
//     {
//         // out out bounds
//         return 0.0;
//     }
//     else if (i == 0 && l1 == 0 && l2 == 0)
//     {
//         // base case
//         return 1.0;
//     }
//
//     double p = exponentA + exponentB;
//     double q = (exponentA * exponentB) / p;
//
//     double result = 0.0;
//     if (l2 == 0)
//     {
//         // decrement l1
//         result = ((1.0 / (2.0 * p)) * E(i - 1, l1 - 1, l2, Q, exponentA, exponentB))
//                - ((q * Q / exponentA) * E(i, l1 - 1, l2, Q, exponentA, exponentB))
//                + ((i + 1) * E(i + 1, l1 - 1, l2, Q, exponentA, exponentB));
//     }
//     else
//     {
//         // decrement l2
//         result = ((1.0 / (2.0 * p)) * E(i - 1, l1, l2 - 1, Q, exponentA, exponentB))
//                + ((q * Q / exponentB) * E(i, l1, l2 - 1, Q, exponentA, exponentB))
//                + ((i + 1) * E(i + 1, l1, l2 - 1, Q, exponentA, exponentB));
//     }
//
//     return result;
// }

// double R(int t, int u, int v, int n, double p, const Vec3& PC, double T)
// {
//     if (t < 0 || u < 0 || v < 0)
//     {
//         return 0.0;
//     }
//
//     double result = 0.0;
//     if (t == 0 && u == 0 && v == 0)
//     {
//         // base case
//         result = std::pow(-2.0 * p, n) * boys(n, T);
//     }
//     else if (t > 0)
//     {
//         result = ((t - 1) * R(t - 2, u, v, n + 1, p, PC, T)) + (PC.x() * R(t - 1, u, v, n + 1, p, PC, T));
//     }
//     else if (u > 0)
//     {
//         result = ((u - 1) * R(t, u - 2, v, n + 1, p, PC, T)) + (PC.y() * R(t, u - 1, v, n + 1, p, PC, T));
//     }
//     else // v > 0
//     {
//         result = ((v - 1) * R(t, u, v - 2, n + 1, p, PC, T)) + (PC.z() * R(t, u, v - 1, n + 1, p, PC, T));
//     }
//
//     return result;
// }

// } // namespace


void IntegralEngine::computeAuxiliaryIntegrals(
    unsigned t_max, unsigned u_max, unsigned v_max, double p, const Vec3& PC, std::span<double> F, RBuffer& R_buffer
)
{
    // This function implements the recurrence relations for R integrals:
    // R_{t,u,v}^{n} = t-1 * R_{t-2,u,v}^{n+1} + (PC)_x * R_{t-1,u,v}^{n+1} (if t>0)
    // R_{t,u,v}^{n} = u-1 * R_{t,u-2,v}^{n+1} + (PC)_y * R_{t,u-1,v}^{n+1} (if u>0)
    // R_{t,u,v}^{n} = v-1 * R_{t,u,v-2}^{n+1} + (PC)_z * R_{t,u,v-1}^{n+1} (if v>0)
    // Base case: R_{0,0,0}^{n} = (-2p)^n * F[n]

    const unsigned n_max = t_max + u_max + v_max;

    for (unsigned n = 0; n <= n_max; ++n) { R_buffer(0, 0, 0, n) = std::pow(-2.0 * p, n) * F[n]; }

    unsigned n = n_max;
    while (n-- > 0)
    {
        // Build up t dimension from t=1 up to max_t
        for (unsigned t = 1; t <= t_max; ++t)
        {
            R_buffer(t, 0, 0, n) = PC.x() * R_buffer(t - 1, 0, 0, n + 1);
            if (t >= 2)
                R_buffer(t, 0, 0, n) += (t - 1) * R_buffer(t - 2, 0, 0, n + 1);
        }

        // Build up u dimension
        for (unsigned u = 1; u <= u_max; ++u)
        {
            for (unsigned t = 0; t <= t_max; ++t)
            {
                R_buffer(t, u, 0, n) = PC.y() * R_buffer(t, u - 1, 0, n + 1);
                if (u >= 2)
                    R_buffer(t, u, 0, n) += (u - 1) * R_buffer(t, u - 2, 0, n + 1);
            }
        }

        // Build up v dimension
        for (unsigned v = 1; v <= v_max; ++v)
        {
            for (unsigned u = 0; u <= u_max; ++u)
            {
                for (unsigned t = 0; t <= t_max; ++t)
                {
                    R_buffer(t, u, v, n) = PC.z() * R_buffer(t, u, v - 1, n + 1);
                    if (v >= 2)
                        R_buffer(t, u, v, n) += (v - 1) * R_buffer(t, u, v - 2, n + 1);
                }
            }
        }
    }
}
