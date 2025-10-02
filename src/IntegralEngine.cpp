#include "IntegralEngine.hpp"

#include "Boys.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

void IntegralEngine::overlap(const Basis& basis, const Shell& shellA, const Shell& shellB, Eigen::MatrixXd& S)
{
    const auto& exps   = basis.getExponents();
    const auto& coeffs = basis.getCoefficients();
    const auto& lx     = basis.getLx();
    const auto& ly     = basis.getLy();
    const auto& lz     = basis.getLz();
    const auto& cx     = basis.getCx();
    const auto& cy     = basis.getCy();
    const auto& cz     = basis.getCz();

    const size_t nprimA = shellA.nprim;
    const size_t nprimB = shellB.nprim;
    const size_t naoA   = shellA.nao;
    const size_t naoB   = shellB.nao;

    // Allocate scratch space for Hermite coefficients once per shell pair.
    EBuffer Ex(shellA.l, shellB.l);
    EBuffer Ey(shellA.l, shellB.l);
    EBuffer Ez(shellA.l, shellB.l);

    // Loop over primitive pairs
    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double alphaA = exps[shellA.primOffset + pA];
        const double cxA    = cx[shellA.primOffset + pA];
        const double cyA    = cy[shellA.primOffset + pA];
        const double czA    = cz[shellA.primOffset + pA];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double alphaB = exps[shellB.primOffset + pB];
            const double cxB    = cx[shellB.primOffset + pB];
            const double cyB    = cy[shellB.primOffset + pB];
            const double czB    = cz[shellB.primOffset + pB];
            const Eigen::Vector3d centerA(cxA, cyA, czA);

            // Compute data for this primitive pair
            const PrimitivePairData primPair = computePrimitivePairData(alphaA, cxA, cyA, czA, alphaB, cxB, cyB, czB);
            double prefactor                 = std::pow(M_PI / primPair.p, 1.5) * primPair.K;

            computeHermiteCoeffs(shellA.l, shellB.l, primPair.p, primPair.PAx, primPair.PBx, Ex);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair.p, primPair.PAy, primPair.PBy, Ey);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair.p, primPair.PAz, primPair.PBz, Ez);

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
    const auto& cx     = basis.getCx();
    const auto& cy     = basis.getCy();
    const auto& cz     = basis.getCz();

    const size_t nprimA = shellA.nprim;
    const size_t nprimB = shellB.nprim;
    const size_t naoA   = shellA.nao;
    const size_t naoB   = shellB.nao;

    // Allocate scratch space. Max angular momentum needed is L+2 for kinetic integrals.
    EBuffer Ex(shellA.l, shellB.l + 2);
    EBuffer Ey(shellA.l, shellB.l + 2);
    EBuffer Ez(shellA.l, shellB.l + 2);

    // Loop over primitive pairs
    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double alphaA = exps[shellA.primOffset + pA];
        const double cxA    = cx[shellA.primOffset + pA];
        const double cyA    = cy[shellA.primOffset + pA];
        const double czA    = cz[shellA.primOffset + pA];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double alphaB = exps[shellB.primOffset + pB];
            const double cxB    = cx[shellB.primOffset + pB];
            const double cyB    = cy[shellB.primOffset + pB];
            const double czB    = cz[shellB.primOffset + pB];

            const PrimitivePairData primPair = computePrimitivePairData(alphaA, cxA, cyA, czA, alphaB, cxB, cyB, czB);
            double prefactor                 = std::pow(M_PI / primPair.p, 1.5) * primPair.K;

            computeHermiteCoeffs(shellA.l, shellB.l + 2, primPair.p, primPair.PAx, primPair.PBx, Ex);
            computeHermiteCoeffs(shellA.l, shellB.l + 2, primPair.p, primPair.PAy, primPair.PBy, Ey);
            computeHermiteCoeffs(shellA.l, shellB.l + 2, primPair.p, primPair.PAz, primPair.PBz, Ez);

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
    const auto& cx     = basis.getCx();
    const auto& cy     = basis.getCy();
    const auto& cz     = basis.getCz();

    const size_t nprimA = shellA.nprim;
    const size_t nprimB = shellB.nprim;
    const size_t naoA   = shellA.nao;
    const size_t naoB   = shellB.nao;

    const unsigned max_L_total = shellA.l + shellB.l;

    // Allocate scratch space
    RBuffer R_buffer(max_L_total, max_L_total, max_L_total);
    std::vector<double> F(R_buffer.max_n + 1);
    EBuffer Ex(shellA.l, shellB.l);
    EBuffer Ey(shellA.l, shellB.l);
    EBuffer Ez(shellA.l, shellB.l);

    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double alphaA = exps[shellA.primOffset + pA];
        const double cxA    = cx[shellA.primOffset + pA];
        const double cyA    = cy[shellA.primOffset + pA];
        const double czA    = cz[shellA.primOffset + pA];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double alphaB = exps[shellB.primOffset + pB];
            const double cxB    = cx[shellB.primOffset + pB];
            const double cyB    = cy[shellB.primOffset + pB];
            const double czB    = cz[shellB.primOffset + pB];

            const PrimitivePairData primPair = computePrimitivePairData(alphaA, cxA, cyA, czA, alphaB, cxB, cyB, czB);
            double prefactor                 = (2.0 * M_PI / primPair.p) * primPair.K;

            // Hermite coefficients are independent of nuclear centers, compute once per AO pair
            computeHermiteCoeffs(shellA.l, shellB.l, primPair.p, primPair.PAx, primPair.PBx, Ex);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair.p, primPair.PAy, primPair.PBy, Ey);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair.p, primPair.PAz, primPair.PBz, Ez);

            for (const auto& atom : geometry)
            {
                const Eigen::Vector3d& C = atom.coords;
                const auto Z             = static_cast<double>(atom.atomicNumber);

                const Eigen::Vector3d PC = Vec3(primPair.Px, primPair.Py, primPair.Pz) - C;
                const double T           = primPair.p * PC.squaredNorm();

                Boys::calculateBoys(R_buffer.max_n, T, F);
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
                        for (unsigned v = 0; v <= n1 + n2; ++v)
                        {
                            double sum_u = 0.0;
                            for (unsigned u = 0; u <= m1 + m2; ++u)
                            {
                                double sum_t = 0.0;
                                for (unsigned t = 0; t <= l1 + l2; ++t)
                                {
                                    sum_t += Ex(l1, l2, t) * R_buffer(t, u, v, 0);
                                }
                                sum_u += Ey(m1, m2, u) * sum_t;
                            }
                            sum += Ez(n1, n2, v) * sum_u;
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
    const auto& cx     = basis.getCx();
    const auto& cy     = basis.getCy();
    const auto& cz     = basis.getCz();

    const size_t nprimA = shellA.nprim, naoA = shellA.nao;
    const size_t nprimB = shellB.nprim, naoB = shellB.nao;
    const size_t nprimC = shellC.nprim, naoC = shellC.nao;
    const size_t nprimD = shellD.nprim, naoD = shellD.nao;

    // Pre-calculate and cache all (cd) pair data
    const size_t size = nprimC * nprimD;
    std::vector<EBuffer> Ex_cd_vec(size, EBuffer(shellC.l, shellD.l));
    std::vector<EBuffer> Ey_cd_vec(size, EBuffer(shellC.l, shellD.l));
    std::vector<EBuffer> Ez_cd_vec(size, EBuffer(shellC.l, shellD.l));
    std::vector<PrimitivePairData> primPairs_cd;
    std::vector<size_t> pC_indices, pD_indices;
    primPairs_cd.reserve(size);
    pC_indices.reserve(size);
    pD_indices.reserve(size);

    const double* alphaC = &exps[shellC.primOffset];
    const double* alphaD = &exps[shellD.primOffset];
    const double* cxC    = &cx[shellC.primOffset];
    const double* cxD    = &cx[shellD.primOffset];
    const double* cyC    = &cy[shellC.primOffset];
    const double* cyD    = &cy[shellD.primOffset];
    const double* czC    = &cz[shellC.primOffset];
    const double* czD    = &cz[shellD.primOffset];

    for (size_t pC = 0; pC < nprimC; ++pC)
    {
        for (size_t pD = 0; pD < nprimD; ++pD)
        {
            pC_indices.emplace_back(pC);
            pD_indices.emplace_back(pD);
            primPairs_cd.emplace_back(
                computePrimitivePairData(alphaC[pC], cxC[pC], cyC[pC], czC[pC], alphaD[pD], cxD[pD], cyD[pD], czD[pD])
            );
            const size_t idx = (pC * nprimD) + pD;

            computeHermiteCoeffs(
                shellC.l, shellD.l, primPairs_cd[idx].p, primPairs_cd[idx].PAx, primPairs_cd[idx].PBx, Ex_cd_vec[idx]
            );
            computeHermiteCoeffs(
                shellC.l, shellD.l, primPairs_cd[idx].p, primPairs_cd[idx].PAy, primPairs_cd[idx].PBy, Ey_cd_vec[idx]
            );
            computeHermiteCoeffs(
                shellC.l, shellD.l, primPairs_cd[idx].p, primPairs_cd[idx].PAz, primPairs_cd[idx].PBz, Ez_cd_vec[idx]
            );
        }
    }

    // Allocate scratch space for (ab) pair and other intermediates
    EBuffer Ex_ab(shellA.l, shellB.l), Ey_ab(shellA.l, shellB.l), Ez_ab(shellA.l, shellB.l);
    const unsigned max_L_total = shellA.l + shellB.l + shellC.l + shellD.l;
    RBuffer R_buffer(max_L_total, max_L_total, max_L_total);
    std::vector<double> F(R_buffer.max_n + 1);

    const double* alphaA = &exps[shellA.primOffset];
    const double* alphaB = &exps[shellB.primOffset];
    const double* cxA    = &cx[shellA.primOffset];
    const double* cxB    = &cx[shellB.primOffset];
    const double* cyA    = &cy[shellA.primOffset];
    const double* cyB    = &cy[shellB.primOffset];
    const double* czA    = &cz[shellA.primOffset];
    const double* czB    = &cz[shellB.primOffset];

    const unsigned* lxA = &lx[shellA.aoOffset];
    const unsigned* lxB = &lx[shellB.aoOffset];
    const unsigned* lxC = &lx[shellC.aoOffset];
    const unsigned* lxD = &lx[shellD.aoOffset];
    const unsigned* lyA = &ly[shellA.aoOffset];
    const unsigned* lyB = &ly[shellB.aoOffset];
    const unsigned* lyC = &ly[shellC.aoOffset];
    const unsigned* lyD = &ly[shellD.aoOffset];
    const unsigned* lzA = &lz[shellA.aoOffset];
    const unsigned* lzB = &lz[shellB.aoOffset];
    const unsigned* lzC = &lz[shellC.aoOffset];
    const unsigned* lzD = &lz[shellD.aoOffset];

    // Create a temporary tensor to store all primitive integrals for this quartet
    std::vector<double> prim_integrals(naoA * naoB * naoC * naoD);

    // Loop over (ab) primitive pairs
    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double* coeffsA = &coeffs[(shellA.coeffOffset + pA * naoA)];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double* coeffsB = &coeffs[(shellB.coeffOffset + pB * naoB)];

            const PrimitivePairData primPair_ab = computePrimitivePairData(
                alphaA[pA], cxA[pA], cyA[pA], czA[pA], alphaB[pB], cxB[pB], cyB[pB], czB[pB]
            );

            computeHermiteCoeffs(shellA.l, shellB.l, primPair_ab.p, primPair_ab.PAx, primPair_ab.PBx, Ex_ab);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair_ab.p, primPair_ab.PAy, primPair_ab.PBy, Ey_ab);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair_ab.p, primPair_ab.PAz, primPair_ab.PBz, Ez_ab);

            // Loop over cached (cd) pairs
            for (size_t cd_idx = 0; cd_idx < primPairs_cd.size(); ++cd_idx)
            {
                const auto& primPair_cd = primPairs_cd[cd_idx];
                const size_t pC         = pC_indices[cd_idx];
                const size_t pD         = pD_indices[cd_idx];
                const EBuffer& Ex_cd    = Ex_cd_vec[cd_idx];
                const EBuffer& Ey_cd    = Ey_cd_vec[cd_idx];
                const EBuffer& Ez_cd    = Ez_cd_vec[cd_idx];

                const double* coeffsC = &coeffs[(shellC.coeffOffset + pC * naoC)];
                const double* coeffsD = &coeffs[(shellD.coeffOffset + pD * naoD)];

                const double prefactor = (2.0 * std::pow(M_PI, 2.5)
                                          / (primPair_ab.p * primPair_cd.p * std::sqrt(primPair_ab.p + primPair_cd.p)))
                                       * primPair_ab.K * primPair_cd.K;

                const double delta = (primPair_ab.p * primPair_cd.p) / (primPair_ab.p + primPair_cd.p);
                const auto P_ab    = Vec3(primPair_ab.Px, primPair_ab.Py, primPair_ab.Pz);
                const auto P_cd    = Vec3(primPair_cd.Px, primPair_cd.Py, primPair_cd.Pz);
                const double T     = delta * (P_ab - P_cd).squaredNorm();

                Boys::calculateBoys(R_buffer.max_n, T, F);
                computeAuxiliaryIntegrals(max_L_total, max_L_total, max_L_total, delta, P_ab - P_cd, F, R_buffer);

                const unsigned naoCD  = naoC * naoD;
                const unsigned naoBCD = naoB * naoCD;

                for (size_t a = 0; a < naoA; ++a)
                {
                    const size_t i          = shellA.aoOffset + a;
                    const size_t big_I_base = i * (i + 1) / 2;
                    const unsigned l1       = lxA[a];
                    const unsigned m1       = lyA[a];
                    const unsigned n1       = lzA[a];

                    for (size_t b = 0; b < naoB; ++b)
                    {
                        const size_t j = shellB.aoOffset + b;
                        if (j > i) // Exploit permutational symmetry
                            continue;
                        const size_t big_I = big_I_base + j;

                        const unsigned l2 = lxB[b];
                        const unsigned m2 = lyB[b];
                        const unsigned n2 = lzB[b];

                        const Eigen::Map<const Eigen::ArrayXd> Ex_ab_array(&Ex_ab(l1, l2, 0), l1 + l2 + 1);
                        const Eigen::Map<const Eigen::ArrayXd> Ey_ab_array(&Ey_ab(m1, m2, 0), m1 + m2 + 1);
                        const Eigen::Map<const Eigen::ArrayXd> Ez_ab_array(&Ez_ab(n1, n2, 0), n1 + n2 + 1);

                        for (size_t c = 0; c < naoC; ++c)
                        {
                            const size_t k          = shellC.aoOffset + c;
                            const size_t big_K_base = k * (k + 1) / 2;
                            const unsigned l3       = lxC[c];
                            const unsigned m3       = lyC[c];
                            const unsigned n3       = lzC[c];

                            for (size_t d = 0; d < naoD; ++d)
                            {
                                const size_t l = shellD.aoOffset + d;

                                if (l > k)
                                    continue;
                                const size_t big_K = big_K_base + l;
                                if ((&shellA == &shellC && &shellB == &shellD) && big_I < big_K)
                                    continue;

                                const unsigned l4 = lxD[d];
                                const unsigned m4 = lyD[d];
                                const unsigned n4 = lzD[d];

                                const Eigen::Map<const Eigen::ArrayXd> Ex_cd_array(&Ex_cd(l3, l4, 0), l3 + l4 + 1);
                                const Eigen::Map<const Eigen::ArrayXd> Ey_cd_array(&Ey_cd(m3, m4, 0), m3 + m4 + 1);
                                const Eigen::Map<const Eigen::ArrayXd> Ez_cd_array(&Ez_cd(n3, n4, 0), n3 + n4 + 1);

                                const unsigned l12 = l1 + l2;
                                const unsigned l34 = l3 + l4;
                                const unsigned m12 = m1 + m2;
                                const unsigned m34 = m3 + m4;
                                const unsigned n12 = n1 + n2;
                                const unsigned n34 = n3 + n4;
                                double prim_int    = 0.0;
                                for (unsigned v1 = 0; v1 <= n12; ++v1)
                                {
                                    double sum_v2 = 0.0;
                                    for (unsigned v2 = 0; v2 <= n34; ++v2)
                                    {
                                        const size_t vStride = (v1 + v2) * R_buffer.stride_v;
                                        double sum_u1        = 0.0;
                                        for (unsigned u1 = 0; u1 <= m12; ++u1)
                                        {
                                            double sum_u2 = 0.0;
                                            for (unsigned u2 = 0; u2 <= m34; ++u2)
                                            {
                                                const size_t uvStride = (u1 + u2) * R_buffer.stride_u + vStride;
                                                const Eigen::Map<const Eigen::ArrayXd> Ruv_array(
                                                    &R_buffer.data[uvStride], l12 + l34 + 1
                                                );

                                                const unsigned u2v2    = u2 + v2;
                                                const auto signFactors = Eigen::ArrayXd::NullaryExpr(
                                                    l34 + 1, [&u2v2](Eigen::Index i) { return 1 - 2 * ((i + u2v2) & 1); }
                                                );

                                                double sum_t1 = 0.0;
                                                for (unsigned t1 = 0; t1 <= l12; ++t1)
                                                {
                                                    // Vectorized inner loop over t2
                                                    sum_t1 += Ex_ab_array[t1]
                                                            * (signFactors * Ex_cd_array * Ruv_array.segment(t1, l34 + 1))
                                                                  .sum();
                                                }
                                                sum_u2 += Ey_cd_array[u2] * sum_t1;
                                            }
                                            sum_u1 += Ey_ab_array[u1] * sum_u2;
                                        }
                                        sum_v2 += Ez_cd_array[v2] * sum_u1;
                                    }
                                    prim_int += Ez_ab_array[v1] * sum_v2;
                                }
                                prim_integrals[a * naoBCD + b * naoCD + c * naoD + d] = prim_int;
                            }
                        }
                    }
                }

                // Contract Integrals
                for (size_t a = 0; a < naoA; ++a)
                {
                    const size_t i          = shellA.aoOffset + a;
                    const size_t big_I_base = i * (i + 1) / 2;
                    const double cA         = coeffsA[a];

                    for (size_t b = 0; b < naoB; ++b)
                    {
                        const size_t j = shellB.aoOffset + b;
                        if (j > i)
                            continue;
                        const size_t big_I = big_I_base + j;

                        const double cB = coeffsB[b];

                        for (size_t c = 0; c < naoC; ++c)
                        {
                            const size_t k          = shellC.aoOffset + c;
                            const size_t big_K_base = k * (k + 1) / 2;

                            const double cC = coeffsC[c];

                            for (size_t d = 0; d < naoD; ++d)
                            {
                                const size_t l = shellD.aoOffset + d;
                                if (l > k)
                                    continue;
                                const size_t big_K = big_K_base + l;

                                if ((&shellA == &shellC && &shellB == &shellD) && big_I < big_K)
                                    continue;

                                const double cD = coeffsD[d];

                                const size_t index = big_I >= big_K ? (big_I * (big_I + 1) / 2 + big_K)
                                                                    : (big_K * (big_K + 1) / 2 + big_I);

                                G[index] += cA * cB * cC * cD * prefactor
                                          * prim_integrals[a * naoBCD + b * naoCD + c * naoD + d];
                            }
                        }
                    }
                }
            }
        }
    }
}

void IntegralEngine::electronRepulsion(
    const Basis& basis, const Shell& shellA, const Shell& shellB, const Shell& shellC, const Shell& shellD, std::vector<double>& G
)
{
    const auto& exps   = basis.getExponents();
    const auto& coeffs = basis.getCoefficients();
    const auto& lx     = basis.getLx();
    const auto& ly     = basis.getLy();
    const auto& lz     = basis.getLz();
    const auto& cx     = basis.getCx();
    const auto& cy     = basis.getCy();
    const auto& cz     = basis.getCz();

    const size_t nprimA = shellA.nprim, naoA = shellA.nao;
    const size_t nprimB = shellB.nprim, naoB = shellB.nao;
    const size_t nprimC = shellC.nprim, naoC = shellC.nao;
    const size_t nprimD = shellD.nprim, naoD = shellD.nao;

    // Pre-calculate and cache all (cd) pair data
    const size_t size = nprimC * nprimD;
    std::vector<EBuffer> Ex_cd_vec(size, EBuffer(shellC.l, shellD.l));
    std::vector<EBuffer> Ey_cd_vec(size, EBuffer(shellC.l, shellD.l));
    std::vector<EBuffer> Ez_cd_vec(size, EBuffer(shellC.l, shellD.l));
    std::vector<PrimitivePairData> primPairs_cd;
    std::vector<size_t> pC_indices, pD_indices;
    primPairs_cd.reserve(size);
    pC_indices.reserve(size);
    pD_indices.reserve(size);

    const double* alphaC = &exps[shellC.primOffset];
    const double* alphaD = &exps[shellD.primOffset];
    const double* cxC    = &cx[shellC.primOffset];
    const double* cxD    = &cx[shellD.primOffset];
    const double* cyC    = &cy[shellC.primOffset];
    const double* cyD    = &cy[shellD.primOffset];
    const double* czC    = &cz[shellC.primOffset];
    const double* czD    = &cz[shellD.primOffset];

    for (size_t pC = 0; pC < nprimC; ++pC)
    {
        for (size_t pD = 0; pD < nprimD; ++pD)
        {
            pC_indices.emplace_back(pC);
            pD_indices.emplace_back(pD);
            primPairs_cd.emplace_back(
                computePrimitivePairData(alphaC[pC], cxC[pC], cyC[pC], czC[pC], alphaD[pD], cxD[pD], cyD[pD], czD[pD])
            );
            const size_t idx = (pC * nprimD) + pD;

            computeHermiteCoeffs(
                shellC.l, shellD.l, primPairs_cd[idx].p, primPairs_cd[idx].PAx, primPairs_cd[idx].PBx, Ex_cd_vec[idx]
            );
            computeHermiteCoeffs(
                shellC.l, shellD.l, primPairs_cd[idx].p, primPairs_cd[idx].PAy, primPairs_cd[idx].PBy, Ey_cd_vec[idx]
            );
            computeHermiteCoeffs(
                shellC.l, shellD.l, primPairs_cd[idx].p, primPairs_cd[idx].PAz, primPairs_cd[idx].PBz, Ez_cd_vec[idx]
            );
        }
    }

    // Allocate scratch space for (ab) pair and other intermediates
    EBuffer Ex_ab(shellA.l, shellB.l), Ey_ab(shellA.l, shellB.l), Ez_ab(shellA.l, shellB.l);
    const unsigned max_L_total = shellA.l + shellB.l + shellC.l + shellD.l;
    RBuffer R_buffer(max_L_total, max_L_total, max_L_total);
    std::vector<double> F(R_buffer.max_n + 1);

    const double* alphaA = &exps[shellA.primOffset];
    const double* alphaB = &exps[shellB.primOffset];
    const double* cxA    = &cx[shellA.primOffset];
    const double* cxB    = &cx[shellB.primOffset];
    const double* cyA    = &cy[shellA.primOffset];
    const double* cyB    = &cy[shellB.primOffset];
    const double* czA    = &cz[shellA.primOffset];
    const double* czB    = &cz[shellB.primOffset];

    const unsigned* lxA = &lx[shellA.aoOffset];
    const unsigned* lxB = &lx[shellB.aoOffset];
    const unsigned* lxC = &lx[shellC.aoOffset];
    const unsigned* lxD = &lx[shellD.aoOffset];
    const unsigned* lyA = &ly[shellA.aoOffset];
    const unsigned* lyB = &ly[shellB.aoOffset];
    const unsigned* lyC = &ly[shellC.aoOffset];
    const unsigned* lyD = &ly[shellD.aoOffset];
    const unsigned* lzA = &lz[shellA.aoOffset];
    const unsigned* lzB = &lz[shellB.aoOffset];
    const unsigned* lzC = &lz[shellC.aoOffset];
    const unsigned* lzD = &lz[shellD.aoOffset];

    // Create a temporary tensor to store all primitive integrals for this quartet
    std::vector<double> prim_integrals(naoA * naoB * naoC * naoD);

    // Loop over (ab) primitive pairs
    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double* coeffsA = &coeffs[(shellA.coeffOffset + pA * naoA)];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double* coeffsB = &coeffs[(shellB.coeffOffset + pB * naoB)];

            const PrimitivePairData primPair_ab = computePrimitivePairData(
                alphaA[pA], cxA[pA], cyA[pA], czA[pA], alphaB[pB], cxB[pB], cyB[pB], czB[pB]
            );

            computeHermiteCoeffs(shellA.l, shellB.l, primPair_ab.p, primPair_ab.PAx, primPair_ab.PBx, Ex_ab);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair_ab.p, primPair_ab.PAy, primPair_ab.PBy, Ey_ab);
            computeHermiteCoeffs(shellA.l, shellB.l, primPair_ab.p, primPair_ab.PAz, primPair_ab.PBz, Ez_ab);

            // Loop over cached (cd) pairs
            for (size_t cd_idx = 0; cd_idx < primPairs_cd.size(); ++cd_idx)
            {
                const auto& primPair_cd = primPairs_cd[cd_idx];
                const size_t pC         = pC_indices[cd_idx];
                const size_t pD         = pD_indices[cd_idx];
                const EBuffer& Ex_cd    = Ex_cd_vec[cd_idx];
                const EBuffer& Ey_cd    = Ey_cd_vec[cd_idx];
                const EBuffer& Ez_cd    = Ez_cd_vec[cd_idx];

                const double* coeffsC = &coeffs[(shellC.coeffOffset + pC * naoC)];
                const double* coeffsD = &coeffs[(shellD.coeffOffset + pD * naoD)];

                const double prefactor = (2.0 * std::pow(M_PI, 2.5)
                                          / (primPair_ab.p * primPair_cd.p * std::sqrt(primPair_ab.p + primPair_cd.p)))
                                       * primPair_ab.K * primPair_cd.K;

                const double delta = (primPair_ab.p * primPair_cd.p) / (primPair_ab.p + primPair_cd.p);
                const auto P_ab    = Vec3(primPair_ab.Px, primPair_ab.Py, primPair_ab.Pz);
                const auto P_cd    = Vec3(primPair_cd.Px, primPair_cd.Py, primPair_cd.Pz);
                const double T     = delta * (P_ab - P_cd).squaredNorm();

                Boys::calculateBoys(R_buffer.max_n, T, F);
                computeAuxiliaryIntegrals(max_L_total, max_L_total, max_L_total, delta, P_ab - P_cd, F, R_buffer);

                const unsigned naoCD  = naoC * naoD;
                const unsigned naoBCD = naoB * naoCD;

                for (size_t a = 0; a < naoA; ++a)
                {
                    const size_t i          = shellA.aoOffset + a;
                    const size_t big_I_base = i * (i + 1) / 2;
                    const unsigned l1       = lxA[a];
                    const unsigned m1       = lyA[a];
                    const unsigned n1       = lzA[a];

                    for (size_t b = 0; b < naoB; ++b)
                    {
                        const size_t j = shellB.aoOffset + b;
                        if (j > i) // Exploit permutational symmetry
                            continue;
                        const size_t big_I = big_I_base + j;

                        const unsigned l2 = lxB[b];
                        const unsigned m2 = lyB[b];
                        const unsigned n2 = lzB[b];

                        const Eigen::Map<const Eigen::ArrayXd> Ex_ab_array(&Ex_ab(l1, l2, 0), l1 + l2 + 1);
                        const Eigen::Map<const Eigen::ArrayXd> Ey_ab_array(&Ey_ab(m1, m2, 0), m1 + m2 + 1);
                        const Eigen::Map<const Eigen::ArrayXd> Ez_ab_array(&Ez_ab(n1, n2, 0), n1 + n2 + 1);

                        for (size_t c = 0; c < naoC; ++c)
                        {
                            const size_t k          = shellC.aoOffset + c;
                            const size_t big_K_base = k * (k + 1) / 2;
                            const unsigned l3       = lxC[c];
                            const unsigned m3       = lyC[c];
                            const unsigned n3       = lzC[c];

                            for (size_t d = 0; d < naoD; ++d)
                            {
                                const size_t l = shellD.aoOffset + d;

                                if (l > k)
                                    continue;
                                const size_t big_K = big_K_base + l;
                                if ((&shellA == &shellC && &shellB == &shellD) && big_I < big_K)
                                    continue;

                                const unsigned l4 = lxD[d];
                                const unsigned m4 = lyD[d];
                                const unsigned n4 = lzD[d];

                                const Eigen::Map<const Eigen::ArrayXd> Ex_cd_array(&Ex_cd(l3, l4, 0), l3 + l4 + 1);
                                const Eigen::Map<const Eigen::ArrayXd> Ey_cd_array(&Ey_cd(m3, m4, 0), m3 + m4 + 1);
                                const Eigen::Map<const Eigen::ArrayXd> Ez_cd_array(&Ez_cd(n3, n4, 0), n3 + n4 + 1);

                                const unsigned l12 = l1 + l2;
                                const unsigned l34 = l3 + l4;
                                const unsigned m12 = m1 + m2;
                                const unsigned m34 = m3 + m4;
                                const unsigned n12 = n1 + n2;
                                const unsigned n34 = n3 + n4;
                                double prim_int    = 0.0;
                                for (unsigned v1 = 0; v1 <= n12; ++v1)
                                {
                                    double sum_v2 = 0.0;
                                    for (unsigned v2 = 0; v2 <= n34; ++v2)
                                    {
                                        const size_t vStride = (v1 + v2) * R_buffer.stride_v;
                                        double sum_u1        = 0.0;
                                        for (unsigned u1 = 0; u1 <= m12; ++u1)
                                        {
                                            double sum_u2 = 0.0;
                                            for (unsigned u2 = 0; u2 <= m34; ++u2)
                                            {
                                                const size_t uvStride = (u1 + u2) * R_buffer.stride_u + vStride;
                                                const Eigen::Map<const Eigen::ArrayXd> Ruv_array(
                                                    &R_buffer.data[uvStride], l12 + l34 + 1
                                                );

                                                const unsigned u2v2    = u2 + v2;
                                                const auto signFactors = Eigen::ArrayXd::NullaryExpr(
                                                    l34 + 1, [&u2v2](Eigen::Index i) { return 1 - 2 * ((i + u2v2) & 1); }
                                                );

                                                double sum_t1 = 0.0;
                                                for (unsigned t1 = 0; t1 <= l12; ++t1)
                                                {
                                                    // Vectorized inner loop over t2
                                                    sum_t1 += Ex_ab_array[t1]
                                                            * (signFactors * Ex_cd_array * Ruv_array.segment(t1, l34 + 1))
                                                                  .sum();
                                                }
                                                sum_u2 += Ey_cd_array[u2] * sum_t1;
                                            }
                                            sum_u1 += Ey_ab_array[u1] * sum_u2;
                                        }
                                        sum_v2 += Ez_cd_array[v2] * sum_u1;
                                    }
                                    prim_int += Ez_ab_array[v1] * sum_v2;
                                }
                                prim_integrals[a * naoBCD + b * naoCD + c * naoD + d] = prim_int;
                            }
                        }
                    }
                }

                // Contract Integrals
                for (size_t a = 0; a < naoA; ++a)
                {
                    const size_t i          = shellA.aoOffset + a;
                    const size_t big_I_base = i * (i + 1) / 2;
                    const double cA         = coeffsA[a];

                    for (size_t b = 0; b < naoB; ++b)
                    {
                        const size_t j = shellB.aoOffset + b;
                        if (j > i)
                            continue;
                        const size_t big_I = big_I_base + j;

                        const double cB = coeffsB[b];

                        for (size_t c = 0; c < naoC; ++c)
                        {
                            const size_t k          = shellC.aoOffset + c;
                            const size_t big_K_base = k * (k + 1) / 2;

                            const double cC = coeffsC[c];

                            for (size_t d = 0; d < naoD; ++d)
                            {
                                const size_t l = shellD.aoOffset + d;
                                if (l > k)
                                    continue;
                                const size_t big_K = big_K_base + l;

                                if ((&shellA == &shellC && &shellB == &shellD) && big_I < big_K)
                                    continue;

                                const double cD = coeffsD[d];

                                G[a * naoBCD + b * naoCD + c * naoD + d] += cA * cB * cC * cD * prefactor
                                                                          * prim_integrals[a * naoBCD + b * naoCD + c * naoD + d];
                            }
                        }
                    }
                }
            }
        }
    }
}

IntegralEngine::PrimitivePairData IntegralEngine::computePrimitivePairData(
    double alpha1, double xa, double ya, double za, double alpha2, double xb, double yb, double zb
)
{
    double p    = alpha1 + alpha2;
    double oo_p = 1.0 / p; // one over p

    double Px = (alpha1 * xa + alpha2 * xb) * oo_p;
    double Py = (alpha1 * ya + alpha2 * yb) * oo_p;
    double Pz = (alpha1 * za + alpha2 * zb) * oo_p;

    double PAx = Px - xa;
    double PAy = Py - ya;
    double PAz = Pz - za;

    double PBx = Px - xb;
    double PBy = Py - yb;
    double PBz = Pz - zb;

    double AmB2 = (xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb);
    double K    = std::exp(-alpha1 * alpha2 * oo_p * AmB2);

    return {.p = p, .K = K, .Px = Px, .Py = Py, .Pz = Pz, .PAx = PAx, .PAy = PAy, .PAz = PAz, .PBx = PBx, .PBy = PBy, .PBz = PBz};
}

void IntegralEngine::computeHermiteCoeffs(unsigned l1, unsigned l2, double p, double PA, double PB, EBuffer& E_buffer)
{
    // Base case: E(0,0,0) = 1.0 (K factor is applied by the caller)
    E_buffer(0, 0, 0) = 1.0;

    const double oo_2p = 0.5 / p;

    for (unsigned i = 1; i <= l1; ++i)
    {
        const size_t offset_im1       = E_buffer.getIOffset(i - 1);
        const size_t offset_i         = E_buffer.getIOffset(i);
        const double* const E_im1_ptr = &E_buffer.data[offset_im1];
        double* const E_i_ptr         = &E_buffer.data[offset_i];

        size_t t = i + 1;
        while (t-- > 0)
        {
            double val = 0.0;
            if (t <= i - 1)
                val = PA * E_im1_ptr[t];
            if (t >= 1)
                val += oo_2p * E_im1_ptr[t - 1];
            if (t + 1 <= i - 1)
                val += (t + 1) * E_im1_ptr[t + 1];

            E_i_ptr[t] = val;
        }
    }

    for (unsigned j = 1; j <= l2; ++j)
    {
        for (unsigned i = 0; i <= l1; ++i)
        {
            const size_t offset_i         = E_buffer.getIOffset(i);
            const size_t offset_jm1       = E_buffer.getJOffset(i, j - 1);
            const size_t offset_j         = E_buffer.getJOffset(i, j);
            const double* const E_jm1_ptr = &E_buffer.data[offset_i + offset_jm1];
            double* const E_ij_ptr        = &E_buffer.data[offset_i + offset_j];

            size_t t = i + j + 1;
            while (t-- > 0)
            {
                double val = 0.0;
                if (t <= i + j - 1)
                    val = PB * E_jm1_ptr[t];
                if (t >= 1)
                    val += oo_2p * E_jm1_ptr[t - 1];
                if (t + 1 <= i + j - 1)
                    val += (t + 1) * E_jm1_ptr[t + 1];

                // E_buffer(i, j, t) = val;
                E_ij_ptr[t] = val;
            }
        }
    }
}

void IntegralEngine::computeAuxiliaryIntegrals(
    unsigned t_max, unsigned u_max, unsigned v_max, double p, const Vec3& PC, std::span<double> F, RBuffer& R_buffer
)
{
    // This function implements the recurrence relations for R integrals:
    // R_{t,u,v}^{n} = t-1 * R_{t-2,u,v}^{n+1} + (PC)_x * R_{t-1,u,v}^{n+1} (if t>0)
    // R_{t,u,v}^{n} = u-1 * R_{t,u-2,v}^{n+1} + (PC)_y * R_{t,u-1,v}^{n+1} (if u>0)
    // R_{t,u,v}^{n} = v-1 * R_{t,u,v-2}^{n+1} + (PC)_z * R_{t,u,v-1}^{n+1} (if v>0)
    // Base case: R_{0,0,0}^{n} = (-2p)^n * F[n]
    //
    const unsigned n_max = t_max + u_max + v_max;

    const double neg_2p      = -2.0 * p;
    double neg_2p_pow_n      = 1.0;
    double* const R_000n_ptr = &R_buffer(0, 0, 0, 0);
    size_t index             = 0;
    for (unsigned n = 0; n <= n_max; ++n)
    {
        R_000n_ptr[index] = neg_2p_pow_n * F[n];
        index += R_buffer.stride_n;
        neg_2p_pow_n *= neg_2p;
    }

    unsigned n = n_max;
    while (n-- > 0)
    {
        // Build up t dimension from t=1 up to max_t
        const double* const R_np1_t_ptr = &R_buffer(0, 0, 0, n + 1);
        double* const R_n_t_ptr         = &R_buffer(0, 0, 0, n);
        for (unsigned t = 1; t <= t_max; ++t)
        {
            const double val = PC.x() * R_np1_t_ptr[t - 1];
            if (t >= 2)
                R_n_t_ptr[t] = val + (t - 1) * R_np1_t_ptr[t - 2];
            else
                R_n_t_ptr[t] = val;
        }

        // Build up u dimension
        for (unsigned u = 1; u <= u_max; ++u)
        {
            double* const R_n_ptr             = &R_buffer(0, u, 0, n);
            const double* const R_np1_um1_ptr = &R_buffer(0, u - 1, 0, n + 1);
            if (u >= 2)
            {
                const double* R_np1_um2_ptr = &R_buffer(0, u - 2, 0, n + 1);
                for (unsigned t = 0; t <= t_max; ++t)
                {
                    R_n_ptr[t] = PC.y() * R_np1_um1_ptr[t] + (u - 1) * R_np1_um2_ptr[t];
                }
            }
            else
            { // u is exactly 1
                for (unsigned t = 0; t <= t_max; ++t) { R_n_ptr[t] = PC.y() * R_np1_um1_ptr[t]; }
            }
        }

        // Build up v dimension
        for (unsigned v = 1; v <= v_max; ++v)
        {
            for (unsigned u = 0; u <= u_max; ++u)
            {
                // Get pointers to the start of the relevant rows for efficiency
                double* const R_n_ptr             = &R_buffer(0, u, v, n);
                const double* const R_np1_vm1_ptr = &R_buffer(0, u, v - 1, n + 1);
                if (v >= 2)
                {
                    const double* const R_np1_vm2_ptr = &R_buffer(0, u, v - 2, n + 1);
                    for (unsigned t = 0; t <= t_max; ++t)
                    {
                        R_n_ptr[t] = PC.z() * R_np1_vm1_ptr[t] + (v - 1) * R_np1_vm2_ptr[t];
                    }
                }
                else
                { // v is exactly 1
                    for (unsigned t = 0; t <= t_max; ++t) { R_n_ptr[t] = PC.z() * R_np1_vm1_ptr[t]; }
                }
            }
        }
    }
}
