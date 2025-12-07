#include "IntegralEngine.hpp"

#include "Boys.hpp"
#include "Eigen/Core"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <fmt/core.h>
#include <numeric>
#include <span>
#include <vector>

IntegralEngine::IntegralEngine(const BasisSet& unsortedBasis)
{
    // Sort shells by angular momentum and contraction degree.
    BasisSet basis = unsortedBasis;
    std::ranges::sort(
        basis.shells,
        [](const Shell& a, const Shell& b)
        {
            if (a.l != b.l)
                return a.l < b.l;
            return a.nPrimitives < b.nPrimitives;
        }
    );

    // Set Basis data
    this->nShells   = basis.nShells;
    this->nAOsTotal = basis.nAOs;

    this->lShell.resize(basis.nShells);
    this->primitiveOffsets.resize(basis.nShells);
    this->aoOffsets.resize(basis.nShells);
    this->normFactorOffsets.resize(basis.nShells);
    this->nPrimitives.resize(basis.nShells);
    this->nAOs.resize(basis.nShells);
    this->centers.resize(3, basis.nShells);
    this->alpha.resize(basis.nPrimitives);
    this->coeff.resize(basis.nPrimitives);
    this->normalizationFactors.resize(basis.nNormFactors);
    this->angularMomentum.resize(3, basis.nAOs);

    const auto& shells      = basis.shells;
    size_t primOffset       = 0;
    size_t aoOffset         = 0;
    size_t normFactorOffset = 0;
    for (size_t i = 0; i < basis.nShells; ++i)
    {
        lShell[i]            = shells[i].l;
        primitiveOffsets[i]  = primOffset;
        aoOffsets[i]         = aoOffset;
        normFactorOffsets[i] = normFactorOffset;
        nPrimitives[i]       = shells[i].nPrimitives;
        nAOs[i]              = shells[i].nAOs;
        centers.col(i)       = shells[i].center;

        for (size_t j = 0; j < shells[i].nPrimitives; ++j)
        {
            const auto& primitive = shells[i].primitives[j];
            alpha[primOffset + j] = primitive.exponent;
            coeff[primOffset + j] = primitive.coefficient;
        }

        normalizationFactors.segment(normFactorOffset, shells[i].nNormFactors) = shells[i].normalizationFactors;

        for (size_t j = 0; j < shells[i].nAOs; ++j)
        {
            angularMomentum.col(aoOffset + j) = shells[i].angularMomentum.col(j);
        }

        primOffset += shells[i].nPrimitives;
        aoOffset += shells[i].nAOs;
        normFactorOffset += shells[i].nNormFactors;
    }
}


void IntegralEngine::overlap(size_t shellA_idx, size_t shellB_idx, Eigen::MatrixXd& S) const
{
    const size_t nprimA = nPrimitives[shellA_idx];
    const size_t nprimB = nPrimitives[shellB_idx];
    const size_t naoA   = nAOs[shellA_idx];
    const size_t naoB   = nAOs[shellB_idx];

    const unsigned lA = lShell[shellA_idx];
    const unsigned lB = lShell[shellB_idx];

    // Allocate scratch space for Hermite coefficients once per shell pair.
    EBuffer Ex(lA, lB);
    EBuffer Ey(lA, lB);
    EBuffer Ez(lA, lB);

    const Eigen::Vector3d& cA = centers.col(shellA_idx);
    const Eigen::Vector3d& cB = centers.col(shellB_idx);

    // Loop over primitive pairs
    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        // const double alphaA = exps[shellA.primOffset + pA];
        const double alphaA = alpha[primitiveOffsets[shellA_idx] + pA];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            // const double alphaB = exps[shellB.primOffset + pB];
            const double alphaB = alpha[primitiveOffsets[shellB_idx] + pB];


            // Compute data for this primitive pair
            const PrimitivePairData primPair = computePrimitivePairData(alphaA, cA, alphaB, cB);
            double prefactor                 = std::pow(M_PI / primPair.p, 1.5) * primPair.K;

            computeHermiteCoeffs(lShell[shellA_idx], lShell[shellB_idx], primPair.p, primPair.PAx, primPair.PBx, Ex);
            computeHermiteCoeffs(lShell[shellA_idx], lShell[shellB_idx], primPair.p, primPair.PAy, primPair.PBy, Ey);
            computeHermiteCoeffs(lShell[shellA_idx], lShell[shellB_idx], primPair.p, primPair.PAz, primPair.PBz, Ez);

            const double coeffA    = coeff[primitiveOffsets[shellA_idx] + pA];
            const double coeffb    = coeff[primitiveOffsets[shellB_idx] + pB];
            const double coeffProd = coeffA * coeffb;

            // This primitive pair contributes to all AO pairs in the shell pair
            for (size_t a = 0; a < naoA; ++a)
            {
                const double normA = normalizationFactors[normFactorOffsets[shellA_idx] + pA * naoA + a];

                const Eigen::Vector3i& am1 = angularMomentum.col(aoOffsets[shellA_idx] + a);
                const unsigned l1 = am1.x(), m1 = am1.y(), n1 = am1.z();

                for (size_t b = 0; b < naoB; ++b)
                {
                    const double normB = normalizationFactors[normFactorOffsets[shellB_idx] + pB * naoB + b];

                    const Eigen::Vector3i& am2 = angularMomentum.col(aoOffsets[shellB_idx] + b);
                    const unsigned l2 = am2.x(), m2 = am2.y(), n2 = am2.z();

                    S(aoOffsets[shellA_idx] + a,
                      aoOffsets[shellB_idx]
                          + b) += normA * normB * coeffProd * prefactor * Ex(l1, l2, 0) * Ey(m1, m2, 0) * Ez(n1, n2, 0);
                }
            }
        }
    }
}


Eigen::MatrixXd IntegralEngine::overlapMatrix() const
{
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(nAOsTotal, nAOsTotal);

#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < nShells; ++i)
    {
        // Only compute upper triangle as S is symmetric.
        for (size_t j = i; j < nShells; ++j) { overlap(i, j, S); }
    }

    S = S.selfadjointView<Eigen::Upper>();

    return S;
}


void IntegralEngine::kinetic(size_t shellA_idx, size_t shellB_idx, Eigen::MatrixXd& T) const
{
    const size_t nprimA = nPrimitives[shellA_idx];
    const size_t nprimB = nPrimitives[shellB_idx];
    const size_t naoA   = nAOs[shellA_idx];
    const size_t naoB   = nAOs[shellB_idx];

    const size_t primOffsetA = primitiveOffsets[shellA_idx];
    const size_t primOffsetB = primitiveOffsets[shellB_idx];
    const size_t aoOffsetA   = aoOffsets[shellA_idx];
    const size_t aoOffsetB   = aoOffsets[shellB_idx];

    const unsigned lA = lShell[shellA_idx];
    const unsigned lB = lShell[shellB_idx];

    // Allocate scratch space. Max angular momentum needed is L+2 for kinetic integrals.
    EBuffer Ex(lA, lB + 2);
    EBuffer Ey(lA, lB + 2);
    EBuffer Ez(lA, lB + 2);

    const Eigen::Vector3d& cA = centers.col(shellA_idx);
    const Eigen::Vector3d& cB = centers.col(shellB_idx);

    // Loop over primitive pairs
    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double alphaA = alpha[primOffsetA + pA];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double alphaB = alpha[primOffsetB + pB];

            const PrimitivePairData primPair = computePrimitivePairData(alphaA, cA, alphaB, cB);
            double prefactor                 = std::pow(M_PI / primPair.p, 1.5) * primPair.K;

            computeHermiteCoeffs(lA, lB + 2, primPair.p, primPair.PAx, primPair.PBx, Ex);
            computeHermiteCoeffs(lA, lB + 2, primPair.p, primPair.PAy, primPair.PBy, Ey);
            computeHermiteCoeffs(lA, lB + 2, primPair.p, primPair.PAz, primPair.PBz, Ez);

            const double coeffA    = coeff[primOffsetA + pA];
            const double coeffB    = coeff[primOffsetB + pB];
            const double coeffProd = coeffA * coeffB;

            for (size_t a = 0; a < naoA; ++a)
            {
                const Eigen::Vector3i& am1 = angularMomentum.col(aoOffsetA + a);
                const unsigned l1 = am1.x(), m1 = am1.y(), n1 = am1.z();
                const double normA = normalizationFactors[normFactorOffsets[shellA_idx] + pA * naoA + a];

                for (size_t b = 0; b < naoB; ++b)
                {
                    const Eigen::Vector3i& am2 = angularMomentum.col(aoOffsetB + b);
                    const unsigned l2 = am2.x(), m2 = am2.y(), n2 = am2.z();
                    const double normB = normalizationFactors[normFactorOffsets[shellB_idx] + pB * naoB + b];

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

                    T(aoOffsetA + a, aoOffsetB + b) += normA * normB * coeffProd * prefactor * (term1 + term2 + term3);
                }
            }
        }
    }
}


Eigen::MatrixXd IntegralEngine::kineticMatrix() const
{

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(nAOsTotal, nAOsTotal);

#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < nShells; ++i)
    {
        // only compute upper triangle as S is symmetric
        for (size_t j = i; j < nShells; ++j) { IntegralEngine::kinetic(i, j, T); }
    }

    T = T.selfadjointView<Eigen::Upper>();

    return T;
}


void IntegralEngine::nuclearAttraction(
    const std::vector<Atom>& geometry, size_t shellA_idx, size_t shellB_idx, Eigen::MatrixXd& V
) const
{
    const size_t nprimA = nPrimitives[shellA_idx];
    const size_t nprimB = nPrimitives[shellB_idx];
    const size_t naoA   = nAOs[shellA_idx];
    const size_t naoB   = nAOs[shellB_idx];

    const size_t primOffsetA = primitiveOffsets[shellA_idx];
    const size_t primOffsetB = primitiveOffsets[shellB_idx];
    const size_t aoOffsetA   = aoOffsets[shellA_idx];
    const size_t aoOffsetB   = aoOffsets[shellB_idx];

    const unsigned lA = lShell[shellA_idx];
    const unsigned lB = lShell[shellB_idx];

    const unsigned max_L_total = lA + lB;

    // Allocate scratch space
    RBuffer R_buffer(max_L_total);
    std::vector<double> F(R_buffer.max_l + 1);
    EBuffer Ex(lA, lB);
    EBuffer Ey(lA, lB);
    EBuffer Ez(lA, lB);

    const Vec3& cA = centers.col(shellA_idx);
    const Vec3& cB = centers.col(shellB_idx);

    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double alphaA = alpha[primOffsetA + pA];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double alphaB = alpha[primOffsetB + pB];

            const PrimitivePairData primPair = computePrimitivePairData(alphaA, cA, alphaB, cB);
            double prefactor                 = (2.0 * M_PI / primPair.p) * primPair.K;

            computeHermiteCoeffs(lA, lB, primPair.p, primPair.PAx, primPair.PBx, Ex);
            computeHermiteCoeffs(lA, lB, primPair.p, primPair.PAy, primPair.PBy, Ey);
            computeHermiteCoeffs(lA, lB, primPair.p, primPair.PAz, primPair.PBz, Ez);

            const double coeffA    = coeff[primOffsetA + pA];
            const double coeffB    = coeff[primOffsetB + pB];
            const double coeffProd = coeffA * coeffB;

            for (const auto& atom : geometry)
            {
                const Eigen::Vector3d& C = atom.coords;
                const auto Z             = static_cast<double>(atom.atomicNumber);

                const Eigen::Vector3d PC = Vec3(primPair.Px, primPair.Py, primPair.Pz) - C;
                const double T           = primPair.p * PC.squaredNorm();

                Boys::calculateBoys(R_buffer.max_l, T, F);
                computeAuxiliaryIntegrals(primPair.p, PC, F, R_buffer);

                for (size_t a = 0; a < naoA; ++a)
                {
                    const Eigen::Vector3i& am1 = angularMomentum.col(aoOffsetA + a);
                    const unsigned l1 = am1.x(), m1 = am1.y(), n1 = am1.z();
                    const double normA = normalizationFactors[normFactorOffsets[shellA_idx] + pA * naoA + a];

                    for (size_t b = 0; b < naoB; ++b)
                    {
                        const Eigen::Vector3i& am2 = angularMomentum.col(aoOffsetB + b);
                        const unsigned l2 = am2.x(), m2 = am2.y(), n2 = am2.z();
                        const double normB = normalizationFactors[normFactorOffsets[shellB_idx] + pB * naoB + b];

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

                        V(aoOffsetA + a, aoOffsetB + b) -= normA * normB * coeffProd * prefactor * Z * sum;
                    }
                }
            }
        }
    }
}


Eigen::MatrixXd IntegralEngine::nuclearAttractionMatrix(const std::vector<Atom>& geometry) const
{
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(nAOsTotal, nAOsTotal);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < nShells; ++i)
    {
        // only compute upper triangle as S is symmetric
        for (size_t j = i; j < nShells; ++j) { IntegralEngine::nuclearAttraction(geometry, i, j, V); }
    }

    V = V.selfadjointView<Eigen::Upper>();

    return V;
}


void IntegralEngine::electronRepulsion(
    size_t shellA_idx, size_t shellB_idx, size_t shellC_idx, size_t shellD_idx, ElectronRepulsionTensor& G
) const
{
    const size_t nprimA = nPrimitives[shellA_idx], naoA = nAOs[shellA_idx];
    const size_t nprimB = nPrimitives[shellB_idx], naoB = nAOs[shellB_idx];
    const size_t nprimC = nPrimitives[shellC_idx], naoC = nAOs[shellC_idx];
    const size_t nprimD = nPrimitives[shellD_idx], naoD = nAOs[shellD_idx];

    const size_t primOffsetA = primitiveOffsets[shellA_idx];
    const size_t primOffsetB = primitiveOffsets[shellB_idx];
    const size_t primOffsetC = primitiveOffsets[shellC_idx];
    const size_t primOffsetD = primitiveOffsets[shellD_idx];

    const size_t aoOffsetA = aoOffsets[shellA_idx];
    const size_t aoOffsetB = aoOffsets[shellB_idx];
    const size_t aoOffsetC = aoOffsets[shellC_idx];
    const size_t aoOffsetD = aoOffsets[shellD_idx];

    const unsigned lA = lShell[shellA_idx];
    const unsigned lB = lShell[shellB_idx];
    const unsigned lC = lShell[shellC_idx];
    const unsigned lD = lShell[shellD_idx];

    // Pre-calculate and cache all (cd) pair data
    const size_t size = nprimC * nprimD;
    std::vector<EBuffer> Ex_cd_vec(size, EBuffer(lC, lD));
    std::vector<EBuffer> Ey_cd_vec(size, EBuffer(lC, lD));
    std::vector<EBuffer> Ez_cd_vec(size, EBuffer(lC, lD));
    std::vector<PrimitivePairData> primPairs_cd;
    std::vector<size_t> pC_indices, pD_indices;
    primPairs_cd.reserve(size);
    pC_indices.reserve(size);
    pD_indices.reserve(size);

    const double* alphaC = &alpha[primOffsetC];
    const double* alphaD = &alpha[primOffsetD];

    const Vec3& cC = centers.col(shellC_idx);
    const Vec3& cD = centers.col(shellD_idx);

    for (size_t pC = 0; pC < nprimC; ++pC)
    {
        for (size_t pD = 0; pD < nprimD; ++pD)
        {
            pC_indices.emplace_back(pC);
            pD_indices.emplace_back(pD);
            primPairs_cd.emplace_back(computePrimitivePairData(alphaC[pC], cC, alphaD[pD], cD));
            const size_t idx = (pC * nprimD) + pD;

            const auto& primPair_cd = primPairs_cd[idx];
            computeHermiteCoeffs(lC, lD, primPair_cd.p, primPair_cd.PAx, primPair_cd.PBx, Ex_cd_vec[idx]);
            computeHermiteCoeffs(lC, lD, primPair_cd.p, primPair_cd.PAy, primPair_cd.PBy, Ey_cd_vec[idx]);
            computeHermiteCoeffs(lC, lD, primPair_cd.p, primPair_cd.PAz, primPair_cd.PBz, Ez_cd_vec[idx]);
        }
    }

    // Allocate scratch space for (ab) pair and other intermediates
    EBuffer Ex_ab(lA, lB), Ey_ab(lA, lB), Ez_ab(lA, lB);
    const unsigned max_L_total = lA + lB + lC + lD;
    RBuffer R_buffer(max_L_total);
    std::vector<double> F(R_buffer.max_l + 1);

    const double* alphaA = &alpha[primOffsetA];
    const double* alphaB = &alpha[primOffsetB];

    const Vec3& cA = centers.col(shellA_idx);
    const Vec3& cB = centers.col(shellB_idx);

    std::vector<double> prim_integrals(naoA * naoB * naoC * naoD);

    // Loop over (ab) primitive pairs
    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double coeffA = coeff[primOffsetA + pA];
        const double* normA = &normalizationFactors[normFactorOffsets[shellA_idx] + pA * naoA];

        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double coeffB = coeff[primOffsetB + pB];
            const double* normB = &normalizationFactors[normFactorOffsets[shellB_idx] + pB * naoB];

            const PrimitivePairData primPair_ab = computePrimitivePairData(alphaA[pA], cA, alphaB[pB], cB);

            computeHermiteCoeffs(lA, lB, primPair_ab.p, primPair_ab.PAx, primPair_ab.PBx, Ex_ab);
            computeHermiteCoeffs(lA, lB, primPair_ab.p, primPair_ab.PAy, primPair_ab.PBy, Ey_ab);
            computeHermiteCoeffs(lA, lB, primPair_ab.p, primPair_ab.PAz, primPair_ab.PBz, Ez_ab);

            // Loop over cached (cd) pairs
            for (size_t cd_idx = 0; cd_idx < primPairs_cd.size(); ++cd_idx)
            {
                const auto& primPair_cd = primPairs_cd[cd_idx];
                const size_t pC         = pC_indices[cd_idx];
                const size_t pD         = pD_indices[cd_idx];
                const EBuffer& Ex_cd    = Ex_cd_vec[cd_idx];
                const EBuffer& Ey_cd    = Ey_cd_vec[cd_idx];
                const EBuffer& Ez_cd    = Ez_cd_vec[cd_idx];

                const double prefactor = (2.0 * std::pow(M_PI, 2.5)
                                          / (primPair_ab.p * primPair_cd.p * std::sqrt(primPair_ab.p + primPair_cd.p)))
                                       * primPair_ab.K * primPair_cd.K;

                const double delta = (primPair_ab.p * primPair_cd.p) / (primPair_ab.p + primPair_cd.p);
                const auto P_ab    = Vec3(primPair_ab.Px, primPair_ab.Py, primPair_ab.Pz);
                const auto P_cd    = Vec3(primPair_cd.Px, primPair_cd.Py, primPair_cd.Pz);
                const double T     = delta * (P_ab - P_cd).squaredNorm();

                Boys::calculateBoys(R_buffer.max_l, T, F);
                computeAuxiliaryIntegrals(delta, P_ab - P_cd, F, R_buffer);

                const unsigned naoCD  = naoC * naoD;
                const unsigned naoBCD = naoB * naoCD;

                for (size_t a = 0; a < naoA; ++a)
                {
                    const size_t i             = aoOffsetA + a;
                    const size_t big_I_base    = i * (i + 1) / 2;
                    const Eigen::Vector3i& am1 = angularMomentum.col(aoOffsetA + a);
                    const unsigned l1 = am1.x(), m1 = am1.y(), n1 = am1.z();

                    for (size_t b = 0; b < naoB; ++b)
                    {
                        const size_t j = aoOffsetB + b;
                        if (j > i) // Exploit permutational symmetry
                            continue;
                        const size_t big_I         = big_I_base + j;
                        const Eigen::Vector3i& am2 = angularMomentum.col(aoOffsetB + b);
                        const unsigned l2 = am2.x(), m2 = am2.y(), n2 = am2.z();

                        const Eigen::Map<const Eigen::ArrayXd> Ex_ab_array(&Ex_ab(l1, l2, 0), l1 + l2 + 1);
                        const Eigen::Map<const Eigen::ArrayXd> Ey_ab_array(&Ey_ab(m1, m2, 0), m1 + m2 + 1);
                        const Eigen::Map<const Eigen::ArrayXd> Ez_ab_array(&Ez_ab(n1, n2, 0), n1 + n2 + 1);

                        for (size_t c = 0; c < naoC; ++c)
                        {
                            const size_t k             = aoOffsetC + c;
                            const size_t big_K_base    = k * (k + 1) / 2;
                            const Eigen::Vector3i& am3 = angularMomentum.col(aoOffsetC + c);
                            const unsigned l3 = am3.x(), m3 = am3.y(), n3 = am3.z();

                            for (size_t d = 0; d < naoD; ++d)
                            {
                                const size_t l = aoOffsetD + d;

                                if (l > k)
                                    continue;
                                const size_t big_K = big_K_base + l;
                                if ((shellA_idx == shellC_idx && shellB_idx == shellD_idx) && big_I < big_K)
                                    continue;
                                const Eigen::Vector3i& am4 = angularMomentum.col(aoOffsetD + d);
                                const unsigned l4 = am4.x(), m4 = am4.y(), n4 = am4.z();

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

                const double coeffC = coeff[primOffsetC + pC];
                const double coeffD = coeff[primOffsetD + pD];
                const double* normC = &normalizationFactors[normFactorOffsets[shellC_idx] + pC * naoC];
                const double* normD = &normalizationFactors[normFactorOffsets[shellD_idx] + pD * naoD];

                const double coeffProd = coeffA * coeffB * coeffC * coeffD;

                // Contract Integrals
                for (size_t a = 0; a < naoA; ++a)
                {
                    const size_t i          = aoOffsetA + a;
                    const size_t big_I_base = i * (i + 1) / 2;
                    const double nA         = normA[a];

                    for (size_t b = 0; b < naoB; ++b)
                    {
                        const size_t j = aoOffsetB + b;
                        if (j > i)
                            continue;
                        const size_t big_I = big_I_base + j;

                        const double nB = normB[b];

                        for (size_t c = 0; c < naoC; ++c)
                        {
                            const size_t k          = aoOffsetC + c;
                            const size_t big_K_base = k * (k + 1) / 2;

                            const double nC = normC[c];

                            for (size_t d = 0; d < naoD; ++d)
                            {
                                const size_t l = aoOffsetD + d;
                                if (l > k)
                                    continue;
                                const size_t big_K = big_K_base + l;

                                if ((shellA_idx == shellC_idx && shellB_idx == shellD_idx) && big_I < big_K)
                                    continue;

                                const double nD = normD[d];

                                const size_t index = big_I >= big_K ? (big_I * (big_I + 1) / 2 + big_K)
                                                                    : (big_K * (big_K + 1) / 2 + big_I);

                                G[index] += nA * nB * nC * nD * coeffProd * prefactor
                                          * prim_integrals[a * naoBCD + b * naoCD + c * naoD + d];
                            }
                        }
                    }
                }
            }
        }
    }
}


ElectronRepulsionTensor IntegralEngine::electronRepulsionTensor(double threshold) const
{
    Eigen::MatrixXd Q(nShells, nShells);
    ElectronRepulsionTensor Vee(nAOsTotal);

#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < nShells; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            IntegralEngine::electronRepulsion(i, j, i, j, Vee);

            double maxVal = 0;
            for (int a = 0; a < nAOs[i]; ++a)
            {
                for (int b = 0; b < nAOs[j]; ++b)
                {
                    double val = Vee(aoOffsets[i] + a, aoOffsets[j] + b, aoOffsets[i] + a, aoOffsets[j] + b);
                    maxVal     = std::max(maxVal, val);
                }
            }

            Q(j, i) = maxVal;
        }
    }

    Q = Q.selfadjointView<Eigen::Upper>();
    Q = Q.cwiseSqrt();

#pragma omp parallel for collapse(4) schedule(dynamic, 64)
    for (size_t i = 0; i < nShells; ++i)
    {
        for (size_t j = 0; j < nShells; ++j)
        {
            for (size_t k = 0; k < nShells; ++k)
            {
                for (size_t l = 0; l < nShells; ++l)
                {
                    // no need to recalculate identical elements of the tensor (enforce quartet symmetry)
                    if (j > i || l > k || (i * (i + 1) / 2 + j) < (k * (k + 1) / 2 + l))
                        continue;

                    // make sure we did not already calculate this integral in Schwartz screening loop
                    if ((i == k && j == l) || (i == l && j == k))
                        continue;

                    // apply Schwartz screening at the shell level
                    if (Q(i, j) * Q(k, l) < threshold)
                        continue;

                    IntegralEngine::electronRepulsion(i, j, k, l, Vee);
                }
            }
        }
    }

    return Vee;
}


Eigen::MatrixXd IntegralEngine::schwartzScreeningMatrix() const
{
    Eigen::MatrixXd Q(nShells, nShells);

#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < nShells; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            std::vector<double> ERIs(nAOs[i] * nAOs[j] * nAOs[i] * nAOs[j], 0.0);
            IntegralEngine::electronRepulsion(i, j, i, j, ERIs);

            const size_t aOffset = nAOs[j] * nAOs[i] * nAOs[j] + nAOs[j];
            const size_t bOffset = nAOs[i] * nAOs[j] + 1;

            double maxVal = 0;
            for (int a = 0; a < nAOs[i]; ++a)
            {
                for (int b = 0; b < nAOs[j]; ++b)
                {
                    double val = ERIs[a * aOffset + b * bOffset];
                    maxVal     = std::max(maxVal, val);
                }
            }

            Q(j, i) = maxVal;
        }
    }

    Q = Q.selfadjointView<Eigen::Upper>();
    Q = Q.cwiseSqrt();

    return Q;
}


void IntegralEngine::electronRepulsion(
    size_t shellA_idx, size_t shellB_idx, size_t shellC_idx, size_t shellD_idx, std::vector<double>& G
) const
{
    const size_t nprimA = nPrimitives[shellA_idx], naoA = nAOs[shellA_idx];
    const size_t nprimB = nPrimitives[shellB_idx], naoB = nAOs[shellB_idx];
    const size_t nprimC = nPrimitives[shellC_idx], naoC = nAOs[shellC_idx];
    const size_t nprimD = nPrimitives[shellD_idx], naoD = nAOs[shellD_idx];

    const size_t primOffsetA = primitiveOffsets[shellA_idx];
    const size_t primOffsetB = primitiveOffsets[shellB_idx];
    const size_t primOffsetC = primitiveOffsets[shellC_idx];
    const size_t primOffsetD = primitiveOffsets[shellD_idx];

    const size_t aoOffsetA = aoOffsets[shellA_idx];
    const size_t aoOffsetB = aoOffsets[shellB_idx];
    const size_t aoOffsetC = aoOffsets[shellC_idx];
    const size_t aoOffsetD = aoOffsets[shellD_idx];

    const unsigned lA = lShell[shellA_idx];
    const unsigned lB = lShell[shellB_idx];
    const unsigned lC = lShell[shellC_idx];
    const unsigned lD = lShell[shellD_idx];

    // Pre-calculate and cache all (cd) pair data
    const size_t size = nprimC * nprimD;
    std::vector<EBuffer> Ex_cd_vec(size, EBuffer(lC, lD));
    std::vector<EBuffer> Ey_cd_vec(size, EBuffer(lC, lD));
    std::vector<EBuffer> Ez_cd_vec(size, EBuffer(lC, lD));
    std::vector<PrimitivePairData> primPairs_cd;
    std::vector<size_t> pC_indices, pD_indices;
    primPairs_cd.reserve(size);
    pC_indices.reserve(size);
    pD_indices.reserve(size);

    const double* alphaC = &alpha[primOffsetC];
    const double* alphaD = &alpha[primOffsetD];

    const Vec3& cC = centers.col(shellC_idx);
    const Vec3& cD = centers.col(shellD_idx);

    for (size_t pC = 0; pC < nprimC; ++pC)
    {
        for (size_t pD = 0; pD < nprimD; ++pD)
        {
            pC_indices.emplace_back(pC);
            pD_indices.emplace_back(pD);
            primPairs_cd.emplace_back(computePrimitivePairData(alphaC[pC], cC, alphaD[pD], cD));
            const size_t idx = (pC * nprimD) + pD;

            const auto& primPair_cd = primPairs_cd[idx];
            computeHermiteCoeffs(lC, lD, primPair_cd.p, primPair_cd.PAx, primPair_cd.PBx, Ex_cd_vec[idx]);
            computeHermiteCoeffs(lC, lD, primPair_cd.p, primPair_cd.PAy, primPair_cd.PBy, Ey_cd_vec[idx]);
            computeHermiteCoeffs(lC, lD, primPair_cd.p, primPair_cd.PAz, primPair_cd.PBz, Ez_cd_vec[idx]);
        }
    }

    // Allocate scratch space for (ab) pair and other intermediates
    EBuffer Ex_ab(lA, lB), Ey_ab(lA, lB), Ez_ab(lA, lB);
    const unsigned max_L_total = lA + lB + lC + lD;
    RBuffer R_buffer(max_L_total);
    std::vector<double> F(R_buffer.max_l + 1);

    const double* alphaA = &alpha[primOffsetA];
    const double* alphaB = &alpha[primOffsetB];

    const Vec3& cA = centers.col(shellA_idx);
    const Vec3& cB = centers.col(shellB_idx);


    // Create a temporary tensor to store all primitive integrals for this quartet
    std::vector<double> prim_integrals(naoA * naoB * naoC * naoD);

    // Loop over (ab) primitive pairs
    for (size_t pA = 0; pA < nprimA; ++pA)
    {
        const double coeffA = coeff[primOffsetA + pA];
        const double* normA = &normalizationFactors[normFactorOffsets[shellA_idx] + pA * naoA];


        for (size_t pB = 0; pB < nprimB; ++pB)
        {
            const double coeffB = coeff[primOffsetB + pB];
            const double* normB = &normalizationFactors[normFactorOffsets[shellB_idx] + pB * naoB];

            const PrimitivePairData primPair_ab = computePrimitivePairData(alphaA[pA], cA, alphaB[pB], cB);

            computeHermiteCoeffs(lA, lB, primPair_ab.p, primPair_ab.PAx, primPair_ab.PBx, Ex_ab);
            computeHermiteCoeffs(lA, lB, primPair_ab.p, primPair_ab.PAy, primPair_ab.PBy, Ey_ab);
            computeHermiteCoeffs(lA, lB, primPair_ab.p, primPair_ab.PAz, primPair_ab.PBz, Ez_ab);

            // Loop over cached (cd) pairs
            for (size_t cd_idx = 0; cd_idx < primPairs_cd.size(); ++cd_idx)
            {
                const auto& primPair_cd = primPairs_cd[cd_idx];
                const size_t pC         = pC_indices[cd_idx];
                const size_t pD         = pD_indices[cd_idx];
                const EBuffer& Ex_cd    = Ex_cd_vec[cd_idx];
                const EBuffer& Ey_cd    = Ey_cd_vec[cd_idx];
                const EBuffer& Ez_cd    = Ez_cd_vec[cd_idx];

                const double prefactor = (2.0 * std::pow(M_PI, 2.5)
                                          / (primPair_ab.p * primPair_cd.p * std::sqrt(primPair_ab.p + primPair_cd.p)))
                                       * primPair_ab.K * primPair_cd.K;

                const double delta = (primPair_ab.p * primPair_cd.p) / (primPair_ab.p + primPair_cd.p);
                const auto P_ab    = Vec3(primPair_ab.Px, primPair_ab.Py, primPair_ab.Pz);
                const auto P_cd    = Vec3(primPair_cd.Px, primPair_cd.Py, primPair_cd.Pz);
                const double T     = delta * (P_ab - P_cd).squaredNorm();

                Boys::calculateBoys(R_buffer.max_l, T, F);
                computeAuxiliaryIntegrals(delta, P_ab - P_cd, F, R_buffer);

                const unsigned naoCD  = naoC * naoD;
                const unsigned naoBCD = naoB * naoCD;

                for (size_t a = 0; a < naoA; ++a)
                {
                    const size_t i             = aoOffsetA + a;
                    const size_t big_I_base    = i * (i + 1) / 2;
                    const Eigen::Vector3i& am1 = angularMomentum.col(aoOffsetA + a);
                    const unsigned l1 = am1.x(), m1 = am1.y(), n1 = am1.z();

                    for (size_t b = 0; b < naoB; ++b)
                    {
                        const size_t j = aoOffsetB + b;
                        if (j > i) // Exploit permutational symmetry
                            continue;
                        const size_t big_I         = big_I_base + j;
                        const Eigen::Vector3i& am2 = angularMomentum.col(aoOffsetB + b);
                        const unsigned l2 = am2.x(), m2 = am2.y(), n2 = am2.z();

                        const Eigen::Map<const Eigen::ArrayXd> Ex_ab_array(&Ex_ab(l1, l2, 0), l1 + l2 + 1);
                        const Eigen::Map<const Eigen::ArrayXd> Ey_ab_array(&Ey_ab(m1, m2, 0), m1 + m2 + 1);
                        const Eigen::Map<const Eigen::ArrayXd> Ez_ab_array(&Ez_ab(n1, n2, 0), n1 + n2 + 1);

                        for (size_t c = 0; c < naoC; ++c)
                        {
                            const size_t k             = aoOffsetC + c;
                            const size_t big_K_base    = k * (k + 1) / 2;
                            const Eigen::Vector3i& am3 = angularMomentum.col(aoOffsetC + c);
                            const unsigned l3 = am3.x(), m3 = am3.y(), n3 = am3.z();

                            for (size_t d = 0; d < naoD; ++d)
                            {
                                const size_t l = aoOffsetD + d;

                                if (l > k)
                                    continue;
                                const size_t big_K = big_K_base + l;
                                if ((shellA_idx == shellC_idx && shellB_idx == shellD_idx) && big_I < big_K)
                                    continue;
                                const Eigen::Vector3i& am4 = angularMomentum.col(aoOffsetD + d);
                                const unsigned l4 = am4.x(), m4 = am4.y(), n4 = am4.z();

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

                const double coeffC = coeff[primOffsetC + pC];
                const double coeffD = coeff[primOffsetD + pD];
                const double* normC = &normalizationFactors[normFactorOffsets[shellC_idx] + pC * naoC];
                const double* normD = &normalizationFactors[normFactorOffsets[shellD_idx] + pD * naoD];

                const double coeffProd = coeffA * coeffB * coeffC * coeffD;

                // Contract Integrals
                for (size_t a = 0; a < naoA; ++a)
                {
                    const size_t i          = aoOffsetA + a;
                    const size_t big_I_base = i * (i + 1) / 2;
                    const double nA         = normA[a];

                    for (size_t b = 0; b < naoB; ++b)
                    {
                        const size_t j = aoOffsetB + b;
                        if (j > i)
                            continue;
                        const size_t big_I = big_I_base + j;

                        const double nB = normB[b];

                        for (size_t c = 0; c < naoC; ++c)
                        {
                            const size_t k          = aoOffsetC + c;
                            const size_t big_K_base = k * (k + 1) / 2;

                            const double nC = normC[c];

                            for (size_t d = 0; d < naoD; ++d)
                            {
                                const size_t l = aoOffsetD + d;
                                if (l > k)
                                    continue;
                                const size_t big_K = big_K_base + l;

                                if ((shellA_idx == shellC_idx && shellB_idx == shellD_idx) && big_I < big_K)
                                    continue;

                                const double nD = normD[d];

                                G[a * naoBCD + b * naoCD + c * naoD + d] += nA * nB * nC * nD * coeffProd * prefactor
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
    double alphaA, const Vec3& centerA, double alphaB, const Vec3& centerB
)
{
    double p    = alphaA + alphaB;
    double oo_p = 1.0 / p; // one over p

    const Vec3 P  = (alphaA * centerA + alphaB * centerB) * oo_p;
    const Vec3 PA = P - centerA;
    const Vec3 PB = P - centerB;

    double AmB2 = (centerA - centerB).squaredNorm();
    double K    = std::exp(-alphaA * alphaB * oo_p * AmB2);

    return {
        .p   = p,
        .K   = K,
        .Px  = P.x(),
        .Py  = P.y(),
        .Pz  = P.z(),
        .PAx = PA.x(),
        .PAy = PA.y(),
        .PAz = PA.z(),
        .PBx = PB.x(),
        .PBy = PB.y(),
        .PBz = PB.z()
    };
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

                E_ij_ptr[t] = val;
            }
        }
    }
}


void IntegralEngine::computeAuxiliaryIntegrals(double p, const Vec3& PC, std::span<double> F, RBuffer& R_buffer)
{
    // This function implements the recurrence relations for R integrals:
    // R_{t,u,v}^{n} = t-1 * R_{t-2,u,v}^{n+1} + (PC)_x * R_{t-1,u,v}^{n+1} (if t>0)
    // R_{t,u,v}^{n} = u-1 * R_{t,u-2,v}^{n+1} + (PC)_y * R_{t,u-1,v}^{n+1} (if u>0)
    // R_{t,u,v}^{n} = v-1 * R_{t,u,v-2}^{n+1} + (PC)_z * R_{t,u,v-1}^{n+1} (if v>0)
    // Base case: R_{0,0,0}^{n} = (-2p)^n * F[n]

    const unsigned max_l = R_buffer.max_l;

    const double neg_2p      = -2.0 * p;
    double neg_2p_pow_n      = 1.0;
    double* const R_000n_ptr = &R_buffer(0, 0, 0, 0);
    size_t index             = 0;
    for (unsigned n = 0; n <= max_l; ++n)
    {
        R_000n_ptr[index] = neg_2p_pow_n * F[n];
        index += R_buffer.stride_n;
        neg_2p_pow_n *= neg_2p;
    }

    unsigned n = max_l;
    while (n-- > 0)
    {
        // Build up t dimension from t=1 up to max_t
        const double* const R_np1_t_ptr = &R_buffer(0, 0, 0, n + 1);
        double* const R_n_t_ptr         = &R_buffer(0, 0, 0, n);
        for (unsigned t = 1; t <= max_l - n; ++t)
        {
            const double val = PC.x() * R_np1_t_ptr[t - 1];
            if (t >= 2)
                R_n_t_ptr[t] = val + (t - 1) * R_np1_t_ptr[t - 2];
            else
                R_n_t_ptr[t] = val;
        }

        // Build up u dimension
        for (unsigned u = 1; u <= max_l - n; ++u)
        {
            double* const R_n_ptr             = &R_buffer(0, u, 0, n);
            const double* const R_np1_um1_ptr = &R_buffer(0, u - 1, 0, n + 1);

            if (u >= 2)
            {
                const double* R_np1_um2_ptr = &R_buffer(0, u - 2, 0, n + 1);
                for (unsigned t = 0; t <= max_l - n - u; ++t)
                {
                    R_n_ptr[t] = PC.y() * R_np1_um1_ptr[t] + (u - 1) * R_np1_um2_ptr[t];
                }
            }
            else
            { // u is exactly 1
                for (unsigned t = 0; t <= max_l - n - 1; ++t) { R_n_ptr[t] = PC.y() * R_np1_um1_ptr[t]; }
            }
        }

        // Build up v dimension
        for (unsigned v = 1; v <= max_l - n; ++v)
        {
            for (unsigned u = 0; u <= max_l - n - v; ++u)
            {
                // Get pointers to the start of the relevant rows for efficiency
                double* const R_n_ptr             = &R_buffer(0, u, v, n);
                const double* const R_np1_vm1_ptr = &R_buffer(0, u, v - 1, n + 1);

                if (v >= 2)
                {
                    const double* const R_np1_vm2_ptr = &R_buffer(0, u, v - 2, n + 1);
                    for (unsigned t = 0; t <= max_l - n - u - v; ++t)
                    {
                        R_n_ptr[t] = PC.z() * R_np1_vm1_ptr[t] + (v - 1) * R_np1_vm2_ptr[t];
                    }
                }
                else
                { // v is exactly 1
                    for (unsigned t = 0; t <= max_l - n - u - 1; ++t) { R_n_ptr[t] = PC.z() * R_np1_vm1_ptr[t]; }
                }
            }
        }
    }
}
