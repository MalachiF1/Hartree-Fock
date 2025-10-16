#include "PRISMEngine.hpp"

#include "Basis.hpp"
#include "Boys.hpp"

#include <algorithm>
#include <fmt/core.h>
#include <ranges>

PRISMEngine::PRISMEngine(const BasisSet& basis) : shellPairs(computeShellPairData(basis)) {}

void PRISMEngine::sortShellPairs(std::vector<RawShellPair>& shellPairs) const
{
    std::ranges::sort(
        shellPairs,
        [](const RawShellPair& a, const RawShellPair& b)
        {
            if (a.shellA.l != b.shellA.l)
                return a.shellA.l < b.shellA.l;

            if (a.shellB.l != b.shellB.l)
                return a.shellB.l < b.shellB.l;

            const unsigned kAB = a.significantPrimPairs.count();
            const unsigned kCD = b.significantPrimPairs.count();
            return kAB < kCD;
        }
    );
}


void PRISMEngine::assignShellAndPrimPairSignificance(std::vector<RawShellPair>& shellPairs, double cutoff) const
{
    const double C = std::sqrt(8.0 * std::pow(M_PI, 2.5));

    for (auto& sp : shellPairs)
    {
        bool isShellPairSignificant = false;
        const double ABsq           = (sp.shellA.center - sp.shellB.center).squaredNorm();
        size_t primPairIndex        = 0;
        for (const auto& primA : sp.shellA.primitives)
        {
            for (const auto& primB : sp.shellB.primitives)
            {
                const double alpha = primA.exponent;
                const double beta  = primB.exponent;
                const double D_a   = primA.coefficient;
                const double D_b   = primB.coefficient;
                const double G     = std::exp(-ABsq * alpha * beta / (alpha + beta));

                const double prefactor = std::abs(D_a * D_b * C * G);

                if (prefactor < cutoff)
                    sp.significantPrimPairs(primPairIndex) = false;
                else
                    isShellPairSignificant = true;

                primPairIndex++;
            }
        }

        sp.significant = isShellPairSignificant;
    }
}


ShellPairs PRISMEngine::computeShellPairData(const BasisSet& basis) const
{

    const size_t nShells = basis.nShells;
    const auto& shells   = basis.shells;

    // Create raw shell pairs
    std::vector<RawShellPair> rawShellPairs;
    rawShellPairs.reserve(nShells * (nShells + 1) / 2);

    for (size_t i = 0; i < nShells; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            // Ensure that in each pair, the shell with the lower angular momentum comes first
            Shell shellA  = shells[i];
            Shell shellB  = shells[j];
            size_t indexA = i;
            size_t indexB = j;
            if (shellA.l > shellB.l)
            {
                std::swap(shellA, shellB);
                std::swap(indexA, indexB);
            }

            rawShellPairs.emplace_back(shellA, shellB, indexA, indexB, shells[i].nPrimitives * shells[j].nPrimitives);
        }
    }

    // 1e-16 cutoff is temporary, probably should be 0.01 * integral threshold.
    assignShellAndPrimPairSignificance(rawShellPairs, 1e-16);
    sortShellPairs(rawShellPairs);

    ShellPairs shellPairs;
    // Set metadata
    shellPairs.nShellPairs     = 0;
    shellPairs.nPrimitivePairs = 0;

    size_t nNormFactors = 0;
    for (const auto& sp : rawShellPairs)
    {
        if (!sp.significant)
            continue;

        shellPairs.nShellPairs++;

        const size_t nPrimPairs = sp.significantPrimPairs.count();
        shellPairs.nPrimitivePairs += nPrimPairs;
        nNormFactors += nPrimPairs * sp.shellA.nAOs * sp.shellB.nAOs;
    }

    // Allocate per shell-pair data
    shellPairs.indices.resize(2, shellPairs.nShellPairs);
    shellPairs.l.resize(2, shellPairs.nShellPairs);
    shellPairs.A.resize(3, shellPairs.nShellPairs);
    shellPairs.B.resize(3, shellPairs.nShellPairs);
    shellPairs.AB.resize(3, shellPairs.nShellPairs);
    shellPairs.K.resize(shellPairs.nShellPairs);
    shellPairs.primPairOffsets.resize(shellPairs.nShellPairs);
    shellPairs.nAOsA.resize(shellPairs.nShellPairs);
    shellPairs.nAOsB.resize(shellPairs.nShellPairs);

    // Allocate per primitive-pair data
    shellPairs.alpha.resize(2, shellPairs.nPrimitivePairs);
    shellPairs.twoAlpha.resize(2, shellPairs.nPrimitivePairs);
    shellPairs.P.resize(3, shellPairs.nPrimitivePairs);
    shellPairs.invTwoZeta.resize(shellPairs.nPrimitivePairs);
    shellPairs.Up.resize(shellPairs.nPrimitivePairs);

    // Fill data
    size_t primPairIndex  = 0;
    size_t shellPairIndex = 0;
    for (const auto& sp : rawShellPairs)
    {
        if (!sp.significant)
            continue;

        const auto& shellA = sp.shellA;
        const auto& shellB = sp.shellB;

        shellPairs.indices.col(shellPairIndex) = Eigen::Vector2i(sp.indexA, sp.indexB);

        shellPairs.l.col(shellPairIndex)           = Eigen::Vector2i(shellA.l, shellB.l);
        shellPairs.A.col(shellPairIndex)           = shellA.center;
        shellPairs.B.col(shellPairIndex)           = shellB.center;
        shellPairs.AB.col(shellPairIndex)          = shellA.center - shellB.center;
        shellPairs.K(shellPairIndex)               = sp.significantPrimPairs.count();
        shellPairs.primPairOffsets(shellPairIndex) = primPairIndex;
        shellPairs.nAOsA(shellPairIndex)           = shellA.nAOs;
        shellPairs.nAOsB(shellPairIndex)           = shellB.nAOs;

        const double ABsq = shellPairs.AB.col(shellPairIndex).squaredNorm();

        size_t localPrimPairIndex = 0;
        for (size_t j = 0; j < shellA.nPrimitives; ++j)
        {
            for (size_t k = 0; k < shellB.nPrimitives; ++k)
            {
                if (!sp.significantPrimPairs(localPrimPairIndex))
                {
                    localPrimPairIndex++;
                    continue;
                }
                localPrimPairIndex++;

                const Primitive& primA = shellA.primitives[j];
                const Primitive& primB = shellB.primitives[k];

                const double& alphaA    = primA.exponent;
                const double& alphaB    = primB.exponent;
                const double zeta       = alphaA + alphaB;
                const double invZeta    = 1.0 / zeta;
                const double invTwoZeta = 0.5 * invZeta;
                const double G          = std::exp(-alphaA * alphaB * invZeta * ABsq);

                shellPairs.alpha.col(primPairIndex)    = Eigen::Vector2d(alphaA, alphaB);
                shellPairs.twoAlpha.col(primPairIndex) = 2 * Eigen::Vector2d(alphaA, alphaB);
                shellPairs.P.col(primPairIndex)        = (alphaA * shellA.center + alphaB * shellB.center) / zeta;
                shellPairs.invTwoZeta(primPairIndex)   = invTwoZeta;
                shellPairs.Up(primPairIndex) = primA.coefficient * primB.coefficient * G * std::pow(M_PI * invZeta, 1.5)
                                             * std::pow(invTwoZeta, shellA.l + shellB.l);

                primPairIndex++;
            }
        }

        shellPairIndex++;
    }

    // Allocate and compute normalization factors
    shellPairs.normFactorValues.resize(nNormFactors);
    shellPairs.normFactors.reserve(shellPairs.nShellPairs);

    size_t normFactorOffset = 0;
    for (size_t i = 0; i < shellPairs.nShellPairs; ++i)
    {
        const size_t nPrimsPairs = shellPairs.K(i);
        const unsigned lA        = shellPairs.l(0, i);
        const unsigned lB        = shellPairs.l(1, i);

        shellPairs.normFactors.emplace_back(
            &shellPairs.normFactorValues(normFactorOffset), shellPairs.nAOsA(i) * shellPairs.nAOsB(i), nPrimsPairs
        );

        for (size_t j = 0; j < nPrimsPairs; ++j)
        {
            const auto& exponents = shellPairs.alpha.col(shellPairs.primPairOffsets[i] + j);
            const double alpha    = exponents[0];
            const double beta     = exponents[1];

            for (int lAx = lA; lAx >= 0; --lAx)
            {
                for (int lAy = lA - lAx; lAy >= 0; --lAy)
                {
                    int lAz = lA - lAx - lAy;

                    const double prefactor1 = std::pow(2.0 * alpha / M_PI, 0.75);
                    const double prefactor2 = std::pow(4.0 * alpha, (lAx + lAy + lAz) / 2.0);
                    const double denominator = std::sqrt(dfact((2 * lAx) - 1) * dfact((2 * lAy) - 1) * dfact((2 * lAz) - 1));
                    const double normFactorA = prefactor1 * prefactor2 / denominator;

                    for (int lBx = lB; lBx >= 0; --lBx)
                    {
                        for (int lBy = lB - lBx; lBy >= 0; --lBy)
                        {
                            int lBz = lB - lBx - lBy;

                            const double prefactor1  = std::pow(2.0 * beta / M_PI, 0.75);
                            const double prefactor2  = std::pow(4.0 * beta, (lBx + lBy + lBz) / 2.0);
                            const double denominator = std::sqrt(
                                dfact((2 * lBx) - 1) * dfact((2 * lBy) - 1) * dfact((2 * lBz) - 1)
                            );
                            const double normFactorB = prefactor1 * prefactor2 / denominator;

                            const double normProduct = normFactorA * normFactorB;

                            shellPairs.normFactorValues(normFactorOffset) = normProduct;

                            normFactorOffset++;
                        }
                    }
                }
            }
        }
    }

    return shellPairs;
}


Eigen::MatrixXd PRISMEngine::T0(size_t AB, size_t CD) const
{
    const Eigen::Vector3d R = shellPairs.AB.col(AB) - shellPairs.AB.col(CD);
    const double Rsq        = R.squaredNorm();

    const size_t Kab = shellPairs.K(AB);
    const size_t Kcd = shellPairs.K(CD);

    const size_t primOffsetAB = shellPairs.primPairOffsets(AB);
    const size_t primOffsetCD = shellPairs.primPairOffsets(CD);

    const unsigned lAB  = shellPairs.l.col(AB).sum();
    const unsigned lCD  = shellPairs.l.col(CD).sum();
    const unsigned maxL = lAB + lCD;

    Eigen::MatrixXd baseIntegrals = Eigen::MatrixXd(maxL + 1, Kab * Kcd);

    for (size_t i = 0; i < Kab; ++i)
    {
        for (size_t j = 0; j < Kcd; ++j)
        {
            const double twoNuSq = 1 / (shellPairs.invTwoZeta(primOffsetAB + i) + shellPairs.invTwoZeta(primOffsetCD + j));
            const double T = 0.5 * twoNuSq * Rsq;
            const double U = shellPairs.Up(primOffsetAB) * shellPairs.Up(primOffsetCD);

            double* const colPtr = baseIntegrals.col(i * Kcd + j).data();

            Boys::calculateBoys(maxL, T, std::span(colPtr, maxL + 1));

            double val = U * std::sqrt(twoNuSq);
            for (size_t m = 0; m <= maxL; ++m)
            {
                colPtr[m] *= val;
                val *= twoNuSq;
            }
        }
    }

    return baseIntegrals;
}
