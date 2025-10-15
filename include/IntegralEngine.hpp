#pragma once
#include "Basis.hpp"
#include "Molecule.hpp"

#include <cassert>
#include <span>

class IntegralEngine
{
  public:
    IntegralEngine(const BasisSet& basis);

    void overlap(size_t shellA_idx, size_t shellB_idx, Eigen::MatrixXd& S) const;
    Eigen::MatrixXd overlapMatrix() const;

    void kinetic(size_t shellA_idx, size_t shellB_idx, Eigen::MatrixXd& T) const;
    Eigen::MatrixXd kineticMatrix() const;

    void nuclearAttraction(const std::vector<Atom>& geometry, size_t shellA_idx, size_t shellB_idx, Eigen::MatrixXd& V) const;
    Eigen::MatrixXd nuclearAttractionMatrix(const std::vector<Atom>& geometry) const;

    void electronRepulsion(
        size_t shellA_idx, size_t shellB_idx, size_t shellC_idx, size_t shellD_idx, ElectronRepulsionTensor& G
    ) const;

    void electronRepulsion(
        size_t shellA_idx, size_t shellB_idx, size_t shellC_idx, size_t shellD_idx, std::vector<double>& G
    ) const;

    ElectronRepulsionTensor electronRepulsionTensor(double threshold) const;

    Eigen::MatrixXd schwartzScreeningMatrix() const;

    const Eigen::VectorXi& getNAOsPerShell() const { return nAOs; }
    const Eigen::VectorXi& getAOOffsetsPerShell() const { return aoOffsets; }

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
    static PrimitivePairData computePrimitivePairData(double alphaA, const Vec3& centerA, double alphaB, const Vec3& centerB);


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

        size_t getIOffset(unsigned i) const { return i * (l2 + 1) * (i + l2 + 1) / 2; }
        size_t getJOffset(unsigned i, unsigned j) const { return (j * (i + 1)) + j * (j - 1) / 2; }

        double& operator()(unsigned i, unsigned j, unsigned t)
        {
            const size_t offset_i = getIOffset(i);
            const size_t offset_j = getJOffset(i, j);
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
    static void computeHermiteCoeffsTest(unsigned l1, unsigned l2, double invTwozeta, double PA, double PB, EBuffer& E_buffer);


    struct RBuffer
    {
        size_t max_l;
        size_t stride_u, stride_v, stride_n;

        std::vector<double> data;
        std::vector<bool> needed;

        RBuffer() = default;
        RBuffer(size_t max_l) :
            max_l(max_l),
            stride_u(max_l + 1),
            stride_v((max_l + 1) * (max_l + 1)),
            stride_n((max_l + 1) * (max_l + 1) * (max_l + 1)),
            needed(stride_n * (max_l + 1), true)
        {
            data.resize(stride_n * (max_l + 1));
        }

        double& operator()(size_t t, size_t u, size_t v, size_t n) noexcept
        {
            const size_t index = t + u * stride_u + v * stride_v + n * stride_n;
            assert(index < data.size());
            return data[index];
        }

        const double& operator()(size_t t, size_t u, size_t v, size_t n) const noexcept
        {
            const size_t index = t + u * stride_u + v * stride_v + n * stride_n;
            assert(index < data.size());
            return data[index];
        }
    };

    static void computeAuxiliaryIntegrals(double p, const Vec3& PC, std::span<double> F, RBuffer& R_buffer);

    // Per basis
    size_t nShells;
    size_t nAOsTotal;
    size_t nShellPairs;
    size_t nPrimitivePairsTot;

    // Per shell
    Eigen::VectorXi lShell;
    Eigen::VectorXi primitiveOffsets;
    Eigen::VectorXi aoOffsets;
    Eigen::VectorXi normFactorOffsets;
    Eigen::VectorXi nPrimitives;
    Eigen::VectorXi nAOs;
    Eigen::Matrix<double, 3, Eigen::Dynamic> centers;

    // per primitive
    Eigen::VectorXd alpha;
    Eigen::VectorXd coeff;

    // per AO
    Eigen::Matrix<int, 3, Eigen::Dynamic> angularMomentum;

    // per primitive per AO
    Eigen::VectorXd normalizationFactors;
};
