#pragma once

#include "Basis.hpp"
#include "Utils.hpp"

#include <cstddef>

/**
 * Struct to hold (significant) shell-pair information and precomputed values for PRISM algorithm.
 */
struct ShellPairs
{
    // Metadata
    size_t nShellPairs;
    size_t nPrimitivePairs;

    /*** Per shell-pair data ***/
    Eigen::Matrix<int, 2, Eigen::Dynamic> indices; // Indices of shell A and shell B for each shell pair.

    Eigen::Matrix<int, 2, Eigen::Dynamic> l; // Angular momentum. Rows: shell A or shell B, cols: index of shell pair.

    Eigen::Matrix<double, 3, Eigen::Dynamic> A;  // Center of shell A
    Eigen::Matrix<double, 3, Eigen::Dynamic> B;  // Center of shell B
    Eigen::Matrix<double, 3, Eigen::Dynamic> AB; // A - B

    Eigen::VectorXi K; // Degree of contraction (k_a * k_b) for each shell pair (number of primitive pairs).
    Eigen::VectorXi primPairOffsets; // Offset into the per primitive-pair data arrays for each shell pair.
    Eigen::VectorXi nAOsA;           // Number of AOs in shell A for each shell pair.
    Eigen::VectorXi nAOsB;           // Number of AOs in shell B for each shell pair.

    /*** Per primitive-pair data ***/

    Eigen::Matrix<double, 2, Eigen::Dynamic> alpha; // Exponents.Rows: primitive A or primitive B, cols: index of shell pair.
    Eigen::Matrix<double, 2, Eigen::Dynamic> twoAlpha; // 2 * alpha
    Eigen::Matrix<double, 3, Eigen::Dynamic> P;        // alpha * A + beta * B / (alpha + beta)
    Eigen::VectorXd invTwoZeta;                        // 1 / (2 * (alpha + beta))
    Eigen::VectorXd Up;                                // D_a * D_b * G_AB * (pi / zeta)^(3/2) * invTwoZeta^(l_a + l_b)

    /*** Normalization Factors ***/

    // Normalization factors depend on the specific angular momentum of each AO in the shell (total l is not
    // sufficient). Therefore, for each shell-pair we store a AO-pairs x primitive-pairs matrix of normalization
    // factors. The values are stored in a flattened array for better memory locality. Longest stride is the shell-pair
    // dimension, followed by the primitive-pair dimension, followed by the AO-pair dimension. For easier access, we
    // store a vector of Eigen::Map objects that map each shell-pair's normalization factors to a matrix view.

    Eigen::VectorXd normFactorValues;
    std::vector<Eigen::Map<Eigen::MatrixXd>> normFactors;
};

/**
 * @brief Class implementing the PRISM algorithm for efficient computation of electron repulsion integrals.
 *
 * This implementation is based on the work of Gill and Pople:
 *  - Gill, P. M., & Pople, J. A. (1991). The prism algorithm for two‐electron integrals. International journal of
 * quantum chemistry, 40(6), 753-772.
 *  - Gill, P. M. (1994). Molecular integrals over Gaussian basis functions. In Advances in quantum chemistry (Vol. 25,
 * pp. 141-205). Academic Press.
 */
class PRISMEngine
{
  public:
    /**
     * Constructs a PRISMEngine for a given basis set.
     * The constructor precomputes all necessary (significant) shell-pair data for the PRISM algorithm.
     *
     * @param basis The basis set to use for the PRISM algorithm.
     */
    PRISMEngine(const BasisSet& basis);

    private:

    /**
     * Represents a shell pair (A, B) before significance processing.
     */
    struct RawShellPair
    {
        RawShellPair() = default;
        RawShellPair(const Shell& shellA, const Shell& shellB, size_t indexA, size_t indexB, size_t nPrimPairs) :
            shellA(shellA),
            shellB(shellB),
            indexA(indexA),
            indexB(indexB),
            significant(true),
            significantPrimPairs(Eigen::Vector<bool, Eigen::Dynamic>::Constant(nPrimPairs, true))
        {
        }

        Shell shellA;
        Shell shellB;
        size_t indexA;
        size_t indexB;
        bool significant;
        Eigen::Vector<bool, Eigen::Dynamic> significantPrimPairs;
    };

    /**
     * Sort shell-pairs by increasing angular momentum of shellA, followed by shellB, and finally by increasing
     * (effective) degree of contraction (k_a * k_b - (insignficant primitive pairs)).
     *
     * @param shellPairs Vector of RawShellPair objects to sort.
     */
    void sortShellPairs(std::vector<RawShellPair>& shellPairs) const;

    /**
     * Precomputes shell-pair data required for the PRISM algorithm.
     *
     * @param basis The basis set to compute shell-pair data for.
     * @return A ShellPairs struct containing all precomputed data.
     */
    ShellPairs computeShellPairData(const BasisSet& basis) const;

    /**
     * Assigns which shell-pairs and primitive-pairs are significant.
     * A primitive-pair is considered neglibible if the distance between the primitives is large, relative to their
     * diffuseness, such that their overlap is negligible. Mathematically, a shell pair (A, B) is negligible if D_a *
     * D_b * C * G_AB < cutoff. A shell-pair is considered negligible if all its primitive-pairs are negligible.
     *
     * This implementation is based on: Neese, F. (2023). The SHARK integral generation and digestion system. Journal of
     * Computational Chemistry, 44(3), 381-396.
     *
     * @param shellPairs Vector of RawShellPairs to assign significance.
     * @param cutoff The cutoff value below which a shell pair is considered negligible.
     */
    void assignShellAndPrimPairSignificance(std::vector<RawShellPair>& shellPairs, double cutoff) const;

    // Precomputed shell-pair data for significant shell pairs only.
    const ShellPairs shellPairs;
};
