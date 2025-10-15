#include "SCF.hpp"

#include "DIIS.hpp"
#include "IntegralEngine.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fmt/core.h>
#include <ios>
#include <iostream>

SCF::SCF(const Molecule& molecule, const BasisSet& basis, const SCFOptions& options, std::shared_ptr<Output> output) :
    output(output),
    molecule(std::make_unique<Molecule>(molecule)),
    basis(basis),
    integralEngine(std::make_unique<IntegralEngine>(basis)),
    options(options)
{
}


void SCF::run()
{
    // Initialize the SCF calculation. (This should be called only once.)
    this->initialize(this->options.schwartzThreshold, this->options.direct);

    // Get the initial guess density.
    this->computeInitialGuessDensity();

    if (this->options.useDIIS) // Initialize the DIIS handler with the specified maximum size.
        diis_handler = std::make_unique<DIIS>(this->options.DIISmaxSize);

    // Print header for the iteration table.
    size_t width = 45;
    if (this->options.useDIIS)
        width += 15;
    if (this->options.damp > 0)
        width += 15;
    if (this->options.levelShift > 0)
        width += 15;
    std::stringstream ss;
    ss << fmt::format("\n{:-<{}}\n", "", width);
    ss << fmt::format("{:^15}{:^15}{:^15}", "Iteration", "ΔE", "ΔD");
    if (this->options.useDIIS)
        ss << fmt::format("{:^15}", "DIIS Error");
    ss << fmt::format("\n{:-<{}}\n", "", width);
    this->output->write(ss.str());

    // SCF cycles
    for (this->iteration = 0; this->iteration < this->options.maxIter; ++this->iteration)
    {
        // Store the energy and density from the previous iteration to check for convergence.
        double lastElectronicEnergy = this->electronicEnergy;

        // Build the new Fock matrix from the current density.
        if (this->options.direct)
            this->buildFockMatrixDirect();
        else
            this->buildFockMatrix();

        // Store the error vector and fock matrix for DIIS if applicable.
        if (this->options.useDIIS)
        {
            // Calculate the commutator of the fock and density matrices as the error vector.
            Eigen::MatrixXd e_alpha = this->F_alpha * this->D_alpha * this->S - this->S * this->D_alpha * this->F_alpha;
            Eigen::MatrixXd e_beta  = this->F_beta * this->D_beta * this->S - this->S * this->D_beta * this->F_beta;
            Eigen::MatrixXd errorVector = e_alpha + e_beta;
            diis_handler->update(this->F_alpha, this->F_beta, errorVector);
        }

        // Diagonalize the Fock matrix and compute the new density and coefficients.
        if (this->options.useDIIS)
        {
            // Extrapolate the orthogonal Fock matrix using DIIS.
            std::tie(this->F_alpha, this->F_beta) = diis_handler->extrapolate(this->F_alpha, this->F_beta);
        }
        this->diagonalizeAndUpdate();

        // Calculate the change in energy and density.
        this->deltaE = std::abs(this->electronicEnergy - lastElectronicEnergy);
        this->deltaD = (this->D_tot - D_tot_prev).norm();
        if (this->options.useDIIS)
            this->DIISError = diis_handler->getErrorNorm();

        // Print the current iteration results.
        this->printIteration();

        // Check for convergence.
        if (this->deltaE < this->options.energyTol && this->deltaD < this->options.densityTol
            && (!this->options.useDIIS || diis_handler->getErrorNorm() < this->options.DIISErrorTol))
        {
            // Remove the level shift from the virtual eigenvalues
            if (this->useLevelShiftingAlpha)
            {
                for (size_t i = this->occupiedCountAlpha; i < this->basisCount; ++i)
                {
                    this->eigenvalues_alpha(i) -= this->options.levelShift;
                }
            }
            if (this->useLevelShiftingBeta)
            {
                for (size_t i = this->occupiedCountBeta; i < this->basisCount; ++i)
                {
                    this->eigenvalues_beta(i) -= this->options.levelShift;
                }
            }

            // basis If converged, print the final results and exit.
            this->output->writeSeperator('-', width);
            this->printFinalResults(true);
            return;
        }
    }

    // If the loop finishes, the calculation did not converge.

    // Remove the level shift from the virtual eigenvalues
    if (this->useLevelShiftingAlpha)
    {
        for (size_t i = this->occupiedCountAlpha; i < this->basisCount; ++i)
        {
            this->eigenvalues_alpha(i) -= this->options.levelShift;
        }
    }
    if (this->useLevelShiftingBeta)
    {
        for (size_t i = this->occupiedCountBeta; i < this->basisCount; ++i)
        {
            this->eigenvalues_beta(i) -= this->options.levelShift;
        }
    }

    this->output->writeSeperator('-', width);
    this->printFinalResults(false);
}

void SCF::initialize(double schwartzThreshold, bool direct)
{
    // Set constant parameters.
    this->basisCount         = basis.nAOs;
    this->occupiedCountAlpha = (molecule->getElectronCount() + molecule->getMultiplicity() - 1) / 2;
    this->occupiedCountBeta  = (molecule->getElectronCount() - molecule->getMultiplicity() + 1) / 2;

    this->output->write(printJobSpec());

    // Calculate one-electron integrals and nuclear repulsion.
    this->nuclearEnergy = molecule->nuclearRepulsion();
    this->T             = integralEngine->kineticMatrix();
    this->V             = integralEngine->nuclearAttractionMatrix(molecule->getGeometry());
    this->h             = this->T + this->V;

    // Calculate the overlap and orthogonalization matrices.
    this->S = integralEngine->overlapMatrix();
    this->X = inverseSqrtMatrix(this->S);

    // Calculate two-electron integrals (only if not using direct method).
    if (!direct)
        this->Vee = integralEngine->electronRepulsionTensor(schwartzThreshold);
    else // If using direct method, only pre-calculate the Schwartz screening matrix.
        this->Q = integralEngine->schwartzScreeningMatrix();
}


void SCF::computeInitialGuessDensity()
{
    // The initial Fock matrix is just the core Hamiltonian.
    this->F_alpha      = this->h;
    this->F_alpha_prev = this->F_alpha;

    // Diagonalize the initial Fock matrix to get guess orbitals.
    Eigen::MatrixXd F_alpha_prime = this->X.transpose() * this->F_alpha * this->X;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> alpha_solver(F_alpha_prime);
    this->C_alpha               = this->X * alpha_solver.eigenvectors();
    this->eigenvalues_alpha     = alpha_solver.eigenvalues();
    Eigen::MatrixXd C_alpha_occ = this->C_alpha.leftCols(this->occupiedCountAlpha);

    // Form the density matrix from the occupied orbitals.
    this->D_alpha_prev = Eigen::MatrixXd::Zero(this->basisCount, this->basisCount);
    this->D_alpha      = C_alpha_occ * C_alpha_occ.transpose();

    this->D_beta_prev = Eigen::MatrixXd::Zero(this->basisCount, this->basisCount);
    if (this->options.unrestricted)
    {
        // HOMO-LUMO mixing for breaking symmetry in UHF.
        if (this->options.guessMix > 0)
        {
            double k = static_cast<double>(this->options.guessMix) / 10.0;

            Eigen::VectorXd homoAlpha = C_alpha.col(this->occupiedCountAlpha - 1);
            Eigen::VectorXd lumoAlpha = C_alpha.col(this->occupiedCountAlpha);

            Eigen::VectorXd mixedHomo = (1 / std::sqrt(1 + k * k)) * (homoAlpha + k * lumoAlpha);
            Eigen::VectorXd mixedLumo = (1 / std::sqrt(1 + k * k)) * (lumoAlpha - k * homoAlpha);

            C_alpha.col(this->occupiedCountAlpha - 1) = mixedHomo;
            C_alpha.col(this->occupiedCountAlpha)     = mixedLumo;
            Eigen::MatrixXd C_alpha_occ               = this->C_alpha.leftCols(this->occupiedCountAlpha);
            this->D_alpha                             = C_alpha_occ * C_alpha_occ.transpose();
        }

        this->F_beta                 = this->h;
        this->F_beta_prev            = this->F_beta;
        Eigen::MatrixXd F_beta_prime = this->X.transpose() * this->F_beta * this->X;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> beta_solver(F_beta_prime);
        this->C_beta               = this->X * beta_solver.eigenvectors();
        this->eigenvalues_beta     = beta_solver.eigenvalues();
        Eigen::MatrixXd C_beta_occ = this->C_beta.leftCols(this->occupiedCountBeta);
        this->D_beta               = C_beta_occ * C_beta_occ.transpose();
    }
    else
    {
        this->F_beta           = this->F_alpha;
        this->F_beta_prev      = this->F_alpha_prev;
        this->C_beta           = this->C_alpha;
        this->eigenvalues_beta = this->eigenvalues_alpha;
        this->D_beta           = this->D_alpha;
    }

    this->D_tot_prev = Eigen::MatrixXd::Zero(this->basisCount, this->basisCount);
    this->D_tot      = this->D_alpha + this->D_beta;

    // Calculate the initial electronic energy.
    this->electronicEnergy = 0.5 * (this->D_alpha * (this->h + this->F_alpha)).trace()
                           + 0.5 * (this->D_beta * (this->h + this->F_beta)).trace();
}

void SCF::buildFockMatrix()
{
    const size_t N_ao       = this->basisCount;
    Eigen::MatrixXd J       = Eigen::MatrixXd::Zero(N_ao, N_ao);
    Eigen::MatrixXd K_alpha = Eigen::MatrixXd::Zero(N_ao, N_ao);
    Eigen::MatrixXd K_beta;
    if (this->options.unrestricted)
        K_beta = Eigen::MatrixXd::Zero(N_ao, N_ao);

#pragma omp parallel
    {
        // Thread-local matrices to void data races
        const double* const Dtot_ptr   = this->D_tot.data();
        const double* const Dalpha_ptr = this->D_alpha.data();
        const double* const Dbeta_ptr  = this->D_beta.data();
        Eigen::MatrixXd J_p            = Eigen::MatrixXd::Zero(N_ao, N_ao);
        double* const J_ptr            = J_p.data();
        Eigen::MatrixXd K_alpha_p      = Eigen::MatrixXd::Zero(N_ao, N_ao);
        double* const K_alpha_ptr      = K_alpha_p.data();
        Eigen::MatrixXd K_beta_p;
        double* K_beta_ptr = nullptr;
        if (this->options.unrestricted)
        {
            K_beta_p   = Eigen::MatrixXd::Zero(N_ao, N_ao);
            K_beta_ptr = K_beta_p.data();
        }

#pragma omp for schedule(guided, 1) nowait
        for (size_t i = 0; i < N_ao; ++i) // (ii|ii) contributions
        {
            const size_t big_I = i * (i + 1) / 2 + i;
            const size_t index = big_I * (big_I + 1) / 2 + big_I;
            const double eri   = this->Vee[index];

            if (std::abs(eri) < this->options.schwartzThreshold)
                continue;

            const size_t ii_idx = i + i * N_ao;

            J_ptr[ii_idx] += Dtot_ptr[ii_idx] * eri;
            K_alpha_ptr[ii_idx] += Dalpha_ptr[ii_idx] * eri;

            if (this->options.unrestricted)
            {
                K_beta_ptr[ii_idx] += Dbeta_ptr[ii_idx] * eri;
            }
        }

#pragma omp for schedule(guided, 1) nowait
        for (size_t i = 0; i < N_ao; ++i) // (ii|kk) contributions
        {
            const size_t iN_ao      = i * N_ao;
            const size_t ii_idx     = i + iN_ao;
            const size_t big_I      = i * (i + 1) / 2 + i;
            const size_t index_base = big_I * (big_I + 1) / 2;

            for (size_t k = 0; k < i; ++k)
            {
                const size_t big_K = k * (k + 1) / 2 + k;
                const size_t index = index_base + big_K;
                const double eri   = this->Vee[index];

                if (std::abs(eri) < this->options.schwartzThreshold)
                    continue;

                const size_t kN_ao  = k * N_ao;
                const size_t ki_idx = k + iN_ao;
                const size_t kk_idx = k + kN_ao;

                J_ptr[ii_idx] += Dtot_ptr[kk_idx] * eri;
                J_ptr[kk_idx] += Dtot_ptr[ii_idx] * eri;
                K_alpha_ptr[ki_idx] += Dalpha_ptr[ki_idx] * eri;

                if (this->options.unrestricted)
                {
                    K_beta_ptr[ki_idx] += Dbeta_ptr[ki_idx] * eri;
                }
            }
        }

#pragma omp for schedule(guided, 1) nowait
        for (size_t i = 0; i < N_ao; ++i) // (ij|ij) contributions
        {
            const size_t iN_ao      = i * N_ao;
            const size_t ii_idx     = i + iN_ao;
            const size_t big_I_base = i * (i + 1) / 2;

            for (size_t j = 0; j < i; ++j)
            {
                const size_t big_I = big_I_base + j;
                const size_t index = big_I * (big_I + 1) / 2 + big_I;
                const double eri   = this->Vee[index];

                if (std::abs(eri) < this->options.schwartzThreshold)
                    continue;

                const size_t jN_ao  = j * N_ao;
                const size_t ji_idx = j + iN_ao;
                const size_t jj_idx = j + jN_ao;

                J_ptr[ji_idx] += 2 * Dtot_ptr[ji_idx] * eri;
                K_alpha_ptr[ji_idx] += Dalpha_ptr[ji_idx] * eri;
                K_alpha_ptr[ii_idx] += Dalpha_ptr[jj_idx] * eri;
                K_alpha_ptr[jj_idx] += Dalpha_ptr[ii_idx] * eri;

                if (this->options.unrestricted)
                {
                    K_beta_ptr[ji_idx] += Dbeta_ptr[ji_idx] * eri;
                    K_beta_ptr[ii_idx] += Dbeta_ptr[jj_idx] * eri;
                    K_beta_ptr[jj_idx] += Dbeta_ptr[ii_idx] * eri;
                }
            }
        }

#pragma omp for schedule(guided, 1) nowait
        for (size_t i = 0; i < N_ao; ++i) // (ii|kl) contributions
        {
            const size_t iN_ao        = i * N_ao;
            const size_t ii_idx       = i + iN_ao;
            const size_t big_I        = i * (i + 1) / 2 + i;
            const size_t index_base_i = big_I * (big_I + 1) / 2;

            for (size_t k = 0; k <= i; ++k)
            {
                const size_t kN_ao      = k * N_ao;
                const size_t ki_idx     = k + iN_ao;
                const size_t big_K_base = k * (k + 1) / 2;
                const size_t index_base = index_base_i + big_K_base;

                for (size_t l = 0; l < k; ++l)
                {
                    const size_t index = index_base + l;
                    const double eri   = this->Vee[index];

                    if (std::abs(eri) < this->options.schwartzThreshold)
                        continue;

                    const size_t li_idx = l + iN_ao;
                    const size_t lk_idx = l + kN_ao;

                    const unsigned symFactor = i == k ? 2 : 1;

                    J_ptr[ii_idx] += 2 * Dtot_ptr[lk_idx] * eri;
                    J_ptr[lk_idx] += Dtot_ptr[ii_idx] * eri;
                    K_alpha_ptr[li_idx] += Dalpha_ptr[ki_idx] * eri;
                    K_alpha_ptr[ki_idx] += symFactor * Dalpha_ptr[li_idx] * eri;

                    if (this->options.unrestricted)
                    {
                        K_beta_ptr[li_idx] += Dbeta_ptr[ki_idx] * eri;
                        K_beta_ptr[ki_idx] += symFactor * Dbeta_ptr[li_idx] * eri;
                    }
                }
            }
        }

#pragma omp for schedule(guided, 1) nowait
        for (size_t i = 0; i < N_ao; ++i) // (ij|kk) contributions
        {
            const size_t iN_ao      = i * N_ao;
            const size_t big_I_base = i * (i + 1) / 2;
            for (size_t j = 0; j < i; ++j)
            {
                const size_t jN_ao      = j * N_ao;
                const size_t ji_idx     = j + iN_ao;
                const size_t big_I      = big_I_base + j;
                const size_t index_base = big_I * (big_I + 1) / 2;

                for (size_t k = 0; k < i; ++k)
                {
                    const size_t big_K = k * (k + 1) / 2 + k;
                    const size_t index = index_base + big_K;
                    const double eri   = this->Vee[index];

                    if (std::abs(eri) < this->options.schwartzThreshold)
                        continue;

                    const size_t kN_ao  = k * N_ao;
                    const size_t kk_idx = k + kN_ao;
                    const size_t ki_idx = k + iN_ao;
                    const size_t kj_idx = k + jN_ao;

                    J_ptr[ji_idx] += Dtot_ptr[kk_idx] * eri;
                    J_ptr[kk_idx] += 2 * Dtot_ptr[ji_idx] * eri;
                    K_alpha_ptr[ki_idx] += Dalpha_ptr[kj_idx] * eri;

                    size_t idx               = j < k ? j + kN_ao : kj_idx;
                    const unsigned symFactor = j == k ? 2 : 1;
                    K_alpha_ptr[idx] += symFactor * Dalpha_ptr[ki_idx] * eri;

                    if (this->options.unrestricted)
                    {
                        K_beta_ptr[ki_idx] += Dbeta_ptr[kj_idx] * eri;
                        K_beta_ptr[idx] += symFactor * Dbeta_ptr[ki_idx] * eri;
                    }
                }
            }
        }

#pragma omp for schedule(guided, 1) nowait
        for (size_t i = 0; i < N_ao; ++i) // (ij|kl) contributions
        {
            const size_t iN_ao      = i * N_ao;
            const size_t big_I_base = i * (i + 1) / 2;
            for (size_t j = 0; j < i; ++j)
            {
                const size_t jN_ao        = j * N_ao;
                const size_t ji_idx       = j + iN_ao;
                const size_t big_I        = big_I_base + j;
                const size_t index_base_I = big_I * (big_I + 1) / 2;

                for (size_t k = 0; k <= i; ++k)
                {
                    const size_t kN_ao        = k * N_ao;
                    const size_t ki_idx       = k + iN_ao;
                    const size_t kj_idx       = k + jN_ao;
                    const size_t big_K_base   = k * (k + 1) / 2;
                    const size_t index_base_K = index_base_I + big_K_base;

                    const size_t l_max = (i == k) ? j : k;
                    for (size_t l = 0; l < l_max; ++l)
                    {
                        if (i == k && j == l) // exclude (ij|ij)
                            continue;

                        const size_t index = index_base_K + l;
                        const double eri   = this->Vee[index];

                        if (std::abs(eri) < this->options.schwartzThreshold)
                            continue;

                        const size_t li_idx = l + iN_ao;
                        const size_t lj_idx = l + jN_ao;
                        const size_t lk_idx = l + kN_ao;

                        const unsigned symFactor_ik = i == k ? 2 : 1;
                        const unsigned symFactor_jl = j == l ? 2 : 1;
                        const unsigned symFactor_jk = j == k ? 2 : 1;

                        J_ptr[ji_idx] += 2 * Dtot_ptr[lk_idx] * eri;
                        J_ptr[lk_idx] += 2 * Dtot_ptr[ji_idx] * eri;
                        K_alpha_ptr[li_idx] += Dalpha_ptr[kj_idx] * eri;

                        K_alpha_ptr[ki_idx] += symFactor_ik * Dalpha_ptr[lj_idx] * eri;

                        size_t idx1 = j < l ? j + l * N_ao : lj_idx;
                        K_alpha_ptr[idx1] += symFactor_jl * Dalpha_ptr[ki_idx] * eri;

                        size_t idx2 = j < k ? j + kN_ao : kj_idx;
                        K_alpha_ptr[idx2] += symFactor_jk * Dalpha_ptr[li_idx] * eri;

                        if (this->options.unrestricted)
                        {
                            K_beta_ptr[li_idx] += Dbeta_ptr[kj_idx] * eri;
                            K_beta_ptr[idx1] += symFactor_jl * Dbeta_ptr[ki_idx] * eri;
                            K_beta_ptr[ki_idx] += symFactor_ik * Dbeta_ptr[lj_idx] * eri;
                            K_beta_ptr[idx2] += symFactor_jk * Dbeta_ptr[li_idx] * eri;
                        }
                    }
                }
            }
        }

#pragma omp critical(build_fock_conventional)
        {
            // Accumulate the thread-local matrices into the global ones.
            J += J_p;
            K_alpha += K_alpha_p;
            if (this->options.unrestricted)
                K_beta += K_beta_p;
        }
    }

    // Symmetrize J and K
    J       = J.selfadjointView<Eigen::Upper>();
    K_alpha = K_alpha.selfadjointView<Eigen::Upper>();

    // Build Fock matrix
    this->F_alpha = this->h + J - K_alpha;
    if (this->options.unrestricted)
    {
        K_beta       = K_beta.selfadjointView<Eigen::Upper>();
        this->F_beta = this->h + J - K_beta;
    }
    else
    {
        this->F_beta = this->F_alpha;
    }
}

void SCF::buildFockMatrixDirect()
{
    const size_t N_ao     = this->basisCount;
    const auto& nAOs      = this->integralEngine->getNAOsPerShell();
    const auto& aoOffsets = this->integralEngine->getAOOffsetsPerShell();
    const size_t nShells  = this->basis.nShells;

    Eigen::MatrixXd J       = Eigen::MatrixXd::Zero(N_ao, N_ao);
    Eigen::MatrixXd K_alpha = Eigen::MatrixXd::Zero(N_ao, N_ao);
    Eigen::MatrixXd K_beta;
    if (this->options.unrestricted)
        K_beta = Eigen::MatrixXd::Zero(N_ao, N_ao);

    const Eigen::MatrixXd deltaD_tot    = this->D_tot - this->D_tot_prev;
    const double* const deltaDtot_ptr   = deltaD_tot.data();
    const Eigen::MatrixXd deltaD_alpha  = this->D_alpha - this->D_alpha_prev;
    const double* const deltaDalpha_ptr = deltaD_alpha.data();
    const Eigen::MatrixXd deltaD_beta   = this->D_beta - this->D_beta_prev;
    const double* const deltaDbeta_ptr  = this->options.unrestricted ? deltaD_beta.data() : nullptr;

    // Pre-calculate the maximum density elements for each shell pair for density screening.
    Eigen::MatrixXd deltaDtotMax   = Eigen::MatrixXd::Zero(nShells, nShells);
    Eigen::MatrixXd deltaDalphaMax = Eigen::MatrixXd::Zero(nShells, nShells);
    Eigen::MatrixXd deltaDbetaMax  = Eigen::MatrixXd::Zero(nShells, nShells);

#pragma omp parallel for collapse(2)
    for (size_t s1 = 0; s1 < nShells; ++s1)
    {
        for (size_t s2 = 0; s2 <= s1; ++s2)
        {
            double maxValTot = 0.0, maxValAlpha = 0.0, maxValBeta = 0.0;
            for (int a = aoOffsets[s1]; a < aoOffsets[s1] + nAOs[s1]; ++a)
            {
                for (int b = aoOffsets[s2]; a < aoOffsets[s2] + nAOs[s2]; ++b)
                {
                    maxValTot   = std::max(maxValTot, std::abs(deltaD_tot(a, b)));
                    maxValAlpha = std::max(maxValAlpha, std::abs(deltaD_alpha(a, b)));
                    if (this->options.unrestricted)
                        maxValBeta = std::max(maxValBeta, std::abs(deltaD_beta(a, b)));
                }
            }
            deltaDtotMax(s2, s1)   = maxValTot;
            deltaDalphaMax(s2, s1) = maxValAlpha;
            if (this->options.unrestricted)
                deltaDbetaMax(s2, s1) = maxValBeta;
        }
    }
    deltaDtotMax   = deltaDtotMax.selfadjointView<Eigen::Upper>();
    deltaDalphaMax = deltaDalphaMax.selfadjointView<Eigen::Upper>();
    if (this->options.unrestricted)
        deltaDbetaMax = deltaDbetaMax.selfadjointView<Eigen::Upper>();

#pragma omp parallel
    {
        // Thread-local matrices to void data races
        Eigen::MatrixXd J_p       = Eigen::MatrixXd::Zero(N_ao, N_ao);
        double* const J_ptr       = J_p.data();
        Eigen::MatrixXd K_alpha_p = Eigen::MatrixXd::Zero(N_ao, N_ao);
        double* const K_alpha_ptr = K_alpha_p.data();
        Eigen::MatrixXd K_beta_p;
        double* K_beta_ptr = nullptr;
        if (this->options.unrestricted)
        {
            K_beta_p   = Eigen::MatrixXd::Zero(N_ao, N_ao);
            K_beta_ptr = K_beta_p.data();
        }

#pragma omp for collapse(4) schedule(dynamic, 32) nowait
        for (size_t a = 0; a < nShells; ++a)
        {
            for (size_t b = 0; b < nShells; ++b)
            {
                for (size_t c = 0; c < nShells; ++c)
                {
                    for (size_t d = 0; d < nShells; ++d)
                    {
                        // Only calculate unique integrals
                        if (b > a || d > c || (a * (a + 1) / 2 + b) < (c * (c + 1) / 2 + d))
                            continue;

                        // Schwartz screening
                        if (Q(a, b) * Q(c, d) < this->options.schwartzThreshold)
                            continue;

                        // Density screening
                        double maxDensity = std::max(
                            {2 * deltaDtotMax(a, b),
                             2 * deltaDtotMax(c, d),
                             deltaDalphaMax(a, c),
                             deltaDalphaMax(a, d),
                             deltaDalphaMax(b, c),
                             deltaDalphaMax(b, d)}
                        );
                        if (this->options.unrestricted)
                            maxDensity = std::max(
                                {maxDensity, deltaDbetaMax(a, c), deltaDbetaMax(a, d), deltaDbetaMax(b, c), deltaDbetaMax(b, d)}
                            );
                        if (maxDensity < this->options.densityThreshold)
                            continue;

                        const int nAOa = nAOs[a];
                        const int nAOb = nAOs[b];
                        const int nAOc = nAOs[c];
                        const int nAOd = nAOs[d];

                        std::vector<double> ERIs(nAOa * nAOb * nAOc * nAOd, 0.0);
                        integralEngine->electronRepulsion(a, b, c, d, ERIs);

                        const size_t naoCD  = nAOc * nAOd;
                        const size_t naoBCD = nAOb * naoCD;
                        const size_t naoD   = nAOd;

                        for (int A = 0; A < nAOa; ++A)
                        {
                            const size_t i           = A + aoOffsets[a];
                            const size_t iN_ao       = i * N_ao;
                            const size_t index_baseA = A * naoBCD;

                            for (int B = 0; B < nAOb; ++B)
                            {
                                const size_t j = B + aoOffsets[b];
                                if (j > i)
                                    continue;

                                const size_t jN_ao       = j * N_ao;
                                const size_t ji_idx      = j + iN_ao;
                                const size_t index_baseB = index_baseA + B * naoCD;

                                for (int C = 0; C < nAOc; ++C)
                                {
                                    const size_t k           = C + aoOffsets[c];
                                    const size_t kN_ao       = k * N_ao;
                                    const size_t ki_idx      = k + iN_ao;
                                    const size_t kj_idx      = k + jN_ao;
                                    const size_t index_baseC = index_baseB + C * naoD;

                                    const size_t l_max = i == k ? j : k;
                                    for (int D = 0; D < nAOd; ++D)
                                    {
                                        const size_t l = D + aoOffsets[d];
                                        if (l > l_max)
                                            continue;

                                        const size_t index = index_baseC + D;
                                        const double eri   = ERIs[index];
                                        if (std::abs(eri) < this->options.schwartzThreshold)
                                            continue;

                                        const size_t li_idx = l + iN_ao;
                                        const size_t lj_idx = l + jN_ao;
                                        const size_t lk_idx = l + kN_ao;

                                        const unsigned symFactor_kl = k == l ? 1 : 2;
                                        const unsigned symFactor_ij = i == j ? 1 : 2;

                                        J_ptr[ji_idx] += symFactor_kl * deltaDtot_ptr[lk_idx] * eri;
                                        J_ptr[lk_idx] += (l != j || k != i) * symFactor_ij * deltaDtot_ptr[ji_idx] * eri;

                                        const unsigned symFactor_jk = j == k ? 2 : 1;
                                        const unsigned symFactor_jl = j == l && i != k ? 2 : 1;
                                        const unsigned symFactor_ik = i == k && j != l ? 2 : 1;

                                        const size_t idx1 = i < k ? i + k * N_ao : ki_idx;
                                        const size_t idx2 = j < l ? j + l * N_ao : lj_idx;
                                        const size_t idx3 = j < k ? j + k * N_ao : kj_idx;

                                        K_alpha_ptr[li_idx] += deltaDalpha_ptr[kj_idx] * eri;
                                        K_alpha_ptr[idx1] += (idx1 != li_idx) * symFactor_ik * deltaDalpha_ptr[lj_idx] * eri;
                                        K_alpha_ptr[idx2] += (idx2 != li_idx && idx2 != idx1) * symFactor_jl
                                                           * deltaDalpha_ptr[ki_idx] * eri;
                                        K_alpha_ptr[idx3] += (idx3 != li_idx && idx3 != idx1 && idx3 != idx2)
                                                           * symFactor_jk * deltaDalpha_ptr[li_idx] * eri;

                                        if (this->options.unrestricted)
                                        {
                                            K_beta_ptr[li_idx] += deltaDbeta_ptr[kj_idx] * eri;
                                            K_beta_ptr[idx1] += (idx1 != li_idx) * symFactor_ik * deltaDbeta_ptr[lj_idx]
                                                              * eri;
                                            K_beta_ptr[idx2] += (idx2 != li_idx && idx2 != idx1) * symFactor_jl
                                                              * deltaDbeta_ptr[ki_idx] * eri;
                                            K_beta_ptr[idx3] += (idx3 != li_idx && idx3 != idx1 && idx3 != idx2)
                                                              * symFactor_jk * deltaDbeta_ptr[li_idx] * eri;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

#pragma omp critical(build_fock_direct)
        {
            // Accumulate the thread-local matrices into the global ones.
            J += J_p;
            K_alpha += K_alpha_p;
            if (this->options.unrestricted)
                K_beta += K_beta_p;
        }
    }

    J       = J.selfadjointView<Eigen::Upper>();
    K_alpha = K_alpha.selfadjointView<Eigen::Upper>();

    // F_m = F_{m-1} + I(D_m - D_{m-1})
    this->F_alpha      = this->F_alpha_prev + J - K_alpha;
    this->F_alpha_prev = this->F_alpha;
    if (this->options.unrestricted)
    {
        K_beta            = K_beta.selfadjointView<Eigen::Upper>();
        this->F_beta      = this->F_beta_prev + J - K_beta;
        this->F_beta_prev = this->F_beta;
    }
    else
    {
        this->F_beta      = this->F_alpha;
        this->F_beta_prev = this->F_alpha_prev;
    }
}


void SCF::diagonalizeAndUpdate()
{
    // Check whether to use damping.
    if (this->options.damp > 0 && this->options.maxDampIter != 0 && this->options.maxDampIter > this->iteration
        && (this->deltaE > this->options.stopDampThresh || this->iteration == 0))
    {
        this->dampCoeff = static_cast<double>(this->options.damp) / 100;
    }
    else
    {
        this->dampCoeff = 0;
    }

    // Check whether to use level shifting.
    if (this->options.levelShift > 0 && this->iteration < this->options.maxLshiftIter
        && (this->deltaE > this->options.stopLshiftThresh || this->iteration == 0))
    {
        double alphaHomoLumoGap = this->eigenvalues_alpha(this->occupiedCountAlpha)
                                - this->eigenvalues_alpha(this->occupiedCountAlpha - 1);
        double betaHomoLumoGap = this->eigenvalues_beta(this->occupiedCountBeta)
                               - this->eigenvalues_beta(this->occupiedCountBeta - 1);
        this->useLevelShiftingAlpha = alphaHomoLumoGap < this->options.lshiftGapTol;
        this->useLevelShiftingBeta  = betaHomoLumoGap < this->options.lshiftGapTol;
    }
    else
    {
        this->useLevelShiftingAlpha = false;
        this->useLevelShiftingBeta  = false;
    }

    // Transform the Fock matrix to the orthogonal basis.
    Eigen::MatrixXd F_alpha_prime = this->X.transpose() * this->F_alpha * this->X;

    // Apply level shifting
    if (this->useLevelShiftingAlpha)
    {
        Eigen::MatrixXd shift = Eigen::MatrixXd::Zero(this->basisCount, this->basisCount);
        for (size_t i = this->occupiedCountAlpha; i < this->basisCount; ++i) { shift(i, i) = this->options.levelShift; }
        F_alpha_prime += this->X * this->S * this->C_alpha * shift * (this->X * this->S * this->C_alpha).transpose();
    }

    // Diagonalize the Fock matrix.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> alpha_solver(F_alpha_prime);
    this->C_alpha           = this->X * alpha_solver.eigenvectors();
    this->eigenvalues_alpha = alpha_solver.eigenvalues();

    // Form the new density matrix from the occupied orbitals.
    Eigen::MatrixXd C_alpha_occ = this->C_alpha.leftCols(this->occupiedCountAlpha);

    // Form the new density matrix with damping (the damping coefficient can be zero).
    this->D_alpha_prev = this->D_alpha;
    this->D_alpha = (1 - this->dampCoeff) * C_alpha_occ * C_alpha_occ.transpose() + this->dampCoeff * this->D_alpha;

    // Deal with beta spin if unrestricted.
    this->D_beta_prev = this->D_beta;
    if (this->options.unrestricted)
    {
        Eigen::MatrixXd F_beta_prime = this->X.transpose() * this->F_beta * this->X;

        if (this->useLevelShiftingBeta)
        {
            Eigen::MatrixXd shift = Eigen::MatrixXd::Zero(this->basisCount, this->basisCount);
            for (size_t i = this->occupiedCountBeta; i < this->basisCount; ++i)
            {
                shift(i, i) = this->options.levelShift;
            }
            F_beta_prime += this->X * this->S * this->C_beta * shift * (this->X * this->S * this->C_beta).transpose();
        }

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> beta_solver(F_beta_prime);
        this->C_beta               = this->X * beta_solver.eigenvectors();
        this->eigenvalues_beta     = beta_solver.eigenvalues();
        Eigen::MatrixXd C_beta_occ = this->C_beta.leftCols(this->occupiedCountBeta);
        this->D_beta = (1 - this->dampCoeff) * C_beta_occ * C_beta_occ.transpose() + this->dampCoeff * this->D_beta;
    }
    else
    {
        this->C_beta           = this->C_alpha;
        this->eigenvalues_beta = this->eigenvalues_alpha;
        this->D_beta           = this->D_alpha;
    }

    this->D_tot_prev = this->D_tot;
    this->D_tot      = this->D_alpha + this->D_beta;

    // Calculate the new electronic energy.
    this->electronicEnergy = 0.5 * (this->D_alpha * (this->h + this->F_alpha)).trace()
                           + 0.5 * (this->D_beta * (this->h + this->F_beta)).trace();
}

double SCF::computeSpinSquared() const
{
    // <S^2> = S(S + 1) + N_beta - sum_ij |<psi_i_alpha | psi_j_beta>|^2
    double Sz               = 0.5 * (this->occupiedCountAlpha - this->occupiedCountBeta);
    double alphaBetaOverlap = (this->C_alpha.leftCols(this->occupiedCountAlpha).transpose() * this->S
                               * this->C_beta.leftCols(this->occupiedCountBeta))
                                  .squaredNorm();
    return (Sz * (Sz + 1)) + this->occupiedCountBeta - alphaBetaOverlap;
}


std::vector<std::string> SCF::getAOLabels() const
{
    std::vector<std::string> aoLabels;
    aoLabels.reserve(basisCount);
    const auto& basis = this->basis;

    const auto& shells = basis.shells;

    const auto& geometry = this->molecule->getGeometry();

    for (const auto& shell : shells)
    {
        const auto& angularMomentum = shell.angularMomentum;
        const auto& atom            = geometry[shell.atomIndex];

        // Get the number of atoms of same atomic number before this atom.
        size_t atomCount = std::ranges::count_if(
            geometry.begin(),
            geometry.begin() + shell.atomIndex,
            [&atom](const Atom& a) { return a.atomicNumber == atom.atomicNumber; }
        );

        for (size_t i = 0; i < shell.nAOs; ++i)
        {
            std::stringstream label;
            label << fmt::format(
                "{:<3} {:<2} {:<2} ", aoLabels.size() + 1, Utils::atomicNumberToName.at(atom.atomicNumber), atomCount + 1
            );
            const Eigen::Vector3i& am = angularMomentum.row(i);
            const unsigned l = am.x(), m = am.y(), n = am.z();

            switch (shell.l)
            {
                case 0: label << "s"; break; // s-orbital
                case 1:
                {
                    std::string sublabel;
                    if (l == 1)
                        sublabel = "x";
                    else if (m == 1)
                        sublabel = "y";
                    else if (n == 1)
                        sublabel = "z";
                    label << fmt::format("p{:<4}", sublabel);
                    break;
                }
                case 2:
                {
                    std::string sublabel;
                    if (l == 2)
                        sublabel = "x2";
                    else if (m == 2)
                        sublabel = "y2";
                    else if (n == 2)
                        sublabel = "z2";
                    else if (l == 1 && m == 1)
                        sublabel = "xy";
                    else if (l == 1 && n == 1)
                        sublabel = "xz";
                    else if (m == 1 && n == 1)
                        sublabel = "yz";
                    label << fmt::format("d{:<4}", sublabel);
                    break;
                }
                case 3:
                {
                    std::string sublabel;
                    if (l == 3)
                        sublabel = "x3";
                    else if (m == 3)
                        sublabel = "y3";
                    else if (n == 3)
                        sublabel = "z3";
                    else if (l == 2 && m == 1)
                        label << "x2y";
                    else if (l == 2 && n == 1)
                        label << "x2z";
                    else if (m == 2 && l == 1)
                        label << "y2x";
                    else if (m == 2 && n == 1)
                        label << "y2z";
                    else if (n == 2 && l == 1)
                        label << "z2x";
                    else if (n == 2 && m == 1)
                        label << "z2y";
                    else if (l == 1 && m == 1 && n == 1)
                        sublabel = "xyz";
                    label << fmt::format("f{:<4}", sublabel);
                    break;
                }
                case 4:
                {
                    std::string sublabel;
                    if (l == 4)
                        sublabel = "x4";
                    else if (m == 4)
                        sublabel = "y4";
                    else if (n == 4)
                        sublabel = "z4";
                    else if (l == 3 && m == 1)
                        sublabel = "x3y";
                    else if (l == 3 && n == 1)
                        sublabel = "x3z";
                    else if (m == 3 && l == 1)
                        sublabel = "y3x";
                    else if (m == 3 && n == 1)
                        sublabel = "y3z";
                    else if (n == 3 && l == 1)
                        sublabel = "z3x";
                    else if (n == 3 && m == 1)
                        sublabel = "z3y";
                    else if (l == 2 && m == 2)
                        sublabel = "x2y2";
                    else if (l == 2 && n == 2)
                        sublabel = "x2z2";
                    else if (m == 2 && n == 2)
                        sublabel = "y2z2";
                    else if (l == 1 && m == 1 && n == 2)
                        sublabel = "xyz2";
                    else if (l == 1 && n == 1 && m == 2)
                        sublabel = "xy2z";
                    else if (m == 1 && n == 1 && l == 2)
                        sublabel = "x2yz";
                    label << fmt::format("g{:<4}", sublabel);
                    break;
                }
                default:
                {
                    label << fmt::format("{:<4}", fmt::format("({}{}{})", l, m, n));
                    break;
                }
            }

            aoLabels.emplace_back(label.str());
        }
    }
    return aoLabels;
}


void SCF::printIteration() const
{
    std::stringstream ss;

    ss << fmt::format("{:^15}{:^15.5e}{:^15.5e}", this->iteration + 1, this->deltaE, this->deltaD);
    if (this->options.useDIIS)
        ss << fmt::format("{:^15.5e}", this->DIISError);
    if (this->dampCoeff > 0)
        ss << fmt::format("{:^15}", fmt::format("damp: {:.2f}", static_cast<double>(this->options.damp) / 100));
    if (this->useLevelShiftingAlpha || this->useLevelShiftingBeta)
        ss << fmt::format("{:^15}", fmt::format("lshift: {:.4f}", this->options.levelShift));
    ss << "\n";

    this->output->write(ss.str());
}


void SCF::printFinalResults(bool converged) const
{
    std::stringstream ss;
    ss << std::left;
    if (converged)
        ss << "\nSCF converged!\n";
    else
        ss << "\nWarning: Maximum iterations reached. SCF did not converge!\n";

    double totalEnergy = this->electronicEnergy + this->nuclearEnergy;
    ss << fmt::format(
        "\n--- SCF Results ---\n"
        "{0}"
        "Electronic Energy: {1:>15.10f}\n"
        "Nuclear Repulsion: {2:>15.10f}\n"
        "Total SCF Energy:  {3:>15.10f}\n",
        this->options.unrestricted ? fmt::format("<S^2>: {:>15.10f}\n", computeSpinSquared()) : "",
        this->electronicEnergy,
        this->nuclearEnergy,
        totalEnergy
    );

    std::vector<std::string> aoLabels = getAOLabels();

    ss << fmt::format("\n{:-<{}}\n", "", 99);
    std::string title = this->options.unrestricted ? "Alpha Orbital Energies (a.u.)" : "Orbital Energies (a.u.)";
    ss << fmt::format("{:^99}\n", title);
    ss << fmt::format("{:-<{}}\n", "", 99);

    ss << "-- Occupied --\n";
    ss << printShortMOs(this->eigenvalues_alpha.head(this->occupiedCountAlpha), 5, 7);
    ss << "-- Virtual --\n";
    ss << printShortMOs(this->eigenvalues_alpha.tail(this->basisCount - this->occupiedCountAlpha), 5, 7);
    ss << fmt::format("{:-<{}}\n", "", 99);

    if (this->options.unrestricted)
    {
        ss << fmt::format("\n{:-<{}}\n", "", 99);
        ss << fmt::format("{:^99}\n", "Beta Orbital Energies (a.u.)");
        ss << fmt::format("{:-<{}}\n", "", 99);

        ss << "-- Occupied --\n";
        ss << printShortMOs(this->eigenvalues_beta.head(this->occupiedCountBeta), 5, 7);
        ss << "-- Virtual --\n";
        ss << printShortMOs(this->eigenvalues_beta.tail(this->basisCount - this->occupiedCountBeta), 5, 7);
        for (size_t i = 0; i < 99; ++i) { ss << "-"; }
        ss << "\n";
    }

    if (this->options.printFullMOs)
    {
        ss << fmt::format("\n{:-<{}}\n", "", 99);
        std::string title = this->options.unrestricted ? "Alpha Molecular Orbital Coefficients"
                                                       : "Molecular Orbital Coefficients";
        ss << fmt::format("{:^99}\n", title);
        ss << fmt::format("{:-<{}}\n", "", 99);

        ss << "\nOccupied orbitals:\n";
        ss << "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n";
        ss << printFullMOs(
            this->C_alpha.leftCols(this->occupiedCountAlpha),
            this->eigenvalues_alpha.head(this->occupiedCountAlpha),
            aoLabels,
            5,
            5
        );

        ss << "\nVirtual orbitals:\n";
        ss << "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n";
        ss << printFullMOs(
            this->C_alpha.rightCols(this->basisCount - this->occupiedCountAlpha),
            this->eigenvalues_alpha.tail(this->basisCount - this->occupiedCountAlpha),
            aoLabels,
            5,
            5
        );
        ss << fmt::format("{:-<{}}\n", "", 99);

        if (this->options.unrestricted)
        {
            ss << fmt::format("\n{:-<{}}\n", "", 99);
            ss << fmt::format("{:^99}\n", "Beta Molecular Orbital Coefficients");
            ss << fmt::format("{:-<{}}\n", "", 99);

            ss << "\nOccupied orbitals:\n";
            ss << "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n";
            ss << printFullMOs(
                this->C_beta.leftCols(this->occupiedCountBeta),
                this->eigenvalues_beta.head(this->occupiedCountBeta),
                aoLabels,
                5,
                5
            );

            ss << "\nVirtual orbitals:\n";
            ss << "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n";
            ss << printFullMOs(
                this->C_beta.rightCols(this->basisCount - this->occupiedCountBeta),
                this->eigenvalues_beta.tail(this->basisCount - this->occupiedCountBeta),
                aoLabels,
                5,
                5
            );
            ss << fmt::format("{:-<{}}\n", "", 99);
        }
    }

    this->output->write(ss.str());
}

std::string SCF::printShortMOs(const Eigen::VectorXd& eigenvalues, size_t precision, size_t MOsPerRow) const
{
    const int width = 6 + precision;

    auto format = [&precision, &width](double num) -> std::string
    { return fmt::format("{}{:<{}.{}e} ", std::signbit(num) ? "-" : " ", std::abs(num), width - 1, precision); };

    std::stringstream ss;
    ss << std::left;

    size_t rowsLimit = std::ceil(static_cast<double>(eigenvalues.size()) / static_cast<double>(MOsPerRow));
    for (size_t k = 0; k < rowsLimit; ++k)
    {
        size_t startIndex = k * MOsPerRow;
        size_t numMOs     = eigenvalues.size();
        size_t endIndex   = std::min(startIndex + MOsPerRow, numMOs);

        for (size_t i = startIndex; i < endIndex; ++i) { ss << format(eigenvalues(i)); }
        ss << "\n";
    }

    return ss.str();
}

std::string SCF::printFullMOs(
    const Eigen::MatrixXd& MOs,
    const Eigen::VectorXd& eigenvalues,
    const std::vector<std::string>& aoLabels,
    size_t precision,
    size_t MOsPerRow
) const
{
    const int width = 10 + precision;

    auto format = [&precision, &width](double num) -> std::string
    { return fmt::format("{}{:<{}.{}e}", std::signbit(num) ? "-" : " ", std::abs(num), width - 1, precision); };

    std::stringstream ss;
    ss << std::left;

    size_t rowsLimit = std::ceil(static_cast<double>(MOs.cols()) / static_cast<double>(MOsPerRow));
    for (size_t k = 0; k < rowsLimit; ++k)
    {
        size_t startIndex = k * MOsPerRow;
        size_t numMOs     = MOs.cols();
        size_t endIndex   = std::min(startIndex + MOsPerRow, numMOs);
        ss << fmt::format("{:<{}}\t", " ", width);
        for (size_t i = startIndex; i < endIndex; ++i) { ss << fmt::format("{:^{}}\t", i + 1, width - 1, 0); }
        ss << "\n";

        ss << fmt::format("{:<{}}\t", "eigenvalues:", width);
        for (size_t i = startIndex; i < endIndex; ++i) { ss << format(eigenvalues(i)) << "\t"; }
        ss << "\n";

        for (size_t i = 0; i < this->basisCount; ++i)
        {
            ss << fmt::format("{:<{}}\t", aoLabels[i], width);
            for (size_t j = startIndex; j < endIndex; ++j) { ss << format(MOs(i, j)) << "\t"; }
            ss << "\n";
        }
        ss << "\n";
    }

    return ss.str();
}

std::string SCF::printJobSpec() const
{
    std::stringstream ss;

    ss << fmt::format(
        "Starting {} calculation.\n\n"
        "Maximum Iterations: {}.\n"
        "Energy Convergence Threshold: {:.6e}.\n"
        "Density Convergence Threshold: {:.6e}.\n"
        "Schwartz Screening Threshold: {:.6e}.\n\n",
        this->options.unrestricted ? "UHF" : "RHF",
        this->options.maxIter,
        this->options.energyTol,
        this->options.densityTol,
        this->options.schwartzThreshold
    );

    if (this->options.direct)
    {
        ss << fmt::format("Direct SCF enabled.\nDensity screening threshold: {:.6e}.\n\n", this->options.densityThreshold);
    }
    if (this->options.useDIIS)
    {
        ss << fmt::format(
            "DIIS enabled.\nSize of DIIS history: {}.\nDIIS error threshold: {:.6e}.\n\n",
            this->options.DIISmaxSize,
            this->options.DIISErrorTol
        );
    }
    if (this->options.guessMix > 0)
    {
        ss << fmt::format("Requested {}% mixing of HOMO and LUMO orbitals in initial guess.\n\n", this->options.guessMix * 10);
    }
    if (this->options.damp > 0)
    {
        ss << fmt::format("Damping enabled. Mixing coefficient set to: {:.2f}.\n", static_cast<double>(this->options.damp) / 100);
        if (this->options.maxDampIter < this->options.maxIter)
            ss << fmt::format("Damping will stop after iteration {}.\n", this->options.maxDampIter);
        if (this->options.stopDampThresh != 0)
            ss << fmt::format("Damping will stop once ΔE < {:.6e}.\n", this->options.stopDampThresh);
        ss << std::endl;
    }
    if (options.levelShift > 0)
    {
        ss << fmt::format("Level-shifting enabled. Level-shift set to {:.4f}", this->options.levelShift);
        if (this->options.maxLshiftIter < this->options.maxLshiftIter)
            ss << fmt::format("Level-shifting will stop after iteration {}.\n", this->options.maxLshiftIter);
        if (this->options.stopLshiftThresh > 0)
            ss << fmt::format("Level-shifting will stop once ΔE < {:.6e}.\n", this->options.stopLshiftThresh);
        if (!std::isinf(this->options.lshiftGapTol))
            ss << fmt::format("Level-shifting will be applied when HOMO-LUMO gap is less than {:.4f}.\n", this->options.lshiftGapTol);
        ss << std::endl;
    }

    ss << fmt::format(
        "Molecule with {} alpha and {} beta electrons.\n"
        "Number of basis functions: {}.\n"
        "Geometry:\n",
        this->occupiedCountAlpha,
        this->occupiedCountBeta,
        this->basisCount
    );

    for (const auto& atom : this->molecule->getGeometry())
    {
        ss << fmt::format(
            "\t{: <2} {: >13.8f} {: >13.8f} {: >13.8f}\n",
            Utils::atomicNumberToName.at(atom.atomicNumber),
            atom.coords.x(),
            atom.coords.y(),
            atom.coords.z()
        );
    }

    return ss.str();
}
