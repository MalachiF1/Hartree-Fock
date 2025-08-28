#include "SCF.hpp"

#include "DIIS.hpp"
#include "Utils.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <ios>
#include <iostream>

SCF::SCF(const Molecule& molecule, const SCFOptions& options, std::shared_ptr<Output> output) :
    output(output), molecule(std::make_unique<Molecule>(molecule)), options(options)
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

    this->output->write("\n--- Starting SCF Iterations ---\n");

    // SCF cycles
    for (size_t iteration = 0; iteration < this->options.maxIter; ++iteration)
    {
        // Store the energy and density from the previous iteration to check for convergence.
        double lastElectronicEnergy = this->electronicEnergy;
        Eigen::MatrixXd D_last      = this->D_tot;

        // Build the new Fock matrix from the current density.
        if (this->options.direct)
            this->buildFockMatrix(this->options.schwartzThreshold, this->options.densityThreshold);
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

        if (this->options.useDIIS && iteration >= this->options.DIISstart)
        {
            // Extrapolate the orthogonal Fock matrix using DIIS.
            std::tie(this->F_alpha, this->F_beta) = diis_handler->extrapolate(this->F_alpha, this->F_beta);
        }

        this->diagonalizeAndUpdate();

        // Calculate the change in energy and density.
        double dE = std::abs(this->electronicEnergy - lastElectronicEnergy);
        double dD = (this->D_tot - D_last).norm();

        // Print the current iteration results.
        if (this->options.useDIIS)
            this->printIteration(iteration + 1, dE, dD, diis_handler->getErrorNorm());
        else
            this->printIteration(iteration + 1, dE, dD);

        // Check for convergence.
        if (dE < this->options.energyTol && dD < this->options.densityTol
            && (!this->options.useDIIS || diis_handler->getErrorNorm() < this->options.DIISErrorTol))
        {
            // basis If converged, print the final results and exit.
            this->printFinalResults(true);
            return;
        }
    }

    // If the loop finishes, the calculation did not converge.
    this->printFinalResults(false);
}

void SCF::initialize(double schwartzThreshold, bool direct)
{
    // Set constant parameters.
    this->basisCount         = molecule->getBasisFunctionCount();
    this->occupiedCountAlpha = (molecule->getElectronCount() + molecule->getMultiplicity() - 1) / 2;
    this->occupiedCountBeta  = (molecule->getElectronCount() - molecule->getMultiplicity() + 1) / 2;

    // Calculate one-electron integrals and nuclear repulsion.
    this->nuclearEnergy = molecule->nuclearRepulsion();
    this->T             = molecule->kineticMatrix();
    this->V             = molecule->nuclearAttractionMatrix();
    this->h             = this->T + this->V;

    // Calculate the overlap and orthogonalization matrices.
    this->S = molecule->overlapMatrix();
    this->X = inverseSqrtMatrix(this->S);

    // Calculate two-electron integrals (only if not using direct method).
    if (!direct)
        this->Vee = molecule->electronRepulsionTensor(schwartzThreshold);
    else
        // If using direct method, only pre-calculate the Schwartz screening matrix.
        this->Q = molecule->schwartzScreeningMatrix();

    this->output->write(printJobSpec());
}


void SCF::computeInitialGuessDensity()
{
    // The initial Fock matrix is just the core Hamiltonian.
    this->F_alpha = this->h;

    // Diagonalize the initial Fock matrix to get guess orbitals.
    Eigen::MatrixXd F_alpha_prime = this->X.transpose() * this->F_alpha * this->X;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> alpha_solver(F_alpha_prime);
    this->C_alpha               = this->X * alpha_solver.eigenvectors();
    this->eigenvalues_alpha     = alpha_solver.eigenvalues();
    Eigen::MatrixXd C_alpha_occ = this->C_alpha.leftCols(this->occupiedCountAlpha);

    // Form the density matrix from the occupied orbitals.
    this->D_alpha = C_alpha_occ * C_alpha_occ.transpose();


    if (this->options.unrestricted)
    {
        // HOMO-LUMO mixing for breaking symmetry in UHF.
        if (this->options.guessMix > 0)
        {
            double k = static_cast<double>(this->options.guessMix) / 10.0;

            Eigen::VectorXd homoAlpha = C_alpha.col(this->occupiedCountAlpha - 1);
            Eigen::VectorXd lumoAlpha = C_alpha.col(this->occupiedCountAlpha);

            Eigen::VectorXd mixedHomo = (1 / (std::sqrt(1 + (k * k)))) * (homoAlpha + k * lumoAlpha);
            Eigen::VectorXd mixedLumo = (1 / (std::sqrt(1 + (k * k)))) * (lumoAlpha - k * homoAlpha);

            C_alpha.col(this->occupiedCountAlpha - 1) = mixedHomo;
            C_alpha.col(this->occupiedCountAlpha)     = mixedLumo;
            Eigen::MatrixXd C_alpha_occ               = this->C_alpha.leftCols(this->occupiedCountAlpha);
            this->D_alpha                             = C_alpha_occ * C_alpha_occ.transpose();
        }

        this->F_beta                 = this->h;
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
        this->C_beta           = this->C_alpha;
        this->eigenvalues_beta = this->eigenvalues_alpha;
        this->D_beta           = this->D_alpha;
    }

    this->D_tot = this->D_alpha + this->D_beta;

    // Calculate the initial electronic energy.
    this->electronicEnergy = 0.5 * (this->D_alpha * (this->h + this->F_alpha)).trace()
                           + 0.5 * (this->D_beta * (this->h + this->F_beta)).trace();
}


void SCF::buildFockMatrix()
{
    size_t N_ao             = this->basisCount;
    Eigen::MatrixXd G_alpha = Eigen::MatrixXd::Zero(N_ao, N_ao);
    Eigen::MatrixXd G_beta;
    if (this->options.unrestricted)
        G_beta = Eigen::MatrixXd::Zero(N_ao, N_ao);

// Calculate the G matrix (the two-electron part of the Fock matrix).
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < N_ao; ++i)
    {
        for (size_t j = 0; j < N_ao; ++j)
        {
            if (j <= i)
            {
                double G_alpha_ij = 0.0;
                double G_beta_ij  = 0.0;
                for (size_t k = 0; k < N_ao; ++k)
                {
                    for (size_t l = 0; l < N_ao; ++l)
                    {
                        double J_kl = this->Vee(i, j, k, l);
                        double K_kl = this->Vee(i, k, j, l);
                        G_alpha_ij += (this->D_alpha(k, l) + this->D_beta(k, l)) * J_kl - D_alpha(k, l) * K_kl;

                        if (this->options.unrestricted)
                            G_beta_ij += (this->D_alpha(k, l) + this->D_beta(k, l)) * J_kl - D_beta(k, l) * K_kl;
                    }
                }
                G_alpha(i, j) = G_alpha(j, i) = G_alpha_ij;

                if (this->options.unrestricted)
                    G_beta(i, j) = G_beta(j, i) = G_beta_ij;
            }
        }
    }
    this->F_alpha = this->h + G_alpha;

    if (this->options.unrestricted)
        this->F_beta = this->h + G_beta;
    else
        this->F_beta = this->F_alpha;
}

void SCF::buildFockMatrix(double schwartzThreshold, double densityThreshold)
{
    const size_t& N_ao = this->basisCount;
    const auto& AOs    = this->molecule->getAtomicOrbitals();

    Eigen::MatrixXd J       = Eigen::MatrixXd::Zero(N_ao, N_ao);
    Eigen::MatrixXd K_alpha = Eigen::MatrixXd::Zero(N_ao, N_ao);
    Eigen::MatrixXd K_beta;
    if (this->options.unrestricted)
        K_beta = Eigen::MatrixXd::Zero(N_ao, N_ao);


#pragma omp parallel
    {
        // Thread-local matrices to void data races
        Eigen::MatrixXd J_p       = Eigen::MatrixXd::Zero(N_ao, N_ao);
        Eigen::MatrixXd K_alpha_p = Eigen::MatrixXd::Zero(N_ao, N_ao);
        Eigen::MatrixXd K_beta_p;
        if (this->options.unrestricted)
            K_beta_p = Eigen::MatrixXd::Zero(N_ao, N_ao);


#pragma omp for collapse(4) schedule(dynamic, 1)
        for (size_t i = 0; i < N_ao; ++i)
        {
            for (size_t j = 0; j < N_ao; ++j)
            {
                for (size_t k = 0; k < N_ao; ++k)
                {
                    for (size_t l = 0; l < N_ao; ++l)
                    {
                        // Only calculate unique integrals
                        if (j > i || l > k || ((i * (i + 1) / 2 + j) < (k * (k + 1) / 2 + l)))
                            continue;

                        // Schwartz screening
                        if (Q(i, j) * Q(k, l) < schwartzThreshold)
                            continue;

                        // Density screening

                        if (std::abs(this->D_tot(k, l)) < densityThreshold && std::abs(this->D_tot(i, j)) < densityThreshold
                            && std::abs(this->D_tot(i, k)) < densityThreshold
                            && std::abs(this->D_tot(i, l)) < densityThreshold && std::abs(this->D_tot(j, k)) < densityThreshold
                            && std::abs(this->D_tot(j, l)) < densityThreshold)
                            continue;


                        double eri = AtomicOrbital::electronRepulsion(AOs[i], AOs[j], AOs[k], AOs[l]);

                        if (i == j && j == k && k == l) // (ii|ii)
                        {
                            J_p(i, i) += D_tot(i, i) * eri;
                            K_alpha_p(i, i) += D_alpha(i, i) * eri;

                            if (this->options.unrestricted)
                            {
                                K_beta_p(i, i) += D_beta(i, i) * eri;
                            }
                        }
                        else if (i == j && j == k) // (ii|il)
                        {
                            J_p(i, i) += 2 * this->D_tot(i, l) * eri;
                            J_p(l, i) += this->D_tot(i, i) * eri;
                            J_p(i, l) += this->D_tot(i, i) * eri;
                            K_alpha_p(i, i) += 2 * this->D_alpha(i, l) * eri;
                            K_alpha_p(l, i) += this->D_alpha(i, i) * eri;
                            K_alpha_p(i, l) += this->D_alpha(i, i) * eri;

                            if (this->options.unrestricted)
                            {
                                K_beta_p(i, i) += 2 * this->D_beta(i, l) * eri;
                                K_beta_p(l, i) += this->D_beta(i, i) * eri;
                                K_beta_p(i, l) += this->D_beta(i, i) * eri;
                            }
                        }
                        else if (i == j && k == l) // (ii|kk)
                        {
                            J_p(i, i) += this->D_tot(k, k) * eri;
                            J_p(k, k) += this->D_tot(i, i) * eri;
                            K_alpha_p(i, k) += this->D_alpha(i, k) * eri;
                            K_alpha_p(k, i) += this->D_alpha(i, k) * eri;

                            if (this->options.unrestricted)
                            {
                                K_beta_p(i, k) += this->D_beta(i, k) * eri;
                                K_beta_p(k, i) += this->D_beta(i, k) * eri;
                            }
                        }
                        else if (i == k && j == l) // (ij|ij)
                        {
                            // (ij|ij) = (ji|ij) = (ij|ji) = (ji|ji)
                            J_p(i, j) += 2 * this->D_tot(i, j) * eri;
                            J_p(j, i) += 2 * this->D_tot(i, j) * eri;
                            K_alpha_p(i, j) += this->D_alpha(j, i) * eri;
                            K_alpha_p(j, i) += this->D_alpha(i, j) * eri;
                            K_alpha_p(i, i) += this->D_alpha(j, j) * eri;
                            K_alpha_p(j, j) += this->D_alpha(i, i) * eri;

                            if (this->options.unrestricted)
                            {
                                K_beta_p(i, j) += this->D_beta(j, i) * eri;
                                K_beta_p(j, i) += this->D_beta(i, j) * eri;
                                K_beta_p(i, i) += this->D_beta(j, j) * eri;
                                K_beta_p(j, j) += this->D_beta(i, i) * eri;
                            }
                        }
                        else if (i == j) // (ii|kl)
                        {
                            J_p(i, i) += 2 * this->D_tot(k, l) * eri;
                            J_p(k, l) += this->D_tot(i, i) * eri;
                            J_p(l, k) += this->D_tot(i, i) * eri;
                            K_alpha_p(i, l) += this->D_alpha(i, k) * eri;
                            K_alpha_p(i, k) += this->D_alpha(i, l) * eri;
                            K_alpha_p(k, i) += this->D_alpha(l, i) * eri;
                            K_alpha_p(l, i) += this->D_alpha(k, i) * eri;

                            if (this->options.unrestricted)
                            {
                                K_beta_p(i, l) += this->D_beta(i, k) * eri;
                                K_beta_p(i, k) += this->D_beta(i, l) * eri;
                                K_beta_p(k, i) += this->D_beta(l, i) * eri;
                                K_beta_p(l, i) += this->D_beta(k, i) * eri;
                            }
                        }
                        else if (k == l) // (ij|kk)
                        {
                            J_p(i, j) += this->D_tot(k, k) * eri;
                            J_p(j, i) += this->D_tot(k, k) * eri;
                            J_p(k, k) += 2 * this->D_tot(i, j) * eri;
                            K_alpha_p(i, k) += this->D_alpha(j, k) * eri;
                            K_alpha_p(j, k) += this->D_alpha(i, k) * eri;
                            K_alpha_p(k, j) += this->D_alpha(k, i) * eri;
                            K_alpha_p(k, i) += this->D_alpha(k, j) * eri;
                            if (this->options.unrestricted)
                            {
                                K_beta_p(i, k) += this->D_beta(j, k) * eri;
                                K_beta_p(j, k) += this->D_beta(i, k) * eri;
                                K_beta_p(k, j) += this->D_beta(k, i) * eri;
                                K_beta_p(k, i) += this->D_beta(k, j) * eri;
                            }
                        }
                        else // (ij|kl)
                        {
                            J_p(i, j) += 2 * this->D_tot(k, l) * eri;
                            J_p(j, i) += 2 * this->D_tot(k, l) * eri;
                            J_p(k, l) += 2 * this->D_tot(i, j) * eri;
                            J_p(l, k) += 2 * this->D_tot(i, j) * eri;
                            K_alpha_p(i, l) += this->D_alpha(j, k) * eri;
                            K_alpha_p(j, l) += this->D_alpha(i, k) * eri;
                            K_alpha_p(i, k) += this->D_alpha(j, l) * eri;
                            K_alpha_p(j, k) += this->D_alpha(i, l) * eri;
                            K_alpha_p(k, j) += this->D_alpha(l, i) * eri;
                            K_alpha_p(l, j) += this->D_alpha(k, i) * eri;
                            K_alpha_p(k, i) += this->D_alpha(l, j) * eri;
                            K_alpha_p(l, i) += this->D_alpha(k, j) * eri;

                            if (this->options.unrestricted)
                            {
                                K_beta_p(i, l) += this->D_beta(j, k) * eri;
                                K_beta_p(j, l) += this->D_beta(i, k) * eri;
                                K_beta_p(i, k) += this->D_beta(j, l) * eri;
                                K_beta_p(j, k) += this->D_beta(i, l) * eri;
                                K_beta_p(k, j) += this->D_beta(l, i) * eri;
                                K_beta_p(l, j) += this->D_beta(k, i) * eri;
                                K_beta_p(k, i) += this->D_beta(l, j) * eri;
                                K_beta_p(l, i) += this->D_beta(k, j) * eri;
                            }
                        }
                    }
                }
            }
        }
#pragma omp critical
        {
            // Accumulate the thread-local G matrix into the global G matrix.
            J += J_p;
            K_alpha += K_alpha_p;
            if (this->options.unrestricted)
                K_beta += K_beta_p;
        }
    }

    this->F_alpha = this->h + J - K_alpha;
    if (this->options.unrestricted)
        this->F_beta = this->h + J - K_beta;
    else
        this->F_beta = this->F_alpha;
}


void SCF::diagonalizeAndUpdate()
{
    // Transform the Fock matrix to the orthogonal basis.
    Eigen::MatrixXd F_alpha_prime = this->X.transpose() * this->F_alpha * this->X;

    // Diagonalize the Fock matrix.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> alpha_solver(F_alpha_prime);
    this->C_alpha           = this->X * alpha_solver.eigenvectors();
    this->eigenvalues_alpha = alpha_solver.eigenvalues();

    // Form the new density matrix from the occupied orbitals.
    Eigen::MatrixXd C_alpha_occ = this->C_alpha.leftCols(this->occupiedCountAlpha);
    this->D_alpha               = C_alpha_occ * C_alpha_occ.transpose();

    if (this->options.unrestricted)
    {
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
        this->C_beta           = this->C_alpha;
        this->eigenvalues_beta = this->eigenvalues_alpha;
        this->D_beta           = this->D_alpha;
    }

    this->D_tot = this->D_alpha + this->D_beta;


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

void SCF::printIteration(int iter, double dE, double dD) const
{
    std::stringstream ss;
    ss << "Iteration " << std::setw(3) << iter << ": "
       << " | ΔE = " << std::scientific << dE << " | ΔD = " << dD << std::endl;
    this->output->write(ss.str());
}


void SCF::printIteration(int iter, double dE, double dD, double DIISError) const
{
    std::stringstream ss;
    ss << "Iteration " << std::setw(3) << iter << ": "
       << " | ΔE = " << std::scientific << dE << " | ΔD = " << dD << " | DIIS error = " << DIISError << std::endl;
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
    ss << "\n--- SCF Results ---\n";
    ss << std::fixed << std::setprecision(10);
    if (this->options.unrestricted)
        ss << "<S^2>: " << computeSpinSquared() << "\n";
    ss << "Electronic Energy: " << this->electronicEnergy << "\n";
    ss << "Nuclear Repulsion: " << this->nuclearEnergy << "\n";
    ss << "Total SCF Energy: " << totalEnergy << "\n";

    std::vector<std::string> aoLabels = this->molecule->getAOLabels();

    ss << "\n";
    for (size_t i = 0; i < 99; ++i) { ss << "-"; }
    ss << "\n";
    std::string title      = this->options.unrestricted ? "Alpha Orbital Energies (a.u.)" : "Orbital Energies (a.u.)";
    size_t whitespaceWidth = (99 - title.length()) / 2;
    ss << std::string(whitespaceWidth, ' ') << title << std::string(whitespaceWidth, ' ') << "\n";
    for (size_t i = 0; i < 99; ++i) { ss << "-"; }
    ss << "\n";

    ss << "-- Occupied --\n";
    ss << printShortMOs(this->eigenvalues_alpha.head(this->occupiedCountAlpha), 5, 7);
    ss << "-- Virtual --\n";
    ss << printShortMOs(this->eigenvalues_alpha.tail(this->basisCount - this->occupiedCountAlpha), 5, 7);
    for (size_t i = 0; i < 99; ++i) { ss << "-"; }
    ss << "\n";

    if (this->options.unrestricted)
    {
        ss << "\n";
        for (size_t i = 0; i < 99; ++i) { ss << "-"; }
        ss << "\n";
        whitespaceWidth = (99 - std::string("Beta Orbital Energies (a.u.)").length()) / 2;
        ss << std::string(whitespaceWidth, ' ') << "Beta Orbital Energies (a.u.)" << std::string(whitespaceWidth, ' ')
           << "\n";
        for (size_t i = 0; i < 99; ++i) { ss << "-"; }
        ss << "\n";

        ss << "-- Occupied --\n";
        ss << printShortMOs(this->eigenvalues_beta.head(this->occupiedCountBeta), 5, 7);
        ss << "-- Virtual --\n";
        ss << printShortMOs(this->eigenvalues_beta.tail(this->basisCount - this->occupiedCountBeta), 5, 7);
        for (size_t i = 0; i < 99; ++i) { ss << "-"; }
        ss << "\n";
    }

    if (this->options.printFullMOs)
    {
        ss << "\n";
        for (size_t i = 0; i < 99; ++i) { ss << "-"; }
        ss << "\n";
        std::string title = this->options.unrestricted ? "Alpha Molecular Orbital Coefficients"
                                                       : "Molecular Orbital Coefficients";
        whitespaceWidth   = (99 - title.length()) / 2;
        ss << std::string(whitespaceWidth, ' ') << title << std::string(whitespaceWidth, ' ') << "\n";
        for (size_t i = 0; i < 99; ++i) { ss << "-"; }
        ss << "\n";

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
        for (size_t i = 0; i < 99; ++i) { ss << "-"; }
        ss << "\n";

        if (this->options.unrestricted)
        {
            ss << "\n";
            for (size_t i = 0; i < 99; ++i) { ss << "-"; }
            ss << "\n";
            whitespaceWidth = (99 - std::string("Beta Molecular Orbital Coefficients").length()) / 2;
            ss << std::string(whitespaceWidth, ' ') << "Beta Molecular Orbital Coefficients"
               << std::string(whitespaceWidth, ' ') << "\n";
            for (size_t i = 0; i < 99; ++i) { ss << "-"; }
            ss << "\n";

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
            for (size_t i = 0; i < 99; ++i) { ss << "-"; }
        }
    }

    this->output->write(ss.str());
}

std::string SCF::printShortMOs(const Eigen::VectorXd& eigenvalues, size_t precision, size_t MOsPerRow) const
{
    const int width = 6 + precision;

    auto format = [&precision, &width](double num) -> std::string
    {
        std::stringstream ss;

        if (std::signbit(num))
            ss << "-";
        else
            ss << " ";

        ss << std::scientific << std::left << std::setw(width - 1) << std::setprecision(precision) << std::abs(num);

        return ss.str();
    };

    std::stringstream ss;
    ss << std::left;

    size_t rowsLimit = std::ceil(static_cast<double>(eigenvalues.size()) / static_cast<double>(MOsPerRow));
    for (size_t k = 0; k < rowsLimit; ++k)
    {
        size_t startIndex = k * MOsPerRow;
        size_t numMOs     = eigenvalues.size();
        size_t endIndex   = std::min(startIndex + MOsPerRow, numMOs);

        for (size_t i = startIndex; i < endIndex; ++i) { ss << format(eigenvalues(i)) << std::setw(2) << " "; }
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
    {
        std::stringstream ss;

        if (std::signbit(num))
            ss << "-";
        else
            ss << " ";

        ss << std::scientific << std::left << std::setw(width - 1) << std::setprecision(precision) << std::abs(num);

        return ss.str();
    };

    std::stringstream ss;
    ss << std::left;

    size_t rowsLimit = std::ceil(static_cast<double>(MOs.cols()) / static_cast<double>(MOsPerRow));
    for (size_t k = 0; k < rowsLimit; ++k)
    {
        size_t startIndex = k * MOsPerRow;
        size_t numMOs     = MOs.cols();
        size_t endIndex   = std::min(startIndex + MOsPerRow, numMOs);
        ss << std::setw(width) << " " << "\t";
        for (size_t i = startIndex; i < endIndex; ++i)
        {
            double whitespaceWidth = precision + 8;
            ss << std::string(std::floor((whitespaceWidth - 1) / 2), ' ');
            ss << i + 1;
            ss << std::string(std::ceil((whitespaceWidth - 1) / 2), ' ') << "\t";
        }
        ss << "\n";

        ss << std::setw(width) << "eigenvalues:" << "\t";
        for (size_t i = startIndex; i < endIndex; ++i) { ss << format(eigenvalues(i)) << "\t"; }
        ss << "\n";

        for (size_t i = 0; i < this->basisCount; ++i)
        {
            ss << std::setw(width) << aoLabels[i] << "\t";
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
    ss << "Starting " << (this->options.unrestricted ? "UHF" : "RHF") << " calculation.\n";
    ss << "Maximum iterations: " << this->options.maxIter << "\n";
    ss << "Energy convergence threshold: " << this->options.energyTol << "\n";
    ss << "Density convergence threshold: " << this->options.densityTol << "\n";
    if (this->options.useDIIS)
    {
        ss << "DIIS enabled. Size of DIIS history: " << this->options.DIISmaxSize << ", start from iteration: "
           << this->options.DIISstart << ", error threshold: " << this->options.DIISErrorTol << "\n";
    }
    ss << "Schwartz screening threshold: " << this->options.schwartzThreshold << "\n";
    if (this->options.direct)
    {
        ss << "Using direct SCF.\n";
        ss << "Density screening threshold: " << this->options.densityThreshold << "\n";
    }
    if (this->options.guessMix > 0)
    {
        ss << "Requested " << this->options.guessMix * 10 << "% mixing of HOMO and LUMO orbitals in initial guess.\n";
    }
    ss << "Molecule with " << this->occupiedCountAlpha << " alpha and " << this->occupiedCountBeta << " beta electrons.\n";
    ss << "Number of basis functions: " << this->basisCount << "\n";
    ss << "Geometry:\n";

    auto format = [](double num) -> std::string
    {
        std::stringstream ss;

        if (std::signbit(num))
            ss << "-";
        else
            ss << " ";

        ss << std::left << std::fixed << std::setprecision(8) << std::abs(num);

        return ss.str();
    };

    for (const auto& atom : this->molecule->getGeometry())
    {
        ss << "\t" << Utils::atomicNumberToName.at(atom.atomicNumber) << "\t" << format(atom.coords.x()) << "\t"
           << format(atom.coords.y()) << "\t" << format(atom.coords.z()) << "\n";
    }
    if (this->options.useSymmetry)
        ss << "Point Group: " << this->molecule->getPointGroup().toString() << "\n";

    return ss.str();
}
