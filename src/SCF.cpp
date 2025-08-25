#include "SCF.hpp"

#include "DIIS.hpp"
#include "Utils.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <ios>
#include <iostream>

SCF::SCF(const Molecule& molecule) : molecule(molecule) {}

void SCF::run(
    size_t maxIter,
    double energyTol,
    double densityTol,
    double schwartzThreshold,
    bool useDIIS,
    size_t DIISmaxSize,
    unsigned DIISstart,
    double DIISErrorTol,
    bool direct,
    double densityThreshold
)
{
    // Initialize the SCF calculation. (This should be called only once.)
    this->initialize(schwartzThreshold, direct);

    // Get the initial guess density.
    this->computeInitialGuessDensity();

    if (useDIIS) // Initialize the DIIS handler with the specified maximum size.
        diis_handler = std::make_unique<DIIS>(DIISmaxSize);

    std::cout << "\n--- Starting SCF Iterations ---\n" << std::endl;

    // SCF cycles
    for (size_t iteration = 0; iteration < maxIter; ++iteration)
    {
        // Store the energy and density from the previous iteration to check for convergence.
        double lastElectronicEnergy = this->electronicEnergy;
        Eigen::MatrixXd D_last      = this->D;

        // Build the new Fock matrix from the current density.
        if (direct)
            this->buildFockMatrix(schwartzThreshold, densityThreshold);
        else
            this->buildFockMatrix();

        // Store the error vector and fock matrix for DIIS if applicable.
        if (useDIIS)
        {
            // Calculate the commutator of the fock and density matrices as the error vector.
            Eigen::MatrixXd errorVector = this->F * this->D * this->S - this->S * this->D * this->F;
            diis_handler->update(this->F, errorVector);
        }

        // Diagonalize the Fock matrix and compute the new density and coefficients.

        if (useDIIS && iteration >= DIISstart) // Extrapolate the orthogonal Fock matrix using DIIS.
            this->F = diis_handler->extrapolate(this->F);

        this->diagonalizeAndUpdate();

        // Calculate the change in energy and density.
        double dE = std::abs(this->electronicEnergy - lastElectronicEnergy);
        double dD = (this->D - D_last).norm();

        // Print the current iteration results.
        if (useDIIS)
            this->printIteration(iteration + 1, dE, dD, diis_handler->getErrorNorm());
        else
            this->printIteration(iteration + 1, dE, dD);

        // Check for convergence.
        if (dE < energyTol && dD < densityTol && (!useDIIS || diis_handler->getErrorNorm() < DIISErrorTol))
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
    std::cout << "Starting SCF calculation..." << std::endl;
    std::cout << molecule.toString() << std::endl;

    // Set constant parameters.
    this->basisCount    = molecule.getBasisFunctionCount();
    this->occupiedCount = molecule.getElectronCount() / 2; // Assuming closed-shell system (RHF)

    // Calculate one-electron integrals and nuclear repulsion.
    this->nuclearEnergy = molecule.nuclearRepulsion();
    this->T             = molecule.kineticMatrix();
    this->V             = molecule.nuclearAttractionMatrix();
    this->h             = this->T + this->V;

    // Calculate the overlap and orthogonalization matrices.
    this->S = molecule.overlapMatrix();
    this->X = inverseSqrtMatrix(this->S);

    // Calculate two-electron integrals (only if not using direct method).
    if (!direct)
        this->Vee = molecule.electronRepulsionTensor(schwartzThreshold);
    else
        // If using direct method, only pre-calculate the Schwartz screening matrix.
        this->Q = molecule.schwartzScreeningMatrix();
}


void SCF::computeInitialGuessDensity()
{
    // The initial Fock matrix is just the core Hamiltonian.
    this->F = this->h;

    // Diagonalize the initial Fock matrix to get guess orbitals.
    Eigen::MatrixXd F_prime = this->X.transpose() * this->F * this->X;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(F_prime);
    this->C           = this->X * solver.eigenvectors();
    this->eigenvalues = solver.eigenvalues();

    // Form the density matrix from the occupied orbitals.
    Eigen::MatrixXd C_occ = this->C.leftCols(this->occupiedCount);
    this->D               = 2 * C_occ * C_occ.transpose();

    // Calculate the initial electronic energy.
    this->electronicEnergy = (this->D.array() * (this->h + this->F).array()).sum() * 0.5;
}


void SCF::buildFockMatrix()
{
    size_t N_ao       = this->basisCount;
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(N_ao, N_ao);

// Calculate the G matrix (the two-electron part of the Fock matrix).
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < N_ao; ++i)
    {
        for (size_t j = 0; j < N_ao; ++j)
        {
            if (j <= i)
            {
                double G_ij = 0.0;
                for (size_t k = 0; k < N_ao; ++k)
                {
                    for (size_t l = 0; l < N_ao; ++l)
                    {
                        double J_kl = this->Vee(i, j, k, l);
                        double K_kl = this->Vee(i, k, j, l);
                        G_ij += this->D(k, l) * (J_kl - 0.5 * K_kl);
                    }
                }
                G(i, j) = G(j, i) = G_ij;
            }
        }
    }
    this->F = this->h + G;
}

void SCF::buildFockMatrix(double schwartzThreshold, double densityThreshold)
{
    const size_t& N_ao = this->basisCount;
    Eigen::MatrixXd G  = Eigen::MatrixXd::Zero(N_ao, N_ao);
    const auto& AOs    = this->molecule.getAtomicOrbitals();

#pragma omp parallel
    {
        Eigen::MatrixXd G_private = Eigen::MatrixXd::Zero(N_ao, N_ao); // Thread-local G matrix to void data races
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
                        if (std::abs(this->D(k, l)) < densityThreshold && std::abs(this->D(i, j)) < densityThreshold
                            && std::abs(this->D(i, k)) < densityThreshold && std::abs(this->D(i, l)) < densityThreshold
                            && std::abs(this->D(j, k)) < densityThreshold && std::abs(this->D(j, l)) < densityThreshold)
                            continue;


                        double eri = AtomicOrbital::electronRepulsion(AOs[i], AOs[j], AOs[k], AOs[l]);

                        if (i == j && j == k && k == l) // (ii|ii)
                        {
                            G_private(i, i) += 0.5 * this->D(i, i) * eri;
                        }
                        else if (i == j && j == k) // (ii|il)
                        {
                            G_private(i, i) += this->D(i, l) * eri;
                            G_private(i, l) += 0.5 * this->D(i, i) * eri;
                            G_private(l, i) += 0.5 * this->D(i, i) * eri;
                        }
                        else if (i == j && k == l) // (ii|kk)
                        {
                            G_private(i, i) += this->D(k, k) * eri;
                            G_private(k, k) += this->D(i, i) * eri;
                            G_private(i, k) -= 0.5 * this->D(k, i) * eri;
                            G_private(k, i) -= 0.5 * this->D(k, i) * eri;
                        }
                        else if (i == k && j == l) // (ij|ij)
                        {
                            G_private(i, j) += 1.5 * this->D(i, j) * eri;
                            G_private(j, i) += 1.5 * this->D(i, j) * eri;
                            G_private(j, j) -= 0.5 * this->D(i, i) * eri;
                            G_private(i, i) -= 0.5 * this->D(j, j) * eri;
                        }
                        else if (i == j) // (ii|kl)
                        {
                            G_private(i, i) += 2 * this->D(k, l) * eri;
                            G_private(k, l) += this->D(i, i) * eri;
                            G_private(l, k) += this->D(i, i) * eri;
                            G_private(i, l) -= 0.5 * this->D(i, k) * eri;
                            G_private(l, i) -= 0.5 * this->D(i, k) * eri;
                            G_private(i, k) -= 0.5 * this->D(i, l) * eri;
                            G_private(k, i) -= 0.5 * this->D(i, l) * eri;
                        }
                        else if (k == l) // (ij|kk)
                        {
                            G_private(k, k) += 2 * this->D(i, j) * eri;
                            G_private(i, j) += this->D(k, k) * eri;
                            G_private(j, i) += this->D(k, k) * eri;
                            G_private(i, k) -= 0.5 * this->D(k, j) * eri;
                            G_private(k, i) -= 0.5 * this->D(k, j) * eri;
                            G_private(j, k) -= 0.5 * this->D(i, k) * eri;
                            G_private(k, j) -= 0.5 * this->D(i, k) * eri;
                        }
                        else // (ij|kl)
                        {
                            G_private(i, j) += 2 * this->D(k, l) * eri;
                            G_private(j, i) += 2 * this->D(k, l) * eri;
                            G_private(k, l) += 2 * this->D(i, j) * eri;
                            G_private(l, k) += 2 * this->D(i, j) * eri;
                            G_private(i, l) -= 0.5 * this->D(k, j) * eri;
                            G_private(l, i) -= 0.5 * this->D(k, j) * eri;
                            G_private(j, l) -= 0.5 * this->D(i, k) * eri;
                            G_private(l, j) -= 0.5 * this->D(i, k) * eri;
                            G_private(i, k) -= 0.5 * this->D(j, l) * eri;
                            G_private(k, i) -= 0.5 * this->D(j, l) * eri;
                            G_private(j, k) -= 0.5 * this->D(i, l) * eri;
                            G_private(k, j) -= 0.5 * this->D(i, l) * eri;
                        }
                    }
                }
            }
        }
#pragma omp critical
        {
            G += G_private; // Accumulate the thread-local G matrix into the global G matrix.
        }
    }
    this->F = this->h + G;
}


void SCF::diagonalizeAndUpdate()
{
    // Transform the Fock matrix to the orthogonal basis.
    Eigen::MatrixXd F_prime = this->X.transpose() * this->F * this->X;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(F_prime);
    this->C               = this->X * solver.eigenvectors();
    this->eigenvalues     = solver.eigenvalues();
    Eigen::MatrixXd C_occ = this->C.leftCols(this->occupiedCount);
    this->D               = 2 * C_occ * C_occ.transpose();


    // Calculate the new electronic energy.
    this->electronicEnergy = (this->D.array() * (this->h + this->F).array()).sum() * 0.5;
}


void SCF::printIteration(int iter, double dE, double dD) const
{
    std::cout << "Iteration " << std::setw(3) << iter << ": "
              << " | ΔE = " << std::scientific << dE << " | ΔD = " << dD << std::endl;
}


void SCF::printIteration(int iter, double dE, double dD, double DIISError) const
{
    std::cout << "Iteration " << std::setw(3) << iter << ": "
              << " | ΔE = " << std::scientific << dE << " | ΔD = " << dD << " | DIIS error = " << DIISError << std::endl;
}


void SCF::printFinalResults(bool converged) const
{
    std::stringstream ss;
    ss << std::left;
    if (converged)
        ss << "\nSCF converged!\n";
    else
        ss << "\nWarning: SCF did not converge!\n";

    double totalEnergy = this->electronicEnergy + this->nuclearEnergy;
    ss << "\n--- SCF Results ---\n";
    ss << std::fixed << std::setprecision(10);
    ss << "Electronic Energy: " << this->electronicEnergy << " Hartree\n";
    ss << "Nuclear Repulsion: " << this->nuclearEnergy << " Hartree\n";
    ss << "Total SCF Energy:  " << totalEnergy << " Hartree\n";

    std::vector<std::string> aoLabels = this->molecule.getAOLabels();

    ss << "\n";
    for (size_t i = 0; i < 99; ++i) { ss << "-"; }
    ss << "\n";
    size_t whitespaceWidth = (99 - std::string("Orbital Energies (a.u.)").length()) / 2;
    ss << std::string(whitespaceWidth, ' ') << "Orbital Energies (a.u.)" << std::string(whitespaceWidth, ' ') << "\n";
    for (size_t i = 0; i < 99; ++i) { ss << "-"; }
    ss << "\n";


    ss << "\n-- Occupied --\n";
    ss << printShortMOs(this->eigenvalues.head(this->occupiedCount), 5, 8);
    ss << "-- Virtual --\n";
    ss << printShortMOs(this->eigenvalues.tail(this->basisCount - this->occupiedCount), 5, 8);
    for (size_t i = 0; i < 99; ++i) { ss << "-"; }
    ss << "\n";

    ss << "\n";
    for (size_t i = 0; i < 99; ++i) { ss << "-"; }
    ss << "\n";
    whitespaceWidth = (99 - std::string("Molecular Orbital Coefficients").length()) / 2;
    ss << std::string(whitespaceWidth, ' ') << "Molecular Orbital Coefficients" << std::string(whitespaceWidth, ' ')
       << "\n";
    for (size_t i = 0; i < 99; ++i) { ss << "-"; }
    ss << "\n";

    ss << "Occupied orbitals:\n";
    ss << "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n";
    ss << printFullMOs(this->C.leftCols(this->occupiedCount), this->eigenvalues.head(this->occupiedCount), aoLabels, 5, 5);

    ss << "\nVirtual orbitals:\n";
    ss << "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n";
    ss << printFullMOs(
        this->C.rightCols(this->basisCount - this->occupiedCount),
        this->eigenvalues.tail(this->basisCount - this->occupiedCount),
        aoLabels,
        5,
        5
    );
    for (size_t i = 0; i < 99; ++i) { ss << "-"; }

    std::cout << ss.str() << std::endl;
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
