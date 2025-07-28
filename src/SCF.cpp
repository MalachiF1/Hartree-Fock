#include "SCF.hpp"

#include "DIIS.hpp"
#include "Utils.hpp"

#include <iomanip>
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
    double DIISErrorTol
)
{
    // Initialize the SCF calculation. (This should be called only once.)
    this->initialize(schwartzThreshold);

    // Get the initial guess density.
    this->computeInitialGuessDensity();

    if (useDIIS)
    {
        // Initialize the DIIS handler with the specified maximum size.
        diis_handler = std::make_unique<DIIS>(DIISmaxSize);
    }

    std::cout << "\n--- Starting SCF Iterations ---\n" << std::endl;

    // SCF cycles
    for (size_t iteration = 0; iteration < maxIter; ++iteration)
    {
        // Store the energy and density from the previous iteration to check for convergence.
        double lastElectronicEnergy = this->electronicEnergy;
        Eigen::MatrixXd D_last      = this->D;

        // Build the new Fock matrix from the current density.
        this->buildFockMatrix();

        // Store the error vector and fock matrix for DIIS if applicable.
        if (useDIIS)
        {
            diis_handler->update(this->F, this->D, this->S, this->X);
        }

        // Diagonalize the Fock matrix and compute the new density and coefficients.
        if (useDIIS && iteration >= DIISstart)
        {
            // Extrapolate the orthogonal Fock matrix using DIIS.
            Eigen::MatrixXd F_prime = diis_handler->extrapolate(this->F);

            this->diagonalizeAndUpdate(F_prime);
        }
        else
        {
            this->diagonalizeAndUpdate();
        }

        // Calculate the change in energy and density.
        double dE = std::abs(this->electronicEnergy - lastElectronicEnergy);
        double dD = (this->D - D_last).norm();

        // Print the current iteration results.
        if (useDIIS)
        {
            this->printIteration(iteration + 1, dE, dD, diis_handler->getErrorNorm());
        }
        else
        {
            this->printIteration(iteration + 1, dE, dD);
        }

        // Check for convergence.
        if (dE < energyTol && dD < densityTol && (!useDIIS || diis_handler->getErrorNorm() < DIISErrorTol))
        {
            // If converged, print the final results and exit.
            this->printFinalResults(true);
            return;
        }
    }

    // If the loop finishes, the calculation did not converge.
    this->printFinalResults(false);
}

void SCF::initialize(double schwartzThreshold)
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

    // Calculate two-electron integrals.
    this->Vee = molecule.electronRepulsionTensor(schwartzThreshold);

    // Calculate the overlap and orthogonalization matrices.
    this->S = molecule.overlapMatrix();
    this->X = inverseSqrtMatrix(this->S);
}


void SCF::computeInitialGuessDensity()
{
    // The initial Fock matrix is just the core Hamiltonian.
    this->F = this->h;

    // Diagonalize the initial Fock matrix to get guess orbitals.
    Eigen::MatrixXd F_prime = this->X.transpose() * this->F * this->X;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(F_prime);
    this->C = this->X * solver.eigenvectors();

    // Form the density matrix from the occupied orbitals.
    Eigen::MatrixXd C_occ = this->C.leftCols(this->occupiedCount);
    this->D               = 2 * C_occ * C_occ.transpose();

    // Calculate the initial electronic energy.
    this->electronicEnergy = (this->D.array() * (this->h + this->F).array()).sum() * 0.5;
}


void SCF::buildFockMatrix()
{
    size_t N_ao = this->basisCount;

    // Calculate the G matrix (the two-electron part of the Fock matrix).
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(N_ao, N_ao);
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

    // The full Fock matrix is the core Hamiltonian plus the G matrix.
    this->F = this->h + G;
}


void SCF::diagonalizeAndUpdate()
{
    // Transform the Fock matrix to the orthogonal basis.
    Eigen::MatrixXd F_prime = this->X.transpose() * this->F * this->X;

    // Find the eigenvalues and eigenvectors.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(F_prime);

    // Transform the eigenvectors back to the original basis.
    this->C = this->X * solver.eigenvectors();

    // Form the new density matrix from the occupied orbitals.
    Eigen::MatrixXd C_occ = this->C.leftCols(this->occupiedCount);
    this->D               = 2 * C_occ * C_occ.transpose();

    // Calculate the new electronic energy.
    this->electronicEnergy = (this->D.array() * (this->h + this->F).array()).sum() * 0.5;
}


void SCF::diagonalizeAndUpdate(const Eigen::MatrixXd& F_prime)
{
    // Find the eigenvalues and eigenvectors.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(F_prime);

    // Transform the eigenvectors back to the original basis.
    this->C = this->X * solver.eigenvectors();

    // Form the new density matrix from the occupied orbitals.
    Eigen::MatrixXd C_occ = this->C.leftCols(this->occupiedCount);
    this->D               = 2 * C_occ * C_occ.transpose();

    // Calculate the new electronic energy.
    this->electronicEnergy = (this->D.array() * (this->h + this->F).array()).sum() * 0.5;
}


void SCF::printIteration(int iter, double dE, double dD) const
{
    std::cout << "Iteration " << std::setw(3) << iter << ": "
              << " | dE = " << std::scientific << dE << " | dD = " << dD << std::endl;
}


void SCF::printIteration(int iter, double dE, double dD, double DIISError) const
{
    std::cout << "Iteration " << std::setw(3) << iter << ": "
              << " | dE = " << std::scientific << dE << " | dD = " << dD << " | DIIS error = " << DIISError << std::endl;
}


void SCF::printFinalResults(bool converged) const
{
    if (converged)
    {
        std::cout << "\nSCF converged!" << std::endl;
    }
    else
    {
        std::cout << "\nWarning: SCF did not converge!" << std::endl;
    }

    double totalEnergy = this->electronicEnergy + this->nuclearEnergy;
    std::cout << "\n--- SCF Results ---" << std::endl;
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Electronic Energy: " << this->electronicEnergy << " Hartree" << std::endl;
    std::cout << "Nuclear Repulsion: " << this->nuclearEnergy << " Hartree" << std::endl;
    std::cout << "Total SCF Energy:  " << totalEnergy << " Hartree" << std::endl;
}
