#include "SCF.hpp"

#include "Utils.hpp"

#include <iomanip>
#include <iostream>

SCF::SCF(const Molecule& molecule) : molecule(molecule) {}

void SCF::run(size_t maxIter, double energyTol, double densityTol, double schwartzThreshold)
{
    // Initialize the SCF calculation. (This should be called only once.)
    this->initialize(schwartzThreshold);

    // Get the initial guess density.
    this->computeInitialGuessDensity();

    std::cout << "\n--- Starting SCF Iterations ---\n" << std::endl;

    // SCF cycles
    for (size_t iteration = 0; iteration < maxIter; ++iteration)
    {
        // Store the energy and density from the previous iteration to check for convergence.
        double lastElectronicEnergy = this->electronicEnergy;
        Eigen::MatrixXd D_last      = this->D;

        // Build the new Fock matrix from the current density.
        this->buildFockMatrix();

        // Diagonalize the Fock matrix and compute the new density and coefficients.
        this->diagonalizeAndUpdate();

        // Calculate the change in energy and density.
        double dE = std::abs(this->electronicEnergy - lastElectronicEnergy);
        double dD = (this->D - D_last).norm();

        this->printIteration(iteration + 1, dE, dD);

        // Check for convergence.
        if (dE < energyTol && dD < densityTol)
        {
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
    std::cout << "Number of electrons: " << molecule.getElectronCount() << std::endl;
    std::cout << "Number of basis functions: " << molecule.getBasisFunctionCount() << std::endl;
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

// Computes the initial guess for the density matrix using the core Hamiltonian.
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

// Builds the Fock matrix G from the two-electron integrals and the density matrix.
void SCF::buildFockMatrix()
{
    size_t N_ao       = this->basisCount;
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(N_ao, N_ao);

    for (size_t i = 0; i < N_ao; ++i)
    {
        for (size_t j = 0; j < N_ao; ++j)
        {
            for (size_t k = 0; k < N_ao; ++k)
            {
                for (size_t l = 0; l < N_ao; ++l)
                {
                    double J_ij = this->Vee[(i * N_ao * N_ao * N_ao) + (j * N_ao * N_ao) + (k * N_ao) + l];
                    double K_ik = this->Vee[(i * N_ao * N_ao * N_ao) + (k * N_ao * N_ao) + (j * N_ao) + l];
                    G(i, j) += this->D(k, l) * (J_ij - 0.5 * K_ik);
                }
            }
        }
    }

    // The full Fock matrix is the core Hamiltonian plus the G matrix.
    this->F = this->h + G;
}

// Solves the eigenvalue problem F'C' = E'C' and computes the new density matrix.
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

// Prints the status for a single SCF iteration.
void SCF::printIteration(int iter, double dE, double dD) const
{
    std::cout << "Iteration " << std::setw(3) << iter << ": "
              << " | dE = " << std::scientific << dE << " | dD = " << dD << std::endl;
}

// Prints the final results of the SCF calculation.
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
