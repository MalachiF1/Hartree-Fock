#include "Molecule.hpp"

#include "Basis.hpp"
#include "Symmetry/Symmetry.hpp"
#include "Utils.hpp"

#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>

Molecule::Molecule(
    int charge, int multiplicity, const std::string& basisName, const std::vector<Atom>& geometry, double symmetryTolerance
) :
    charge(charge), multiplicity(multiplicity)
{

    auto newGeometry                       = Symmetry::translateCOMToOrigin(geometry);
    auto [principalMoments, principalAxes] = Symmetry::diagonalizeInertiaTensor(newGeometry);
    PointGroup pointGroup = Symmetry::classifyPointGroup(newGeometry, principalMoments, symmetryTolerance);
    Symmetry::alignCanonically(newGeometry, pointGroup, principalMoments, principalAxes, symmetryTolerance);
    this->geometry   = newGeometry;
    this->pointGroup = pointGroup;

    buildBasis(basisName);
    this->electronCount      = countElectrons();
    this->basisFunctionCount = atomicOrbitals.size();
}

std::string Molecule::toString() const
{
    std::stringstream ss;
    ss << "Molecule with charge " << charge << " and multiplicity " << multiplicity << ":\n";
    ss << "Number of electrons: " << electronCount << "\n";
    ss << "Number of basis functions: " << basisFunctionCount << "\n";
    ss << "Geometry (atom, x, y, z):\n";
    ss << std::fixed << std::setprecision(8);
    for (const auto& atom : geometry)
    {
        ss << "\t" << Utils::atomicNumberToName.at(atom.atomicNumber) << "\t" << atom.coords.x() << "\t"
           << atom.coords.y() << "\t" << atom.coords.z() << "\n";
    }
    ss << "Point Group: " << pointGroup.toString() << "\n";
    return ss.str();
}

size_t Molecule::countElectrons() const
{
    size_t totalElectrons = 0;
    for (const auto& atom : geometry) { totalElectrons += atom.atomicNumber; }
    return totalElectrons - charge;
}


void Molecule::buildBasis(const std::string& basisName)
{
    std::vector<int> elements;
    for (const auto& atom : geometry)
    {
        // Ensure that no duplicate atomic numbers are added to the elements vector.
        if (std::find(elements.begin(), elements.end(), atom.atomicNumber) == elements.end())
            elements.push_back(atom.atomicNumber);
    }

    auto basisData = Basis::getBasis(basisName, elements);

    for (const auto& atom : geometry)
    {
        const auto& shellsForAtom = basisData.at(atom.atomicNumber);
        for (const auto& shell : shellsForAtom)
        {
            // For a given angular momentum L, generate all components (e.g., L=1 -> px, py, pz)
            // This loop handles the generation of Cartesian Gaussians.
            for (int i = shell.angularMomentum; i >= 0; --i)
            {
                for (int j = shell.angularMomentum - i; j >= 0; --j)
                {
                    int k = shell.angularMomentum - i - j;

                    std::vector<PrimitiveGaussian> primitives;
                    for (size_t p = 0; p < shell.exponents.size(); ++p)
                    {
                        primitives.emplace_back(
                            shell.exponents[p], shell.coefficients[p], atom.coords, Eigen::Vector3i(i, j, k)
                        );
                    }
                    atomicOrbitals.emplace_back(atom.coords, primitives);
                }
            }
        }
    }
}


double Molecule::nuclearRepulsion() const
{
    // The nuclear repulsion energy is calculated using the formula:
    // E_NN = sum_{i < j} (Z_i * Z_j) / r_ij
    double energy = 0.0;
    for (size_t i = 0; i < geometry.size(); ++i)
    {
        for (size_t j = i + 1; j < geometry.size(); ++j)
        {
            double dist = (geometry[i].coords - geometry[j].coords).norm();
            energy += (geometry[i].atomicNumber * geometry[j].atomicNumber) / dist;
        }
    }
    return energy;
}


Eigen::MatrixXd Molecule::overlapMatrix() const
{
    Eigen::MatrixXd S(basisFunctionCount, basisFunctionCount);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < basisFunctionCount; ++i)
    {
        for (size_t j = i; j < basisFunctionCount; ++j) // only compute upper triangle as S is symmetric
        {
            S(i, j) = S(j, i) = AtomicOrbital::overlap(atomicOrbitals[i], atomicOrbitals[j]);
        }
    }
    return S;
}


Eigen::MatrixXd Molecule::kineticMatrix() const
{
    Eigen::MatrixXd T(basisFunctionCount, basisFunctionCount);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < basisFunctionCount; ++i)
    {
        for (size_t j = i; j < basisFunctionCount; ++j) // only compute upper triangle as T is symmetric
        {
            T(i, j) = T(j, i) = AtomicOrbital::kinetic(atomicOrbitals[i], atomicOrbitals[j]);
        }
    }
    return T;
}


Eigen::MatrixXd Molecule::nuclearAttractionMatrix() const
{
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(basisFunctionCount, basisFunctionCount);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < basisFunctionCount; ++i)
    {
        for (size_t j = i; j < basisFunctionCount; ++j) // only compute upper triangle as V is symmetric
        {
            double v_ij = 0.0;
            for (const auto& atom : geometry)
            {
                v_ij += atom.atomicNumber
                      * AtomicOrbital::nuclearAttraction(atomicOrbitals[i], atomicOrbitals[j], atom.coords);
            }
            V(i, j) = V(j, i) = v_ij;
        }
    }
    return V;
}


ElectronRepulsionTensor Molecule::electronRepulsionTensor(double threshold) const
{
    size_t N_ao = basisFunctionCount;
    ElectronRepulsionTensor Vee(N_ao);

    // Pre-calculate Schwartz screening matrix
    Eigen::MatrixXd Q(N_ao, N_ao);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < N_ao; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            double integral = AtomicOrbital::electronRepulsion(
                atomicOrbitals[i], atomicOrbitals[j], atomicOrbitals[i], atomicOrbitals[j]
            );
            Q(i, j) = Q(j, i) = std::sqrt(integral);

            // We can set these elements in the tensor instead of calculating them again in the next loop.
            // The ElectronRepulsionTensor class handles the symmetry
            Vee(i, j, i, j) = integral;
        }
    }

#pragma omp parallel for collapse(4) schedule(dynamic, 1)
    for (size_t i = 0; i < N_ao; ++i)
    {
        for (size_t j = 0; j < N_ao; ++j)
        {
            for (size_t k = 0; k < N_ao; ++k)
            {
                for (size_t l = 0; l < N_ao; ++l)
                {
                    // no need to recalculate identical elements of the tensor (enforce quartet symmetry)
                    if (j > i || l > k || (i * (i + 1) / 2 + j) < (k * (k + 1) / 2 + l))
                        continue;

                    // make sure we did not already calculate this integral in Schwartz screening loop
                    if ((i == k && j == l) || (i == l && j == k))
                        continue;

                    // apply Schwartz screening at the AO level
                    if (Q(i, j) * Q(k, l) < threshold)
                        continue;

                    double integral = AtomicOrbital::electronRepulsion(
                        atomicOrbitals[i], atomicOrbitals[j], atomicOrbitals[k], atomicOrbitals[l]
                    );

                    // Store with 8-fold symmetry
                    Vee(i, j, k, l) = integral;
                }
            }
        }
    }
    return Vee;
}

Eigen::MatrixXd Molecule::schwartzScreeningMatrix() const
{
    size_t N_ao = basisFunctionCount;

    // Pre-calculate Schwartz screening matrix
    Eigen::MatrixXd Q(N_ao, N_ao);
#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < N_ao; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            double integral = AtomicOrbital::electronRepulsion(
                atomicOrbitals[i], atomicOrbitals[j], atomicOrbitals[i], atomicOrbitals[j]
            );
            Q(i, j) = Q(j, i) = std::sqrt(integral);
        }
    }
    return Q;
}
