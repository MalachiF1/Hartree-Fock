#include "Molecule.hpp"

#include "Eigen/Core"
#include "Symmetry/Symmetry.hpp"
#include "Utils.hpp"

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <ranges>

Molecule::Molecule(int charge, int multiplicity, const std::vector<Atom>& geometry, bool detectSymmetry, double symmetryTolerance) :
    charge(charge), multiplicity(multiplicity), symmetryTolerance(symmetryTolerance)
{

    if (detectSymmetry)
    {
        this->geometry = geometry;
        Symmetry::translateOrigin(this->geometry, Symmetry::findCOM(this->geometry));

        auto SEAs                   = Symmetry::findSEAs(this->geometry, symmetryTolerance);
        auto [symbol, paxis, saxis] = Symmetry::findPointGroup(this->geometry, symmetryTolerance);
        Symmetry::rotateToNewAxes(this->geometry, saxis, paxis.cross(saxis), paxis);
        this->pointGroupName = symbol;
    }
    else
    {
        this->geometry = geometry;
    }

    this->electronCount = countElectrons();
}

std::string Molecule::toString() const
{
    return fmt::format(
        "Molecule with charge {} and multiplicity {}:\nNumber of electrons: {}\nNumber of basis functions: "
        "{}\nGeometry (atom, x, y, z):\n{}",
        charge,
        multiplicity,
        electronCount,
        basisFunctionCount,
        fmt::join(
            geometry
                | std::views::transform(
                    [](const Atom& atom)
                    {
                        return fmt::format(
                            "\t{:<2} {:>3} {:>13.8f} {:>13.8f} {:>13.8f}\n",
                            Utils::atomicNumberToName.at(atom.atomicNumber),
                            atom.atomicNumber,
                            atom.coords.x(),
                            atom.coords.y(),
                            atom.coords.z()
                        );
                    }
                ),
            ""
        )
    );
}


size_t Molecule::countElectrons() const
{
    size_t totalElectrons = 0;
    for (const auto& atom : geometry) { totalElectrons += atom.atomicNumber; }
    return totalElectrons - charge;
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
