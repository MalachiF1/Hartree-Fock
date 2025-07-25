#include "Molecule.hpp"
#include "Utils.hpp"

#include <iostream>

int main()
{
    // Water molecule (H2O) example with STO-3G basis set
    Molecule H2O(
        0,        // charge
        1,        // multiplicity
        "STO-3G", // basis set name
        {
            {1, Vec3(0.0, 1.43233673, -0.96104039)},  // H1
            {1, Vec3(0.0, -1.43233673, -0.96104039)}, // H2
            {8, Vec3(0.0, 0.0, 0.24026010)}           // O
        }
    );

    Molecule H2(
        0,        // charge
        1,        // multiplicity
        "STO-3G", // basis set name
        {
            {1, Vec3(0.0, 0.0, 0.0)}, // H1
            {1, Vec3(0.0, 0.0, 1.4)}  // H2
        }
    );

    Molecule HF(
        0,        // charge
        1,        // multiplicity
        "STO-3G", // basis set name
        {
            {1, Vec3(0.0, 0.0, -0.1527810236)}, // H
            {9, Vec3(0.0, 0.0, 1.6527810236)}   // F
        }
    );

    auto S   = H2O.overlapMatrix();
    auto T   = H2O.kineticMatrix();
    auto Vne = H2O.nuclearAttractionMatrix();
    auto Vee = H2O.electronRepulsionTensor();

    std::cout << "H2O: " << H2O.toString() << "\n";
    std::cout << "S:\n" << S << "\n\n";
    std::cout << "T:\n" << T << "\n\n";
    std::cout << "Vne:\n" << Vne << "\n\n";
    std::cout << "Vee (first 10 elements):\n";
    for (size_t i = 0; i < 10 && i < Vee.size(); ++i) { std::cout << Vee[i] << " "; }

    return 0;
}
