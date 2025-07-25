#include "Molecule.hpp"
#include "SCF.hpp"
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

    Molecule propane(
        0,
        1,
        "6-31g",
        {
            {6, Vec3(-10.7589077582, 6.2565685212, -2.5241907065)},
            {6, Vec3(-8.6512753667, 7.4733819239, -0.9589621263)},
            {1, Vec3(-10.1054011389, 5.8121286194, -4.4169695293)},
            {1, Vec3(-11.4014728952, 4.5108119925, -1.6601017898)},
            {1, Vec3(-12.3764854323, 7.5036136768, -2.6960181430)},
            {6, Vec3(-6.3364934421, 5.7599542046, -0.6783909163)},
            {1, Vec3(-8.0875192406, 9.2432998819, -1.8334766919)},
            {1, Vec3(-9.3740209847, 7.9533035327, 0.9018145467)},
            {1, Vec3(-5.5367252480, 5.2996231945, -2.5102398991)},
            {1, Vec3(-4.8715538317, 6.6605988574, 0.4369528611)},
            {1, Vec3(-6.8344536616, 3.9976985952, 0.2455653944)},
        }
    );

    // SCF scf_h2o(H2O);
    // scf_h2o.run(50, 1e-6, 1e-6, 1e-6, true, 8, 1, 1e-6);

    SCF scf_propane(propane);
    scf_propane.run(50, 1e-6, 1e-6, 1e-6, true, 8, 1, 1e-6);
    // SCF scf_h2(H2);
    // scf_h2.run(50, 1e-10, 1e-10);
    // SCF scf_hf(HF);
    // scf_hf.run(50, 1e-10, 1e-10);

    return 0;
}
