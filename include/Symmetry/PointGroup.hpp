#pragma once

#include "SymmetryOperation.hpp"

#include <cstddef>
#include <string>

/**
 * @brief Represents a point group classification of a molecule.
 */
struct PointGroup
{
    enum PointGroupName
    {
        Dinfh,
        Cinfh,
        I,
        Ih,
        O,
        Oh,
        T,
        Td,
        Th,
        Dn,
        Dnh,
        Dnd,
        Cn,
        Cnh,
        Cnv,
        Cs,
        Ci,
        S2n,
        C1,
    };

    PointGroup();
    PointGroup(PointGroupName name, size_t n, const std::vector<SymmetryOperation>& operations);

    std::string toString() const;

    PointGroupName name;
    size_t n; // only relavent for Cn, Cnv, Cnh, Dn, Dnd, Dnh, S2n point groups, 1 otherwise.
    size_t order;
    std::vector<SymmetryOperation> operations;
};
