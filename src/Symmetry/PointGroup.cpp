#include "Symmetry/PointGroup.hpp"

PointGroup::PointGroup() : name(PointGroupName::C1), n(1), order(1), operations({SymmetryOperation::Identity()}) {}

PointGroup::PointGroup(PointGroupName name, size_t n, const std::vector<SymmetryOperation>& operations) :
    name(name), n(n), operations(operations)
{
    // make sure identiy operation is present
    bool hasIdentity = false;
    for (const auto& op : operations)
    {
        if (op.name == SymmetryOperation::E)
        {
            hasIdentity = true;
            break;
        }
    }
    if (!hasIdentity)
    {
        this->operations.emplace_back(SymmetryOperation::E, Eigen::Matrix3d::Identity(), Vec3(0, 0, 0));
    }

    order = operations.size();
}

std::string PointGroup::toString() const
{
    std::string n_str = std::to_string(n);
    switch (name)
    {
        case PointGroupName::Dinfh: return "D∞h";
        case PointGroupName::Cinfh: return "C∞h";
        case PointGroupName::I: return "I";
        case PointGroupName::Ih: return "Ih";
        case PointGroupName::O: return "O";
        case PointGroupName::Oh: return "Oh";
        case PointGroupName::T: return "T";
        case PointGroupName::Td: return "Td";
        case PointGroupName::Th: return "Th";
        case PointGroupName::Dn: return "D" + n_str;
        case PointGroupName::Dnh: return "D" + n_str + "h";
        case PointGroupName::Dnd: return "D" + n_str + "d";
        case PointGroupName::Cn: return "C" + n_str;
        case PointGroupName::Cnh: return "C" + n_str + "h";
        case PointGroupName::Cnv: return "C" + n_str + "v";
        case PointGroupName::Cs: return "Cs";
        case PointGroupName::Ci: return "Ci";
        case PointGroupName::S2n: return "S" + std::to_string(this->n * 2);
        case PointGroupName::C1: return "C1";
        default: return "Unknown point group";
    }
}
