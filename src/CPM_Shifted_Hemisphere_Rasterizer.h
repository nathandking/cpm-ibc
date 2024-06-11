#pragma once

#include <nova/Dynamics/Hierarchy/Rasterizers/Hierarchical_Rasterizer.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/Tools/Log/Log.h>

namespace cpm
{
template<class Struct_type,class T,int d>
class CPM_Shifted_Hemisphere_Rasterizer: public Hierarchical_Rasterizer<Struct_type,T,d>
{
    using Base                  = Hierarchical_Rasterizer<Struct_type, T, d>;
    using Td                    = double;
    using TV                    = Vector<Td, d>;
    using T_INDEX               = Vector<int,d>;
    using T_NODE                = std::pair<unsigned,T_INDEX>;
    using Hierarchy             = Grid_Hierarchy<Struct_type,T,d>;

  public:
    using Base::hierarchy;
    TV min_corner;
    double tube_radius, dx;
    double R;

    CPM_Shifted_Hemisphere_Rasterizer(Hierarchy& hierarchy_input, double R_input, double tube_radius_input, double dx_input, TV min_corner_input)
        :Base(hierarchy_input), R(R_input), tube_radius(tube_radius_input), dx(dx_input), min_corner(min_corner_input)
        {}

    bool Consume(const T_NODE& node) override
    {
        const unsigned level = node.first;
        const T_INDEX& index = node.second;

        TV X;
        for(int axis = 0; axis < 3; ++axis)
        {
            X(axis) = min_corner(axis) + dx * (index(axis) - 1);
        }

        TV cpx;
        Td theta = std::atan2(X(1), X(0));
        Td dist = (Td)0.;
        if(X(2) > (Td)0.) 
        {
            dist = std::abs(X.Norm() - R);
        }
        else
        {
            cpx(0) = R * std::cos(theta);
            cpx(1) = R * std::sin(theta);
            cpx(2) = (Td)0.;

            dist = (X - cpx).Norm();
        }

        if(level==0 && dist <= tube_radius)
        {
            hierarchy.Activate_Cell(level, index, Node_Active);
            return false;
        }
        return (level > 0);
    }
};
}
