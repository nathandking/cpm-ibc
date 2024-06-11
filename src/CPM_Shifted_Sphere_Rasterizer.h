#pragma once

#include <nova/Dynamics/Hierarchy/Rasterizers/Hierarchical_Rasterizer.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/Tools/Log/Log.h>

namespace cpm
{
template<class Struct_type,class T,int d>
class CPM_Shifted_Sphere_Rasterizer: public Hierarchical_Rasterizer<Struct_type,T,d>
{
    using Base                  = Hierarchical_Rasterizer<Struct_type, T, d>;
    using TV                    = Vector<double, d>;
    using T_INDEX               = Vector<int,d>;
    using T_NODE                = std::pair<unsigned,T_INDEX>;
    using Hierarchy             = Grid_Hierarchy<Struct_type,T,d>;

  public:
    using Base::hierarchy;
    TV min_corner;
    double tube_radius, dx;
    double R;

    CPM_Shifted_Sphere_Rasterizer(Hierarchy& hierarchy_input, double R_input, double tube_radius_input, double dx_input, TV min_corner_input)
        :Base(hierarchy_input), R(R_input), tube_radius(tube_radius_input), dx(dx_input), min_corner(min_corner_input)
        {}

    bool Consume(const T_NODE& node) override
    {
        const unsigned level = node.first;
        const T_INDEX& index = node.second;

        TV Xd;
        for(int axis = 0; axis < 3; ++axis)
        {
            Xd(axis) = min_corner(axis) + dx * (index(axis) - 1);
        }
        double dist = Xd.Norm() - R;

        if(level==0 && fabs(dist) <= tube_radius)
        {
            hierarchy.Activate_Cell(level, index, Node_Active);
            return false;
        }
        return (level > 0);
    }
};
}
