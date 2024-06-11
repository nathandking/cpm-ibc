#pragma once

#include "Defines.h"
#include "Scalar.h"
#include "Surface.h"

#include <nova/Dynamics/Hierarchy/Rasterizers/Hierarchical_Rasterizer.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
#include <nova/Tools/Log/Log.h>

namespace cpm
{
template<class Struct_type,class T,int d>
class CPM_Rasterizer: public Hierarchical_Rasterizer<Struct_type,T,d>
{
    using Base                  = Hierarchical_Rasterizer<Struct_type, T, d>;
    using TV                    = Vector<Scalar, d>;
    using T_INDEX               = Vector<int, d>;
    using T_NODE                = std::pair<unsigned, T_INDEX>;
    using Hierarchy             = Grid_Hierarchy<Struct_type,T,d>;

  public:
    using Base::hierarchy;
    Surface& surface;
    TV min_corner;
    Scalar tube_radius, dx;

    CPM_Rasterizer(Hierarchy& hierarchy_input, Surface& surface_input, Scalar tube_radius_input, Scalar dx_input, TV min_corner_input)
        :Base(hierarchy_input), surface(surface_input), tube_radius(tube_radius_input), dx(dx_input), min_corner(min_corner_input)
        {}

    bool Consume(const T_NODE& node) override
    {
        const unsigned level = node.first;
        const T_INDEX& index = node.second;

        std::vector<Scalar> Xd(d, (Scalar)0.);
        for(int axis = 0; axis < d; ++axis)
        {
            Xd[axis] = min_corner(axis) + dx * (index(axis) - 1);
        }

        Scalar dist;
        std::vector<Scalar> cpx(d, (Scalar)0.);
        size_t bdy;
        surface.ClosestPoint(Xd, cpx, dist, bdy);

        if(level==0 && fabs(dist) <= tube_radius)
        {
            hierarchy.Activate_Cell(level, index, Node_Active);
            return false;
        }
        return (level > 0);
    }
};
}
