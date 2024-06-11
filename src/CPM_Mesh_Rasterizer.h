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
class CPM_Mesh_Rasterizer: public Hierarchical_Rasterizer<Struct_type,T,d>
{
    using Base                  = Hierarchical_Rasterizer<Struct_type, T, d>;
    using TV                    = Vector<Scalar, d>;
    using T_INDEX               = Vector<int, d>;
    using T_NODE                = std::pair<unsigned, T_INDEX>;
    using Hierarchy             = Grid_Hierarchy<Struct_type,T,d>;

  public:
    using Base::hierarchy;
    const Surface& surface;
    TV min_corner;
    Scalar tube_radius, dx;

    CPM_Mesh_Rasterizer(Hierarchy& hierarchy_input, Surface& surface_input, Scalar tube_radius_input, Scalar dx_input, TV min_corner_input)
        :Base(hierarchy_input), surface(surface_input), tube_radius(tube_radius_input), dx(dx_input), min_corner(min_corner_input)
        {}

    bool Consume(const T_NODE& node) override
    {
        const unsigned level = node.first;
        const T_INDEX& index = node.second;

        std::vector<Scalar> Xd(3, (Scalar)0.);
        for(int axis = 0; axis < 3; ++axis)
        {
            Xd[axis] = min_corner(axis) + dx * (index(axis) - 1);
        }
        
        Scalar dist = surface.GetDistanceToMesh(Xd);

        if(level==0 && fabs(dist) <= tube_radius)
        {
            hierarchy.Activate_Cell(level, index, Node_Active);
            return false;
        }
        return (level > 0);
    }
};
}
