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
class CPM_Rasterizer_Test: public Hierarchical_Rasterizer<Struct_type,T,d>
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
    std::vector<std::vector<int>> index_offsets;

    CPM_Rasterizer_Test(Hierarchy& hierarchy_input, Surface& surface_input, Scalar tube_radius_input, Scalar dx_input, TV min_corner_input, std::vector<std::vector<int>>& index_offsets_input)
        :Base(hierarchy_input), surface(surface_input), tube_radius(tube_radius_input), dx(dx_input), min_corner(min_corner_input), index_offsets(index_offsets_input)
        {}

    bool Consume(const T_NODE& node) override
    {
        const unsigned level = node.first;
        const T_INDEX& index = node.second;

        const Scalar cur_dx = dx * std::pow((Scalar)2., level);
        std::vector<Scalar> Xd(d, (Scalar)0.);
        for(int axis = 0; axis < d; ++axis)
        {
            Xd[axis] = min_corner(axis) + cur_dx * (index(axis) - 1);
        }
        if(level > 0)
        {
            // check dist for its 4/8 neighbors
            // only when they all fall within the range will the children nodes be further considered
            bool to_discard = true;
            for(int id = 0; id < index_offsets.size(); ++id)
            {
                std::vector<Scalar> cur_Xd = Xd;
                for(int axis=0;axis<d;++axis)
                {
                    cur_Xd[axis] += (Scalar).5 * cur_dx * (index_offsets[id][axis]==0?(Scalar)-1:(Scalar)1.);
                }
                Scalar dist;
                std::vector<Scalar> cpx(d, (Scalar)0.);
                size_t bdy;
                surface.ClosestPoint(cur_Xd, cpx, dist, bdy);

                if(fabs(dist) <= tube_radius)
                {
                    to_discard = false;
                    break;
                }
            }
            if(to_discard)
            {
                // std::cout << "skip" << std::endl;
            }
            else 
            {
                std::cout << "reserve" << std::endl;
            }
            return !to_discard;
        }
        else 
        {
            Scalar dist;
            std::vector<Scalar> cpx(d, (Scalar)0.);
            size_t bdy;
            surface.ClosestPoint(Xd, cpx, dist, bdy);
            if(fabs(dist) <= tube_radius)
            {
                hierarchy.Activate_Cell(level, index, Node_Active);
            }
            return false;
        }
    }
};
}
