#pragma once

#include <vector>
#include <iostream>
#include "Scalar.h"

namespace cpm {

struct IBCMetaData
{
    bool is_oriented = false; // required
    std::vector<Scalar> point_on_side1; // required if is_oriented = true

    std::vector<size_t> boundary_type{0, 0}; // required if solving PDE, default is Dirichlet-Dirichlet interior BC
    size_t boundary_order = 1; // required if solving PDE
    
    Scalar on_ibc_tolerance; // required if solving PDE

    bool use_cp_diff_directions = false; // use cp_diff as directions or cp_diff with normals and tangential components projected out
    bool use_Jcp_surface_normals = true;
    bool visualize = false;

    Scalar ibc_DOFSubset_radius_scale = 1.0;
};

} // namespace cpm