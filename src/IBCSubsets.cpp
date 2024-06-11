#include "IBCSubsets.h" 


namespace cpm {

IBCSubsets::IBCSubsets()
{

}

void IBCSubsets::Initialize(Tube &surface_tube, vector<SurfaceSpecifier> ibc_surface_specs)
{
#ifdef TRACK_WHERE
    std::cout << "IBCSubsets::Initialize(0)" << std::endl;
#endif
    vector<Scalar> bounding_box; //NOTE: bounding box is not set, this may cause issues if used improperly later
    m_num_ibcs = ibc_surface_specs.size();
    for(size_t c = 0; c < m_num_ibcs; ++c)
    {
        Surface temp_surface;
        temp_surface.SetSurfaceSpecs(ibc_surface_specs[c]);
        m_ibc_surface.push_back(std::move(temp_surface));
    }

    Initialize(surface_tube);
}


void IBCSubsets::Initialize(Tube &surface_tube)
{
#ifdef TRACK_WHERE
    std::cout << "IBCSubsets::Initialize(2)" << std::endl;
#endif
    ComputeSubsets(surface_tube);

    IBCStartingColumns(surface_tube);

    TotalDOFs();

    for(size_t c = 0; c < m_num_ibcs; ++c)
    {
        if(meta(c).visualize)
        {
            if(surface_tube.dim() == 3)
            {
                VisualizeSperp(surface_tube, c);
            }
            VisualizeSubsets(surface_tube, c);
        }
    }
}


void IBCSubsets::ComputeSubsets(Tube& surface_tube)
{
#ifdef TRACK_WHERE
    std::cout << "IBCSubsets::ComputeSubsets" << std::endl;
#endif
    for(size_t c = 0; c < m_num_ibcs; ++c)
    {
        IBCSubset temp_subset;
        temp_subset.ComputeSubset(m_ibc_surface[c], surface_tube);
        m_ibc_subset.push_back(std::move(temp_subset));
    }
    for(size_t c = 0; c < m_num_ibcs; ++c)
    {
        IBCDOFSubset temp_subset;
        temp_subset.ComputeDOFSubset(m_ibc_surface[c], m_ibc_subset[c], surface_tube);
        m_ibc_DOF_subset.push_back(std::move(temp_subset));
    }
}


void IBCSubsets::IBCStartingColumns(Tube& surface_tube)
{
#ifdef TRACK_WHERE
    std::cout << "IBCSubsets::IBCStartingColumns" << std::endl;
#endif
    m_ibc_starting_column.resize(m_num_ibcs);
    m_dirichlet_ibc_starting_column.resize(m_num_ibcs);
    if(m_ibc_starting_column.size() > 0)
    {
        m_ibc_starting_column[0] = surface_tube.nNodes();
        m_dirichlet_ibc_starting_column[0] = surface_tube.nNodes();
    }
    for(size_t c = 1; c < m_num_ibcs; ++c)
    {
        m_ibc_starting_column[c] = m_ibc_starting_column[c-1] + m_ibc_DOF_subset[c-1].size();
        if( (meta(c).boundary_type[0] == 0 || meta(c).boundary_type[1] == 0) && meta(c).boundary_order == 2)
        {
            m_dirichlet_ibc_starting_column[c] = m_dirichlet_ibc_starting_column[c-1] + m_ibc_DOF_subset[c-1].DirichletSetIndexFromSubsetIndex().size();
        }
    }
}


void IBCSubsets::TotalDOFs()
{
    m_identity_rows_start = 0;
    m_num_2nd_order_dirichlet_rows = 0;
    m_total_ibc_DOFs = 0;
    for(size_t c = 0; c < m_num_ibcs; ++c)
    {
        if((meta(c).boundary_type[0] == 0 || meta(c).boundary_type[1] == 0) && meta(c).boundary_order == 2) // if either side of BC is 2nd order Dirichlet 
        {
            m_identity_rows_start += 1;
            m_num_2nd_order_dirichlet_rows += m_ibc_DOF_subset[c].DirichletSetIndexFromSubsetIndex().size();
        }
        m_total_ibc_DOFs += m_ibc_DOF_subset[c].size();
    }

    if(m_identity_rows_start != 0)
    {
        m_identity_rows_start = m_total_ibc_DOFs;
    }
}


void IBCSubsets::VisualizeSperp(Tube &surface_tube, size_t ibc_index)
{
#ifdef POLYSCOPE
    glm::vec3 green{0.4660, 0.6740, 0.1880};

    Geometry geom;
    vector<vector<Scalar>> surface_normals = geom.UnorientedSurfaceNormals(surface_tube);

    // interpolate surface normals onto ibc polyline
    vector<vector<Scalar>> surface_normals_on_ibc_polyline = geom.InterpolateUnorientedVectors(surface_tube, surface_normals, m_ibc_surface[ibc_index].xp(), true /*normalize*/);

    // flip all normals to point in the same direction
    for(size_t i = 1; i < m_ibc_surface[ibc_index].xp().size(); ++i)
    {
        if(DotProduct(surface_normals_on_ibc_polyline[i], surface_normals_on_ibc_polyline[i-1]) < 0)
        {
            surface_normals_on_ibc_polyline[i] *= -1;
        }
    }

    polyscope::registerCurveNetworkLine("IBC Curve " + to_string(ibc_index), m_ibc_surface[ibc_index].xp())->setColor(glm::vec3(1,1,1))->addNodeVectorQuantity("Surface Normals", surface_normals_on_ibc_polyline);       

    // create normal surface faces
    vector<vector<size_t>> normal_surface_faces(2 * (m_ibc_surface[ibc_index].xp().size() - 1), vector<size_t>(4));
    vector<vector<Scalar>> normal_surface_vertices(6 * (m_ibc_surface[ibc_index].xp().size() - 1), vector<Scalar>(surface_tube.dim()));
    for(size_t i = 1; i < m_ibc_surface[ibc_index].xp().size(); ++i)
    {
        normal_surface_faces[2 * (i - 1)][0] = 6 * (i - 1);
        normal_surface_faces[2 * (i - 1)][1] = 6 * (i - 1) + 1;
        normal_surface_faces[2 * (i - 1)][2] = 6 * (i - 1) + 4;
        normal_surface_faces[2 * (i - 1)][3] = 6 * (i - 1) + 2;
        normal_surface_faces[2 * (i - 1) + 1][0] = 6 * (i - 1);
        normal_surface_faces[2 * (i - 1) + 1][1] = 6 * (i - 1) + 1;
        normal_surface_faces[2 * (i - 1) + 1][2] = 6 * (i - 1) + 5;
        normal_surface_faces[2 * (i - 1) + 1][3] = 6 * (i - 1) + 3;
        
        normal_surface_vertices[6 * (i - 1)] = m_ibc_surface[ibc_index].xp()[i-1];
        normal_surface_vertices[6 * (i - 1) + 1] = m_ibc_surface[ibc_index].xp()[i];
        normal_surface_vertices[6 * (i - 1) + 2] = normal_surface_vertices[6 * (i - 1)] + surface_tube.TubeRadius() * surface_normals_on_ibc_polyline[i-1];
        normal_surface_vertices[6 * (i - 1) + 3] = normal_surface_vertices[6 * (i - 1)] - surface_tube.TubeRadius() * surface_normals_on_ibc_polyline[i-1];
        normal_surface_vertices[6 * (i - 1) + 4] = normal_surface_vertices[6 * (i - 1) + 1] + surface_tube.TubeRadius() * surface_normals_on_ibc_polyline[i];
        normal_surface_vertices[6 * (i - 1) + 5] = normal_surface_vertices[6 * (i - 1) + 1] - surface_tube.TubeRadius() * surface_normals_on_ibc_polyline[i];
    }

    polyscope::registerSurfaceMesh("Normal Surface IBC " + to_string(ibc_index), normal_surface_vertices, normal_surface_faces)->setSurfaceColor(green);
#endif
}


void IBCSubsets::VisualizeSubsets(Tube &surface_tube, size_t ibc_index)
{
    VisualizeSubset(surface_tube, TaggingSubset()[ibc_index], "Tagging Subset ", ibc_index);

    VisualizeSubset(surface_tube, DOFSubset()[ibc_index], "DOF Subset ", ibc_index);
    VisualizeTubeSubsetAddScalarQuantity(DOFSubset()[ibc_index].which_side(), "DOF Subset " + to_string(ibc_index), "which_side");
}


template <typename T>
void IBCSubsets::VisualizeSubset(Tube &surface_tube, const T &subset, string subset_name, size_t ibc_index)
{
    VisualizeTubeSubset(surface_tube, subset.SetIndexFromSubsetIndex(), subset.cpx(), subset_name + to_string(ibc_index));
    VisualizeTubeSubsetAddVectorQuantity(subset.cp_diff(), subset_name + to_string(ibc_index), "cp_diff");
    VisualizeTubeSubsetAddVectorQuantity(subset.directions(), subset_name + to_string(ibc_index), "directions");

    if(!meta(ibc_index).use_cp_diff_directions)
    {
        VisualizeTubeSubsetAddVectorQuantity(subset.surface_normals(), subset_name + to_string(ibc_index), "surface normals");
        if(surface_tube.dim() == 3)
        {
            VisualizeTubeSubsetAddVectorQuantity(subset.tangents(), subset_name + to_string(ibc_index), "tangents");
        }
        VisualizeTubeSubsetAddVectorQuantity(subset.binormals(), subset_name + to_string(ibc_index), "binormals");
    }
}


template void IBCSubsets::VisualizeSubset<IBCSubset>(Tube &surface_tube, const IBCSubset &subset, string subset_name, size_t ibc_index);
template void IBCSubsets::VisualizeSubset<IBCDOFSubset>(Tube &surface_tube, const IBCDOFSubset &subset, string subset_name, size_t ibc_index);


} // namespace cpm