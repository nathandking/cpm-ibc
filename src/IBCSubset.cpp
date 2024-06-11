#include "IBCSubset.h" 

// #define IBC_TESTS

namespace cpm {

IBCSubset::IBCSubset()
{

}

void IBCSubset::ComputeSubset(Surface& ibc_surface, Tube &surface_tube)
{
#ifdef TRACK_WHERE
    std::cout << "IBCSubset::ComputeSubset" << std::endl;
#endif
    m_ibc_meta = ibc_surface.SurfaceSpecs().IBCMeta();

    ComputeSubsetNearIBC(ibc_surface, surface_tube);   
    ComputeInverseMapping(surface_tube.nNodes(), m_set_idx_from_subset_idx, m_subset_idx_from_set_idx);

    ComputeCPdiff(ibc_surface, surface_tube);

    if(m_ibc_meta.use_cp_diff_directions)
    {
        m_directions.resize(size(), vector<Scalar>(surface_tube.dim()));
        for(size_t i = 0; i < size(); ++i)
        {
            m_directions[i] = m_cp_diff[i];
        }
    }
    else
    {
        ComputeFrameBasedDirections(surface_tube);
    }
}


void IBCSubset::ComputeSubsetNearIBC(Surface& ibc_surface, Tube& surface_tube)
{
#ifdef TRACK_WHERE
    std::cout << "IBCSubset::ComputeSubsetNearIBC" << std::endl;
#endif
    Scalar radius = 2.0 * surface_tube.TubeRadius();

    size_t start_idx = FindStartingPoint(ibc_surface, surface_tube);

    vector<Scalar> cpx_c(surface_tube.dim());
    Scalar dist_c;
    size_t bdy;

    // do a breadth-first traversal to construct the tube
    // grow outward by moving to neighbouring grid nodes within the tube
    list<size_t> queue;
    queue.push_back(start_idx);
    vector<bool> is_visited(surface_tube.nNodes(), false);

    while(!queue.empty())
    {
        for(size_t nbr = 0; nbr < 2 * surface_tube.dim(); ++nbr)
        {
            if(surface_tube.TubeNode(queue.front()).neighbour(nbr) != nullptr)
            {
                size_t nbr_tube_idx = surface_tube.TubeNode(queue.front()).neighbour(nbr)->TubeIndex();
                
                if(!is_visited[nbr_tube_idx])
                {
                    is_visited[nbr_tube_idx] = true;
                    ibc_surface.ClosestPoint(surface_tube.x()[nbr_tube_idx], cpx_c, dist_c, bdy);

                    if(dist_c < radius)
                    {
                        m_cpx_c.push_back(cpx_c);
                        m_dist_c.push_back(dist_c);
                        m_bdy.push_back(bdy);

                        m_set_idx_from_subset_idx.push_back(nbr_tube_idx);
                        queue.push_back(nbr_tube_idx);
                    }
                }
            }
        }
        queue.pop_front();
    }

#ifdef IBC_TESTS
    // test if we have all the points
    vector<size_t> brute_idx = ComputeSubsetNearIBCBruteForce(ibc_surface, surface_tube, radius);

    if(m_set_idx_from_subset_idx.size() != brute_idx.size())
    {
        cout << "The breadth first search and brute force approaches give different points" << endl;
        cout << "breadth first search size = " << m_set_idx_from_subset_idx.size() << endl;
        cout << "brute force size = " << brute_idx.size() << endl;
    }
    else
    {
        bool is_equal = isSortedVectorsEqual(m_set_idx_from_subset_idx, brute_idx);
        if(is_equal == false)
        {
            cout << "We are missing some points or have too many points" << endl;
        }
        else
        {
            cout << "The sets of grid points match for both ibc tube constructions :)" << endl;
        }
    }
#endif
}


size_t IBCSubset::FindStartingPoint(Surface& ibc_surface, Tube &surface_tube)
{
    vector<Scalar> point_on_ibc = ibc_surface.xp()[0];
    return surface_tube.NearestTubeIndex(point_on_ibc);
}


void IBCSubset::ComputeInverseMapping(size_t set_size, const vector<size_t> &set_from_subset, vector<ssize_t> &subset_from_set)
{
    // compute the inverse mapping from an index in the full set to the subset index
    subset_from_set.assign(set_size, -1); // return -1 if the index in the full set is not in the subset
    for(size_t i = 0; i < set_from_subset.size(); ++i)
    {
        subset_from_set[set_from_subset[i]] = i;
    }
}


vector<size_t> IBCSubset::ComputeSubsetNearIBCBruteForce(Surface& ibc_surface, Tube &surface_tube)
{
    Scalar radius = 2.0 * surface_tube.TubeRadius();

    vector<size_t> set_idx_from_subset_idx;

    // compute the indices in the full set that are within the radius distance from the ibc
    vector<Scalar> cpx_c(surface_tube.dim());
    Scalar dist_c;
    size_t bdy;

    for(size_t i = 0; i < surface_tube.nNodes(); ++i)
    {
        // compute the closest point cpx_c on the ibc from the point x and check the distance is less than the radius
        ibc_surface.ClosestPoint(surface_tube.x()[i], cpx_c, dist_c, bdy);

        if(dist_c <= radius)
        {
            set_idx_from_subset_idx.push_back(i);
        }
    }

    return set_idx_from_subset_idx;
}


void IBCSubset::ComputeCPdiff(Surface &ibc_surface, Tube &surface_tube)
{
    m_cp_diff.resize(size(), vector<Scalar>(surface_tube.dim()));
    for(size_t i = 0; i < size(); ++i)
    {
        if(m_bdy[i] > 0 && meta().use_cp_diff_directions)
        {
            vector<Scalar> projection = 2.0 * cpx()[i] - surface_tube.x()[SetIndexFromSubsetIndex()[i]];

            vector<Scalar> cpBar_s(surface_tube.dim()); // not quite cpBar_s, it is actually cpBar_c below, but here it is cp_s of the projection point onto the end of the ibc
            vector<Scalar> cpBar_c(surface_tube.dim());
            
            surface_tube.ClosestPoint(projection, cpBar_s);
            ibc_surface.ClosestPoint(projection, cpBar_c);
            
            m_cp_diff[i] = cpBar_c - cpBar_s; // -ve cp_diff of projection point
        }
        else
        {
            m_cp_diff[i] = surface_tube.cpx()[SetIndexFromSubsetIndex()[i]] - cpx()[i]; 
        }
    }

    // if point is on Sperp, we get a cp_diff by interpolating from its neighbours. The direction points to the side of the first vector used in the interpolation. This controls which side the point belongs to
    for(size_t i = 0; i < size(); ++i)
    {
        if(Norm(m_cp_diff[i]) < m_ibc_meta.on_ibc_tolerance)
        {
            m_geom.InterpolateUnorientedVectorFromNeighbours(surface_tube, m_cp_diff, SetIndexFromSubsetIndex(), SubsetIndexFromSetIndex(), i, m_ibc_meta.on_ibc_tolerance);
        }
    }
}


void IBCSubset::ComputeFrameBasedDirections(Tube &surface_tube)
{
    if(meta().use_Jcp_surface_normals) 
    {
        if(surface_tube.dim() == 2)
        {
            vector<Matrix2> eigenvectors;
            vector<Vector2> eigenvalues;
            m_geom.EigenDecompositionOfJcp(surface_tube, eigenvectors, eigenvalues);

            m_surface_normals_on_tube.resize(surface_tube.nNodes(), vector<Scalar>(surface_tube.dim()));
            for(size_t i = 0; i < surface_tube.nNodes(); ++i)
            {
                for(size_t d = 0; d < surface_tube.dim(); ++d)
                {
                    m_surface_normals_on_tube[i][d] = eigenvectors[i](d, 0);
                }
            }
        }
        else
        {
            vector<Matrix3> eigenvectors;
            vector<Vector3> eigenvalues;
            m_geom.EigenDecompositionOfJcp(surface_tube, eigenvectors, eigenvalues);

            m_surface_normals_on_tube.resize(surface_tube.nNodes(), vector<Scalar>(surface_tube.dim()));
            for(size_t i = 0; i < surface_tube.nNodes(); ++i)
            {
                for(size_t d = 0; d < surface_tube.dim(); ++d)
                {
                    m_surface_normals_on_tube[i][d] = eigenvectors[i](d, 0);
                }
            }
        }
    }
    else
    {
        m_surface_normals_on_tube = m_geom.UnorientedSurfaceNormals(surface_tube);
    }
    m_surface_normals = m_geom.InterpolateUnorientedVectors(surface_tube, m_surface_normals_on_tube, cpx());    

    m_binormals.resize(size(), vector<Scalar>(surface_tube.dim()));
    m_directions.resize(size(), vector<Scalar>(surface_tube.dim()));
    if(surface_tube.dim() == 2)
    {
        // "binormal" in 2D is the tangent and is a 90 degree rotation of the surface normal vector
        for(size_t i = 0; i < size(); ++i)
        {
            m_binormals[i][0] = m_surface_normals[i][1];
            m_binormals[i][1] = -m_surface_normals[i][0];
            m_directions[i] = DotProduct(m_binormals[i], m_cp_diff[i]) * m_binormals[i];
        }
    }
    else if(surface_tube.dim() == 3)
    {
        vector<Matrix3> eigenvectors;
        vector<Vector3> eigenvalues;
        m_geom.EigenDecompositionOfJcp(surface_tube, cpx(), SetIndexFromSubsetIndex(), SubsetIndexFromSetIndex(), eigenvectors, eigenvalues);

        m_tangents.resize(size(), vector<Scalar>(surface_tube.dim(), 0.0));
        for(size_t i = 0; i < size(); ++i)
        {
            if(abs(eigenvalues[i][2]) < 1e-8) // the ibc is a point, there is no defined tangent direction, J_cp_C is the zero matrix
            {
                m_binormals[i] = vector<Scalar>(surface_tube.dim(), 0.0);
                m_directions[i] = m_cp_diff[i] - DotProduct(m_surface_normals[i], m_cp_diff[i]) * m_surface_normals[i];
            }
            else
            {
                for(size_t d = 0; d < surface_tube.dim(); ++d)
                {
                    m_tangents[i][d] = eigenvectors[i](d, 2);
                }
                m_binormals[i] = CrossProduct(m_tangents[i], m_surface_normals[i]);
                m_directions[i] = DotProduct(m_binormals[i], m_cp_diff[i]) * m_binormals[i];
            }
        }
    }
}


} // namespace cpm