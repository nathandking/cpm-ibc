#include "IBCDOFSubset.h" 

namespace cpm {

IBCDOFSubset::IBCDOFSubset()
{

}

void IBCDOFSubset::ComputeDOFSubset(Surface &ibc_surface, IBCSubset &ibc_subset, Tube &surface_tube)
{
#ifdef TRACK_WHERE
    std::cout << "IBCDOFSubset::ComputeDOFSubset" << std::endl;
#endif
    m_ibc_meta = ibc_subset.meta();
    
    // compute the smaller subset of just the new DOFs for the constriant
    for(size_t i = 0; i < ibc_subset.size(); ++i)
    {
        if(ibc_subset.dist()[i] < m_ibc_meta.ibc_DOFSubset_radius_scale * surface_tube.TubeRadius())
        {
            m_set_idx_from_subset_idx.push_back(ibc_subset.SetIndexFromSubsetIndex()[i]);

            m_cpx_c.push_back(ibc_subset.cpx()[i]);
            m_cp_diff.push_back(ibc_subset.cp_diff()[i]);
            m_dist_c.push_back(ibc_subset.dist()[i]);
            m_bdy.push_back(ibc_subset.bdy()[i]);
        }
    }
    ComputeInverseMapping(surface_tube.nNodes(), m_set_idx_from_subset_idx, m_subset_idx_from_set_idx);

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
        for(size_t i = 0; i < ibc_subset.size(); ++i)
        {
            if(ibc_subset.dist()[i] < m_ibc_meta.ibc_DOFSubset_radius_scale * surface_tube.TubeRadius())
            {
                m_surface_normals.push_back(ibc_subset.surface_normals()[i]);
                if(surface_tube.dim() == 3)
                {
                    m_tangents.push_back(ibc_subset.tangents()[i]);
                }
                m_binormals.push_back(ibc_subset.binormals()[i]);
                m_directions.push_back(ibc_subset.directions()[i]);
            }
        }
    }

    m_which_side.assign(size(), 0);
    if(m_ibc_meta.is_oriented)
    {   
        OrientSubset(ibc_surface, surface_tube);

        for(size_t i = 0; i < size(); ++i)
        {
            assert(m_which_side[i] != 0);
        }
    }

    if( (m_ibc_meta.boundary_type[0] == 0 || m_ibc_meta.boundary_type[1] == 0) && m_ibc_meta.boundary_order == 2) // only for 2nd order Dirichlet
    {
        for(size_t i = 0; i < size(); ++i)
        {
            const size_t bc_type_idx = m_ibc_meta.is_oriented ? m_which_side[i] - 1 : 0;
            if(m_ibc_meta.boundary_type[bc_type_idx] == 0)
            {
                m_dirichlet_set_idx_from_subset_idx.push_back(i);
            }
        }
        ComputeInverseMapping(size(), m_dirichlet_set_idx_from_subset_idx, m_dirichlet_subset_idx_from_set_idx);
    }

    if(m_ibc_meta.boundary_order == 2)
    {
        m_cp_c_bar.resize(size(), vector<Scalar>(surface_tube.dim()));
        for(size_t i = 0; i < size(); ++i)
        {
            vector<Scalar> projection = 2.0 * cpx()[i] - surface_tube.x()[SetIndexFromSubsetIndex()[i]];

            surface_tube.ClosestPoint(projection, m_cp_c_bar[i]);
        }
    }

    ComputeExtensionMatrix(surface_tube);
}


void IBCDOFSubset::ComputeInverseMapping(size_t set_size, const vector<size_t> &set_from_subset, vector<ssize_t> &subset_from_set)
{
#ifdef TRACK_WHERE
    std::cout << "IBCDOFSubset::ComputeInverseMapping" << std::endl;
#endif
    // compute the inverse mapping from an index in the full set to the subset index
    subset_from_set.assign(set_size, -1); // return -1 if the index in the full set is not in the subset
    for(size_t i = 0; i < set_from_subset.size(); ++i)
    {
        subset_from_set[set_from_subset[i]] = i;
    }
}


void IBCDOFSubset::ComputeExtensionMatrix(Tube& surface_tube)
{
    vector<vector<Scalar>> xquery;
    if(m_ibc_meta.boundary_order == 1)
    {
        xquery = cpx();
    }
    else if(m_ibc_meta.boundary_order == 2)
    {
        xquery = m_cp_c_bar;
    }

    Interpolation interp(surface_tube, xquery); 
    m_E.resize(size(), surface_tube.nNodes());
    interp.BuildInterpolationMatrix(surface_tube, m_E);
}


void IBCDOFSubset::OrientSubset(Surface& ibc_surface, Tube &surface_tube)
{
    vector<Scalar> point_on_side1 = m_ibc_meta.point_on_side1;
    
    // find starting point that is not on Sperp, i.e. cp_diff > 0
    size_t start_idx;
    for(size_t i = 0; i < size(); ++i)
    {
        if(Norm(m_directions[i]) >= m_ibc_meta.on_ibc_tolerance) // TO DO: should this be >= m_ibc_meta.on_ibc_tolerance?
        {
            start_idx = i;
            break;
        }
    }

    m_which_side[start_idx] = 1;
    list<size_t> queue;
    queue.push_back(start_idx);
    vector<bool> is_visited(size(), false);

    while(!queue.empty())
    {
        if(is_visited[queue.front()] == false)
        {
            size_t set_idx = SetIndexFromSubsetIndex()[queue.front()];
            is_visited[queue.front()] = true;

            for(size_t nbr = 0; nbr < 2 * surface_tube.dim(); ++nbr)
            {
                if(surface_tube.TubeNode(set_idx).neighbour(nbr) != nullptr)
                {
                    ssize_t nbr_subset_idx = SubsetIndexFromSetIndex()[surface_tube.TubeNode(set_idx).neighbour(nbr)->TubeIndex()];
                        
                    if(nbr_subset_idx >= 0)
                    {
                        if(m_which_side[nbr_subset_idx] == 0)
                        {
                            if(DotProduct(m_directions[queue.front()], m_directions[nbr_subset_idx]) < 0)
                            {
                                if(m_which_side[queue.front()] == 1)
                                {
                                    m_which_side[nbr_subset_idx] = 2;
                                }
                                else
                                {
                                    m_which_side[nbr_subset_idx] = 1;
                                }
                            }
                            else
                            {
                                m_which_side[nbr_subset_idx] = m_which_side[queue.front()];
                            }
                        }

                        queue.push_back(nbr_subset_idx);
                    }
                }
            }
        }
        queue.pop_front();
    }

    if(point_on_side1.size() == surface_tube.dim())
    {
        // make sure the which_side tag is such that points near point_on_side1 have opposite tag of which_side = 2
        size_t nearest_subset_idx;
        Scalar nearest_dist = numeric_limits<Scalar>::max(); // find nearest subset index
        for(size_t i = 0; i < size(); ++i)
        {
            Scalar dist = Distance(point_on_side1, surface_tube.x()[SetIndexFromSubsetIndex()[i]]); 
            
            vector<Scalar> cpx_c(surface_tube.dim());
            
            ibc_surface.ClosestPoint(point_on_side1, cpx_c);
            vector<Scalar> cp_diff = point_on_side1 - cpx_c;

            if(dist < nearest_dist && Norm(m_cp_diff[i]) >= 2 * surface_tube.dx() && DotProduct(cp_diff, m_cp_diff[i]) > 0)
            {
                nearest_dist = dist;
                nearest_subset_idx = i; 
            }
        }


        if(m_which_side[nearest_subset_idx] == 1)
        {
            for(size_t i = 0; i < m_which_side.size(); ++i)
            {
                if(m_which_side[i] == 1)
                {
                    m_which_side[i] = 2;
                }
                else if(m_which_side[i] == 2)
                {
                    m_which_side[i] = 1;
                }
            }
        }
    }
}


} // namespace cpm