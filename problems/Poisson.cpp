#include "Poisson.h"

namespace cpm {

Poisson::Poisson(SurfaceSpecifier &surface_specs, Scalar dx, size_t bdy_order, vector<size_t> &which_bdy)
{
    m_surface_tube.ConstructTube(surface_specs, 3 /* interpolation order */, dx);
    Initialize(bdy_order, which_bdy);
}


Poisson::Poisson(SurfaceSpecifier &surface_specs, Scalar dx)
{
    m_surface_tube.ConstructTube(surface_specs, 3 /* interpolation order */, dx);

    vector<size_t> which_bdy_dummy;
    Initialize(0, which_bdy_dummy);
}


VectorX Poisson::GetRHS(function<Scalar(const vector<Scalar>&)> f, function<Scalar(const vector<Scalar>&)> exact)
{
#ifdef TRACK_WHERE
    std::cout << "Poisson::GetRHS(0)" << std::endl;
#endif
    VectorX rhs(m_surface_tube.nNodes());
    for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
    {
        rhs[i] = f(m_surface_tube.cpx()[i]);
    }

    if(m_which_bdy.size() > 0)
    {    
        if(m_bdy_order == 1)
        {
            for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
            {
                for(size_t bdy_idx = 0; bdy_idx < m_which_bdy.size(); ++bdy_idx)
                {
                    if(m_surface_tube.bdy()[i] == m_which_bdy[bdy_idx])
                    {
                        rhs[i] = exact(m_surface_tube.cpx()[i]);
                    }
                }
            }
        }
        else if(m_bdy_order == 2) 
        {
            if(m_surface_tube.dim() == 3)
            {
                cout << "Poisson::Solve - Second order boundary conditions not implemented for 3D yet" << endl;
            }
            else if(m_surface_tube.dim() == 2)
            {
            }
        }
    }
    return rhs;
}

VectorX Poisson::GetRHS(function<Scalar(const Scalar&, const Scalar&)> f, function<Scalar(const Scalar&, const Scalar&)> exact)
{
#ifdef TRACK_WHERE
    std::cout << "Poisson::GetRHS(1)" << std::endl;
#endif
    VectorX rhs(m_surface_tube.nNodes());
    for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
    {
        Scalar theta = atan2(m_surface_tube.cpx()[i][1], m_surface_tube.cpx()[i][0]);
        if(m_surface_tube.dim() == 2)
        {
            rhs[i] = f(theta, 0);
        }
        else if(m_surface_tube.dim() == 3)
        {
            Scalar phi = acos(m_surface_tube.cpx()[i][2] / m_surface_tube.surface().surfaceParams()[0]);
            rhs[i] = f(theta, phi);
        }
    }

    if(m_which_bdy.size() > 0)
    {    
        if(m_bdy_order == 1)
        {
            for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
            {
                for(size_t bdy_idx = 0; bdy_idx < m_which_bdy.size(); ++bdy_idx)
                {
                    if(m_surface_tube.bdy()[i] == m_which_bdy[bdy_idx])
                    {
                        Scalar theta = atan2(m_surface_tube.cpx()[i][1], m_surface_tube.cpx()[i][0]);
                        if(m_surface_tube.dim() == 2)
                        {
                            rhs[i] = exact(theta, 0);
                        }
                        else if(m_surface_tube.dim() == 3)
                        {
                            Scalar phi = acos(m_surface_tube.cpx()[i][2] / m_surface_tube.surface().surfaceParams()[0]);
                            rhs[i] = exact(theta, phi);
                        }
                    }
                }
            }
        }
        else if(m_bdy_order == 2) 
        {
            if(m_surface_tube.dim() == 3)
            {
                cout << "Poisson::Solve - Second order boundary conditions not implemented for 3D yet" << endl;
            }
            else if(m_surface_tube.dim() == 2)
            {
                VectorX rhs_full(m_surface_tube.nNodes() + 2);
                for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
                {
                    rhs_full[i] = rhs[i];
                }
                rhs_full[m_surface_tube.nNodes()] = exact(m_surface_tube.surface().surfaceParams()[1], 0);
                rhs_full[m_surface_tube.nNodes() + 1] = exact(m_surface_tube.surface().surfaceParams()[2], 0);
                return rhs_full;
            }
        }
    }
    return rhs;
}

VectorX Poisson::GetRHS(function<Scalar(const Scalar&, const Scalar&)> f)
{
#ifdef TRACK_WHERE
    std::cout << "Poisson::GetRHS(2)" << std::endl;
#endif
    VectorX rhs(m_surface_tube.nNodes());
    for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
    {
        Scalar theta = atan2(m_surface_tube.cpx()[i][1], m_surface_tube.cpx()[i][0]);
        if(m_surface_tube.dim() == 2)
        {
            rhs[i] = f(theta, 0);
        }
        else if(m_surface_tube.dim() == 3)
        {
            Scalar phi = acos(m_surface_tube.cpx()[i][2] / m_surface_tube.surface().surfaceParams()[0]);
            rhs[i] = f(theta, phi);
        }
    }

    return rhs;
}

SpMat Poisson::GetPlotInterpMatrix()
{
    SpMat plotE;
    m_h.GetPlotInterpMatrix(m_surface_tube, plotE);

    return plotE;
}


////////////////////////////////////////////////////////////////////////////////////
// Private Functions
////////////////////////////////////////////////////////////////////////////////////


void Poisson::Initialize(size_t bdy_order, vector<size_t> &which_bdy)
{
    m_bdy_order = bdy_order;
    m_which_bdy = which_bdy;

    if(m_bdy_order == 2)
    {
        m_h.ReplaceClosestPointWithcpBar(m_surface_tube, m_surface_tube.bdy());
    }
    
    // construct the discrete Laplacian and cp extension matrices
    m_h.SetupMatrices(m_surface_tube, m_E, m_L);

    // negate interpolation weights for second order Dirichlet BCs
    if(m_bdy_order == 2 && m_which_bdy.size() > 0)
    {
        for (size_t k=0; k < m_E.outerSize(); ++k) 
        {
            for (SpMat::InnerIterator it(m_E, k); it; ++it)
            {
                for(size_t bdy_idx = 0; bdy_idx < m_which_bdy.size(); ++bdy_idx)
                {
                    if(m_surface_tube.bdy()[it.row()] == m_which_bdy[bdy_idx])
                    {
                        it.valueRef() *= -1.0;
                    }
                }
            }
        }
    }

    m_identity_rows = vector<bool>(m_surface_tube.nNodes(), false);
    if(m_which_bdy.size() > 0)
    {
        if(m_bdy_order == 1)
        {
            for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
            {
                for(size_t bdy_idx = 0; bdy_idx < m_which_bdy.size(); ++bdy_idx)
                {
                    if(m_surface_tube.bdy()[i] == m_which_bdy[bdy_idx])
                    {
                        m_identity_rows[i] = true;
                    }
                }
            }
            m_h.SetIdentityRows(m_E, m_identity_rows);
            m_Efull = m_E;
            m_Lfull = m_L;
        }
        else if(m_bdy_order == 2)
        {
            m_bdy = vector<size_t>(m_surface_tube.nNodes(), 0);
            for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
            {
                for(size_t bdy_idx = 0; bdy_idx < m_which_bdy.size(); ++bdy_idx)
                {
                    if(m_surface_tube.bdy()[i] == m_which_bdy[bdy_idx]) // Dirichlet BC on the first boundary
                    {
                        m_bdy[i] = m_which_bdy[bdy_idx];
                    }
                }
            }
            m_h.LaplacianSharpWithDirichletConditions2ndOrder(m_L, m_E, m_Lfull, m_Efull, m_bdy);
            m_identity_rows = vector<bool>(m_Lfull.rows(), false);
            m_identity_rows[m_identity_rows.size() - 1] = 1;
            m_identity_rows[m_identity_rows.size() - 2] = 1;
        }
    }
    else
    {
        m_Lfull = m_L;
        m_Efull = m_E;
    }
}
} // namespace cpm