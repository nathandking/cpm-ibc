#include "Heat.h"

namespace cpm {

Heat::Heat(SurfaceSpecifier &surface_specs, Scalar dx, Scalar final_time, bool use_implicit, size_t bdy_order, vector<size_t> &which_bdy)
{
    m_surface_tube.ConstructTube(surface_specs, 3 /* interpolation order */, dx);
    Initialize(final_time, bdy_order, which_bdy, use_implicit);
}


Heat::Heat(SurfaceSpecifier &surface_specs, Scalar dx, Scalar final_time, bool use_implicit)
{
    m_surface_tube.ConstructTube(surface_specs, 3 /* interpolation order */, dx);

    vector<size_t> which_bdy_dummy;
    Initialize(final_time, 0, which_bdy_dummy, use_implicit);
}


// for open surfaces
void Heat::TimeStep(VectorX &u, function<Scalar(const Scalar&, const Scalar&)> exact)
{
    for(size_t t = 0; t < m_num_time_steps; ++t)
    {
        if(!m_use_implicit) // explicit Euler
        {
            u = u + m_dt * (m_L * u);
            u = m_E * u; // closest point extension

            if(m_which_bdy.size() > 0)
            {
                DirichletBC(t+1, u, exact);
            }
        }
        else
        {
            if(m_bdy_order == 2)
            {
                VectorX ufull(m_surface_tube.nNodes() + 2);
                for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
                {
                    ufull[i] = u[i];
                }
                ufull[m_surface_tube.nNodes()] = exact(t * m_dt, m_surface_tube.surface().surfaceParams()[1]);
                ufull[m_surface_tube.nNodes() + 1] = exact(t * m_dt, m_surface_tube.surface().surfaceParams()[2]);
                
                cpm::VectorX res = sv.Solve(ufull);

                for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
                {
                    u[i] = res[i];
                }
            }
            else if(m_bdy_order == 1)
            {
                if(m_which_bdy.size() > 0)
                {
                    DirichletBC(t, u, exact);
                }
                
                VectorX ufull(m_surface_tube.nNodes());
                for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
                {
                    ufull[i] = u[i];
                }

                cpm::VectorX res = sv.Solve(ufull);
                u = res;
            }
        }
    }
}


// for closed surfaces
void Heat::TimeStep(VectorX &u)
{
    for(size_t t = 0; t < m_num_time_steps; ++t)
    {
        if(!m_use_implicit) // forward Euler time stepping
        {
            u = u + m_dt * (m_L * u);
            u = m_E * u; // closest point extension
        }
        else // implicit Euler time stepping
        {
            cpm::VectorX ufull = sv.Solve(u);
            u = ufull;
        }
    }
}


SpMat Heat::GetPlotInterpMatrix()
{
    SpMat plotE;
    m_h.GetPlotInterpMatrix(m_surface_tube, plotE);

    return plotE;
}


////////////////////////////////////////////////////////////////////////////////////
// Private Functions
////////////////////////////////////////////////////////////////////////////////////


void Heat::Initialize(Scalar final_time, size_t bdy_order, vector<size_t> &which_bdy, bool use_implicit)
{
    m_final_time = final_time;
    m_use_implicit = use_implicit;
    m_bdy_order = bdy_order;
    m_which_bdy = which_bdy;
    nD = 0;

    // for time stepping heat flow
    m_dt = 0.2 * m_surface_tube.dx() * m_surface_tube.dx(); // gives second-order with both explicit and implicit Euler, as expected
    m_num_time_steps = ceil(m_final_time / m_dt);
    m_dt = m_final_time / m_num_time_steps; // adjust dt so we have an integer number of time steps
    
    if(m_bdy_order == 2 && m_which_bdy.size() > 0)
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
            for (SpMat::InnerIterator it(m_E,k); it; ++it)
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

    if(m_use_implicit) // use implicit Euler
    {
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
                m_Lfull = m_L;
                m_Efull = m_E;
            }
            else if(m_bdy_order == 2)
            {
                vector<size_t> bdy(m_surface_tube.nNodes(), 0);
                for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
                {
                    for(size_t bdy_idx = 0; bdy_idx < m_which_bdy.size(); ++bdy_idx)
                    {
                        if(m_surface_tube.bdy()[i] == m_which_bdy[bdy_idx]) // Dirichlet BC on the first boundary
                        {
                            bdy[i] = m_which_bdy[bdy_idx];
                        }
                    }
                }
                m_h.LaplacianSharpWithDirichletConditions2ndOrder(m_L, m_E, m_Lfull, m_Efull, bdy);
                m_identity_rows = vector<bool>(m_Lfull.rows(), false);
                nD = 2;
            }
        }
        else
        {
            m_Efull = m_E;
            m_Lfull = m_L;
            m_identity_rows = vector<bool>(m_Lfull.rows(), false);
        }
        
        sv.Init((Scalar)1., -m_dt, m_Lfull, m_Efull, Preconditioner::diagonal, m_identity_rows, nD);
    }
}


void Heat::DirichletBC(Scalar t, VectorX &u, function<Scalar(const Scalar&, const Scalar&)> exact)
{
    for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
    {
        for(size_t bdy_idx = 0; bdy_idx < m_which_bdy.size(); ++bdy_idx)
        {
            if(m_surface_tube.bdy()[i] == m_which_bdy[bdy_idx])
            {
                if(m_bdy_order == 1)
                {
                    u[i] = exact(t * m_dt, m_surface_tube.surface().surfaceParams()[m_surface_tube.bdy()[i]]);
                }
                else if (m_bdy_order == 2)
                {
                    u[i] += 2.0 * exact(t * m_dt, m_surface_tube.surface().surfaceParams()[m_surface_tube.bdy()[i]]);
                }
            }
        }
    }
}


} // namespace cpm