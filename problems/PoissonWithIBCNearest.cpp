#include "PoissonWithIBCNearest.h"

namespace cpm {
    
PoissonWithIBCNearest::PoissonWithIBCNearest(SurfaceSpecifier &surface_specs, vector<SurfaceSpecifier> &ibc_specs, Scalar dx, bool freeze_nearest, size_t interp_deg)
{
    m_freeze_nearest = freeze_nearest;
    m_surface_tube.ConstructTube(surface_specs, interp_deg, dx);

    m_ibc.Initialize(m_surface_tube, ibc_specs);

    m_ibc_idx.resize(m_ibc.NumIBCs());
    if(m_freeze_nearest)
    {
        for(size_t c = 0; c < m_ibc.NumIBCs(); ++c)
        {
            m_ibc_idx[c].resize(m_ibc.surface(c).xp().size());
            for(size_t i = 0; i < m_ibc.surface(c).xp().size(); ++i)
            {
                m_ibc_idx[c][i] = m_surface_tube.NearestTubeIndex(m_ibc.xp(c)[i]);
            }
        }
    }
    else
    {
        for(size_t c = 0; c < m_ibc.NumIBCs(); ++c)
        {
            m_ibc_idx[c] = m_ibc.DOFSubset()[c].SetIndexFromSubsetIndex();
        }
    }

    Initialize();
}


VectorX PoissonWithIBCNearest::GetRHS(function<Scalar(const vector<Scalar>&)> f, function<Scalar(const vector<Scalar>&)> surface_exact)
{
    VectorX rhs(m_surface_tube.nNodes());
    for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
    {
        rhs[i] = f(m_surface_tube.cpx()[i]);
    }

    for(size_t c = 0; c < m_ibc.NumIBCs(); ++c)
    {
        for(size_t i = 0; i < m_ibc_idx[c].size(); ++i)
        {
            if(m_freeze_nearest)
            {
                rhs[m_ibc_idx[c][i]] = surface_exact(m_ibc.xp(c)[i]);
            }
            else
            {
                rhs[m_ibc_idx[c][i]] = surface_exact(m_ibc.DOFSubset()[c].cpx()[i]);
            }
        }
    }

    return rhs;
}

SpMat PoissonWithIBCNearest::GetPlotInterpMatrix()
{
    SpMat plotE;
    m_h.GetPlotInterpMatrix(m_surface_tube, plotE);

    return plotE;
}


////////////////////////////////////////////////////////////////////////////////////
// Private Functions
////////////////////////////////////////////////////////////////////////////////////


void PoissonWithIBCNearest::Initialize()
{   
    // construct the discrete Laplacian and cp extension matrices
    m_h.SetupMatrices(m_surface_tube, m_E, m_L);

    m_identity_rows.assign(m_surface_tube.nNodes(), false);
    for(size_t c = 0; c < m_ibc_idx.size(); ++c)
    {
        for(size_t i = 0; i < m_ibc_idx[c].size(); ++i)
        {
            m_identity_rows[m_ibc_idx[c][i]] = true;
        }
    }
    m_h.SetIdentityRows(m_E, m_identity_rows);
}
} // namespace cpm