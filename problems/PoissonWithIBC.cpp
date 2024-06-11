#include <chrono>
#include "PoissonWithIBC.h"

// #define TIME
using namespace std;
namespace cpm {
    
PoissonWithIBC::PoissonWithIBC(SurfaceSpecifier &surface_specs, vector<SurfaceSpecifier> &ibc_specs, Scalar dx, size_t interp_deg)
{
#ifdef TRACK_WHERE
    std::cout << "PoissonWithIBC::PoissonWithIBC(0)" << std::endl;
#endif
#ifdef TIME
    m_start = chrono::system_clock::now();
#endif

    m_surface_tube.ConstructTube(surface_specs, interp_deg, dx);
#ifdef TIME
    m_end = chrono::system_clock::now();
    m_elapsed_seconds = m_end - m_start;
    cout << "ConstructTube time: " << m_elapsed_seconds.count() << "s\n";
#endif

#ifdef TIME
    m_start = chrono::system_clock::now();
#endif

    m_ibc.Initialize(m_surface_tube, ibc_specs);
    
#ifdef TIME    
    m_end = chrono::system_clock::now();
    m_elapsed_seconds = m_end - m_start;
    cout << "IBCSubsets Initialize time: " << m_elapsed_seconds.count() << "s\n";
#endif

#ifdef TIME
    m_start = chrono::system_clock::now();
#endif
    
    Initialize();

#ifdef TIME
    m_end = chrono::system_clock::now();
    m_elapsed_seconds = m_end - m_start;
    cout << "PoissonWithIBC Initialize time: " << m_elapsed_seconds.count() << "s\n";
#endif
}


VectorX PoissonWithIBC::GetRHS(function<Scalar(const vector<Scalar>&)> f, function<Scalar(const vector<Scalar>&)> surface_exact)
{
    VectorX rhsfull(m_surface_tube.nNodes() + m_ibc.TotalIBCDOFs() + m_ibc.Num2ndOrderDirichletRows());
    for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
    {
        rhsfull[i] = f(m_surface_tube.cpx()[i]);
    }

    for(size_t c = 0; c < m_ibc.NumIBCs(); ++c)
    {
        for(size_t i = 0; i < m_ibc.DOFSubset()[c].size(); ++i)
        {
            const size_t bc_type_idx = m_ibc.meta(c).is_oriented ? m_ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
            if(m_ibc.meta(c).boundary_type[bc_type_idx] == 1 || (m_ibc.meta(c).boundary_type[bc_type_idx] == 0 && m_ibc.meta(c).boundary_order == 2))
            {
                rhsfull[i + m_ibc.IBCColumnStartIndex(c)] = f(m_ibc.DOFSubset()[c].cpx()[i]);
            }
        }
    }

    for(size_t c = 0; c < m_ibc.NumIBCs(); ++c)
    {
        for(size_t i = 0; i < m_ibc.DOFSubset()[c].size(); ++i)
        {
            const size_t bc_type_idx = m_ibc.meta(c).is_oriented ? m_ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
            if(m_ibc.meta(c).boundary_type[bc_type_idx] == 0)
            {
                if(m_ibc.meta(c).boundary_order == 2)
                {
                    rhsfull[m_ibc.IdentityRowsStart() + m_ibc.DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + m_ibc.DirichletIBCColumnStartIndex(c)] = surface_exact(m_ibc.DOFSubset()[c].cpx()[i]);
                }
                else
                {
                    rhsfull[m_ibc.IdentityRowsStart() + i + m_ibc.IBCColumnStartIndex(c)] = surface_exact(m_ibc.DOFSubset()[c].cpx()[i]);
                }
            }
        }
    }

    return rhsfull;
}

VectorX PoissonWithIBC::GetRHS(function<Scalar(const vector<Scalar>&)> f, function<Scalar(size_t &, const size_t &, const vector<Scalar>&)> ibc_exact)
{
    VectorX rhsfull(m_surface_tube.nNodes() + m_ibc.TotalIBCDOFs() + m_ibc.Num2ndOrderDirichletRows());

    for(size_t i = 0; i < m_surface_tube.nNodes(); ++i)
    {
        rhsfull[i] = f(m_surface_tube.cpx()[i]);
    }

    for(size_t c = 0; c < m_ibc.NumIBCs(); ++c)
    {
        for(size_t i = 0; i < m_ibc.DOFSubset()[c].size(); ++i)
        {
            const size_t bc_type_idx = m_ibc.meta(c).is_oriented ? m_ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
            if(m_ibc.meta(c).boundary_type[bc_type_idx] == 1 || (m_ibc.meta(c).boundary_type[bc_type_idx] == 0 && m_ibc.meta(c).boundary_order == 2))
            {
                rhsfull[i + m_ibc.IBCColumnStartIndex(c)] = f(m_ibc.DOFSubset()[c].cpx()[i]);
            }
        }
    }

    for(size_t c = 0; c < m_ibc.NumIBCs(); ++c)
    {
        for(size_t i = 0; i < m_ibc.DOFSubset()[c].size(); ++i)
        {
            const size_t bc_type_idx = m_ibc.meta(c).is_oriented ? m_ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
            if(m_ibc.meta(c).boundary_type[bc_type_idx] == 0)
            {
                if(m_ibc.meta(c).boundary_order == 2)
                {
                    rhsfull[m_ibc.IdentityRowsStart() + m_ibc.DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + m_ibc.DirichletIBCColumnStartIndex(c)] = ibc_exact(c, m_ibc.DOFSubset()[c].which_side()[i], m_ibc.DOFSubset()[c].cpx()[i]);
                }
                else
                {
                    rhsfull[m_ibc.IdentityRowsStart() + i + m_ibc.IBCColumnStartIndex(c)] = ibc_exact(c, m_ibc.DOFSubset()[c].which_side()[i], m_ibc.DOFSubset()[c].cpx()[i]);
                }
            }
        }
    }

    return rhsfull;
}

SpMat PoissonWithIBC::GetPlotInterpMatrix()
{
    SpMat plotE;
    m_h.GetPlotInterpMatrix(m_surface_tube, plotE);

    SpMat plotEfull = m_ibc_mat_manip.AddIBCPlotInterpolation(plotE, m_surface_tube, m_ibc, m_surface_tube.surface().xp());

    return plotEfull;
}


SpMat PoissonWithIBC::GetPlotInterpMatrix(vector<vector<size_t>> &director_set_index_from_subset_index, vector<vector<size_t>> &director_which_side, vector<vector<bool>> &director_is_on_ibc)
{
    SpMat plotE;
    m_h.GetPlotInterpMatrix(m_surface_tube, plotE);

    SpMat plotEfull = m_ibc_mat_manip.AddIBCPlotInterpolation(plotE, m_surface_tube, m_ibc, m_surface_tube.surface().xp(), director_set_index_from_subset_index, director_which_side, director_is_on_ibc);

    return plotEfull;
}


////////////////////////////////////////////////////////////////////////////////////
// Private Functions
////////////////////////////////////////////////////////////////////////////////////


void PoissonWithIBC::Initialize()
{   
#ifdef TRACK_WHERE
    std::cout << "PoissonWithIBC::Initialize" << std::endl;
#endif
    // construct the discrete Laplacian and cp extension matrices
    m_h.SetupMatrices(m_surface_tube, m_E, m_L);

    m_Lfull = m_ibc_mat_manip.AddIBCFiniteDifference(m_L, m_surface_tube, m_ibc, m_surface_tube.x());
    m_Efull = m_ibc_mat_manip.AddIBCExtension(m_E, m_surface_tube, m_ibc, m_surface_tube.cpx());

    m_identity_rows = vector<bool>(m_surface_tube.nNodes() + m_ibc.TotalIBCDOFs() + m_ibc.Num2ndOrderDirichletRows(), false);
    for(size_t c = 0; c < m_ibc.NumIBCs(); ++c)
    {
        for(size_t i = 0; i < m_ibc.DOFSubset()[c].size(); ++i)
        {
            const size_t bc_type_idx = m_ibc.meta(c).is_oriented ? m_ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
            if(m_ibc.meta(c).boundary_type[bc_type_idx] == 0 /* Dirichlet */)
            {
                if(m_ibc.meta(c).boundary_order == 2)
                {
                    m_identity_rows[m_ibc.IdentityRowsStart() + m_ibc.DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + m_ibc.DirichletIBCColumnStartIndex(c)] = true;
                }
                else
                {
                    m_identity_rows[m_ibc.IdentityRowsStart() + i + m_ibc.IBCColumnStartIndex(c)] = true;
                }
            }
        }
    }
}


void PoissonWithIBC::BuildExtensionMatrix()
{
    SpMat E;
    E.resize(m_surface_tube.cpx().size(), m_surface_tube.nNodes());
    m_h.GetExtensionMatrix(m_surface_tube, E);
}

} // namespace cpm