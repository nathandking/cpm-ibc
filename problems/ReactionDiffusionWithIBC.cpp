#include "ReactionDiffusionWithIBC.h"

// #define TIME

namespace cpm {
    
ReactionDiffusionWithIBC::ReactionDiffusionWithIBC(SurfaceSpecifier &surface_specs, vector<SurfaceSpecifier> &ibc_specs, Scalar dx, Scalar dt, size_t interp_deg)
{
    m_dt = dt;

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
    cout << "ReactionDiffusionWithIBC Initialize time: " << m_elapsed_seconds.count() << "s\n";
#endif
}


void ReactionDiffusionWithIBC::SetDirichletValue(VectorX &ufull, function<Scalar(size_t &, const size_t &, const vector<Scalar>&)> ibc_exact)
{
#ifdef TRACK_WHERE
    std::cout << "ReactionDiffusionWithIBC::SetDirichletValue" << std::endl;
#endif
    //set exact solution for tube points extra "virtual" ibc points
#ifdef OPENMP
#pragma omp parallel for
#endif
    for(size_t c = 0; c < m_ibc.NumIBCs(); ++c)
    {
        for(size_t i = 0; i < m_ibc.DOFSubset()[c].size(); ++i)
        {
            const size_t bc_type_idx = m_ibc.meta(c).is_oriented ? m_ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
            if(m_ibc.meta(c).boundary_type[bc_type_idx] == 0)
            {
                if(m_ibc.meta(c).boundary_order == 2)
                {
                    ufull[m_ibc.IdentityRowsStart() + m_ibc.DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + m_ibc.DirichletIBCColumnStartIndex(c)] = ibc_exact(c, m_ibc.DOFSubset()[c].which_side()[i], m_ibc.DOFSubset()[c].cpx()[i]);
                }
                else
                {
                    ufull[m_ibc.IdentityRowsStart() + i + m_ibc.IBCColumnStartIndex(c)] = ibc_exact(c, m_ibc.DOFSubset()[c].which_side()[i], m_ibc.DOFSubset()[c].cpx()[i]);
                }
            }
        }
    }
}


void ReactionDiffusionWithIBC::TimeStep(VectorX &ufull, VectorX &vfull, function<Scalar(size_t &, const size_t &, const vector<Scalar>&)> ibc_exact_u, function<Scalar(size_t &, const size_t &, const vector<Scalar>&)> ibc_exact_v)
{
    Scalar nu_u = 1.0 / pow(3.0 / m_surface_tube.dx(), 2);
    Scalar nu_v = nu_u / 3.0;
    Scalar F = 0.054;
    Scalar k = 0.063;

    VectorX ufull_new = ufull + nu_u * m_dt * (m_Lfull * ufull);
    VectorX vfull_new = vfull + nu_v * m_dt * (m_Lfull * vfull);

    for(size_t i = 0; i < ufull_new.size(); ++i)
    {
        ufull_new[i] += m_dt * (-ufull[i] * vfull[i] * vfull[i] + F * (1.0 - ufull[i]));
        vfull_new[i] += m_dt * (ufull[i] * vfull[i] * vfull[i] - (F + k) * vfull[i]);
    }


    ufull_new = m_Efull * ufull_new; // closest point extension
    vfull_new = m_Efull * vfull_new; // closest point extension

    SetDirichletValue(ufull_new, ibc_exact_u);
    SetDirichletValue(vfull_new, ibc_exact_v);

    ufull = ufull_new;
    vfull = vfull_new;
}

SpMat ReactionDiffusionWithIBC::GetPlotInterpMatrix()
{
    SpMat plotE;
    m_h.GetPlotInterpMatrix(m_surface_tube, plotE);

    SpMat plotEfull = m_ibc_mat_manip.AddIBCPlotInterpolation(plotE, m_surface_tube, m_ibc, m_surface_tube.surface().xp());

    return plotEfull;
}


SpMat ReactionDiffusionWithIBC::GetInterpMatrix(vector<vector<Scalar>> &xp)
{
    SpMat plotE;
    plotE.resize(xp.size(), m_surface_tube.nNodes());

    Interpolation plotInterp(m_surface_tube, xp);
    plotInterp.BuildInterpolationMatrix(m_surface_tube, plotE);


    SpMat plotEfull = m_ibc_mat_manip.AddIBCPlotInterpolation(plotE, m_surface_tube, m_ibc, xp);

    return plotEfull;
}

////////////////////////////////////////////////////////////////////////////////////
// Private Functions
////////////////////////////////////////////////////////////////////////////////////


void ReactionDiffusionWithIBC::Initialize()
{   
#ifdef TRACK_WHERE
    std::cout << "ReactionDiffusionWithIBC::Initialize" << std::endl;
#endif
    // construct the discrete Laplacian and cp extension matrices
    m_h.SetupMatrices(m_surface_tube, m_E, m_L);

    m_Lfull = m_ibc_mat_manip.AddIBCFiniteDifference(m_L, m_surface_tube, m_ibc, m_surface_tube.x());
    m_Efull = m_ibc_mat_manip.AddIBCExtension(m_E, m_surface_tube, m_ibc, m_surface_tube.cpx());

    m_identity_rows = vector<bool>(m_surface_tube.nNodes() + m_ibc.TotalIBCDOFs() + m_ibc.Num2ndOrderDirichletRows(), false);

#ifdef OPENMP
#pragma omp parallel for
#endif
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
    m_h.SetIdentityRows(m_Efull, m_identity_rows);
}


} // namespace cpm