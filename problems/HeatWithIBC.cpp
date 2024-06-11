#include "HeatWithIBC.h"

// #define TIME

namespace cpm {
    
HeatWithIBC::HeatWithIBC(SurfaceSpecifier &surface_specs, vector<SurfaceSpecifier> &ibc_specs, Scalar dx, Scalar dt, size_t time_stepping_method, size_t interp_deg)
{
    m_dt = dt;
    m_time_stepping_method = time_stepping_method;

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
    cout << "HeatWithIBC Initialize time: " << m_elapsed_seconds.count() << "s\n";
#endif
}


void HeatWithIBC::SetDirichletValue(VectorX &ufull, function<Scalar(size_t &, const size_t &, const vector<Scalar>&)> ibc_exact)
{
#ifdef TRACK_WHERE
    std::cout << "HeatWithIBC::SetDirichletValue" << std::endl;
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


VectorX HeatWithIBC::TimeStep(VectorX &ufull, function<Scalar(size_t &, const size_t &, const vector<Scalar>&)> ibc_exact)
{
    if(m_time_stepping_method > 0) // implicit methods
    {
        size_t DOFs = m_Lfull.rows();
        if(m_time_stepping_method == 2) // Crank-Nicolson
        {
            cpm::VectorX tmp = cpm::VectorX::Zero(DOFs);
            sv.Multiply(ufull, tmp, (Scalar)1., 0.5 * m_dt, m_Lfull, m_Efull);
            for(int i=0;i<DOFs;++i) ufull(i) = tmp(i);
        }
        SetDirichletValue(ufull, ibc_exact);
        cpm::VectorX res = sv.Solve(ufull);
        return res;
    }
    else
    {
        ufull = ufull + m_dt * (m_Lfull * ufull);
        ufull = m_Efull * ufull; // closest point extension
        SetDirichletValue(ufull, ibc_exact);
    }

    return ufull;
}

SpMat HeatWithIBC::GetPlotInterpMatrix()
{
    SpMat plotE;
    m_h.GetPlotInterpMatrix(m_surface_tube, plotE);

    SpMat plotEfull = m_ibc_mat_manip.AddIBCPlotInterpolation(plotE, m_surface_tube, m_ibc, m_surface_tube.surface().xp());

    return plotEfull;
}


SpMat HeatWithIBC::GetInterpMatrix(vector<vector<Scalar>> &xp)
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


void HeatWithIBC::Initialize()
{   
#ifdef TRACK_WHERE
    std::cout << "HeatWithIBC::Initialize" << std::endl;
#endif
    // construct the discrete Laplacian and cp extension matrices
    m_h.SetupMatrices(m_surface_tube, m_E, m_L);

    m_Lfull = m_ibc_mat_manip.AddIBCFiniteDifference(m_L, m_surface_tube, m_ibc, m_surface_tube.x());
    m_Efull = m_ibc_mat_manip.AddIBCExtension(m_E, m_surface_tube, m_ibc, m_surface_tube.cpx());

    m_identity_rows = vector<bool>(m_surface_tube.nNodes() + m_ibc.TotalIBCDOFs() + m_ibc.IdentityRowsStart(), false);
#ifdef OPENMP
#pragma omp parallel for
#endif
    for(size_t c = 0; c < m_ibc.NumIBCs(); ++c)
    {
        if(m_ibc.meta(c).boundary_type[0] == 0)
        {
            for(size_t i = 0; i < m_ibc.DOFSubset()[c].size(); ++i)
            {
                m_identity_rows[m_ibc.IdentityRowsStart() + i + m_ibc.IBCColumnStartIndex(c)] = true;
            }
        }
    }
    // m_h.SetIdentityRows(m_Efull, identity_rows); // TO DO: speed this up by only doing it for contraint rows and only if Neumann
    
    if(m_time_stepping_method > 0) // implicit methods
    {
        m_Lfull.makeCompressed();
        m_Efull.makeCompressed();
        Scalar r = m_time_stepping_method == 1 ? 1.0 : 0.5;
        sv.Init((Scalar)1., -r * m_dt, m_Lfull, m_Efull, Preconditioner::diagonal, m_identity_rows);
    }
}


void HeatWithIBC::ComputeSurfaceNormals()
{
    m_surface_normals = m_ibc.TaggingSubset()[0].surface_normals_on_tube();
}


vector<vector<Scalar>> HeatWithIBC::InterpolateUnorientedVectors(const vector<vector<Scalar>> &tube_vectors, const vector<vector<Scalar>> &xq, bool is_normalized)
{
    // build interpolation matrix for interpolating from the tube points onto query points xq
    SpMat E(xq.size(), m_surface_tube.nNodes());
    Interpolation interp(m_surface_tube, xq);
    interp.BuildInterpolationMatrix(m_surface_tube, E);

    // now apply interpolation one weight at a time, checking if there is a flip in the sign of the vector, since they are unoriented
    vector<vector<Scalar>> vq(xq.size(), vector<Scalar>(tube_vectors[0].size(), 0.0));

    vector<bool> first_row_visit(E.rows(), true);
    vector<vector<Scalar>> first_direction(E.rows(), vector<Scalar>(tube_vectors[0].size(), 0.0));
    for (size_t k=0; k < E.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(E,k); it; ++it)
        {   
            if(first_row_visit[it.row()])
            {
                first_direction[it.row()] = tube_vectors[it.col()];
                first_row_visit[it.row()] = false;
            }

            Scalar dot_product = DotProduct(first_direction[it.row()], tube_vectors[it.col()]);
            assert(dot_product != 0.0);

            if(dot_product > 0.0)
            {
                vq[it.row()] += it.valueRef() * tube_vectors[it.col()];
            }
            else // flip vector direction before using it in the interpolation
            {
                vq[it.row()] -= it.valueRef() * tube_vectors[it.col()];
            }
        }
    }

    if(is_normalized)
    {
        // normalize the accumulated vector
        for(size_t i = 0; i < vq.size(); ++i)
        {
            Normalize(vq[i], 0.0);
        }
    }
    else
    {
        // divide by the number of accumlated vectors, i.e., number of points in the interpolation stencil
        size_t num_interp_points = pow(m_surface_tube.p() + 1, m_surface_tube.dim());
        for(size_t i = 0; i < vq.size(); ++i)
        {
            vq[i] /= num_interp_points;
        }
    }

    return vq;
}


} // namespace cpm