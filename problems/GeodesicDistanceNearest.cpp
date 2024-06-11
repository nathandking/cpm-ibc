#include "GeodesicDistanceNearest.h"

namespace cpm {

GeodesicDistanceNearest::GeodesicDistanceNearest(bool visualize)
{
    m_visualize = visualize;
}


VectorX GeodesicDistanceNearest::ComputeDistance(Tube &surface_tube, vector<vector<size_t>> &ibc_idx, vector<vector<size_t>> &inside_ibc_idx, vector<vector<Scalar>> &inside_ibc_idx_dist, Scalar dt)
{
    ComputeSurfaceNormals(surface_tube);
    
    m_h.SetupMatrices(surface_tube, m_E, m_L);

    m_identity_rows.assign(surface_tube.nNodes(), false);
    for(size_t c = 0; c < ibc_idx.size(); ++c)
    {
        for(size_t i = 0; i < ibc_idx[c].size(); ++i)
        {
            m_identity_rows[ibc_idx[c][i]] = true;
        }
    }

    VectorX u = -Step1(surface_tube, ibc_idx, inside_ibc_idx, inside_ibc_idx_dist, dt);
    
    vector<VectorX> gradu = Step2(u, surface_tube);

    VectorX div = ComputeDivergence(gradu, surface_tube);
    
    VectorX phi = Step3(div, surface_tube, ibc_idx);

    for(size_t k = 0; k < 2; ++k)
    {
        // phi = m_Efull * phi;
        for(size_t i = 0; i < surface_tube.nNodes(); ++i)
        {
            u[i] = phi[i];
        }
        u = m_E * u;

        gradu = Step2(u, surface_tube);

        div = ComputeDivergence(gradu, surface_tube);

        phi = Step3(div, surface_tube, ibc_idx);
    }

    SpMat plotE;
    m_h.GetPlotInterpMatrix(surface_tube, plotE);

    return plotE * phi;
}


/////////////////////////////////////////////////////////////////////////
// Step I: Heat Flow From IBC
/////////////////////////////////////////////////////////////////////////
VectorX GeodesicDistanceNearest::Step1(Tube &surface_tube, vector<vector<size_t>> &ibc_idx, vector<vector<size_t>> &inside_ibc_idx, vector<vector<Scalar>> &inside_ibc_idx_dist, Scalar dt)
{   
    Scalar extent = surface_tube.TubeRadius(); // extent of smoothing is 1/2 the ibc tube radius
    Scalar tol = surface_tube.dx();
    Scalar scale = atanh(1 - tol) / extent;

    auto SmoothH = [&scale](Scalar& t)
    {
        return 0.5 * tanh(-scale * t) + 0.5;
    };

    VectorX rhs = VectorX::Zero(surface_tube.nNodes());
    for(size_t c = 0; c < ibc_idx.size(); ++c)
    {
        for(size_t i = 0; i < ibc_idx[c].size(); ++i)
        {
            Scalar t = -Distance(surface_tube.x()[ibc_idx[c][i]], surface_tube.cpx()[ibc_idx[c][i]]);
            rhs[ibc_idx[c][i]] = SmoothH(t);
            // rhs[ibc_idx[c][i]] = 1;
        }
    }

    for(size_t c = 0; c < inside_ibc_idx.size(); ++c)
    {
        for(size_t i = 0; i < inside_ibc_idx[c].size(); ++i)
        {
            Scalar t = -inside_ibc_idx_dist[c][i];
            rhs[inside_ibc_idx[c][i]] = SmoothH(t);
        }
    }

    SpMat Afull = m_h.ImplicitEulerMatrix(m_L, m_E, dt);
    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Ar = Afull.transpose();
    // m_h.SetIdentityRows(Ar, m_identity_rows); // this only needs to be done because of the -ve sign, all the identity rows are given -1, so I just set the identity rows over again
    for (size_t k=0; k < Ar.outerSize(); ++k)
    {
        bool diagonal_entry_exists = false;
        for (Eigen::SparseMatrix<Scalar, Eigen::ColMajor>::InnerIterator it(Ar, k); it; ++it)
        {
            if(m_identity_rows[it.row()])
            {
                if(it.row() == it.col())
                {
                    diagonal_entry_exists = true;
                    it.valueRef() = 1.0;
                }
                else
                {
                    it.valueRef() = 0.0;
                }
            }
        }

        if(!diagonal_entry_exists && m_identity_rows[k]) // if the diagonal entry is zero we must insert a new element into the sparse matrix
        {
            Ar.insert(k,k) = 1.0;
        }
    }
    Ar.prune(0,0);  
    
    Eigen::SparseLU<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver;
    solver.compute(Ar);
    if(solver.info()!=Eigen::Success) {
        cout << "LU decomposition failed" << endl;
    }

    VectorX u = solver.solve(rhs);
    if(solver.info()!=Eigen::Success) {
        cout << "solver failed" << endl;
    }
    // do a cp extension of the heat flow above
    u = m_E * u;  // doing this extension increases the error and makes it not converge as well

#ifdef POLYSCOPE
    if(m_visualize)
    {
        if(surface_tube.dim() == 2)
        {
            polyscope::registerPointCloud2D("surface closest points", surface_tube.cpx())->addScalarQuantity("heat", u);;
        }
        else if(surface_tube.dim() == 3)
        {
            polyscope::registerPointCloud("surface closest points", surface_tube.cpx())->addScalarQuantity("heat", u);;
        }
    }
#endif

    return u;
}


/////////////////////////////////////////////////////////////////////////
// Step II: Compute Negative Normalized Gradient
/////////////////////////////////////////////////////////////////////////
vector<VectorX> GeodesicDistanceNearest::Step2(VectorX &u, Tube &surface_tube)
{
    // build gradient matrices
    FDMatrices FDMat;
    m_Dc.resize(surface_tube.dim());
    for(int d = 0; d < surface_tube.dim(); ++d)
    {
        m_Dc[d].resize(surface_tube.nNodes(), surface_tube.nNodes());
    }
    FDMat.BuildGradientMatrices(surface_tube, m_Dc); 

    //////////////////////////////////////////////////////////////////////
    // compute negative gradient with ibc
    /////////////////////////////////////////////////////////////////////
    vector<VectorX> gradu(surface_tube.dim());
    for(int d = 0; d < surface_tube.dim(); ++d)
    {
        gradu[d] = m_Dc[d] * u;
    }

    for(int d = 0; d < surface_tube.dim(); ++d)
    {
        gradu[d] = m_E * gradu[d];
    }

    // project out normal component of gradient, without this it does not converge
    for(size_t i = 0; i < m_surface_normals.size(); ++i)
    {
        vector<Scalar> grad(surface_tube.dim());
        for(size_t d = 0; d < surface_tube.dim(); ++d)
        {
            grad[d] = gradu[d][i];
        }
        
        grad = grad - DotProduct(grad, m_surface_normals[i]) * m_surface_normals[i];
        
        for(size_t d = 0; d < surface_tube.dim(); ++d)
        {
            gradu[d][i] = grad[d];
        }
    }

    // normalize the negative gradient
    for(int i = 0; i < surface_tube.nNodes(); ++i)
    {
        Scalar gradNorm = 0.0;
        for(int d = 0; d < surface_tube.dim(); ++d)
        {
            gradNorm += pow(gradu[d][i], 2);
        }
        gradNorm = sqrt(gradNorm);

        if(gradNorm != 0.0)
        {
            for(int d = 0; d < surface_tube.dim(); ++d)
            {
                gradu[d][i] /= gradNorm;
            }
        }
    }

#ifdef POLYSCOPE
    if(m_visualize)
    {
        vector<Vector3> gradu_surface(surface_tube.nNodes());
        for(int i = 0; i < surface_tube.nNodes(); ++i)
        {
            for(int d = 0; d < surface_tube.dim(); ++d)
            {
                gradu_surface[i][d] = gradu[d][i];
            }
        }
        if(surface_tube.dim() == 2)
        { 
            polyscope::getPointCloud("surface closest points")->addVectorQuantity2D("grad heat", gradu_surface);
        }
        else
        {
            polyscope::getPointCloud("surface closest points")->addVectorQuantity("grad heat", gradu_surface);
        }
    }
#endif

    return gradu;
}


/////////////////////////////////////////////////////////////////////////
// Step III: Solve Poisson Equation for Distance
/////////////////////////////////////////////////////////////////////////
VectorX GeodesicDistanceNearest::Step3(VectorX &div, Tube &surface_tube, vector<vector<size_t>> &ibc_idx)
{
    // set RHS for ibc point values to zero (for zero distance to ibc point)
    for(size_t c = 0; c < ibc_idx.size(); ++c)
    {
        for(size_t i = 0; i < ibc_idx[c].size(); ++i)
        {
            div[ibc_idx[c][i]] = 0.0;
        }
    }
    Solver sv;
    sv.Init((Scalar)0., (Scalar)1., m_L, m_E, Preconditioner::diagonal, m_identity_rows);
    cpm::VectorX phi=sv.Solve(div);

    // make sure the minimum distance is zero, should it be though if none of the closest points are on the ibc?
    Scalar min_phi = numeric_limits<Scalar>::max();
    for(size_t i = 0; i < phi.size(); ++i)
    {
        if(min_phi > phi[i])
        {
            min_phi = phi[i];
        }
    }

    for(size_t i = 0; i < phi.size(); ++i)
    {
        phi[i] -= min_phi;   
    }

#ifdef POLYSCOPE
    if(m_visualize)
    {
        // needed if computing error on closest points themselves, but no change otherwise, so save computation
        phi = m_E * phi;

        VectorX uexact(surface_tube.nNodes());
        for(int i = 0; i < surface_tube.nNodes(); ++i)
        {
            uexact[i] = acos(DotProduct(surface_tube.cpx()[i], surface_tube.cpx()[ibc_idx[0][0]]));
        }

        polyscope::getPointCloud("surface closest points")->addScalarQuantity("Distance", phi);
        polyscope::getPointCloud("surface closest points")->addScalarQuantity("div", div);

        // VectorX error = uexact - phi_surface;
        // cout << error.lpNorm<Eigen::Infinity>() << endl;
        // cout << error.lpNorm<2>() / surface_tube.nNodes() << endl;
    }
#endif

    return phi;
}


VectorX GeodesicDistanceNearest::ComputeDivergence(vector<VectorX> &gradu, Tube &surface_tube)
{
    // compute divergence of negative normalized gradient of u
    VectorX div = VectorX::Zero(surface_tube.nNodes());
    for(size_t d = 0; d < surface_tube.dim(); ++d)
    { 
        div += m_Dc[d] * gradu[d];
    }

    return div;
}


void GeodesicDistanceNearest::ComputeSurfaceNormals(Tube &surface_tube)
{
    if(surface_tube.dim() == 2)
    {
        vector<Matrix2> eigenvectors;
        vector<Vector2> eigenvalues;
        m_geom.EigenDecompositionOfJcp(surface_tube, eigenvectors, eigenvalues);

        m_surface_normals.resize(surface_tube.nNodes(), vector<Scalar>(surface_tube.dim()));
        for(size_t i = 0; i < surface_tube.nNodes(); ++i)
        {
            for(size_t d = 0; d < surface_tube.dim(); ++d)
            {
                m_surface_normals[i][d] = eigenvectors[i](d, 0);
            }
        }
    }
    else
    {
        vector<Matrix3> eigenvectors;
        vector<Vector3> eigenvalues;
        m_geom.EigenDecompositionOfJcp(surface_tube, eigenvectors, eigenvalues);

        m_surface_normals.resize(surface_tube.nNodes(), vector<Scalar>(surface_tube.dim()));
        for(size_t i = 0; i < surface_tube.nNodes(); ++i)
        {
            for(size_t d = 0; d < surface_tube.dim(); ++d)
            {
                m_surface_normals[i][d] = eigenvectors[i](d, 0);
            }
        }
    }
}


vector<vector<Scalar>> GeodesicDistanceNearest::InterpolateUnorientedVectors(Tube &surface_tube, const vector<vector<Scalar>> &tube_vectors, const vector<vector<Scalar>> &xq, bool is_normalized)
{
    // build interpolation matrix for interpolating from the tube points onto query points xq
    SpMat E(xq.size(), surface_tube.nNodes());
    Interpolation interp(surface_tube, xq);
    interp.BuildInterpolationMatrix(surface_tube, E);

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
        size_t num_interp_points = pow(surface_tube.p() + 1, surface_tube.dim());
        for(size_t i = 0; i < vq.size(); ++i)
        {
            vq[i] /= num_interp_points;
        }
    }

    return vq;
}


} // namespace cpm