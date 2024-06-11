#include "GeodesicDistance.h"

namespace cpm {

GeodesicDistance::GeodesicDistance(bool visualize)
{
    m_visualize = visualize;
}


VectorX GeodesicDistance::ComputeDistance(Tube &surface_tube, IBCSubsets &ibc, Scalar dt)
{
    m_h.SetupMatrices(surface_tube, m_E, m_L);

    m_Lfull = m_ibc_mat_manip.AddIBCFiniteDifference(m_L, surface_tube, ibc, surface_tube.x());
    m_Efull = m_ibc_mat_manip.AddIBCExtension(m_E, surface_tube, ibc, surface_tube.cpx());

    m_identity_rows.assign(surface_tube.nNodes() + ibc.TotalIBCDOFs() + ibc.Num2ndOrderDirichletRows(), false);
    for(size_t c = 0; c < ibc.NumIBCs(); ++c)
    {
        for(size_t i = 0; i < ibc.DOFSubset()[c].size(); ++i)
        {
            const size_t bc_type_idx = ibc.meta(c).is_oriented ? ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
            if(ibc.meta(c).boundary_type[bc_type_idx] == 0 /* Dirichlet */)
            {
                if(ibc.meta(c).boundary_order == 2)
                {
                    m_identity_rows[ibc.IdentityRowsStart() + ibc.DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + ibc.DirichletIBCColumnStartIndex(c)] = true;
                }
                else
                {
                    m_identity_rows[ibc.IdentityRowsStart() + i + ibc.IBCColumnStartIndex(c)] = true;
                }
            }
        }
    }

    VectorX u = -Step1(surface_tube, ibc, dt);
    
    vector<VectorX> gradu = Step2(u, surface_tube, ibc);

    VectorX div = ComputeDivergence(gradu, surface_tube, ibc);
    // div = ExactDivergenceSpherePoint(surface_tube, ibc); // this does not work anymore
    // div = ExactDivergenceCirclePoint(surface_tube, ibc); // this does not work, but the divergence is actually zero, and infinite at the source and cut locus, maybe this is just a really hard problem to solve

    VectorX phi = Step3(div, surface_tube, ibc);

    for(size_t k = 0; k < 2; ++k)
    {
        // phi = m_Efull * phi;
        for(size_t i = 0; i < surface_tube.nNodes(); ++i)
        {
            u[i] = phi[i];
        }
        u = m_E * u;

        gradu = Step2(u, surface_tube, ibc);

        div = ComputeDivergence(gradu, surface_tube, ibc);

        phi = Step3(div, surface_tube, ibc);
    }

    SpMat plotE;
    m_h.GetPlotInterpMatrix(surface_tube, plotE);

    SpMat plotEfull = m_ibc_mat_manip.AddIBCPlotInterpolation(plotE, surface_tube, ibc, surface_tube.surface().xp());
    return plotEfull * phi;
}


/////////////////////////////////////////////////////////////////////////
// Step I: Heat Flow From IBC
/////////////////////////////////////////////////////////////////////////
VectorX GeodesicDistance::Step1(Tube &surface_tube, IBCSubsets &ibc, Scalar dt)
{   
    Scalar extent = surface_tube.TubeRadius(); // extent of smoothing is 1/2 the ibc tube radius
    Scalar tol = surface_tube.dx();
    Scalar scale = atanh(1 - tol) / extent;

    auto SmoothH = [&scale](Scalar& t)
    {
        return 0.5 * tanh(-scale * t) + 0.5;
    };

    VectorX rhsfull = VectorX::Zero(surface_tube.nNodes() + ibc.TotalIBCDOFs() + ibc.Num2ndOrderDirichletRows());
    for(size_t c = 0; c < ibc.NumIBCs(); ++c)
    {
        for(size_t i = 0; i < ibc.DOFSubset()[c].size(); ++i)
        {
            const size_t bc_type_idx = ibc.meta(c).is_oriented ? ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
            if(ibc.meta(c).boundary_type[bc_type_idx] == 0)
            {
                Scalar t = -Norm(ibc.DOFSubset()[c].cp_diff()[i]);
                if(ibc.meta(c).boundary_order == 2)
                {
                    rhsfull[ibc.IdentityRowsStart() + ibc.DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + ibc.DirichletIBCColumnStartIndex(c)] = SmoothH(t);
                }
                else
                {
                    rhsfull[ibc.IdentityRowsStart() + i + ibc.IBCColumnStartIndex(c)] = SmoothH(t);
                }
            }
        }
    }

    // smooth the initial condition
    for(size_t c = 0; c < ibc.NumIBCs(); ++c)
    {
        for(size_t i = 0; i < ibc.DOFSubset()[c].size(); ++i)
        {
            Scalar t = Norm(ibc.DOFSubset()[c].cp_diff()[i]);
            rhsfull[ibc.DOFSubset()[c].SetIndexFromSubsetIndex()[i]] = SmoothH(t);
        }
    }
    
    Solver sv;
    sv.Init((Scalar)1., -dt, m_Lfull, m_Efull, Preconditioner::jacobi, m_identity_rows);
    sv.SetPrecond(2000);
    cpm::VectorX ufull = sv.Solve(rhsfull);

#ifdef POLYSCOPE
    if(m_visualize)
    {
        VectorX u_surface(surface_tube.nNodes());
        for(int i = 0; i < surface_tube.nNodes(); ++i)
        {
            u_surface[i] = ufull[i];
        }

        if(surface_tube.dim() == 2)
        {
            polyscope::registerPointCloud2D("surface closest points", surface_tube.cpx())->addScalarQuantity("heat", u_surface);;
        }
        else if(surface_tube.dim() == 3)
        {
            polyscope::registerPointCloud("surface closest points", surface_tube.cpx())->addScalarQuantity("heat", u_surface);;
        }

        for(size_t c = 0; c < ibc.NumIBCs(); ++c)
        {
            vector<vector<Scalar>> ibc_cpx_s(ibc.DOFSubset()[c].size(), vector<Scalar>(surface_tube.dim()));
            vector<Scalar> u_ibc(ibc.DOFSubset()[c].size());
            for(size_t i = 0; i < ibc.DOFSubset()[c].size(); ++i)
            {
                size_t dof_idx = ibc.IdentityRowsStart() + i + ibc.IBCColumnStartIndex(c);
                u_ibc[i] = ufull[dof_idx];
                for(size_t d = 0; d < surface_tube.dim(); ++d)
                {
                    ibc_cpx_s[i][d] = surface_tube.cpx()[ibc.DOFSubset()[c].SetIndexFromSubsetIndex()[i]][d];
                }
            }
        
            if(surface_tube.dim() == 2)
            {
                polyscope::registerPointCloud2D("ibc closest points " + to_string(c), ibc_cpx_s)->addScalarQuantity("ibc heat", u_ibc);
            }
            else if(surface_tube.dim() == 3)
            {
                polyscope::registerPointCloud("ibc closest points " + to_string(c), ibc_cpx_s)->addScalarQuantity("ibc heat", u_ibc);
            }
        }
    }
#endif

    VectorX u(surface_tube.nNodes());
    for(size_t i = 0; i < surface_tube.nNodes(); ++i)
    {
        u[i] = ufull[i];
    }

    return u;
}


/////////////////////////////////////////////////////////////////////////
// Step II: Compute Negative Normalized Gradient
/////////////////////////////////////////////////////////////////////////
vector<VectorX> GeodesicDistance::Step2(VectorX &u, Tube &surface_tube, IBCSubsets &ibc)
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
    // u = ExactMatchStep1ScalarField(surface_tube, ibc); // FOR TESTING ON SPHERE OR CIRCLE

#ifdef POLYSCOPE
    if(m_visualize)
    {
        if(surface_tube.dim() == 2)
        {
            polyscope::registerPointCloud2D("surface closest points", surface_tube.cpx())->addScalarQuantity("Matched Step 1", u);
        }
        else if(surface_tube.dim() == 3)
        {
            polyscope::registerPointCloud("surface closest points", surface_tube.cpx())->addScalarQuantity("Matched Step 1", u);
        }
    }
#endif

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
    m_surface_normals = ibc.TaggingSubset()[0].surface_normals_on_tube();
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

    // FOR TESTING ON A SPHERE OR CIRCLE
    // vector<VectorX> graduexact = ExactGradientSpherePoint(surface_tube, ibc);
    // gradu = graduexact;
    // vector<Scalar> theta(surface_tube.nNodes());
    // for(size_t i = 0; i < surface_tube.nNodes(); ++i)
    // {
    //     vector<Scalar> grad(surface_tube.dim());
    //     vector<Scalar> gradexact(surface_tube.dim());
    //     for(size_t d = 0; d < surface_tube.dim(); ++d)
    //     {
    //         grad[d] = gradu[d][i];
    //         gradexact[d] = graduexact[d][i];
    //     }

    //     Scalar dot = DotProduct(gradexact, grad);

    //     if(abs(dot) > 1)
    //     {
    //         theta[i] = 0.0;
    //     }
    //     else
    //     {
    //         theta[i] = acos(dot);
    //     }
    // }

#ifdef POLYSCOPE
    if(m_visualize)
    {
        vector<Vector3> gradu_surface(surface_tube.nNodes());
        // vector<Vector3> gradu_surface_exact(surface_tube.nNodes());
        for(int i = 0; i < surface_tube.nNodes(); ++i)
        {
            for(int d = 0; d < surface_tube.dim(); ++d)
            {
                gradu_surface[i][d] = gradu[d][i];
                // gradu_surface_exact[i][d] = graduexact[d][i];
            }
        }
        if(surface_tube.dim() == 2)
        { 
            polyscope::getPointCloud("surface closest points")->addVectorQuantity2D("grad heat", gradu_surface);
            //polyscope::getPointCloud("surface closest points")->addVectorQuantity2D("grad heat exact", gradu_surface_exact);
        }
        else
        {
            polyscope::getPointCloud("surface closest points")->addVectorQuantity("grad heat", gradu_surface);
            //polyscope::getPointCloud("surface closest points")->addVectorQuantity("grad heat exact", gradu_surface_exact);
        }
        //polyscope::getPointCloud("surface closest points")->addScalarQuantity("theta diff", theta);
    }
#endif

    return gradu;
}


/////////////////////////////////////////////////////////////////////////
// Step III: Solve Poisson Equation for Distance
/////////////////////////////////////////////////////////////////////////
VectorX GeodesicDistance::ComputeDivergence(vector<VectorX> &gradu, Tube &surface_tube, IBCSubsets &ibc)
{
    // compute divergence of negative normalized gradient of u
    VectorX div = VectorX::Zero(surface_tube.nNodes());
    for(size_t d = 0; d < surface_tube.dim(); ++d)
    { 
        div += m_Dc[d] * gradu[d];
    }

    return div;
}


VectorX GeodesicDistance::Step3(VectorX &div, Tube &surface_tube, IBCSubsets &ibc)
{
    // set RHS for ibc point values to zero (for zero distance to ibc point)
    VectorX divfull(surface_tube.nNodes() + ibc.TotalIBCDOFs());
    for(size_t i = 0; i < surface_tube.nNodes(); ++i)
    {
        divfull[i] = div[i];
    }

    for(size_t c = 0; c < ibc.NumIBCs(); ++c)
    {
        for(size_t i = 0; i < ibc.DOFSubset()[c].size(); ++i)
        {
            const size_t bc_type_idx = ibc.meta(c).is_oriented ? ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
            if(ibc.meta(c).boundary_type[bc_type_idx] == 0)
            {
                if(ibc.meta(c).boundary_order == 2)
                {
                    divfull[ibc.IdentityRowsStart() + ibc.DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + ibc.DirichletIBCColumnStartIndex(c)] = 0.0;
                }
                else
                {
                    divfull[ibc.IdentityRowsStart() + i + ibc.IBCColumnStartIndex(c)] = 0.0;
                }
            }
        }
    }

    
    Solver sv;
    sv.Init((Scalar)0., (Scalar)1., m_Lfull, m_Efull, Preconditioner::diagonal, m_identity_rows);
    sv.SetPrecond(2000);
    cpm::VectorX phi = sv.Solve(divfull);


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
        phi = m_Efull * phi;

        VectorX phi_surface(surface_tube.nNodes());
        VectorX div_surface(surface_tube.nNodes());
        for(int i = 0; i < surface_tube.nNodes(); ++i)
        {
            phi_surface[i] = phi[i];
            div_surface[i] = div[i];
        }
        polyscope::getPointCloud("surface closest points")->addScalarQuantity("Distance", phi_surface);
        polyscope::getPointCloud("surface closest points")->addScalarQuantity("div", div_surface);

        for(size_t c = 0; c < ibc.NumIBCs(); ++c)
        {
            vector<Scalar> phi_ibc(ibc.DOFSubset()[c].size());
            vector<Scalar> div_ibc(ibc.DOFSubset()[c].size());
            for(size_t i = 0; i < ibc.DOFSubset()[c].size(); ++i)
            {
                size_t dof_idx = ibc.IdentityRowsStart() + i + ibc.IBCColumnStartIndex(c);
                phi_ibc[i] = phi[dof_idx];
                div_ibc[i] = div[dof_idx];
            }
            polyscope::getPointCloud("ibc closest points " + to_string(c))->addScalarQuantity("distance", phi_ibc);
            polyscope::getPointCloud("ibc closest points " + to_string(c))->addScalarQuantity("div", div_ibc);
        }
    }
#endif

    return phi;
}


// For testing using exact gradient on sphere with point source
VectorX GeodesicDistance::ExactMatchStep1ScalarField(Tube &surface_tube, IBCSubsets &ibc)
{
    VectorX u(surface_tube.nNodes());
    for(size_t i = 0; i < surface_tube.nNodes(); ++i)
    {
        Scalar dot = DotProduct(surface_tube.cpx()[i], ibc.xp(0)[0]);
        Scalar theta;
        if(dot > 1)
        {
            theta = 0;
        }
        else if(dot < -1)
        {
            theta = M_PI;
        }
        else
        {
            theta = acos(dot);
        }
        u[i] = theta; //sin(0.5 * theta);
    }

    return u;
}


// For testing using exact gradient on sphere with point source
vector<VectorX> GeodesicDistance::ExactGradientSpherePoint(Tube &surface_tube, IBCSubsets &ibc)
{
    vector<VectorX> gradu(surface_tube.dim(), VectorX(surface_tube.nNodes()));
    for(size_t i = 0; i < surface_tube.nNodes(); ++i)
    {
        vector<Scalar> y = surface_tube.cpx()[i];
        vector<Scalar> x = ibc.xp(0)[0];
        vector<Scalar> t = y-x;
        vector<Scalar> n = y;
        vector<Scalar> X =  t - DotProduct(t, n) * n;
        Normalize(X, 0.0);
        for(size_t d = 0; d < surface_tube.dim(); ++d)
        {
            gradu[d][i] = X[d];
        }
    }

    return gradu;
}


// For testing using exact divergence on sphere with point source
VectorX GeodesicDistance::ExactDivergenceSpherePoint(Tube &surface_tube, IBCSubsets &ibc)
{
    VectorX div(surface_tube.nNodes() + ibc.TotalIBCDOFs());
    for(size_t i = 0; i < surface_tube.nNodes(); ++i)
    {
        vector<Scalar> y = surface_tube.cpx()[i];
        vector<Scalar> x = ibc.xp(0)[0];
        Scalar xdoty = DotProduct(x, y);
        Scalar theta = acos(xdoty);

        if(abs(sin(theta)) > 1e-14)
        {
            div[i] = xdoty / sin(theta);
        }
        else
        {
            div[i] = 40;
        }
    }

    return div;
}


VectorX GeodesicDistance::ExactDivergenceCirclePoint(Tube &surface_tube, IBCSubsets &ibc)
{
    // VectorX div(surface_tube.nNodes() + ibc.TotalIBCDOFs());
    // for(size_t i = 0; i < surface_tube.nNodes(); ++i)
    // {
    //     vector<Scalar> y = surface_tube.cpx()[i];
    //     vector<Scalar> x = ibc.xp(0)[0];
    //     Scalar xdoty = DotProduct(x, y);
    //     Scalar theta = acos(xdoty);

    //     if(abs(sin(theta)) > 1e-14)
    //     {
    //         div[i] = 0;
    //     }
    //     else
    //     {
    //         if(xdoty > 0)
    //         {
    //             div[i] = 40;
    //         }
    //         else
    //         {
    //             div[i] = -40;
    //         }
    //     }
    // }

    // return div;
    return VectorX::Zero(surface_tube.nNodes() + ibc.TotalIBCDOFs());
}

} // namespace cpm