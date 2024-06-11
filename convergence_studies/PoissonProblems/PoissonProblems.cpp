#include <chrono>

#include "Poisson.h"
#include "Defines.h"
#include "PoissonWithIBCNearest.h"
#include "Solver.h"

#ifdef OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace cpm;

Scalar ErrorAndVisualize(Poisson &solver, Scalar dx, cpm::VectorX &u, function<Scalar(const Scalar&, const Scalar&)> exact)
{
    // compute error with exact solution on the surface
    SpMat plotE = solver.GetPlotInterpMatrix();
    cpm::VectorX uplot = plotE * u;

    cpm::VectorX uexact(solver.thetap().size());
    for(size_t i = 0; i < solver.thetap().size(); ++i)
    {
        if(solver.dim() == 2)
        {
            uexact[i] = exact(solver.thetap()[i], 0 /* dummy */);
        }
        else if(solver.dim() == 3)
        {
            uexact[i] = exact(solver.thetap()[i], solver.phip()[i]);
        }
    }

    cpm::VectorX error = uplot - uexact;
    
    // Visualize with polyscope
#ifdef POLYSCOPE
    string surfaceName = "Surface dx = " + to_string(dx);
    if(solver.dim() == 2)
    {  
        polyscope::registerCurveNetworkLine2D(surfaceName, solver.xp())->addNodeScalarQuantity("uplot", uplot)->setEnabled(true);
        polyscope::getCurveNetwork(surfaceName)->addNodeScalarQuantity("uexact", uexact);
        polyscope::getCurveNetwork(surfaceName)->addNodeScalarQuantity("error", error);
    }
    else if(solver.dim() == 3)
    {
        polyscope::registerSurfaceMesh(surfaceName, solver.xp(), solver.faces())->addVertexScalarQuantity("uplot", uplot)->setEnabled(true);
        polyscope::getSurfaceMesh(surfaceName)->addVertexScalarQuantity("uexact", uexact);
        polyscope::getSurfaceMesh(surfaceName)->addVertexScalarQuantity("error", error);
    }
#endif

    return error.lpNorm<Eigen::Infinity>();
}


Scalar ArcDirichletBCs(Scalar dx, size_t bdy_order, vector<size_t>& which_bdy)
{
    vector<Scalar> bounding_box{-1.0, 1.0};
    vector<Scalar> surface_params{1.0, -3 * M_PI_4, M_PI_4};
    SurfaceSpecifier surface_specs("Arc", bounding_box, surface_params, 2 /* embedding dim */);

    Poisson solver(surface_specs, dx, bdy_order, which_bdy);
    
    auto f = [](const Scalar& th, const Scalar &dummy){ return sin(th); };
    auto uexact = [&f](const Scalar& th, const Scalar &dummy){ return -f(th, 0); };

    Solver sv;
    sv.Init((Scalar)0., (Scalar)1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
    cpm::VectorX rhs = solver.GetRHS(f, uexact);
    cpm::VectorX u = sv.Solve(rhs);
    
    return ErrorAndVisualize(solver, dx, u, uexact);
}


Scalar ArcNeumannAndDirichletBCs(Scalar dx, size_t bdy_order, vector<size_t>& which_bdy)
{
    vector<Scalar> bounding_box{-1.0, 1.0};
    vector<Scalar> surface_params{1.0, 0.0, M_PI};
    SurfaceSpecifier surface_specs("Arc", bounding_box, surface_params, 2 /* embedding dim */);

    Poisson solver(surface_specs, dx, bdy_order, which_bdy);
    
    auto f = [](const Scalar& th, const Scalar &dummy){ return cos(th); };
    auto uexact = [&f](const Scalar& th, const Scalar &dummy){ return -f(th, 0); };
    
    Solver sv;
    sv.Init((Scalar)0., (Scalar)1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
    cpm::VectorX rhs = solver.GetRHS(f, uexact);
    cpm::VectorX u = sv.Solve(rhs);

    return ErrorAndVisualize(solver, dx, u, uexact);
}


Scalar ShiftedPoissonCircle(Scalar dx)
{
    SurfaceSpecifier surface_specs("Circle");
    Poisson solver(surface_specs, dx);
    
    auto uexact = [](const Scalar& th, const Scalar& dummy){ return sin(th) + sin(12 * th); };
    auto f = [](const Scalar& th, const Scalar& dummy){ return 2 * sin(th) + 145 * sin(12 * th); };

    Solver sv;
#ifdef EBCG
    sv.Init((Scalar)1., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
#else
    sv.Init((Scalar)1., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
#endif
    cpm::VectorX rhs = solver.GetRHS(f);
    cpm::VectorX u = sv.Solve(rhs);

    return ErrorAndVisualize(solver, dx, u, uexact);
}



Scalar ShiftedPoissonSphere(Scalar dx)
{
    SurfaceSpecifier surface_specs("Sphere");
    Poisson solver(surface_specs, dx);
    
    auto uexact = [](const Scalar& th, const Scalar& phi){ return cos(3*th) * pow(sin(phi), 3) * (9 * pow(cos(phi),2) - 1); };
    auto f = [](const Scalar& th, const Scalar& phi){ return ( -29 * (cos(3*th) * pow(sin(phi), 3) * (9 * pow(cos(phi),2) - 1) )); };

    Solver sv;
#ifdef EBCG
    sv.Init((Scalar)1., (Scalar)1., solver.GetLfull(), solver.GetEfull(), Preconditioner::jacobi, solver.GetIdentityRows());
#else
    sv.Init((Scalar)1., (Scalar)1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
#endif
    cpm::VectorX rhs = solver.GetRHS(f);
    cpm::VectorX u = sv.Solve(rhs);

    return ErrorAndVisualize(solver, dx, u, uexact);
}


Scalar PoissonHemisphere(Scalar dx, size_t bdy_order, vector<size_t>& which_bdy)
{
    SurfaceSpecifier surface_specs("Hemisphere");
    Poisson solver(surface_specs, dx, bdy_order, which_bdy);
    
    auto uexact = [](const Scalar& th, const Scalar& phi){ return cos(3*th) * pow(sin(phi), 3) * (9 * pow(cos(phi),2) - 1); };
    auto f = [](const Scalar& th, const Scalar& phi){ return ( -30 * (cos(3*th) * pow(sin(phi), 3) * (9 * pow(cos(phi),2) - 1) )); };

    Solver sv;
    sv.Init((Scalar)0., (Scalar)1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
    cpm::VectorX rhs = solver.GetRHS(f, uexact);
    cpm::VectorX u = sv.Solve(rhs);

    return ErrorAndVisualize(solver, dx, u, uexact);
}


int main(int argc, char** argv)
{
#ifdef OPENMP
    omp_set_dynamic(0);
	omp_set_num_threads(8);
#endif

#ifdef POLYSCOPE
    polyscope::init();
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
#endif

    vector<Scalar> test_avg_conv_order{1.15925, 1.00101, 1.89298, 2.01658, 1.9876, 2.18548, 2.11767};

    bool all_passed = true;
    for(size_t problem = 0; problem < 7; ++problem)
    {
#ifdef POLYSCOPE
        polyscope::removeAllStructures();
#endif 
        cout << endl;
        cout << "Problem " << problem << endl;
        cout << endl;

        chrono::time_point<chrono::system_clock> start, end;
        start = chrono::system_clock::now();

        // make a lambda for convergence study, we need function format Scalar f(Scalar)
        auto Poisson = [&problem](Scalar& dx)
        {
            vector<size_t> which_bdy;
            switch(problem)
            {
                case 0:
                    which_bdy.push_back(1); which_bdy.push_back(2);
                    return ArcDirichletBCs(dx, 1, which_bdy);
                    break;
                case 1:
                    which_bdy.push_back(2);
                    return ArcNeumannAndDirichletBCs(dx, 1, which_bdy);
                    break;
                case 2:
                    which_bdy.push_back(1); which_bdy.push_back(2);
                    return ArcDirichletBCs(dx, 2, which_bdy);
                    break;
                case 3:
                    which_bdy.push_back(2);
                    return ArcNeumannAndDirichletBCs(dx, 2, which_bdy);
                    break;
                case 4:
                    return ShiftedPoissonCircle(dx);
                    break;
                case 5:
                    return ShiftedPoissonSphere(dx);
                    break;
                case 6:
                    which_bdy.push_back(1);
                    return PoissonHemisphere(dx, 1, which_bdy);
                    break;
                default:
                    cout << "Problem type does not exist" << endl;
                    return (cpm::Scalar)0.0;
            } 
        };

        Helpers h;
        Scalar avg_conv_order;
        if(problem <= 4) // 2D tests
        {
            avg_conv_order = h.ConvergenceStudy(0.05, 5, Poisson);
        }
        else // 3D tests
        {
            avg_conv_order = h.ConvergenceStudy(0.2, 3, Poisson);
        }

        if(abs(avg_conv_order - test_avg_conv_order[problem]) < 1e-2)
        {
            cout << "PASSED: average convergence order test" << endl;
        }
        else
        {
            all_passed = false;
            cout << "FAILED: average convergence order test" << endl;
        }

        end = chrono::system_clock::now();
        chrono::duration<Scalar> elapsed_seconds = end - start;
        cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

#ifdef POLYSCOPE
        polyscope::show();
#endif
    }

    if(all_passed)
    {
        cout << "ALL PASSED :)" << endl;
    }
    else
    {
        cout << "AT LEAST ONE FAILED :(" << endl;
    }

#ifdef CUSTOM_SOLVER
    cout << "solved via matrix-free solver" << "\n";
#endif
#ifdef ESL
    cout << "solved via Eigen::SparseLU" << "\n";
#endif
#ifdef EBCG
    cout << "solved via Eigen::BiCGStab" << "\n";
#endif

    return 0;
}
