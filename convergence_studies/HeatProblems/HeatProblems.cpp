#include <chrono>

#include "Heat.h"

#ifdef POLYSCOPE
#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include "polyscope/surface_mesh.h"
#endif

using namespace std;
using namespace cpm;

Helpers h;

Scalar g_final_time = 0.1;

#define USE_IMPLICIT
#ifdef USE_IMPLICIT
    bool g_use_implicit = true;
#else
    bool g_use_implicit = false;
#endif


// for 2D
Scalar ErrorAndVisualize(Heat &solver, Scalar dx, cpm::VectorX &u, function<Scalar(const Scalar &, const Scalar&)> exact)
{
    // compute error with exact solution on the surface
    SpMat plotE = solver.GetPlotInterpMatrix();
    cpm::VectorX uplot = plotE * u;

    cpm::VectorX uexact(solver.thetap().size());
    for(size_t i = 0; i < solver.thetap().size(); ++i)
    {
        uexact[i] = exact(g_final_time, solver.thetap()[i]);
    }

    cpm::VectorX error = uplot - uexact;
    
    // Visualize with polyscope
#ifdef POLYSCOPE
    string surfaceName = "Surface dx = " + to_string(dx);
    polyscope::registerCurveNetworkLine2D(surfaceName, solver.xp())->addNodeScalarQuantity("uplot", uplot)->setEnabled(true);
    polyscope::getCurveNetwork(surfaceName)->addNodeScalarQuantity("uexact", uexact);
    polyscope::getCurveNetwork(surfaceName)->addNodeScalarQuantity("error", error);
#endif

    return error.lpNorm<Eigen::Infinity>();
}


// for 3D
Scalar ErrorAndVisualize3D(Heat &solver, Scalar dx, cpm::VectorX &u, function<Scalar(const Scalar&, const Scalar&, const Scalar&)> exact)
{
    // compute error with exact solution on the surface
    SpMat plotE = solver.GetPlotInterpMatrix();
    cpm::VectorX uplot = plotE * u;

    cpm::VectorX uexact(solver.thetap().size());
    for(size_t i = 0; i < solver.thetap().size(); ++i)
    {
        uexact[i] = exact(g_final_time, solver.thetap()[i], solver.phip()[i]);
    }

    cpm::VectorX error = uplot - uexact;
    
    // Visualize with polyscope
#ifdef POLYSCOPE
    string surfaceName = "Surface dx = " + to_string(dx);
    polyscope::registerSurfaceMesh(surfaceName, solver.xp(), solver.faces())->addVertexScalarQuantity("uplot", uplot)->setEnabled(true);
    polyscope::getSurfaceMesh(surfaceName)->addVertexScalarQuantity("uexact", uexact);
    polyscope::getSurfaceMesh(surfaceName)->addVertexScalarQuantity("error", error);
#endif

    return error.lpNorm<Eigen::Infinity>();
}


Scalar ArcNeumannAndDirichletBCs(Scalar dx, size_t bdy_order, vector<size_t>& which_bdy)
{
    vector<Scalar> bounding_box{-1.0, 1.0};
    vector<Scalar> surface_params{1.0, 0.0, M_PI};
    SurfaceSpecifier surface_specs("Arc", bounding_box, surface_params, 2 /* embedding dimension */);

    Heat solver(surface_specs, dx, g_final_time, g_use_implicit, bdy_order, which_bdy);
    
    auto f = [](const Scalar& th){ return cos(th); };
    auto uexact = [&f](const Scalar& t, const Scalar& th){ return exp(-t) * f(th); };

    cpm::VectorX u(solver.tube().nNodes());
    for(size_t i = 0; i < solver.tube().nNodes(); ++i)
    {
        Scalar theta = atan2(solver.cpx()[i][1], solver.cpx()[i][0]);
        u[i] = f(theta);
    }

    solver.TimeStep(u, uexact);

    return ErrorAndVisualize(solver, dx, u, uexact);
}


Scalar Circle(Scalar dx)
{
    SurfaceSpecifier surface_specs("Circle");
    Heat solver(surface_specs, dx, g_final_time, g_use_implicit);
    
    auto f = [](const Scalar& th){ return sin(th); };
    auto uexact = [&f](const Scalar& t, const Scalar& th){ return exp(-t) * f(th); };

    cpm::VectorX u(solver.tube().nNodes());
    for(size_t i = 0; i < solver.tube().nNodes(); ++i)
    {
        Scalar theta = atan2(solver.cpx()[i][1], solver.cpx()[i][0]);
        u[i] = f(theta);
    }

    solver.TimeStep(u);

    return ErrorAndVisualize(solver, dx, u, uexact);
}


Scalar Sphere(Scalar dx)
{
    SurfaceSpecifier surface_specs("Sphere");
    Heat solver(surface_specs, dx, g_final_time, g_use_implicit);
    
    auto f = [](const Scalar& th, const Scalar& phi){ return cos(phi); };
    auto uexact = [&f](const Scalar& t, const Scalar& th, const Scalar& phi){ return exp(-2*t) * f(th, phi); };


    cpm::VectorX u(solver.tube().nNodes());
    for(size_t i = 0; i < solver.tube().nNodes(); ++i)
    {
        Scalar theta = atan2(solver.cpx()[i][1], solver.cpx()[i][0]);
        Scalar phi = acos(solver.cpx()[i][2] / solver.surfaceParams()[0]);
        u[i] = f(theta, phi);
    }

    solver.TimeStep(u);

    return ErrorAndVisualize3D(solver, dx, u, uexact);
}


int main(int argc, char** argv)
{
#ifdef POLYSCOPE
    polyscope::init();
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
#endif

    #ifdef USE_IMPLICIT
        vector<Scalar> test_avg_conv_order{0.966198, 2.02785, 2.01011, 2.08643};
    #else
        vector<Scalar> test_avg_conv_order{0.965015, 2.29892, 2.02327, 2.39477};
    #endif

    bool all_passed = true;
    for(size_t problem = 0; problem < 4; ++problem)
    {
#ifdef POLYSCOPE
        polyscope::removeAllStructures();
#endif 
        chrono::time_point<chrono::system_clock> start, end;
        start = chrono::system_clock::now();

        // make a lambda for convergence study, we need function format Scalar f(Scalar)
        auto Heat = [&problem](Scalar& dx)
        {
            vector<size_t> which_bdy;
            switch(problem)
            {
                case 0:
                    which_bdy.push_back(2);
                    return ArcNeumannAndDirichletBCs(dx, 1, which_bdy);
                    break;
                case 1:
                    which_bdy.push_back(2);
                    return ArcNeumannAndDirichletBCs(dx, 2, which_bdy);
                    break;
                case 2:
                    return Circle(dx);
                    break;
                case 3:
                    return Sphere(dx);
                    break;
                default:
                    cout << "Problem type does not exist" << endl;
                    return (cpm::Scalar)0.0;
            } 
        };

        Scalar avg_conv_order;
        if(problem <= 2) // 2D tests
        {
            avg_conv_order = h.ConvergenceStudy(0.05, 3, Heat);
        }
        else // 3D tests
        {
            avg_conv_order = h.ConvergenceStudy(0.2, 3, Heat);
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
