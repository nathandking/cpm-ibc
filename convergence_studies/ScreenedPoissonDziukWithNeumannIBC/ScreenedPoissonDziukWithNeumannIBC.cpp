#include <chrono>

#ifdef OPENMP
#include <omp.h>
#endif

#include "PoissonWithIBC.h"
#include "Solver.h"

#ifdef POLYSCOPE
#include "polyscope/point_cloud.h"
#endif

using namespace std;
using namespace cpm;

size_t interp_deg = 3;

auto uexact = [](const std::vector<cpm::Scalar>& x)
{ 
    return x[0] * x[1];
};

auto f = [](const std::vector<cpm::Scalar>& p)
{
   Scalar x = p[0]; Scalar y = p[1]; Scalar z = p[2];
   return p[0] * p[1] + (2*y*(-32*pow(x, 5)*pow(z, 4) - 16*pow(x, 5)*z * z - 32*pow(x, 4)*y * y*pow(z, 4) - 16*pow(x, 4)*y * y*z * z + 16*pow(x, 4)*pow(z, 6) + 48*pow(x, 4)*pow(z, 4) + 8*pow(x, 4)*z * z + 16*pow(x, 3)*y * y*pow(z, 6) - 12*pow(x, 3)*y * y*z * z + 112*pow(x, 3)*pow(z, 8) - 8*pow(x, 3)*pow(z, 6) + 46*pow(x, 3)*pow(z, 4) + 34*pow(x, 3)*z * z + 3*pow(x, 3) - 32*x * x*pow(y, 4)*pow(z, 4) - 16*x * x*pow(y, 4)*z * z + 112*x * x*y * y*pow(z, 8) - 64*x * x*y * y*pow(z, 6) - 12*x * x*y * y*pow(z, 4) + 16*x * x*y * y*z * z - 144*x * x*pow(z, 10) - 216*x * x*pow(z, 8) - 122*x * x*pow(z, 6) - 18*x * x*pow(z, 4) + 9*x * x*z * z + 2*x * x - 48*x*pow(y, 4)*pow(z, 6) - 16*x*pow(y, 4)*pow(z, 4) + 4*x*pow(y, 4)*z * z - 144*x*y * y*pow(z, 10) - 64*x*y * y*pow(z, 8) + 144*x*y * y*pow(z, 6) + 134*x*y * y*pow(z, 4) + 38*x*y * y*z * z + 3*x*y * y + 48*x*pow(z, 12) + 200*x*pow(z, 10) + 230*x*pow(z, 8) + 36*x*pow(z, 6) - 90*x*pow(z, 4) - 46*x*z * z - 6*x + 48*pow(y, 4)*pow(z, 8) + 48*pow(y, 4)*pow(z, 6) + 4*pow(y, 4)*pow(z, 4) + 48*y * y*pow(z, 12) + 80*y * y*pow(z, 10) - 48*y * y*pow(z, 8) - 138*y * y*pow(z, 6) - 42*y * y*pow(z, 4) - 3*y * y*z * z - 40*pow(z, 12) - 106*pow(z, 10) - 56*pow(z, 8) + 64*pow(z, 6) + 38*pow(z, 4) + 5*z * z)) / pow(4*y * y*z * z + 4*x*z * z - 4*z * z - 1, 3);
};


Scalar ErrorAndVisualize(PoissonWithIBC &solver, Scalar dx, cpm::VectorX &ufull, function<Scalar(const vector<Scalar>&)> exact)
{
    // compute error with exact solution on the surface
    SpMat plotEfull = solver.GetPlotInterpMatrix();

    cpm::VectorX uplot = plotEfull * ufull;

    cpm::VectorX uexact(solver.xp().size());
    for(size_t i = 0; i < solver.xp().size(); ++i)
    {
        uexact[i] = exact(solver.xp()[i]);
    }

    cpm::VectorX error = uplot - uexact;

    // Visualize with polyscope
#ifdef POLYSCOPE
    string surfaceName = "Surface dx = " + to_string(dx);
    polyscope::registerSurfaceMesh(surfaceName, solver.xp(), solver.faces())->addVertexScalarQuantity("uplot", uplot)->setEnabled(true);
    polyscope::getSurfaceMesh(surfaceName)->addVertexScalarQuantity("uexact", uexact);
    polyscope::getSurfaceMesh(surfaceName)->addVertexScalarQuantity("error", error);
    for(size_t c = 0; c < solver.NumIBCs(); ++c)
    {
        polyscope::registerCurveNetworkLine("IBC", solver.ibc_xp(c))->setColor(glm::vec3(1,1,1))->setRadius(0.01);
    }
#endif

    return error.lpNorm<Eigen::Infinity>();
}


Scalar PoissonCPMSolver(Scalar dx, size_t ibc_type, size_t problem)
{
#ifdef TRACK_WHERE
    std::cout << "PoissonCPMSolver" << std::endl;
#endif
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<Scalar> elapsed_seconds;

    SurfaceSpecifier surface_specs("Dziuk");

    vector<SurfaceSpecifier> ibc_surface_specs(1);
    if(problem == 0)
    {
        ibc_surface_specs[0].SetSurface("Circle3D");
    }
    else
    {
        vector<Scalar> ibc_params{1, -3*M_PI_4, M_PI_4};
        ibc_surface_specs[0].SetSurface("Arc", surface_specs.boundingBox(), ibc_params, surface_specs.Dim());
    }

    IBCMetaData meta;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    meta.boundary_order = ibc_type + 1;
    meta.boundary_type = {1, 1};
    ibc_surface_specs[0].SetIBCMetaData(meta);

    PoissonWithIBC solver(surface_specs, ibc_surface_specs, dx, interp_deg);

    Solver sv;
    sv.Init((Scalar)1., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
    cpm::VectorX rhs = solver.GetRHS(f, uexact);
    cpm::VectorX u = sv.Solve(rhs);
    
    Scalar error = ErrorAndVisualize(solver, dx, u, uexact);

    return error;
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
    polyscope::view::upDir = polyscope::UpDir::ZUp;
    polyscope::view::lookAt(glm::vec3{-0.1, -2.9, 1.5}, glm::vec3{0, -0.45, 0}); // right
#endif
    
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    vector<string> ibc_name{"Circle", "Arc"};
    vector<string> ibc_type_name{"First-Order IBC", "Second-Order IBC"};


    for(size_t problem = 0; problem < 2; ++problem)
    {
        cout << " " << endl;
        cout << "IBC = " << ibc_name[problem] << endl; 
        cout << " " << endl;
        
        for(size_t ibc_type = 0; ibc_type < 2; ++ibc_type)
        {
            cout << " " << endl;
            cout << "CPM IBC Type = " << ibc_type_name[ibc_type] << endl; 
            cout << " " << endl;

            auto PoissonCPM = [&ibc_type, &problem](Scalar& dx)
            {
                Scalar error = PoissonCPMSolver(dx, ibc_type, problem);
                return error;
            };

            Helpers h;
            Scalar avg_conv_order = h.ConvergenceStudy(0.1, 3, PoissonCPM);

#ifdef POLYSCOPE
            polyscope::show();
#endif
        }
    }
    end = chrono::system_clock::now();
    chrono::duration<Scalar> elapsed_seconds = end - start;
    cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

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
