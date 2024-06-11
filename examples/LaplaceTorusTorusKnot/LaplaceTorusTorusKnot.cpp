#include <chrono>

#ifdef OPENMP
#include <omp.h>
#endif

#include "PoissonWithIBC.h"
#include "Defines.h"
#include "Helpers.h"
#include "Solver.h"

using namespace std;
using namespace cpm;

size_t g_boundary_order = 2;
size_t interp_deg = 3;


void Laplace(Scalar dx)
{
    vector<Scalar> params{3.0, 1.0};
    vector<Scalar> bounding_box{-4.0, 4.0};
    SurfaceSpecifier surface_specs("Torus", bounding_box, params, 3);
    vector<SurfaceSpecifier> ibc_surface_specs(1);
    vector<Scalar> ibc_params{3.0, 3.0, 7.0};
    ibc_surface_specs[0].SetSurface("Torus Knot", surface_specs.boundingBox(), ibc_params, surface_specs.Dim());

    IBCMetaData meta;
    meta.boundary_order = g_boundary_order;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    ibc_surface_specs[0].SetIBCMetaData(meta);

    PoissonWithIBC solver(surface_specs, ibc_surface_specs, dx, interp_deg);

    auto bc = [&solver](const vector<Scalar> &x){ Scalar theta = solver.IBCParameter(0, x); return sin(theta); };
    auto f = [](const vector<Scalar> &x){ return 0; };

    Solver sv;
    sv.Init((Scalar)0., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
    cpm::VectorX rhs = solver.GetRHS(f, bc);
    cpm::VectorX u = sv.Solve(rhs);

    // visualize result
    SpMat plotEfull = solver.GetPlotInterpMatrix();
    cpm::VectorX uplot = plotEfull * u;

    // Visualize with polyscope
    cpm::VectorX u_ibc(solver.ibc_xp(0).size());
    for (size_t i = 0; i < solver.ibc_xp(0).size(); ++i)
    {
        u_ibc[i] = bc(solver.ibc_xp(0)[i]);
    }

#ifdef POLYSCOPE
    string mesh_name = "Plotting Points dx=" + to_string(dx);
    polyscope::registerSurfaceMesh(mesh_name, solver.xp(), solver.faces());
    polyscope::getSurfaceMesh(mesh_name)->addVertexScalarQuantity("uplot", uplot)->setEnabled(true);
    polyscope::registerCurveNetworkLine("IBC Curve", solver.ibc_xp(0))->addNodeScalarQuantity("u", u_ibc)->setEnabled(true);
#endif
}


int main(int argc, char **argv)
{
#ifdef OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(8);
#endif

#ifdef POLYSCOPE
    polyscope::init(); 
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    Scalar dx = 0.05;
    Laplace(dx);
    
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

#ifdef POLYSCOPE
    polyscope::show();
#endif

    return 0;
}
