#include <chrono>

#ifdef OPENMP
#include <omp.h>
#endif

#include "Defines.h"
#include "ReactionDiffusionWithIBC.h"

using namespace std;
using namespace cpm;

size_t interp_deg = 2;

Scalar g_final_time = 10000;

void Visualize(ReactionDiffusionWithIBC &solver, Scalar dx, cpm::VectorX &ufull, cpm::VectorX &vfull)
{
    // compute error with exact solution on the surface
    SpMat plotEfull = solver.GetPlotInterpMatrix();

    cpm::VectorX uplot = plotEfull * ufull;
    cpm::VectorX vplot = plotEfull * vfull;

    // Visualize with polyscope
#ifdef POLYSCOPE
    string Meshname = "Mesh dx=" + to_string(dx);
    polyscope::registerSurfaceMesh(Meshname, solver.xp(), solver.faces())->addVertexScalarQuantity("u", uplot)->setEnabled(true);
    polyscope::getSurfaceMesh(Meshname)->addVertexScalarQuantity("v", vplot);

    for (size_t c = 0; c < solver.NumIBCs(); ++c)
    {
        polyscope::registerCurveNetworkLine("IBC Curve " + to_string(c), solver.ibc_xp(c))->setColor(glm::vec3(1,1,1))->setRadius(0.003);
    }
#endif
}


void ReactionDiffusion(Scalar dx)
{
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Read surface and IBCs
    //////////////////////////////////////////////////////////////////////////////////////////////
    SurfaceSpecifier surface_specs("../assets/fish_high_res_scaled_centered.obj");

    vector<string> ibc_files{"../assets/fish_eye.obj", "../assets/fish_back_fin.obj", "../assets/fish_tail.obj"};        
    vector<vector<size_t>> ibc_type{{0,0}, {0,1}, {0,0}};

    size_t num_ibcs = 3;
    vector<SurfaceSpecifier> ibc_surface_specs(num_ibcs);
    for(size_t c = 0; c < num_ibcs; ++c)
    {
        ibc_surface_specs[c].SetSurface(ibc_files[c]);

        IBCMetaData meta;
        meta.on_ibc_tolerance = 0.1 * dx * dx;
        meta.boundary_type = ibc_type[c];
        meta.is_oriented = true;
        ibc_surface_specs[c].SetIBCMetaData(meta);
    }

    Scalar dt;
    size_t num_time_steps;

    Scalar nu_u = 1.0 / pow(3.0 / dx, 2);
    Scalar nu_v = nu_u / 3.0;
    dt = 0.1 * ( 1 / max(nu_u,nu_v) ) * dx * dx;
    num_time_steps = ceil(g_final_time / dt);
    dt = g_final_time / num_time_steps; // adjust dt so we have an integer number of time steps

    ReactionDiffusionWithIBC solver(surface_specs, ibc_surface_specs, dx, dt, interp_deg);

    auto uv0 = [&solver](const vector<Scalar> &x)
    {   
        Scalar pert = 0.5 * exp( -pow(10.0 * (x[2] - 0.1), 2) ) + 0.5 * (rand() % 1000) / 1000.0;
        vector<Scalar> uv{1.0 - pert, 0.5 * pert};
        return uv;
    };

    auto uexact_ibc = [](size_t &ibc_index, const size_t &which_side, const vector<Scalar> &x)
    { return 0.0; };

    auto vexact_ibc = [](size_t &ibc_index, const size_t &which_side, const vector<Scalar> &x)
    { return 0.0; };

    cpm::VectorX u(solver.nNodes() + solver.ibc().TotalIBCDOFs() + solver.ibc().Num2ndOrderDirichletRows());
    cpm::VectorX v(solver.nNodes() + solver.ibc().TotalIBCDOFs() + solver.ibc().Num2ndOrderDirichletRows());

#ifdef OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < solver.nNodes(); ++i)
    {
        vector<Scalar> uv = uv0(solver.cpx()[i]);
        u[i] = uv[0];
        v[i] = uv[1];
    }
#ifdef OPENMP
#pragma omp parallel for
#endif
    for (size_t c = 0; c < solver.ibc().NumIBCs(); ++c)
    {
        for (size_t i = 0; i < solver.ibc().DOFSubset()[c].size(); ++i)
        {
            if (solver.ibc().meta(c).boundary_type[0] == 1 || (solver.ibc().meta(c).boundary_type[0] == 0 && solver.ibc().meta(c).boundary_order == 2))
            {
                u[i + solver.ibc().IBCColumnStartIndex(c)] = uexact_ibc(c, solver.ibc().DOFSubset()[c].which_side()[i], solver.ibc().DOFSubset()[c].cpx()[i]);
                v[i + solver.ibc().IBCColumnStartIndex(c)] = vexact_ibc(c, solver.ibc().DOFSubset()[c].which_side()[i], solver.ibc().DOFSubset()[c].cpx()[i]);
            }
        }
    }
#ifdef OPENMP
#pragma omp parallel for
#endif
    for (size_t c = 0; c < solver.ibc().NumIBCs(); ++c)
    {
        for (size_t i = 0; i < solver.ibc().DOFSubset()[c].size(); ++i)
        {
            const size_t bc_type_idx = solver.ibc().meta(c).is_oriented ? solver.ibc().DOFSubset()[c].which_side()[i] - 1 : 0;
            if (solver.ibc().meta(c).boundary_type[bc_type_idx] == 0)
            {
                if (solver.ibc().meta(c).boundary_order == 2)
                {
                    u[solver.ibc().IdentityRowsStart() + solver.ibc().DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + solver.ibc().DirichletIBCColumnStartIndex(c)] = uexact_ibc(c, solver.ibc().DOFSubset()[c].which_side()[i], solver.ibc().DOFSubset()[c].cpx()[i]);
                    v[solver.ibc().IdentityRowsStart() + solver.ibc().DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + solver.ibc().DirichletIBCColumnStartIndex(c)] = vexact_ibc(c, solver.ibc().DOFSubset()[c].which_side()[i], solver.ibc().DOFSubset()[c].cpx()[i]);
                }
                else
                {
                    u[solver.ibc().IdentityRowsStart() + i + solver.ibc().IBCColumnStartIndex(c)] = uexact_ibc(c, solver.ibc().DOFSubset()[c].which_side()[i], solver.ibc().DOFSubset()[c].cpx()[i]);
                    v[solver.ibc().IdentityRowsStart() + i + solver.ibc().IBCColumnStartIndex(c)] = vexact_ibc(c, solver.ibc().DOFSubset()[c].which_side()[i], solver.ibc().DOFSubset()[c].cpx()[i]);
                }
            }
        }
    }

    for (size_t t = 0; t < num_time_steps; ++t)
    {
        cout << "current time step = " << t << " with total time steps = " << num_time_steps << endl;
        solver.TimeStep(u, v, uexact_ibc, vexact_ibc);
    }

    Visualize(solver, dx, u, v);
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
    polyscope::view::lookAt(glm::vec3{1.7, 0.75, 0.4}, glm::vec3{0, -0.15, 0}); // right
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    Scalar dx = 0.01;
    ReactionDiffusion(dx);

    end = chrono::system_clock::now();
    chrono::duration<Scalar> elapsed_seconds = end - start;
    
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
