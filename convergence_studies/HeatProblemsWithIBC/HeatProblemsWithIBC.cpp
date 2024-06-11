#include <chrono>

#ifdef OPENMP
#include <omp.h>
#endif

#include "Defines.h"
#include "HeatWithIBC.h"
#include "Helpers.h"
#include "IBCMetaData.h"

using namespace std;
using namespace cpm;

Helpers h;

size_t g_interp_deg = 3;

Scalar g_final_time = 0.1;

size_t g_boundary_type = 0;
size_t g_boundary_order = 2;


Scalar ErrorAndVisualize(HeatWithIBC &solver, Scalar dx, cpm::VectorX &ufull, function<Scalar(const Scalar &, const vector<Scalar> &)> exact)
{
    // compute error with exact solution on the surface
    SpMat plotEfull = solver.GetPlotInterpMatrix();

    cpm::VectorX uplot = plotEfull * ufull;

    cpm::VectorX uexact(solver.xp().size());
    for (size_t i = 0; i < solver.xp().size(); ++i)
    {
        uexact[i] = exact(g_final_time, solver.xp()[i]);
    }

    cpm::VectorX error = uplot - uexact;

    // Visualize with polyscope
#ifdef POLYSCOPE
    string surfaceName = "Surface dx = " + to_string(dx);
    if (solver.dim() == 2)
    {
        polyscope::registerCurveNetworkLine2D(surfaceName, solver.xp())->addNodeScalarQuantity("uplot", uplot)->setEnabled(true);
        polyscope::getCurveNetwork(surfaceName)->addNodeScalarQuantity("uexact", uexact);
        polyscope::getCurveNetwork(surfaceName)->addNodeScalarQuantity("error", error);
        for (size_t c = 0; c < solver.NumIBCs(); ++c)
        {
            polyscope::registerPointCloud2D("IBC Points " + to_string(c), solver.ibc_xp(c))->setPointColor(glm::vec3(1,1,1))->setPointRadius(0.01);
        }
    }
    else if (solver.dim() == 3)
    {
        polyscope::registerSurfaceMesh(surfaceName, solver.xp(), solver.faces())->addVertexScalarQuantity("uplot", uplot)->setEnabled(true);
        polyscope::getSurfaceMesh(surfaceName)->addVertexScalarQuantity("uexact", uexact);
        polyscope::getSurfaceMesh(surfaceName)->addVertexScalarQuantity("error", error);
        for (size_t c = 0; c < solver.NumIBCs(); ++c)
        {
            polyscope::registerCurveNetworkLine("IBC Curve", solver.ibc_xp(c))->setColor(glm::vec3(1,1,1))->setRadius(0.01);
        }
    }
#endif

    return error.lpNorm<Eigen::Infinity>();
}

vector<Scalar> ComputeCircle3DParams(vector<Scalar> &normal, vector<Scalar> &point_on_plane)
{
    Normalize(normal, 0.0);
    // See https://math.stackexchange.com/questions/943383/determine-circle-of-intersection-of-plane-and-sphere
    Scalar rho = -DotProduct(point_on_plane, normal);
    vector<Scalar> centre = rho * normal;
    Scalar circle_radius = sqrt(1 - rho * rho); // assuming a sphere of radius 1 centred at [0,0,0]

    vector<Scalar> ibc_params;
    ibc_params.push_back(circle_radius);
    ibc_params.push_back(centre[0]);
    ibc_params.push_back(centre[1]);
    ibc_params.push_back(centre[2]);
    ibc_params.push_back(normal[0]);
    ibc_params.push_back(normal[1]);
    ibc_params.push_back(normal[2]);

    return ibc_params;
}

Scalar SphereCircleIBC(Scalar dx, size_t time_stepping_method)
{
    SurfaceSpecifier surface_specs("Sphere");

    vector<Scalar> normal{0, 0, 0};
    vector<Scalar> point_on_plane{0, 0, 0};
    if (g_boundary_type == 0)
    {
        normal = vector<Scalar>({1.1, 1.234, 1.321});
        Scalar r = 0.27; // have tried ranges from 0-0.8, must be less than sphere radius of 1
        Scalar theta = 0.132;
        Scalar phi = 0.221;
        point_on_plane = vector<Scalar>({r * cos(phi) * sin(theta), r * sin(phi) * sin(theta), r * cos(theta)});
    }
    else
    {
        normal = vector<Scalar>({0, 1, 0});
        point_on_plane = vector<Scalar>({0, 0, 0});
    }
    
    vector<SurfaceSpecifier> ibc_surface_specs(1);
    vector<Scalar> ibc_params = ComputeCircle3DParams(normal, point_on_plane);
    ibc_surface_specs[0].SetSurface("Circle3D", surface_specs.boundingBox(), ibc_params, surface_specs.Dim());

    IBCMetaData meta;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    meta.boundary_type[0] = g_boundary_type;
    meta.boundary_order = g_boundary_order;
    ibc_surface_specs[0].SetIBCMetaData(meta);

    Scalar dt;
    size_t num_time_steps;
    if (time_stepping_method > 0)
    {
        dt = 0.1 * dx;
    }
    else
    {
        dt = 0.1 * dx * dx;
    }
    num_time_steps = ceil(g_final_time / dt);
    dt = g_final_time / num_time_steps; // adjust dt so we have an integer number of time steps

    HeatWithIBC solver(surface_specs, ibc_surface_specs, dx, dt, time_stepping_method, g_interp_deg);
    
    Scalar current_time;
    auto f = [&solver](const vector<Scalar> &x)
    { Scalar theta = atan2(x[1], x[0]); Scalar phi = acos(x[2] / solver.surfaceParams()[0]); return cos(phi); };
    auto uexact = [&f](const Scalar &t, const vector<Scalar> &x)
    { return exp(-2 * t) * f(x); };
    auto uexact_ibc = [&f, &current_time](size_t &ibc_index, const size_t &which_side, const vector<Scalar> &x)
    { return exp(-2 * current_time) * f(x); };

    cpm::VectorX u(solver.nNodes() + solver.ibc().TotalIBCDOFs() + solver.ibc().Num2ndOrderDirichletRows());
#ifdef OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < solver.nNodes(); ++i)
    {
        u[i] = f(solver.cpx()[i]);
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
                u[i + solver.ibc().IBCColumnStartIndex(c)] = f(solver.ibc().DOFSubset()[c].cpx()[i]);
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
                    u[solver.ibc().IdentityRowsStart() + solver.ibc().DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + solver.ibc().DirichletIBCColumnStartIndex(c)] = uexact(0.0, solver.ibc().DOFSubset()[c].cpx()[i]);
                }
                else
                {
                    u[solver.ibc().IdentityRowsStart() + i + solver.ibc().IBCColumnStartIndex(c)] = uexact(0.0, solver.ibc().DOFSubset()[c].cpx()[i]);
                }
            }
        }
    }
    for (size_t t = 0; t < num_time_steps; ++t)
    {
        current_time = (t + 1) * dt;
        u = solver.TimeStep(u, uexact_ibc);
    }
    return ErrorAndVisualize(solver, dx, u, uexact);
}

Scalar Circle(Scalar ibc_theta, Scalar dx, size_t time_stepping_method)
{
    SurfaceSpecifier surface_specs("Circle");

    vector<SurfaceSpecifier> ibc_surface_specs(1);
    vector<Scalar> ibc_params{cos(ibc_theta), sin(ibc_theta)};
    ibc_surface_specs[0].SetSurface("Point", surface_specs.boundingBox(), ibc_params, surface_specs.Dim());

    IBCMetaData meta;
    meta.boundary_type[0] = g_boundary_type;
    meta.boundary_order = g_boundary_order;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    ibc_surface_specs[0].SetIBCMetaData(meta);

    Scalar dt;
    size_t num_time_steps;
    if (time_stepping_method > 0)
    {
        dt = 0.1 * dx;
    }
    else
    {
        dt = 0.1 * dx * dx;
    }
    num_time_steps = ceil(g_final_time / dt);
    dt = g_final_time / num_time_steps; // adjust dt so we have an integer number of time steps

    HeatWithIBC solver(surface_specs, ibc_surface_specs, dx, dt, time_stepping_method, g_interp_deg);

    Scalar current_time;
    auto f = [&meta](const vector<Scalar> &x)
    {
        Scalar theta = atan2(x[1], x[0]);
        Scalar result;
        if (meta.boundary_type[0] == 0)
        {
            result = sin(theta);
        }
        else
        {
            result = cos(theta); // need cos(theta) for Neumann because we can only implement zero Neumann now, ibc_theta must also stay zero here. If you use cos(theta) for Dirichlet case, you get 2nd order for both 1st and 2nd order implementations. I hypothesize this is because there is a zero Neumann condition at that point also, which helps the cp extension accuracy?
        }
        return result;
    };

    auto uexact = [&f](const Scalar &t, const vector<Scalar> &x)
    { return exp(-t) * f(x); };
    auto uexact_ibc = [&f, &current_time](size_t &ibc_index, const size_t &which_side, const vector<Scalar> &x)
    { return exp(-current_time) * f(x); };

    cpm::VectorX u(solver.nNodes() + solver.ibc().TotalIBCDOFs() + solver.ibc().Num2ndOrderDirichletRows());
    for (size_t i = 0; i < solver.nNodes(); ++i)
    {
        u[i] = f(solver.cpx()[i]);
    }

    for (size_t c = 0; c < solver.ibc().NumIBCs(); ++c)
    {
        for (size_t i = 0; i < solver.ibc().DOFSubset()[c].size(); ++i)
        {
            if (solver.ibc().meta(c).boundary_type[0] == 1 || (solver.ibc().meta(c).boundary_type[0] == 0 && solver.ibc().meta(c).boundary_order == 2))
            {
                u[i + solver.ibc().IBCColumnStartIndex(c)] = f(solver.ibc().DOFSubset()[c].cpx()[i]);
            }
        }
    }

    for (size_t c = 0; c < solver.ibc().NumIBCs(); ++c)
    {
        for (size_t i = 0; i < solver.ibc().DOFSubset()[c].size(); ++i)
        {
            const size_t bc_type_idx = solver.ibc().meta(c).is_oriented ? solver.ibc().DOFSubset()[c].which_side()[i] - 1 : 0;
            if (solver.ibc().meta(c).boundary_type[bc_type_idx] == 0)
            {
                if (solver.ibc().meta(c).boundary_order == 2)
                {
                    u[solver.ibc().IdentityRowsStart() + solver.ibc().DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + solver.ibc().DirichletIBCColumnStartIndex(c)] = uexact(0.0, solver.ibc().DOFSubset()[c].cpx()[i]);
                }
                else
                {
                    u[solver.ibc().IdentityRowsStart() + i + solver.ibc().IBCColumnStartIndex(c)] = uexact(0.0, solver.ibc().DOFSubset()[c].cpx()[i]);
                }
            }
        }
    }

    for (size_t t = 0; t < num_time_steps; ++t)
    {
        current_time = (t + 1) * dt;
        u = solver.TimeStep(u, uexact_ibc);
    }

    return ErrorAndVisualize(solver, dx, u, uexact);
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

    vector<Scalar> test_avg_conv_order;
    if (g_boundary_type == 0)
    {
        if (g_boundary_order == 1)
        {
            test_avg_conv_order = vector<Scalar>({0.970592, 0.970211, 0.970927, 0.969203, 0.96295, 0.970163});
        }
        else if (g_boundary_order == 2)
        {
            test_avg_conv_order = vector<Scalar>({1.9959, 1.07669, 2.02602, 1.98885, 1.15183, 2.01012});
        }
    }
    else if (g_boundary_type == 1)
    {
        if (g_boundary_order == 1)
        {
            test_avg_conv_order = vector<Scalar>({0.98058, 0.976024, 0.978244, 0.946114, 0.93853, 0.944366});
        }
        else if (g_boundary_order == 2)
        {
            test_avg_conv_order = vector<Scalar>({1.99499, 1.07708, 2.0274, 2.07915, 1.15527, 2.04239});
        }
    }

    bool all_passed = true;
    for (size_t problem = 0; problem < 6; ++problem)
    {
#ifdef POLYSCOPE
        polyscope::removeAllStructures();
#endif 
        chrono::time_point<chrono::system_clock> start, end;
        start = chrono::system_clock::now();

        // make a lambda for convergence study, we need function format Scalar f(Scalar)
        auto Heat = [&problem](Scalar &dx)
        {
            switch (problem)
            {
            case 0:
                return Circle(0.0, dx, 0);
                break;
            case 1:
                return Circle(0.0, dx, 1);
                break;
            case 2:
                return Circle(0.0, dx, 2);
                break;
            case 3:
                return SphereCircleIBC(dx, 0);
                break;
            case 4:
                return SphereCircleIBC(dx, 1);
                break;
            case 5:
                return SphereCircleIBC(dx, 2);
                break;
            default:
                cout << "Problem type does not exist" << endl;
                return (cpm::Scalar)0.0;
            }
        };

        Scalar avg_conv_order;
        if (problem <= 2) // 2D tests
        {
            avg_conv_order = h.ConvergenceStudy(0.05, 5, Heat);
        }
        else // 3D tests
        {
            avg_conv_order = h.ConvergenceStudy(0.1, 3, Heat); 
        }

        end = chrono::system_clock::now();
        chrono::duration<Scalar> elapsed_seconds = end - start;
        if (abs(avg_conv_order - test_avg_conv_order[problem]) < 1e-2)
        {
            cout << "PASSED: average convergence order test" << endl;
        }
        else
        {
            all_passed = false;
            cout << "FAILED: average convergence order test" << endl;
        }
        cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
        
#ifdef POLYSCOPE
            polyscope::show();
#endif
    }

    if (all_passed)
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
