#include <chrono>

#ifdef OPENMP
#include <omp.h>
#endif

#include "PoissonWithIBC.h"
#include "Defines.h"
#include "PoissonWithIBCNearest.h"
#include "Solver.h"

using namespace std;
using namespace cpm;

size_t g_boundary_order = 2;
size_t g_interp_deg = 3;


template <typename T>
void VisualizeSolution(T &solver, Scalar dx, cpm::VectorX &uplot, cpm::VectorX &uexact, cpm::VectorX &error)
{
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
}
template void VisualizeSolution<PoissonWithIBC>(PoissonWithIBC &solver, Scalar dx, cpm::VectorX &uplot, cpm::VectorX &uexact, cpm::VectorX &error);
template void VisualizeSolution<PoissonWithIBCNearest>(PoissonWithIBCNearest &solver, Scalar dx, cpm::VectorX &uplot, cpm::VectorX &uexact, cpm::VectorX &error);


template <typename T>
Scalar ErrorAndVisualize(T &solver, Scalar dx, cpm::VectorX &u, function<Scalar(const vector<Scalar> &)> exact)
{
    // compute error with exact solution on the surface
    SpMat plotE = solver.GetPlotInterpMatrix();

    cpm::VectorX uplot = plotE * u;

    cpm::VectorX uexact(solver.xp().size());
    for (size_t i = 0; i < solver.xp().size(); ++i)
    {
        uexact[i] = exact(solver.xp()[i]);
    }

    cpm::VectorX error = uplot - uexact;

    // Visualize with polyscope
#ifdef POLYSCOPE
    VisualizeSolution(solver, dx, uplot, uexact, error);
#endif

    return error.lpNorm<Eigen::Infinity>();
}
template Scalar ErrorAndVisualize<PoissonWithIBC>(PoissonWithIBC &solver, Scalar dx, cpm::VectorX &u, function<Scalar(const vector<Scalar> &)> exact);
template Scalar ErrorAndVisualize<PoissonWithIBCNearest>(PoissonWithIBCNearest &solver, Scalar dx, cpm::VectorX &u, function<Scalar(const vector<Scalar> &)> exact);


Scalar ErrorAndVisualize(PoissonWithIBC &solver, Scalar dx, cpm::VectorX &ufull, function<Scalar(const vector<Scalar> &)> exact, function<Scalar(size_t &, const size_t &, const vector<Scalar> &)> exact_ibc)
{
    // compute error with exact solution on the surface
    vector<vector<size_t>> director_set_index_from_subset_index;
    vector<vector<size_t>> director_which_side;
    vector<vector<bool>> director_is_on_ibc;
    SpMat plotEfull = solver.GetPlotInterpMatrix(director_set_index_from_subset_index, director_which_side, director_is_on_ibc);

    cpm::VectorX uplot = plotEfull * ufull;

    cpm::VectorX uexact(solver.xp().size());
    cpm::VectorX which_side(solver.xp().size());
    for (size_t i = 0; i < solver.xp().size(); ++i)
    {
        uexact[i] = exact(solver.xp()[i]);
        which_side[i] = 0;
    }

    for (size_t c = 0; c < solver.NumIBCs(); ++c)
    {
        for (size_t i = 0; i < director_set_index_from_subset_index[c].size(); ++i)
        {
            which_side[director_set_index_from_subset_index[c][i]] = director_which_side[c][i];
            if (director_is_on_ibc[c][i])
            {
                uexact[director_set_index_from_subset_index[c][i]] = exact_ibc(c, director_which_side[c][i], solver.xp()[director_set_index_from_subset_index[c][i]]);
            }
        }
    }

    cpm::VectorX error = uplot - uexact;

    // Visualize with polyscope
#ifdef POLYSCOPE
    VisualizeSolution(solver, dx, uplot, uexact, error);
#endif

    return error.lpNorm<Eigen::Infinity>();
}

vector<Scalar> ComputeCircle3DParams(vector<Scalar> &normal)
{
    Normalize(normal, 0.0);
    Scalar r = 0.2; // have tried ranges from 0-0.8, must be less than sphere radius of 1
    Scalar theta = 0.1;
    Scalar phi = 0.2;
    vector<Scalar> point_on_plane{r * cos(phi) * sin(theta), r * sin(phi) * sin(theta), r * cos(theta)};
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

Scalar PoissonCirclePointIBC(Scalar ibc_theta, Scalar dx)
{
    SurfaceSpecifier surface_specs("Circle");

    vector<Scalar> ibc_params{cos(ibc_theta), sin(ibc_theta)};
    vector<SurfaceSpecifier> ibc_surface_specs(1);
    ibc_surface_specs[0].SetSurface("Point", surface_specs.boundingBox(), ibc_params, surface_specs.Dim());
    
    IBCMetaData meta;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    meta.boundary_order = g_boundary_order;
    ibc_surface_specs[0].SetIBCMetaData(meta);

    PoissonWithIBC solver(surface_specs, ibc_surface_specs, dx, g_interp_deg);

    auto uexact = [](const vector<Scalar> &x)
    { Scalar theta = atan2(x[1], x[0]); return sin(theta); };
    auto f = [uexact](const vector<Scalar> &x)
    { return uexact(x); };

    Solver sv;
    sv.Init((Scalar)0., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
    cpm::VectorX rhs = solver.GetRHS(f, uexact);
    cpm::VectorX u = sv.Solve(rhs);

    return ErrorAndVisualize(solver, dx, u, uexact);
}

Scalar PoissonCirclePointIBCWithJump(Scalar ibc_theta, Scalar dx)
{
    SurfaceSpecifier surface_specs("Circle");

    vector<Scalar> ibc_params{cos(ibc_theta), sin(ibc_theta)};
    vector<SurfaceSpecifier> ibc_surface_specs(1);
    ibc_surface_specs[0].SetSurface("Point", surface_specs.boundingBox(), ibc_params, surface_specs.Dim());

    Scalar dtheta = 0.001 * 2 * M_PI;
    vector<Scalar> point_on_side1{cos(ibc_theta + dtheta), sin(ibc_theta + dtheta)};

    IBCMetaData meta;
    meta.is_oriented = true;
    meta.point_on_side1 = point_on_side1;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    meta.boundary_order = g_boundary_order;
    ibc_surface_specs[0].SetIBCMetaData(meta);

    PoissonWithIBC solver(surface_specs, ibc_surface_specs, dx, g_interp_deg);

#ifdef POLYSCOPE
    vector<vector<Scalar>> point_on_side11{point_on_side1};
    polyscope::registerPointCloud2D("Point On Side", point_on_side11);
#endif

    Scalar aa = 2;
    Scalar bb = 20;
    auto theta_func = [&ibc_theta](const vector<Scalar> &x)
    {
        Scalar theta = atan2(x[1], x[0]);
        if (theta < 0)
        {
            theta += 2 * M_PI;
        }

        if (theta - ibc_theta < 0)
        {
            return theta - ibc_theta + 2 * M_PI;
        }
        else
        {
            return theta - ibc_theta;
        }
    };
    auto uexact = [&theta_func, &aa, &bb](const vector<Scalar> &x)
    { Scalar theta_minus_theta_c = theta_func(x); return aa * cos(theta_minus_theta_c) + 0.5 * bb * (theta_minus_theta_c) / M_PI; };
    auto f = [&theta_func, &aa, &bb](const vector<Scalar> &x)
    { Scalar theta_minus_theta_c = theta_func(x); return aa * cos(theta_minus_theta_c); };
    vector<Scalar> ibc_exact{aa, aa + bb};
    auto u_ibc = [&solver, &ibc_exact](size_t &ibc_index, const size_t &which_side, const vector<Scalar> &x)
    { return ibc_exact[which_side - 1]; };

    Solver sv;
    sv.Init((Scalar)0., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::jacobi, solver.GetIdentityRows());
    sv.SetPrecond(500);
    cpm::VectorX rhs = solver.GetRHS(f, u_ibc);
    cpm::VectorX u = sv.Solve(rhs);

    return ErrorAndVisualize(solver, dx, u, uexact, u_ibc);
}

Scalar PoissonSphereCircleIBC(vector<Scalar> &normal, Scalar dx, size_t ibc_type)
{
    SurfaceSpecifier surface_specs("Sphere");

    vector<SurfaceSpecifier> ibc_surface_specs(1);
    ibc_surface_specs[0].SetSurface("Circle3D", surface_specs.boundingBox(), ComputeCircle3DParams(normal), surface_specs.Dim());

    IBCMetaData meta;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    meta.boundary_order = g_boundary_order;
    ibc_surface_specs[0].SetIBCMetaData(meta);

    auto uexact = [](const vector<Scalar> &x)
    { Scalar theta = atan2(x[1], x[0]); Scalar phi = acos(x[2] / 1.0); return cos(3*theta) * pow(sin(phi), 3) * (9 * pow(cos(phi),2) - 1); };
    auto f = [uexact](const vector<Scalar> &x)
    { return 30 * uexact(x); };

    if (ibc_type == 0)
    {
        PoissonWithIBC solver(surface_specs, ibc_surface_specs, dx, g_interp_deg);

        Solver sv;
        sv.Init((Scalar)0., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
        cpm::VectorX rhs = solver.GetRHS(f, uexact);
        cpm::VectorX u = sv.Solve(rhs);

        return ErrorAndVisualize(solver, dx, u, uexact);
    }
    else if (ibc_type > 0)
    {
        bool freeze_nearest = (ibc_type == 1);
        PoissonWithIBCNearest solver(surface_specs, ibc_surface_specs, dx, freeze_nearest, g_interp_deg);

        Solver sv;
        sv.Init((Scalar)0., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
        cpm::VectorX rhs = solver.GetRHS(f, uexact);
        cpm::VectorX u = sv.Solve(rhs);

        return ErrorAndVisualize(solver, dx, u, uexact);
    }

    return (cpm::Scalar)0.0;
}

Scalar PoissonSpherePointIBC(Scalar ibc_theta, Scalar ibc_phi, Scalar dx, size_t ibc_type)
{
    SurfaceSpecifier surface_specs("Sphere");

    vector<Scalar> ibc_params{cos(ibc_phi) * cos(ibc_theta), cos(ibc_phi) * sin(ibc_theta), sin(ibc_phi)};
    vector<SurfaceSpecifier> ibc_surface_specs(1);
    ibc_surface_specs[0].SetSurface("Point", surface_specs.boundingBox(), ibc_params, surface_specs.Dim());

    IBCMetaData meta;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    meta.boundary_order = g_boundary_order;
    ibc_surface_specs[0].SetIBCMetaData(meta);

    auto uexact = [](const vector<Scalar> &x)
    { Scalar theta = atan2(x[1], x[0]); Scalar phi = acos(x[2] / 1.0); return cos(3*theta) * pow(sin(phi), 3) * (9 * pow(cos(phi),2) - 1); };
    auto f = [uexact](const vector<Scalar> &x)
    { return 30 * uexact(x); };

    if (ibc_type == 0)
    {
        PoissonWithIBC solver(surface_specs, ibc_surface_specs, dx, g_interp_deg);
        
        Solver sv;
        sv.Init((Scalar)0., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
        cpm::VectorX rhs = solver.GetRHS(f, uexact);
        cpm::VectorX u = sv.Solve(rhs);

        return ErrorAndVisualize(solver, dx, u, uexact);
    }
    else if (ibc_type > 0)
    {
        bool freeze_nearest = (ibc_type == 1);
        PoissonWithIBCNearest solver(surface_specs, ibc_surface_specs, dx, freeze_nearest, g_interp_deg);
        
        Solver sv;
        sv.Init((Scalar)0., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
        cpm::VectorX rhs = solver.GetRHS(f, uexact);
        cpm::VectorX u = sv.Solve(rhs);

        return ErrorAndVisualize(solver, dx, u, uexact);
    }

    return (cpm::Scalar)0.0;
}

Scalar PoissonSphereArcIBC(Scalar theta1, Scalar theta2, Scalar dx)
{
    SurfaceSpecifier surface_specs("Sphere");

    vector<Scalar> ibc_params{1.0, theta1, theta2};
    vector<SurfaceSpecifier> ibc_surface_specs(1);
    ibc_surface_specs[0].SetSurface("Arc", surface_specs.boundingBox(), ibc_params, surface_specs.Dim());

    IBCMetaData meta;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    meta.boundary_order = g_boundary_order;
    ibc_surface_specs[0].SetIBCMetaData(meta);

    PoissonWithIBC solver(surface_specs, ibc_surface_specs, dx, g_interp_deg);

    auto uexact = [&solver](const vector<Scalar> &x)
    { Scalar theta = atan2(x[1], x[0]); Scalar phi = acos(x[2] / solver.surfaceParams()[0]); return cos(3*theta) * pow(sin(phi), 3) * (9 * pow(cos(phi),2) - 1); };
    auto f = [uexact](const vector<Scalar> &x)
    { return 30 * uexact(x); };

    Solver sv;
    sv.Init((Scalar)0., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
    cpm::VectorX rhs = solver.GetRHS(f, uexact);
    cpm::VectorX u = sv.Solve(rhs);

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

    vector<Scalar> circle_point_order;
    vector<Scalar> circle_point_jump_order;
    vector<Scalar> sphere_circle_order;
    vector<Scalar> sphere_point_order;
    vector<Scalar> sphere_arc_order;

    if(g_boundary_order == 1)
    {
        circle_point_order = vector<Scalar>({0.993515, 0.953104, 0.842597, 1.03335, 1.27055, 0.903129, 0.980151, 0.718754, 0.840704, 1.35021});
        circle_point_jump_order = vector<Scalar>({1.00427, 0.964602, 0.809851, 1.00818, 1.28006, 0.90979, 0.981358, 0.874118, 0.865422, 1.34103});
        sphere_circle_order = vector<Scalar>({0.954594, 0.990541, 1.07685, 0.865972, 0.784669, 0.922806, 0.991011, 0.987983, 0.90721, 0.956464});
        sphere_point_order = vector<Scalar>({2.42046, 2.03873, 1.03422, 0.973755, 1.54908, 2.12488, 2.12239, 1.02281, 1.05949, 1.37816});
        sphere_arc_order = vector<Scalar>({1.76665, 1.79881, 1.72063, 1.74305, 1.72427, 1.75405, 1.73974, 1.80824, 1.6929, 1.6667});
    }
    else if (g_boundary_order == 2)
    {
        circle_point_order = vector<Scalar>({2.01419, 1.97624, 2.06569, 2.04698, 2.1271, 2.12353, 1.99891, 1.96309, 1.83585, 2.1617});
        circle_point_jump_order = vector<Scalar>({2.0, 1.90007, 1.89512, 2.10919, 2.05796, 1.91437, 1.95175, 2.00729, 1.92675, 2.0158});
        sphere_circle_order = vector<Scalar>({2.03624, 2.20571, 2.09173, 2.0828, 2.1572, 2.10696, 1.97858, 1.97003, 2.17008, 1.9772});
        sphere_point_order = vector<Scalar>({2.16283, 2.13198, 1.92074, 2.22093, 2.27176, 2.36021, 1.89006, 2.12619, 2.22663, 2.70325});
        sphere_arc_order = vector<Scalar>({1.91054, 1.90547, 1.88648, 1.92365, 1.90556, 1.90494, 1.89661, 1.91465, 1.81442, 1.73931});
    }

    vector<Scalar> sphere_circle_nearest_order{0.945984, 0.82667, 0.98775, 0.703821, 1.18035, 1.12101, 0.867009, 1.09471, 0.962311, 0.786675};
    vector<Scalar> sphere_circle_nearest_bw_order{0.626507, 0.64894, 0.659936, 0.804202, 0.608543, 0.790434, 0.892918, 0.639867, 0.666589, 0.727743};
    vector<Scalar> sphere_point_nearest_order{2.15041, 2.23765, 0.959175, 0.911574, 1.21328, 2.55332, 1.55977, 1.03094, 1.15633, 1.03559};
    vector<Scalar> sphere_point_nearest_bw_order{2.15603, 1.17371, 0.763936, 0.639966, 0.838463, 0.896607, 0.640536, 0.59289, 0.304371, 1.15492};

    vector<vector<Scalar>> test_avg_conv_order(9);
    test_avg_conv_order[0] = circle_point_order;
    test_avg_conv_order[1] = circle_point_jump_order;
    test_avg_conv_order[2] = sphere_circle_order;
    test_avg_conv_order[3] = sphere_circle_nearest_order;
    test_avg_conv_order[4] = sphere_circle_nearest_bw_order;
    test_avg_conv_order[5] = sphere_point_order;
    test_avg_conv_order[6] = sphere_point_nearest_order;
    test_avg_conv_order[7] = sphere_point_nearest_bw_order;
    test_avg_conv_order[8] = sphere_arc_order;

    // 10 different ibc points around the circle
    vector<Scalar> rand_num{0.0, 0.170, 0.297, 0.331, 0.481, 0.511, 0.631, 0.735, 0.813, 0.973};
    vector<Scalar> ibc_theta(rand_num.size());
    for (size_t i = 0; i < rand_num.size(); ++i)
    {
        ibc_theta[i] = 2.0 * M_PI * rand_num[i];
    }

    // add 10 different ibc points for the endpoints on the arc on the sphere
    vector<Scalar> rand_num2{0.25, 0.874, 0.473, 0.973, 0.812, 0.693, 0.873, 0.912, 0.945, 1.0};
    vector<Scalar> ibc_arc_endpoint1(rand_num.size());
    vector<Scalar> ibc_arc_endpoint2(rand_num2.size());
    for (size_t i = 0; i < rand_num.size(); ++i)
    {
        ibc_arc_endpoint1[i] = ibc_theta[i] - M_PI;
        ibc_arc_endpoint2[i] = 2.0 * M_PI * rand_num2[i] - M_PI;
    }

    // add 10 more points to define the normal for the circle on the sphere
    vector<Scalar> rand_num3{0.55, 0.123, -0.743, 0.547, -0.239, -0.693, 0.731, -0.175, 0.453, -0.364};
    vector<vector<Scalar>> circle_normal(rand_num.size(), vector<Scalar>(3));
    for (size_t i = 0; i < rand_num.size(); ++i)
    {
        circle_normal[i][0] = rand_num[i];
        circle_normal[i][1] = rand_num2[i];
        circle_normal[i][2] = rand_num3[i];
    }

    bool all_passed = true;
    for(size_t problem = 0; problem < 9; ++problem)
    {
#ifdef POLYSCOPE
        polyscope::removeAllStructures();
#endif 
        cout << endl;
        cout << "Problem " << problem << endl;
        cout << endl;

        for (size_t i = 0; i < rand_num.size(); ++i)
        {
            // make a lambda for convergence study, we need function format Scalar f(Scalar)
            Scalar theta = ibc_theta[i];
            Scalar endpoint1 = ibc_arc_endpoint1[i];
            Scalar endpoint2 = ibc_arc_endpoint2[i];
            vector<Scalar> normal = circle_normal[i];
            auto Poisson = [&theta, &endpoint1, &endpoint2, &normal, &problem](Scalar &dx)
            {
                switch (problem)
                {
                case 0:
                    return PoissonCirclePointIBC(theta, dx);
                    break;
                case 1:
                    return PoissonCirclePointIBCWithJump(theta, dx);
                    break;
                case 2:
                    return PoissonSphereCircleIBC(normal, dx, 0);
                    break;
                case 3:
                    return PoissonSphereCircleIBC(normal, dx, 1);
                    break;
                case 4:
                    return PoissonSphereCircleIBC(normal, dx, 2);
                    break;
                case 5: //
                    return PoissonSpherePointIBC(endpoint1, endpoint2, dx, 0);
                    break;
                case 6: // pinning nearest point
                    return PoissonSpherePointIBC(endpoint1, endpoint2, dx, 1);
                    break;
                case 7: // Auer
                    return PoissonSpherePointIBC(endpoint1, endpoint2, dx, 2);
                    break;
                case 8:
                    return PoissonSphereArcIBC(endpoint1, endpoint2, dx);
                    break;
                default:
                    cout << "Problem type does not exist" << endl;
                    return (cpm::Scalar)0.0;
                }
            };

            chrono::time_point<chrono::system_clock> start, end;
            start = chrono::system_clock::now();

            Helpers h;
            Scalar avg_conv_order;
            if (problem <= 1) // 2D tests
            {
                avg_conv_order = h.ConvergenceStudy(0.1, 5, Poisson);
            }
            else // 3D tests
            {
                avg_conv_order = h.ConvergenceStudy(0.2, 3, Poisson); 
            }
            if (abs(avg_conv_order - test_avg_conv_order[problem][i]) < 1e-2)
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
