#include <chrono>

#include "PoissonWithIBC.h"
#include "Solver.h"

using namespace std;
using namespace cpm;


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
    Scalar circle_radius = sqrt(1.25 * 1.25 - rho * rho); // assuming a sphere of radius 1.25 centred at [0,0,0] 
 
    vector<Scalar> ibc_params;
    ibc_params.push_back(circle_radius);
    ibc_params.push_back(centre[0]); ibc_params.push_back(centre[1]); ibc_params.push_back(centre[2]);
    ibc_params.push_back(normal[0]); ibc_params.push_back(normal[1]); ibc_params.push_back(normal[2]);

    return ibc_params;
}


void DiffusionCurves(vector<vector<vector<Scalar>>> &ibc_colors, Scalar dx)
{
    SurfaceSpecifier surface_specs("Torus-Line-Sphere");

    vector<string> ibc_surface{"Torus Knot", "Circle3D"};
    vector<Scalar> circle_normal{1.0, 1.0, 0.0};
    vector<Scalar> circle_params = ComputeCircle3DParams(circle_normal);
    vector<vector<Scalar>> ibc_params{{3.0, 3.0, 7.0}, circle_params};

    vector<SurfaceSpecifier> ibc_surface_specs(ibc_surface.size());
    for(size_t c = 0; c < ibc_surface.size(); ++c)
    {
        ibc_surface_specs[c].SetSurface(ibc_surface[c], surface_specs.boundingBox(), ibc_params[c], surface_specs.Dim());

        IBCMetaData meta;
        meta.is_oriented = true;
        meta.on_ibc_tolerance = 0.1 * dx * dx;

        ibc_surface_specs[c].SetIBCMetaData(meta);
    }

    PoissonWithIBC solver(surface_specs, ibc_surface_specs, dx, 2 /* interp deg */);

    auto f = [](const vector<Scalar>& x){ return 0; };

    vector<cpm::VectorX> color(3);
    for(size_t color_channel = 0; color_channel < 3; ++color_channel)
    {
        auto color_ibc = [&solver, &ibc_colors, color_channel](size_t &ibc_index, const size_t &which_side, const vector<Scalar>& x){ return /*sin(M_PI * solver.IBCParameter(ibc_index, x)) * */ ibc_colors[ibc_index][which_side - 1][color_channel]; };
        
        Solver sv;
        sv.Init((Scalar)0., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
        cpm::VectorX rhsfull = solver.GetRHS(f, color_ibc);
        cpm::VectorX ufull = sv.Solve(rhsfull);

        SpMat plotEfull = solver.GetPlotInterpMatrix();

        color[color_channel] = plotEfull * ufull;
    }

    vector<vector<Scalar>> colors(color[0].size(), vector<Scalar>(3));
    for(size_t i = 0; i < color[0].size(); ++i)
    {
        for(size_t color_channel = 0; color_channel < 3; ++color_channel)
        {
            colors[i][color_channel] = color[color_channel][i];
        }
    }

    // extract vertices for just the lines
    size_t Np = 250000;
    size_t Nline = 1000;
    vector<vector<Scalar>> l1(Nline, vector<Scalar>(3));
    vector<vector<Scalar>> l2(Nline, vector<Scalar>(3));
    vector<vector<Scalar>> lcolors1(Nline, vector<Scalar>(3));
    vector<vector<Scalar>> lcolors2(Nline, vector<Scalar>(3));
    for(size_t i = 0; i < Nline; ++i)
    {
        l1[i] = solver.xp()[2 * Np + i];
        l2[i] = solver.xp()[2 * Np + i + Nline];
        for(size_t color_channel = 0; color_channel < 3; ++color_channel)
        {
            lcolors1[i][color_channel] = color[color_channel][i + 2 * Np];
            lcolors2[i][color_channel] = color[color_channel][i + 2 * Np + Nline];
        }
    }

#ifdef POLYSCOPE
    polyscope::registerSurfaceMesh("Surface", solver.xp(), solver.faces());
    polyscope::getSurfaceMesh("Surface")->addVertexColorQuantity("u", colors)->setEnabled(true);

    polyscope::registerCurveNetworkLine("Line 1", l1)->addNodeColorQuantity("Color", lcolors1)->setEnabled(true);
    polyscope::registerCurveNetworkLine("Line 2", l2)->addNodeColorQuantity("Color", lcolors2)->setEnabled(true);

    for(size_t c = 0; c < solver.NumIBCs(); ++c)
    {
        polyscope::registerCurveNetwork("IBC " + to_string(c), solver.ibc_xp(c), solver.ibc_faces(c))->setRadius(0.00262)->setColor(glm::vec3{1, 1, 1});
    }
#endif
}



int main(int argc, char** argv)
{
#ifdef POLYSCOPE
    polyscope::init(); 
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Colors
    //////////////////////////////////////////////////////////////////////////////////////////////
    vector<Scalar> blue{0.0,0.447,0.741};  
    vector<Scalar> red{0.8196, 0.1333, 0.0392};
    vector<Scalar> green{0.655, 0.949, 0.267};
    vector<Scalar> yellow{1, 1, 0.0};
    vector<vector<vector<Scalar>>> ibc_colors{{red, blue}, {yellow, green}};

    Scalar dx = 0.05;
    DiffusionCurves(ibc_colors, dx);
        
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
