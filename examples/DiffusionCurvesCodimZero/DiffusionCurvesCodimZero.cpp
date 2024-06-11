#include <chrono>

#include "PoissonWithIBC.h"
#include "Solver.h"

using namespace std;
using namespace cpm;


void DiffusionCurves(Scalar dx, bool visualize = true)
{
    SurfaceSpecifier surface_specs("Plane");

    vector<Scalar> ibc_params{3.0, 4.0, 0.5};
    vector<SurfaceSpecifier> ibc_surface_specs(1);
    ibc_surface_specs[0].SetSurface("Planar Curve", surface_specs.boundingBox(), ibc_params, surface_specs.Dim());
    
    IBCMetaData meta;
    meta.is_oriented = true;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    meta.use_cp_diff_directions = true;
    ibc_surface_specs[0].SetIBCMetaData(meta);

    PoissonWithIBC solver(surface_specs, ibc_surface_specs, dx, 2 /* interp deg */);

    auto f = [](const vector<Scalar>& x){ return 0; };

    vector<cpm::VectorX> color(3);
    for(size_t color_channel = 0; color_channel < 3; ++color_channel)
    {
        auto color_ibc = [&solver, color_channel](size_t &ibc_index, const size_t &which_side, const vector<Scalar>& x)
        { 
            Scalar s = solver.IBCParameter(ibc_index, x) / (2 * M_PI);
            Scalar t = (1.5 * M_PI + solver.IBCParameter(ibc_index, x)) / (2 * M_PI);
            if(t > 1.0)
            {
                t = t - 1.0;
            }
            vector<vector<Scalar>> c(2, vector<Scalar>(3));
            c[0][0] = 1.0 - s;
            c[0][1] = s;
            c[0][2] = 0.0;
            c[1][0] = 0.0;
            c[1][1] = t;
            c[1][2] = 1.0-t;
            return c[which_side - 1][color_channel]; 
        };

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

#ifdef POLYSCOPE
    polyscope::registerSurfaceMesh2D("Surface", solver.xp(), solver.faces());
    polyscope::getSurfaceMesh("Surface")->addVertexColorQuantity("Colors", colors)->setEnabled(true);
    polyscope::registerCurveNetwork2D("IBC " + to_string(0), solver.ibc_xp(0), solver.ibc_faces(0))->setRadius(0.005)->setColor(glm::vec3{1.0, 1.0, 1.0});
#endif
}



int main(int argc, char** argv)
{
#ifdef POLYSCOPE
    polyscope::init(); 
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;    
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    Scalar dx = 0.005;
    DiffusionCurves(dx);
        
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
