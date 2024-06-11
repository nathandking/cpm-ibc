#include <chrono>

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef POLYSCOPE
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#endif

#include "GeodesicDistanceNearest.h"

using namespace std;
using namespace cpm;


void GeoDist(Scalar dx)
{
    SurfaceSpecifier surface_specs("../assets/beetle.obj", "../assets/beetle_coarse.obj");
    Tube surface_tube(surface_specs, 2 /* interp deg */, dx);

    size_t num_bc = 11;
    vector<string> bc_file(num_bc);
    bc_file[0] = "../assets/beetle_boundary.obj";
    for(size_t c = 1; c < num_bc; ++c)
    {
        bc_file[c] = "../assets/beetle_boundary" + to_string(c) + ".obj";
    }

    vector<SurfaceSpecifier> bc_surface_specs(num_bc);
    vector<vector<size_t>> constraint_idx(num_bc);
    vector<vector<size_t>> inside_constraint_idx(num_bc);
    vector<vector<Scalar>> inside_constraint_idx_dist(num_bc);
    for(size_t c = 0; c < num_bc; ++c)
    {
        bc_surface_specs[c].SetSurface(bc_file[c]);
        Surface bc_surf(bc_surface_specs[c]);

        for(size_t i = 0; i < surface_tube.nNodes(); ++i)
        {   vector<Scalar> cpc(surface_tube.dim());
            bc_surf.ClosestPoint(surface_tube.x()[i], cpc);

            if(Distance(cpc, surface_tube.cpx()[i]) < 1e-6)
            {
                constraint_idx[c].push_back(i);
            }
            else if(Distance(surface_tube.x()[i], cpc) < surface_tube.TubeRadius())
            {
                inside_constraint_idx[c].push_back(i);
                inside_constraint_idx_dist[c].push_back(Distance(surface_tube.x()[i], cpc));
            }
        }
    }

    GeodesicDistanceNearest solver(false);
    Scalar heat_dt_scale = 1.0;
    Scalar dt = heat_dt_scale * dx * dx;
    cpm::VectorX uplot = solver.ComputeDistance(surface_tube, constraint_idx, inside_constraint_idx, inside_constraint_idx_dist, dt);

#ifdef POLYSCOPE
    polyscope::registerSurfaceMesh("Surface", surface_specs.Mesh().vertices, surface_specs.Mesh().faces)->setSmoothShade(true);
    polyscope::getSurfaceMesh("Surface")->addVertexDistanceQuantity("Distance", uplot)->setEnabled(true);

    for(size_t c = 0; c < num_bc; ++c)
    {
        polyscope::registerCurveNetwork("BC"+to_string(c), bc_surface_specs[c].Mesh().vertices, bc_surface_specs[c].Mesh().faces)->setColor(glm::vec3{0,0,0})->setRadius(0.004);
    }
#endif
}


int main(int argc, char** argv)
{
#ifdef OPENMP
    int threads = 8;
    omp_set_dynamic(0);
	omp_set_num_threads(threads);
#endif

#ifdef POLYSCOPE
    polyscope::init(); 
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
#endif
    
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    Scalar dx = 0.003125;
    dx = 0.0125;
    GeoDist(dx);

    end = chrono::system_clock::now();
    chrono::duration<Scalar> elapsed_seconds = end - start;
    cout << "total elapsed time: " << elapsed_seconds.count() << "s\n";

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
