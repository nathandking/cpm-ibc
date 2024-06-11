#include <chrono>

#ifdef OPENMP
#include <omp.h>
#endif

#include "GeodesicDistance.h"

using namespace std;
using namespace cpm;


void GeoDist(Scalar dx)
{
    SurfaceSpecifier surface_specs("../assets/Dziuk_high_high_res.obj");
    
    size_t source_idx = 58588;
    vector<SurfaceSpecifier> ibc_surface_specs(1);
    ibc_surface_specs[0].SetSurface("Point", surface_specs.boundingBox(), surface_specs.Mesh().vertices[source_idx], surface_specs.Dim());

    IBCMetaData meta;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    ibc_surface_specs[0].SetIBCMetaData(meta);

    Tube surface_tube(surface_specs, 2 /* interp deg */, dx);

    IBCSubsets ibc;
    ibc.Initialize(surface_tube, ibc_surface_specs);

    GeodesicDistance solver(false /* visualize internal steps */);
    cpm::VectorX uplot = solver.ComputeDistance(surface_tube, ibc, dx * dx /* dt */);

#ifdef POLYSCOPE
    polyscope::registerSurfaceMesh("Surface", surface_specs.Mesh().vertices, surface_specs.Mesh().faces)->setSmoothShade(true);
    polyscope::getSurfaceMesh("Surface")->addVertexDistanceQuantity("Distance", uplot)->setEnabled(true);
    polyscope::registerPointCloud("IBC", ibc_surface_specs[0].Mesh().vertices)->setPointColor(glm::vec3{0,0,0})->setPointRadius(0.004);
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
    polyscope::view::lookAt(glm::vec3{-0.1, -2.9, 1.5}, glm::vec3{0, -0.45, 0});
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::view::upDir = polyscope::UpDir::ZUp;
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    Scalar dx = 0.0125;
    dx = 0.05;
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
