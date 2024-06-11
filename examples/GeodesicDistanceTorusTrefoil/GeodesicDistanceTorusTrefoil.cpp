#include <chrono>

#ifdef OPENMP
#include <omp.h>
#endif

#include "GeodesicDistance.h"

using namespace std;
using namespace cpm;

void GeoDistTorusTrefoil(Scalar dx)
{
    SurfaceSpecifier surface_specs("Torus");
    Tube surface_tube(surface_specs, 2 /* interp deg */, dx);

    vector<SurfaceSpecifier> ibc_surface_specs(1);
    vector<Scalar> ibc_params{2.0, 2.0, 3.0};
    ibc_surface_specs[0].SetSurface("Torus Knot", surface_specs.boundingBox(), ibc_params, surface_specs.Dim());

    IBCMetaData meta;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    ibc_surface_specs[0].SetIBCMetaData(meta);

    IBCSubsets ibc;
    ibc.Initialize(surface_tube, ibc_surface_specs);

    GeodesicDistance solver(false /* visualize internal steps */);
    cpm::VectorX uplot = solver.ComputeDistance(surface_tube, ibc, dx * dx /* dt */);

#ifdef POLYSCOPE
    polyscope::registerSurfaceMesh("Surface", surface_tube.surface().xp(), surface_tube.surface().faces())->setSmoothShade(true);
    polyscope::getSurfaceMesh("Surface")->addVertexDistanceQuantity("Distance", uplot)->setEnabled(true);
    polyscope::registerCurveNetworkLoop("IBC", ibc.surface(0).xp())->setColor(glm::vec3{0,0,0})->setRadius(0.004);
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
    polyscope::view::lookAt(glm::vec3{0, 0, -8.25}, glm::vec3{0, -0.25, 0});
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    Scalar dx = 0.0125;
    dx =        0.00625;
    dx = 0.025;
    GeoDistTorusTrefoil(dx);

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
