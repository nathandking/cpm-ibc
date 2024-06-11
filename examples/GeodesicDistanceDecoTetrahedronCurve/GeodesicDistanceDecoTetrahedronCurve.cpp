#include <chrono>

#ifdef OPENMP
#include <omp.h>
#endif

#include "GeodesicDistance.h"

using namespace std;
using namespace cpm;

Helpers h;

void GeoDist(Scalar dx)
{
    SurfaceSpecifier surface_specs("DecoTetrahedron");
    Tube surface_tube(surface_specs, 2 /* interp deg */, dx);

    vector<SurfaceSpecifier> ibc_surface_specs(1);
    ibc_surface_specs[0].SetSurface("../assets/Decotetrahedron_curve.obj");
    IBCMetaData meta;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    ibc_surface_specs[0].SetIBCMetaData(meta);

    IBCSubsets ibc;
    ibc.Initialize(surface_tube, ibc_surface_specs);

    // surface_tube.OverideParameterizationPoints(surface_tube.cpx());

    GeodesicDistance geo_dist;
    cpm::VectorX uplot = geo_dist.ComputeDistance(surface_tube, ibc, dx * dx /* dt */);

#ifdef POLYSCOPE
    polyscope::registerPointCloud("Surface", surface_tube.surface().xp())->setPointRadius(0.002);
    polyscope::getPointCloud("Surface")->addScalarQuantity("Distance", uplot)->setColorMap("blues")->setIsolinesEnabled(true)->setEnabled(true);
    polyscope::registerCurveNetworkLine("IBC", ibc.DOFSubset()[0].cpx())->setColor(glm::vec3{0,0,0})->setRadius(0.004)->setEnabled(true);
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
    polyscope::view::lookAt(glm::vec3{-6.25, 6.25, -6.25}, glm::vec3{0, -0.75, 0});
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    Scalar dx = 0.05;
    GeoDist(dx);

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
