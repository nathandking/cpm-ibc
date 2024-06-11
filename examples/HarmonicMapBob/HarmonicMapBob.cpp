#include <chrono>

#include "HeatWithIBC.h"
#include "SimplePolygonMesh.h"

using namespace std;
using namespace cpm;

size_t interp_deg = 2;

vector<vector<Scalar>> red;
vector<vector<Scalar>> green;
vector<vector<Scalar>> blue;


vector<Scalar> GetColor(const vector<Scalar> &uv)
{
    size_t rows = red.size();
    size_t cols = red[0].size();

    size_t period = 2;

    // find pixel uv coordinate is in image
    size_t xpixel = floor(uv[0] * period * cols);
    size_t ypixel = floor(uv[1] * period * rows);
    xpixel %= cols;
    ypixel %= rows;

    return vector<Scalar>{(Scalar)red[xpixel][ypixel] / 255.0, (Scalar)green[xpixel][ypixel] / 255.0, (Scalar)blue[xpixel][ypixel] / 255.0};
}


void TargetClosestPoint(Surface& target_surface, vector<cpm::VectorX> &ufull)
{
    vector<Scalar> cp(3);        
    for(size_t i = 0; i < ufull[0].size(); ++i)
    {
        vector<Scalar> u{ufull[0][i], ufull[1][i], ufull[2][i]};
        
        target_surface.ClosestPoint(u, cp);
        
        ufull[0][i] = cp[0];
        ufull[1][i] = cp[1];
        ufull[2][i] = cp[2];
    }
}


size_t HarmonicMap(vector<vector<Scalar>> &initial_map, vector<vector<Scalar>> &map_colors, SurfaceSpecifier source_specs, Surface& target_surface, vector<SurfaceSpecifier> &source_ibc_specs, vector<SurfaceSpecifier> &target_ibc_specs, Scalar dx)
{
    // get ibc surfaces on target surface
    size_t num_ibcs = target_ibc_specs.size();
    vector<Surface> target_ibc_surface(num_ibcs);
    for(size_t c = 0; c < num_ibcs; ++c)
    {
        target_ibc_surface[c].SetSurfaceSpecs(target_ibc_specs[c]);
    }

    // construct Heat flow solver on source surface
    size_t time_stepping_method = 0; // explicit Euler
    HeatWithIBC solver(source_specs, source_ibc_specs, dx, 0.1 * dx * dx, time_stepping_method, interp_deg);

    // function for Dirichlet values on the ibcs, the Dirichlet values are the coordinates on the target surface ibcs
    size_t d;
    auto ibc_location = [&solver, &d, &target_ibc_surface](size_t &ibc_index, const size_t &which_side, const vector<Scalar>& x){ vector<Scalar> coord = target_ibc_surface[ibc_index].CoordinateFromIBCParameter(solver.IBCParameter(ibc_index, x)); return coord[d]; };

    // Compute harmonic map using CPM harmonic mapping https://arxiv.org/pdf/1710.09655.pdf
    vector<cpm::VectorX> ufull(3, cpm::VectorX::Zero(solver.nNodes() + solver.ibc().TotalIBCDOFs() + solver.ibc().IdentityRowsStart()));
    for(d = 0; d < 3; ++d)
    {
        for(size_t i = 0; i < initial_map.size(); ++i)
        {
            ufull[d][i] = initial_map[i][d];   
        }

        for(size_t c = 0; c < solver.ibc().NumIBCs(); ++c)
        {
            for(size_t i = 0; i < solver.ibc().DOFSubset()[c].size(); ++i)
            {
                if(source_ibc_specs[c].IBCMeta().boundary_type[0] == 0)
                {
                    ufull[d][solver.ibc().IdentityRowsStart() + solver.ibc().IBCColumnStartIndex(c)] = ibc_location(c, solver.ibc().DOFSubset()[c].which_side()[i], solver.ibc().DOFSubset()[c].cpx()[i]);
                }
            }
        }
    }

    TargetClosestPoint(target_surface, ufull);

    size_t num_time_steps = 1500;
    for(size_t t = 0; t < num_time_steps; ++t)
    {
        for(d = 0; d < 3; ++d)
        {
            ufull[d] = solver.TimeStep(ufull[d], ibc_location);
        }

        TargetClosestPoint(target_surface, ufull);
        cout << "Completed time step = " << (t + 1) << endl;
    }

    // visualize solution
    vector<vector<Scalar>> harmonic_map(solver.nNodes(), vector<Scalar>(3));
    for(size_t i = 0; i < solver.nNodes(); ++i)
    {
        for(size_t d = 0; d < 3; ++d)
        {
            harmonic_map[i][d] = ufull[d][i];
        }
    }
#ifdef POLYSCOPE
    polyscope::registerPointCloud("Harmonic Map", harmonic_map)->addColorQuantity("Colors", map_colors)->setEnabled(true);
    polyscope::getPointCloud("Harmonic Map")->setPointRadius(0.00164);
#endif

    return solver.DOFs();
}


// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
// https://ceng2.ktu.edu.tr/~cakir/files/grafikler/Texture_Mapping.pdf
void BarycentricFromPoint(const vector<Scalar> &p, const vector<Scalar> &a, const vector<Scalar> &b, const vector<Scalar> &c, Scalar &u, Scalar &v, Scalar &w)
{
    vector<Scalar> v0 = b - a, v1 = c - a, v2 = p - a;
    Scalar d00 = DotProduct(v0, v0);
    Scalar d01 = DotProduct(v0, v1);
    Scalar d11 = DotProduct(v1, v1);
    Scalar d20 = DotProduct(v2, v0);
    Scalar d21 = DotProduct(v2, v1);
    Scalar denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0f - v - w;
}


void PointFromBarycentric(vector<Scalar> &p, const vector<Scalar> &a, const vector<Scalar> &b, const vector<Scalar> &c, const Scalar &u, const Scalar &v, const Scalar &w)
{
    p = u * a + v * b + w * c;
}


void ValueFromBarycentric(Scalar &value, const Scalar &a, const Scalar &b, const Scalar &c, const Scalar &u, const Scalar &v, const Scalar &w)
{
    value = u * a + v * b + w * c;
}


int main(int argc, char** argv)
{
#ifdef POLYSCOPE
    polyscope::init();
    polyscope::view::lookAt(glm::vec3{-1.5, 1., -1.4}, glm::vec3{-0.5, 0.15, -0.36});
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::options::shadowDarkness = 0.75;
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    Scalar dx = 0.00625;
    dx = 0.025;
    
    ////////////////////////////////////////////////////////////
    // load the source and target meshes as well as ibcs
    ////////////////////////////////////////////////////////////
    SurfaceSpecifier source_specs("../assets/bob_tri_bff.obj");
    SurfaceSpecifier target_specs("../assets/bob_tri_simple_deformed.obj");

    vector<SurfaceSpecifier> source_ibc_specs(1);
    vector<SurfaceSpecifier> target_ibc_specs(1);
    source_ibc_specs[0].SetSurface("../assets/bob_source_loop.obj");
    target_ibc_specs[0].SetSurface("../assets/bob_target_loop.obj");

    IBCMetaData meta;
    meta.on_ibc_tolerance = 0.1 * dx * dx;
    source_ibc_specs[0].SetIBCMetaData(meta);

#ifdef POLYSCOPE
        polyscope::registerCurveNetwork("Source IBC", source_ibc_specs[0].Mesh().vertices, source_ibc_specs[0].Mesh().faces)->setColor(glm::vec3{1.0,1.0,1.0})->setRadius(0.01)->setEnabled(false);
        polyscope::registerCurveNetwork("Target IBC", target_ibc_specs[0].Mesh().vertices, target_ibc_specs[0].Mesh().faces)->setColor(glm::vec3{1.0,1.0,1.0})->setRadius(0.01);
#endif

    ///////////////////////////////////////////////////////////
    // load texture image
    ///////////////////////////////////////////////////////////
    Helpers h;
    h.ReadCSV("../assets/BobTextureRed.csv", red);
    h.ReadCSV("../assets/BobTextureGreen.csv", green);
    h.ReadCSV("../assets/BobTextureBlue.csv", blue);

    ///////////////////////////////////////////////////////////
    // Construct initial Map
    ///////////////////////////////////////////////////////////
    // For each cpx, find its barycentric coordinates in the triangle, then use those coordinates with the vertices on the target mesh. 
    // Since the meshes are 1-1 coorespondence, this will work

    // NOTE: this is inefficient because I am constructing the tube again in the HeatWithIBC solver in the HarmonicMaps function
    Tube source_tube(source_specs, interp_deg, dx);

    // get the face index that each source cpx belongs to
    vector<size_t> cpx_faces(source_tube.nNodes());
    Surface source_surface(source_specs);
    for(size_t i = 0; i < source_tube.nNodes(); ++i)
    {
        cpx_faces[i] = source_surface.cpTri(source_tube.x()[i]);
    }

    // get the barycentric coordinates of each source cpx and the corresponding point on the target face using those same barycentric coordinates
    vector<vector<Scalar>> initial_map(source_tube.nNodes(), vector<Scalar>(3));
    vector<vector<Scalar>> bary_coords(source_tube.nNodes(), vector<Scalar>(3));
    for(size_t i = 0; i < source_tube.nNodes(); ++i)
    {
        vector<size_t> face = source_specs.Mesh().faces[cpx_faces[i]];
        vector<Scalar> a = source_specs.Mesh().vertices[face[0]];
        vector<Scalar> b = source_specs.Mesh().vertices[face[1]];
        vector<Scalar> c = source_specs.Mesh().vertices[face[2]];

        BarycentricFromPoint(source_tube.cpx()[i], a, b, c, bary_coords[i][0], bary_coords[i][1], bary_coords[i][2]);

        face = target_specs.Mesh().faces[cpx_faces[i]];
        a = target_specs.Mesh().vertices[face[0]];
        b = target_specs.Mesh().vertices[face[1]];
        c = target_specs.Mesh().vertices[face[2]];
        PointFromBarycentric(initial_map[i], a, b, c, bary_coords[i][0], bary_coords[i][1], bary_coords[i][2]);
    }

    /////////////////////////////////////////////////////////////////////
    // Use parametrization of source to texture both source and target
    /////////////////////////////////////////////////////////////////////

    // get uv coordinates at source cpx
    vector<vector<Scalar>> cpx_uv(source_tube.nNodes(), vector<Scalar>(2));
    for(size_t i = 0; i < source_tube.nNodes(); ++i)
    {
        for(size_t uv_dim = 0; uv_dim < 2; ++uv_dim)
        {
            vector<vector<Scalar>> uv = source_specs.Mesh().uv[cpx_faces[i]];
            ValueFromBarycentric(cpx_uv[i][uv_dim], uv[0][uv_dim], uv[1][uv_dim], uv[2][uv_dim], bary_coords[i][0], bary_coords[i][1], bary_coords[i][2]);
        }
    }

    // use colors from the texture image to color source cpx and target points the same color
    vector<vector<Scalar>> colors(source_tube.nNodes());
    for(size_t i = 0; i < source_tube.nNodes(); ++i)
    {
        colors[i] = GetColor(cpx_uv[i]);
    }
#ifdef POLYSCOPE
    polyscope::registerPointCloud("Source cpx", source_tube.cpx())->setEnabled(false);
    polyscope::getPointCloud("Source cpx")->addColorQuantity("Colors", colors)->setEnabled(true);
    polyscope::getPointCloud("Source cpx")->setPointRadius(0.00164);
#endif

    // add some noise to the initial map
    Surface target_surface(target_specs);
    vector<vector<Scalar>> rand_noise(source_tube.nNodes(), vector<Scalar>(3));
    for(size_t i = 0; i < source_tube.nNodes(); ++i)
    {
        for(size_t d = 0; d < 3; ++d)
        {
            rand_noise[i][d] = (Scalar)(rand() % 1000 + 1) / 20000.0;
        }
        vector<Scalar> noisy = initial_map[i] + rand_noise[i];
        target_surface.ClosestPoint(noisy, initial_map[i]);
    }
#ifdef POLYSCOPE
    polyscope::registerPointCloud("Noisy Target Initial Map", initial_map)->setEnabled(false);
    polyscope::getPointCloud("Noisy Target Initial Map")->addColorQuantity("Colors", colors)->setEnabled(true);
    polyscope::getPointCloud("Noisy Target Initial Map")->setPointRadius(0.00164);
#endif

    ///////////////////////////////////////////
    // compute harmonic map
    ///////////////////////////////////////////
    HarmonicMap(initial_map, colors, source_specs, target_surface, source_ibc_specs, target_ibc_specs, dx);

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
