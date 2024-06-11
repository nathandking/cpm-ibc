// Note: Eigen solver used for result in paper "A Closest Point Method for Surface PDEs with Interior Boundary Conditions for Geometry Processing"
#include <chrono>

#include "HeatWithIBC.h"

using namespace std;
using namespace cpm;

vector<size_t> g_subsamples;
vector<vector<Scalar>> point_IBCs{ {1.7076,0.472749,-1.11021}, {1.7076,0.0285716,-1.20633}, {1.6553,0.293987,-1.35848}, {1.53124,0.215062,-0.993774}, {1.61702,-0.342926,1.02236}, {1.61702,-0.643437,0.865341}, {1.71824,-0.583514,1.11673}, {1.40528,-0.445278,0.852175}, {0.487487,-0.369731,0.831606}, {0.546848,-0.66938,-0.962247} };
vector<vector<Scalar>> point_IBC_directions{ {0,-0.926129,-0.377207}, {0,0.999221,0.039457}, {0.762689,-0.136799,0.632132}, {0.850219,0.111347,-0.514519}, {0,0.952989,0.303005}, {0,-0.792953,-0.609283}, {-0.0298822,-0.462903,0.885905}, {-0.957571,0.133467,-0.25543}, {0,-0.900479,-0.434899}, {0,0.829821,-0.558029} };


void ProjectOutNormalComponent(vector<cpm::VectorX> &ufull, const vector<vector<Scalar>> &normals)
{
    for(size_t i = 0; i < normals.size(); ++i)
    {
        vector<Scalar> u{ufull[0][i], ufull[1][i], ufull[2][i]};
        
        u = u - DotProduct(u, normals[i]) * normals[i];
        
        ufull[0][i] = u[0];
        ufull[1][i] = u[1];
        ufull[2][i] = u[2];
    }
}


void NormalizeVector(vector<cpm::VectorX> &ufull)
{
    for(size_t i = 0; i < ufull[0].size(); ++i)
    {
        vector<Scalar> u{ufull[0][i], ufull[1][i], ufull[2][i]};
        
        Normalize(u, 0.0);
        
        ufull[0][i] = u[0];
        ufull[1][i] = u[1];
        ufull[2][i] = u[2];
    }
}


void VisualizeSolution(HeatWithIBC &solver, vector<cpm::VectorX> &ufull)
{
    vector<cpm::VectorX> uplot(3);
    SpMat plotEfull = solver.GetPlotInterpMatrix();
    for(size_t d = 0; d < 3; ++d)
    {
        uplot[d] = plotEfull * ufull[d];
    }

    vector<vector<Scalar>> surface_normals_at_xp = solver.InterpolateUnorientedVectors(solver.surface_normals(), solver.xp(), true); // I think using this is not accurate enough in high curvature regions, therefore the extra projection looks worse
    ProjectOutNormalComponent(uplot, surface_normals_at_xp);
    NormalizeVector(uplot);

    vector<vector<Scalar>> uplot_rearrange(uplot[0].size(), vector<Scalar>(3));
    for(size_t i = 0; i < uplot[0].size(); ++i)
    {
        for(size_t d = 0; d < 3; ++d)
        {
            uplot_rearrange[i][d] = uplot[d][i];
        }
        Normalize(uplot_rearrange[i], 0.0);
    }

#ifdef POLYSCOPE
    polyscope::registerSurfaceMesh("Surface", solver.xp(), solver.faces())->setSmoothShade(true);
    for(size_t c = 0; c < solver.NumIBCs(); ++c)
    {
        if(c < 1)
        {
            polyscope::registerCurveNetwork("IBC " + to_string(c), solver.ibc_xp(c), solver.ibc_faces(c))->setColor(glm::vec3{1.0, 1.0, 1.0})->setRadius(0.004);
        }
        else
        {
            vector<vector<Scalar>> direct(1, vector<Scalar>(3));
            direct[0] = point_IBC_directions[c-1];
            polyscope::registerPointCloud("IBC " + to_string(c), solver.ibc_xp(c))->setPointColor(glm::vec3{1.0, 1.0, 1.0})->setPointRadius(0.009)->addVectorQuantity("Direction", direct)->setVectorLengthScale(0.055)->setVectorRadius(0.01)->setVectorColor(glm::vec3{1.0, 1.0, 1.0})->setEnabled(true);
        }
    }

    vector<vector<Scalar>> xp_subset(g_subsamples.size(), vector<Scalar>(3));
    vector<vector<Scalar>> subset_directions(g_subsamples.size(), vector<Scalar>(3));
    for(size_t i = 0; i < g_subsamples.size(); ++i)
    {
        xp_subset[i] = solver.xp()[g_subsamples[i]];
        subset_directions[i] = uplot_rearrange[g_subsamples[i]];
    }
    polyscope::registerPointCloud("xp subset", xp_subset)->setPointRadius(0.0)->addVectorQuantity("Directions", subset_directions)->setVectorColor(glm::vec3{0.0, 0.0, 0.0})->setEnabled(true);
#endif
}


void VectorFieldDesign(Scalar dx)
{
    SurfaceSpecifier surface_specs("Surface Of Revolution");
    
    IBCMetaData meta;
    meta.is_oriented = true;
    meta.on_ibc_tolerance = 0.1 * dx * dx;

    size_t num_ibcs = 1 + point_IBCs.size();
    vector<SurfaceSpecifier> ibc_surface_specs(num_ibcs);
    ibc_surface_specs[0].SetSurface("Surface Of Revolution Curve", surface_specs.boundingBox(), {3.0, 4.0, 0.5}, surface_specs.Dim());
    ibc_surface_specs[0].SetIBCMetaData(meta);

    for(size_t c = 1; c < num_ibcs; ++c)
    {
        ibc_surface_specs[c].SetSurface("Point", surface_specs.boundingBox(), point_IBCs[c-1], surface_specs.Dim());
        ibc_surface_specs[c].SetIBCMetaData(meta);
    }

    bool time_stepping_method = 1;
    Scalar dt = 0.1 * dx;
    HeatWithIBC solver(surface_specs, ibc_surface_specs, dx, dt, time_stepping_method, 2 /* interp deg */);
    solver.ComputeSurfaceNormals();

    size_t d;
    vector<Scalar> flip_tangent{-1.0, 1.0};
    auto ibc_tangent = [&solver, &d, flip_tangent](size_t &ibc_index, const size_t &which_side, const vector<Scalar>& x)
    {
        if(ibc_index < 1)
        {
            return flip_tangent[which_side - 1] * solver.IBCTangent(ibc_index, x)[d]; 
        }
        else
        {
            return point_IBC_directions[ibc_index-1][d];
        }
    };

    vector<cpm::VectorX> ufull(3, cpm::VectorX::Zero(solver.nNodes() + solver.ibc().TotalIBCDOFs() + solver.ibc().IdentityRowsStart()));
    for(d = 0; d < 3; ++d)
    {
        for(size_t c = 0; c < solver.ibc().NumIBCs(); ++c)
        {
            for(size_t i = 0; i < solver.ibc().DOFSubset()[c].size(); ++i)
            {
                if(solver.ibc().meta(c).boundary_type[0] == 0)
                {
                    ufull[d][solver.ibc().IdentityRowsStart() + solver.ibc().IBCColumnStartIndex(c)] = ibc_tangent(c, solver.ibc().DOFSubset()[c].which_side()[i], solver.ibc().DOFSubset()[c].cpx()[i]);
                }
            }
        }
    }
    NormalizeVector(ufull);

    size_t num_time_steps = 10;
    for(size_t t = 0; t < num_time_steps; ++t)
    {
        cout << "time step = " << t << endl;
        for(d = 0; d < 3; ++d)
        {            
            ufull[d] = solver.TimeStep(ufull[d], ibc_tangent);
        }

        ProjectOutNormalComponent(ufull, solver.surface_normals());
        NormalizeVector(ufull);
    }
    
    VisualizeSolution(solver, ufull);
}


int main(int argc, char** argv)
{
#ifdef POLYSCOPE
    polyscope::init(); 
    polyscope::view::lookAt(glm::vec3{5.5, 1., 0}, glm::vec3{1, 0, 0});
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    Helpers h;
    h.ReadCSV("../assets/SurfaceOfRevolutionSubsamples.csv", g_subsamples);

    Scalar dx = 0.025;
    VectorFieldDesign(dx);

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
