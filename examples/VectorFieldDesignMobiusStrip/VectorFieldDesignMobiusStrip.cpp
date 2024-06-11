#include <chrono>

#ifdef OPENMP
#include <omp.h>
#endif

#include "HeatWithIBC.h"

using namespace std;
using namespace cpm;

vector<size_t> g_subsamples;
vector<size_t> c_indices{110949, 291777, 142056, 191816};
vector<size_t> c_neighbour_indices{18166, 151466, 301228, 201958};

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


void VisualizeSolution(HeatWithIBC &solver, vector<cpm::VectorX> &ufull, vector<vector<Scalar>> &dir)
{
    vector<cpm::VectorX> uplot(3);    
    SpMat plotEfull = solver.GetPlotInterpMatrix();
    for(size_t d = 0; d < 3; ++d)
    {
        uplot[d] = plotEfull * ufull[d];
    }

    vector<vector<Scalar>> surface_normals_at_xp = solver.InterpolateUnorientedVectors(solver.surface_normals(), solver.xp(), true);
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
        if(c < 2)
        {
            polyscope::registerCurveNetwork("IBC " + to_string(c), solver.ibc_xp(c), solver.ibc_faces(c))->setColor(glm::vec3{1.0, 1.0, 1.0})->setRadius(0.004);
        }
        else
        {
            vector<vector<Scalar>> direct(1, vector<Scalar>(3));
            direct[0] = dir[c-2];
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


void VectorFieldDesign(SurfaceSpecifier &surface_specs, vector<SurfaceSpecifier> &ibc_surface_specs, Scalar dx)
{
    bool time_stepping_method = 1; // implicit Euler
    Scalar dt = 0.1 * dx;
    HeatWithIBC solver(surface_specs, ibc_surface_specs, dx, dt, time_stepping_method, 2 /* interp deg */);
    solver.ComputeSurfaceNormals();

    vector<vector<Scalar>> dir(4, vector<Scalar>(3));
    for(size_t i = 0; i < 4; ++i)
    {
        dir[i] = surface_specs.Mesh().vertices[c_neighbour_indices[i]] - surface_specs.Mesh().vertices[c_indices[i]];
        Normalize(dir[i], 0.0);
    }

    size_t d;
    auto ibc_tangent = [&solver, &d, dir](size_t &ibc_index, const size_t &which_side, const vector<Scalar>& x)
    { 
        if(ibc_index < 2)
        {
            return solver.IBCTangent(ibc_index, x)[d]; 
        }
        else
        {
            return dir[ibc_index-2][d];
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

    VisualizeSolution(solver, ufull, dir);
}


int main(int argc, char** argv)
{
#ifdef OPENMP
    omp_set_dynamic(0);
	omp_set_num_threads(8);
#endif

#ifdef POLYSCOPE
    polyscope::init(); 
    polyscope::view::lookAt(glm::vec3{-1.5, 1., -1.4}, glm::vec3{-0.5, 0.15, -0.36});
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    Scalar dx = 0.00625;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Read surface and IBCs
    //////////////////////////////////////////////////////////////////////////////////////////////
    SurfaceSpecifier surface_specs("../assets/Mobius_high_res.obj");

    vector<string> ibc_files{"../assets/M_O.obj", "../assets/M_block.obj"};
    vector<vector<size_t>> ibc_type{{0,0}, {1,1}, {0,0}, {0,0}, {0,0}, {0,0}};        

    size_t num_ibcs = 6;
    vector<SurfaceSpecifier> ibc_surface_specs(num_ibcs);
    for(size_t c = 0; c < num_ibcs; ++c)
    {
        if(c < 2)
        {
            ibc_surface_specs[c].SetSurface(ibc_files[c]);
        }
        else
        {
            ibc_surface_specs[c].SetSurface("Point", surface_specs.boundingBox(), surface_specs.Mesh().vertices[c_indices[c-2]], surface_specs.Dim());
        }

        IBCMetaData meta;
        meta.is_oriented = true;
        meta.on_ibc_tolerance = 0.1 * dx * dx;
        meta.boundary_type = ibc_type[c];
        ibc_surface_specs[c].SetIBCMetaData(meta);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////

    Helpers h;
    h.ReadCSV("../assets/MobiusStripSubsamples.csv", g_subsamples);

    VectorFieldDesign(surface_specs, ibc_surface_specs, dx);

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
