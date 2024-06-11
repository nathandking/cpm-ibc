#include <chrono>

#ifdef OPENMP
#include <omp.h>
#endif

#include "HeatWithIBC.h"

using namespace std;
using namespace cpm;

vector<size_t> c_indices{981232, 960537, 956686};
vector<size_t> c_neighbour_indices{728398, 908484, 805184};

SurfaceSpecifier g_surface_specs;
vector<SurfaceSpecifier> g_ibc_surface_specs;

size_t g_interp_deg = 2;

#ifdef POLYSCOPE
glm::vec3 surface_blue{0.1098, 0.3882, 0.8902};
#endif

#define FLOW_LINES

vector<vector<Scalar>> g_current_x;


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
        polyscope::registerPointCloud("Point Cloud with Directions", g_surface_specs.Mesh().vertices)->addVectorQuantity("direction", uplot_rearrange);
        polyscope::getPointCloud("Point Cloud with Directions")->setPointColor(surface_blue)->setPointRadius(0.001)->setEnabled(false);

        for(size_t c = 0; c < solver.NumIBCs(); ++c)
        {
            if(c < 2)
            {
                polyscope::registerCurveNetwork("IBC " + to_string(c), solver.ibc_xp(c), solver.ibc_faces(c))->setColor(glm::vec3{1.0, 1.0, 1.0})->setRadius(0.002);
            }
            else
            {
                vector<vector<Scalar>> direct(1, vector<Scalar>(3));
                direct[0] = dir[c-2];
                polyscope::registerPointCloud("IBC " + to_string(c), solver.ibc_xp(c))->setPointColor(glm::vec3{1.0, 1.0, 1.0})->setPointRadius(0.0045)->addVectorQuantity("Direction", direct)->setVectorLengthScale(0.025)->setVectorRadius(0.005)->setVectorColor(glm::vec3{1.0, 1.0, 1.0})->setEnabled(true);
            }
        }

    #ifdef FLOW_LINES
        // TO DO: make it so I don't have to load mesh again; need a surface specifier that just takes in a mesh...
        SurfaceSpecifier tri_surface_specs(g_surface_specs.Mesh());
        Surface tri_surface(tri_surface_specs);

        size_t num_steps = 200;
        size_t num_seeds = g_current_x.size();
        vector<vector<vector<Scalar>>> polylines(num_seeds, vector<vector<Scalar>>(num_steps, vector<Scalar>(3)));
        for(size_t i = 0; i < num_seeds; ++i)
        {
            polylines[i][0] = g_current_x[i];
        }

        Scalar step_size = 0.001;
        for(size_t t = 1; t < num_steps; ++t)
        {
            cout << "Visualization time step = " << t << endl;
            SpMat interp = solver.GetInterpMatrix(g_current_x);

            vector<cpm::VectorX> directions(3);
            for(size_t d = 0; d < 3; ++d)
            {
                directions[d] = interp * ufull[d];
            }

            vector<vector<Scalar>> directions_rearrange(num_seeds, vector<Scalar>(3));
            for(size_t i = 0; i < num_seeds; ++i)
            {
                for(size_t d = 0; d < 3; ++d)
                {
                    directions_rearrange[i][d] = directions[d][i];
                }
                Normalize(directions_rearrange[i], 0.0);
            }

            for(size_t i = 0; i < num_seeds; ++i)
            {
                vector<Scalar> point(3);
                point = g_current_x[i] + step_size * directions_rearrange[i];

                tri_surface.ClosestPoint(point, polylines[i][t]);

                g_current_x[i] = polylines[i][t];
            }
        }

        vector<vector<Scalar>> all_polylines(num_seeds * num_steps, vector<Scalar>(polylines[0][0].size()));
        for(size_t i = 0; i < num_seeds; ++i)
        {
            for(size_t t = 0; t < num_steps; ++t)
            {
                all_polylines[num_steps * i + t] = polylines[i][t];
            }
        }

        vector<vector<size_t>> edges;
        for(size_t i = 0; i < num_seeds; ++i)
        {
            for(size_t t = 0; t < num_steps - 1; ++t)
            {
                vector<size_t> edge(2);
                edge[0] = num_steps * i + t;
                edge[1] = num_steps * i + t + 1;
                edges.push_back(edge);
            }
        }

        polyscope::registerCurveNetwork("Flow ", all_polylines, edges)->setRadius(0.001)->setColor(glm::vec3(0.0, 0.0, 0.0));
    #endif
#endif
}


void VectorFieldDesign(Scalar dx)
{
    bool time_stepping_method = 1; // implicit Euler
    Scalar dt = 0.1 * dx;
    HeatWithIBC solver(g_surface_specs, g_ibc_surface_specs, dx, dt, time_stepping_method, g_interp_deg);
    solver.ComputeSurfaceNormals();

    vector<vector<Scalar>> dir(3, vector<Scalar>(3));
    for(size_t i = 0; i < 3; ++i)
    {
        dir[i] = g_surface_specs.Mesh().vertices[c_neighbour_indices[i]] - g_surface_specs.Mesh().vertices[c_indices[i]];
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
    polyscope::view::lookAt(glm::vec3{0.5, 1., 2.5}, glm::vec3{0.0, 0.0, 0.0});
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    Scalar dx = 0.0125;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Read surface and IBCs
    //////////////////////////////////////////////////////////////////////////////////////////////
    g_surface_specs.SetIsPointCloud(true);
    g_surface_specs.SetSurface("../assets/lucy_20L_scaled_centered.obj");

#ifdef POLYSCOPE
    polyscope::registerSurfaceMesh("Mesh", g_surface_specs.Mesh().vertices, g_surface_specs.Mesh().faces)->setSurfaceColor(surface_blue);
#endif

    IBCMetaData meta;
    meta.is_oriented = true;
    meta.on_ibc_tolerance = 0.1 * dx * dx;

    vector<string> ibc_files{"../assets/lucy_20L_belt.obj", "../assets/lucy_20L_skirt.obj"};      

    size_t num_ibcs = 5;
    g_ibc_surface_specs.resize(num_ibcs);
    for(size_t c = 0; c < num_ibcs; ++c)
    {
        if(c < 2)
        {
            g_ibc_surface_specs[c].SetSurface(ibc_files[c]);
        }
        else
        {
            g_ibc_surface_specs[c].SetSurface("Point", g_surface_specs.boundingBox(), g_surface_specs.Mesh().vertices[c_indices[c-2]], g_surface_specs.Dim());
        }

        g_ibc_surface_specs[c].SetIBCMetaData(meta);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////

    Helpers h;
    h.ReadCSV("../assets/LucySubsamplesX.csv", g_current_x);

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
