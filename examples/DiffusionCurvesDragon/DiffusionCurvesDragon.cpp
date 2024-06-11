#include <chrono>

#ifdef OPENMP
#include <omp.h>
#endif

#include "PoissonWithIBC.h"
#include "PoissonWithIBCNearest.h"
#include "Solver.h"

using namespace std;
using namespace cpm;

SurfaceSpecifier g_surface_specs;
vector<SurfaceSpecifier> g_ibc_surface_specs;

//////////////////////////////////////////////////////////////////////////////////////////////
// Colors
//////////////////////////////////////////////////////////////////////////////////////////////
#ifdef POLYSCOPE
glm::vec3 surface_blue{0.1098, 0.3882, 0.8902};
#endif

vector<Scalar> surface_blue_vec{0.1098, 0.3882, 0.8902};
vector<Scalar> blue{0.00,0.05,0.16};  
vector<Scalar> shirt_red{0.5294, 0.1451, 0.0902};
vector<Scalar> red{1, 0.05, 0.05};
vector<Scalar> dark_red{0.4, 0.05, 0.05};
vector<Scalar> skin_light{0.6961,0.3412,0.1765};
vector<Scalar> skin_dark{0.5961,0.2412,0.0765};
vector<Scalar> shirt_green{0.176, 0.49, 0.039};
vector<Scalar> green{0.2, 1.0, 0.2};
vector<Scalar> dark_green{0., 0.4, 0.0};
vector<Scalar> white{0.9, 0.9, 0.9};
vector<Scalar> gold{0.5490, 0.5059, 0.2353};
vector<Scalar> yellow{1.0, 1.0, 0.0};
vector<Scalar> neumann{0.0, 0.0, 0.0}; // dummy color for Neumann BCs

#define SOLVE_COLORS

void DiffusionCurves(vector<vector<vector<Scalar>>> &ibc_colors, Scalar dx)
{
#ifdef TRACK_WHERE
    std::cout << "DiffusionCurves" << std::endl;
#endif

    PoissonWithIBC solver(g_surface_specs, g_ibc_surface_specs, dx, 2 /* interp deg */);

    auto f = [](const vector<Scalar>& x){ return 0; };

#ifdef SOLVE_COLORS
    vector<cpm::VectorX> color(3);
    vector<vector<vector<Scalar>>> ibc_color(solver.NumIBCs(), vector<vector<Scalar>>(3));
    for(size_t color_channel = 0; color_channel < 3; ++color_channel)
    {
        std::cout << "color: " << color_channel << std::endl;
        auto color_ibc = [&solver, &ibc_colors, color_channel](size_t &ibc_index, const size_t &which_side, const vector<Scalar>& x){ return /*sin(M_PI * solver.IBCParameter(ibc_index, x)) * */ ibc_colors[ibc_index][which_side - 1][color_channel]; };

        Solver sv;
        sv.Init((Scalar)0., (Scalar)-1., solver.GetLfull(), solver.GetEfull(), Preconditioner::diagonal, solver.GetIdentityRows());
        cpm::VectorX rhs = solver.GetRHS(f, color_ibc);
        cpm::VectorX ufull = sv.Solve(rhs);

        SpMat plotEfull = solver.GetPlotInterpMatrix();

        color[color_channel] = plotEfull * ufull;

        for(size_t c = 0; c < solver.NumIBCs(); ++c)
        {
            ibc_color[c][color_channel].resize(solver.ibc_xp(c).size());
            for(size_t i = 0; i < solver.ibc_xp(c).size(); ++i)
            {
                ibc_color[c][color_channel][i] = color_ibc(c, 1, solver.ibc_xp(c)[i]);
            }
        }
    }

    vector<vector<Scalar>> colors(color[0].size(), vector<Scalar>(3));
    for(size_t i = 0; i < color[0].size(); ++i)
    {
        for(size_t color_channel = 0; color_channel < 3; ++color_channel)
        {
            colors[i][color_channel] = color[color_channel][i];
        }
    }

    vector<vector<vector<Scalar>>> c_colors(solver.NumIBCs());
    for(size_t c = 0; c < solver.NumIBCs(); ++c)
    {
        c_colors[c].resize(ibc_color[c][0].size(), vector<Scalar>(3));
        for(size_t i = 0; i < ibc_color[c][0].size(); ++i)
        {
            for(size_t color_channel = 0; color_channel < 3; ++color_channel)
            {            
                c_colors[c][i][color_channel] = ibc_color[c][color_channel][i];
            }
        }
    }

#endif

#ifdef POLYSCOPE
    polyscope::registerSurfaceMesh("Surface", solver.xp(), solver.faces())->setSurfaceColor(surface_blue);
    #ifdef SOLVE_COLORS
        polyscope::getSurfaceMesh("Surface")->addVertexColorQuantity("Colors", colors)->setEnabled(true);
    #endif

        glm::vec3 curve_color;
        for(size_t c = 0; c < solver.NumIBCs(); ++c)
        {
    #ifdef SOLVE_COLORS
            if(c == 0)
            {
                curve_color[0] = surface_blue_vec[0]; curve_color[1] = surface_blue_vec[1]; curve_color[2] = surface_blue_vec[2];
            }
            else if(c == 1 || c == 3)
            {
                curve_color[0] = red[0]; curve_color[1] = red[1]; curve_color[2] = red[2];
            }
            else if(c == 2)
            {
                curve_color[0] = shirt_green[0]; curve_color[1] = shirt_green[1]; curve_color[2] = shirt_green[2];
            }

            polyscope::registerCurveNetwork("IBC " + to_string(c), solver.ibc_xp(c), solver.ibc_faces(c))->setRadius(0.001)->setColor(curve_color)->setEnabled(true);
            polyscope::getCurveNetwork("IBC " + to_string(c))->addNodeColorQuantity("uexact_ibc", c_colors[c]);
    #else
            polyscope::registerCurveNetwork("IBC " + to_string(c), solver.ibc_xp(c), solver.ibc_faces(c))->setRadius(0.003)->setEnabled(true);
            if(ibc_meta[c].boundary_type[0] == 1 || ibc_meta[c].boundary_type[1] == 1)
            {
                polyscope::getCurveNetwork("IBC " + to_string(c))->setColor(glm::vec3(0,0,0));
            }
            else
            {
                polyscope::getCurveNetwork("IBC " + to_string(c))->setColor(glm::vec3(1,1,1));
            }
    #endif
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
    polyscope::view::lookAt(glm::vec3{0.2, 0.3, -1.6}, glm::vec3{0.0, -0.18, 0.0});
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    g_surface_specs.SetSurface("../assets/dragon_scaled_flipped_centered.obj");
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Diffusion curve specifications
    //////////////////////////////////////////////////////////////////////////////////////////////
    vector<string> ibc_files{"../assets/dragon_neck.obj", "../assets/dragon_body.obj", "../assets/dragon_tail.obj", "../assets/dragon_scale.obj"};

    vector<vector<Scalar>> neck{surface_blue_vec, red};
    vector<vector<Scalar>> body{red, shirt_green};
    vector<vector<Scalar>> tail{shirt_green, yellow};
    vector<vector<Scalar>> scale{red, surface_blue_vec};
    vector<vector<vector<Scalar>>> ibc_colors{neck, body, tail, scale};
    
    vector<Scalar> p_neck{0.369766, -0.062277, 0.014395};
    vector<Scalar> p_body{0.017274, 0.019249, -0.078513};
    vector<Scalar> p_tail{-0.579615, -0.144089, -0.326012};
    vector<Scalar> p_scale{-0.204649, -0.167266, -0.372070};
    vector<vector<Scalar>> point_on_side1{p_neck, p_body, p_tail, p_scale};

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Scalar dx = 0.003125; // takes around 17 hours
    // dx = 0.00625;
    dx = 0.0125;

    g_ibc_surface_specs.resize(ibc_files.size());
    for(size_t c = 0; c < ibc_files.size(); ++c)
    {
        g_ibc_surface_specs[c].SetSurface(ibc_files[c]);

        IBCMetaData meta;
        meta.is_oriented = true;
        meta.on_ibc_tolerance = 0.1 * dx * dx;
        meta.point_on_side1 = point_on_side1[c];

        g_ibc_surface_specs[c].SetIBCMetaData(meta);
    }

    DiffusionCurves(ibc_colors, dx);

#ifdef POLYSCOPE
    // add a cube to see scale of grid resolution
    vector<vector<Scalar>> verts;
    vector<Scalar> lower_left{-0.05,-0.4083,-0.6200}; // near back foot
    size_t num_cubes_1d = 4;
    for(size_t x = 0; x < num_cubes_1d; ++x)
    {
        for(size_t y = 0; y < num_cubes_1d; ++y)
        {
            for(size_t z = 0; z < num_cubes_1d; ++z)
            {
                vector<Scalar> v0 = lower_left;
                v0[0] += x * dx;
                v0[1] += y * dx;
                v0[2] += z * dx;

                vector<Scalar> v1 = v0;
                v1[0] += dx;
                vector<Scalar> v2 = v1;
                v2[1] += dx;
                vector<Scalar> v3 = v2;
                v3[0] -= dx;
                vector<Scalar> v4 = v0;
                vector<Scalar> v5 = v1;
                vector<Scalar> v6 = v2;
                vector<Scalar> v7 = v3;
                v4[2] -= dx;
                v5[2] -= dx;
                v6[2] -= dx;
                v7[2] -= dx;

                verts.push_back(v0);
                verts.push_back(v1);
                verts.push_back(v2);
                verts.push_back(v3);
                verts.push_back(v4);
                verts.push_back(v5);
                verts.push_back(v6);
                verts.push_back(v7);
            }
        }
    }
    
    vector<vector<size_t>> hexes;
    for(size_t i = 0; i < pow(num_cubes_1d, 3); ++i)
    {
        vector<size_t> hex(8);
        for(size_t j = 0; j < 8; ++j)
        {
            hex[j] = i * 8 + j;
        }
        hexes.push_back(hex);
    }
    polyscope::registerHexMesh("Cube", verts, hexes)->setColor(glm::vec3{1,1,1})->setEdgeColor(glm::vec3{0,0,0})->setEnabled(true);
#endif
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
