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
            if(c < 8 && c != 3)
            {
                curve_color[0] = blue[0]; curve_color[1] = blue[1]; curve_color[2] = blue[2];
            }
            else if(c == 3)
            {
                curve_color[0] = skin_dark[0]; curve_color[1] = skin_dark[1]; curve_color[2] = skin_dark[2];
            }
            else if(c == solver.NumIBCs() - 1)
            {
                curve_color[0] = white[0]; curve_color[1] = white[1]; curve_color[2] = white[2];
            }
            else
            {
                if(c % 2 == 0)
                {
                    curve_color[0] = dark_green[0]; curve_color[1] = dark_green[1]; curve_color[2] = dark_green[2];
                }
                else
                {
                    curve_color[0] = dark_red[0]; curve_color[1] = dark_red[1]; curve_color[2] = dark_red[2];
                }
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
    polyscope::view::lookAt(glm::vec3{0.2, 0.4, 3.0}, glm::vec3{0.01, 0.0, 0.66});
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    g_surface_specs.SetSurface("../assets/Nefertiti.obj", "../assets/Nefertiti_coarse.obj");
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Diffusion curve specifications
    //////////////////////////////////////////////////////////////////////////////////////////////
    vector<string> ibc_files{"../assets/N_hat_rim.obj", "../assets/N_hat_top_band.obj", "../assets/N_hat_bottom_band.obj", "../assets/N_neck.obj", "../assets/N_shirt_top_loop.obj", "../assets/N_shirt_bottom_loop.obj", "../assets/N_hat_emblem_bottom.obj", "../assets/N_hat_emblem_top.obj", "../assets/N_band_top_left1.obj", "../assets/N_band_top_left2.obj", "../assets/N_band_top_left3.obj", "../assets/N_band_top_left4.obj", "../assets/N_band_top_left5.obj", "../assets/N_band_top_left6.obj", "../assets/N_band_bottom_left1.obj", "../assets/N_band_bottom_left2.obj", "../assets/N_band_bottom_left3.obj", "../assets/N_band_bottom_left4.obj", "../assets/N_band_top_right1.obj", "../assets/N_band_top_right2.obj", "../assets/N_band_top_right3.obj", "../assets/N_band_top_right4.obj", "../assets/N_band_top_right5.obj", "../assets/N_band_top_right6.obj", "../assets/N_band_bottom_right1.obj", "../assets/N_band_bottom_right2.obj", "../assets/N_band_bottom_right3.obj", "../assets/N_band_bottom_right4.obj", "../assets/N_bottom.obj"};
    
    vector<vector<Scalar>> hat_rim{skin_light, blue};
    vector<vector<Scalar>> band_top{neumann, blue};
    vector<vector<Scalar>> band_bottom{neumann, blue};
    vector<vector<Scalar>> neck{skin_dark, white};
    vector<vector<Scalar>> shirt_top_loop{blue, shirt_red};
    vector<vector<Scalar>> shirt_bottom_loop{shirt_green, blue};
    vector<vector<Scalar>> emblem_bottom{blue, gold};
    vector<vector<Scalar>> emblem_top{blue, gold};
    vector<vector<Scalar>> band_div{red, dark_green};
    vector<vector<Scalar>> band_div2{dark_red, green};
    vector<vector<Scalar>> bottom{white, white};
    vector<vector<Scalar>> shirt{gold, gold};
    vector<vector<vector<Scalar>>> ibc_colors{hat_rim, band_top, band_bottom, neck, shirt_top_loop, shirt_bottom_loop, emblem_bottom, emblem_top, band_div, band_div2, band_div, band_div2, band_div, band_div2, band_div, band_div2, band_div, band_div2, band_div, band_div2, band_div, band_div2, band_div, band_div2, band_div, band_div2, band_div, band_div2, bottom/*, shirt*/};
    
    vector<vector<size_t>> ibc_type{{0,0}, {1,0}, {1,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}/*, {0,0}*/};
    
    vector<Scalar> p_hat_rim{-0.114368, 0.287988, 0.558695};
    vector<Scalar> p_band_top{-0.339168, 0.600229, 0.224288};
    vector<Scalar> p_band_bottom{-0.330943, 0.522861, 0.250379};
    vector<Scalar> p_neck{-0.016643, -0.743718, 0.259507};
    vector<Scalar> p_shirt_top_loop{-0.039600, -0.861542, 0.297695};
    vector<Scalar> p_shirt_bottom_loop{-0.030970, -0.916690, 0.310422};
    vector<Scalar> p_emblem_bottom{-0.038028, 0.497075, 0.493837};
    vector<Scalar> p_emblem_top{-0.046787, 0.804793, 0.350715};
    vector<Scalar> p_band_top_left1{0.164642, 0.638886, 0.389849};
    vector<Scalar> p_band_top_left2{0.274268, 0.605838, 0.314053};
    vector<Scalar> p_band_top_left3{0.350645, 0.563024, 0.190168};
    vector<Scalar> p_band_top_left4{0.387360, 0.486210, 0.019568};
    vector<Scalar> p_band_top_left5{0.370991, 0.411041, -0.137812};
    vector<Scalar> p_band_top_left6{0.308649, 0.371414, -0.274834};
    vector<Scalar> p_band_bottom_left1{0.238515, 0.234822, -0.223521};
    vector<Scalar> p_band_bottom_left2{0.285281, 0.223652, -0.096570};
    vector<Scalar> p_band_bottom_left3{0.309291, 0.215114, -0.004209};
    vector<Scalar> p_band_bottom_left4{0.321091, 0.223768, 0.129298};
    vector<Scalar> p_band_top_right1{-0.173687, 0.638149, 0.388254};
    vector<Scalar> p_band_top_right2{-0.259027, 0.610701, 0.334402};
    vector<Scalar> p_band_top_right3{-0.353421, 0.561688, 0.196354};
    vector<Scalar> p_band_top_right4{-0.387048, 0.503286, 0.050060};
    vector<Scalar> p_band_top_right5{-0.379138, 0.435644, -0.091840};
    vector<Scalar> p_band_top_right6{-0.347024, 0.395856, -0.203488};
    vector<Scalar> p_band_bottom_right1{-0.250955, 0.246216, -0.226057};
    vector<Scalar> p_band_bottom_right2{-0.297929, 0.235448, -0.112715};
    vector<Scalar> p_band_bottom_right3{-0.316184, 0.225537, -0.006552};
    vector<Scalar> p_band_bottom_right4{-0.316519, 0.201636, 0.147822};
    vector<Scalar> p_bottom{-0.197666, -0.928798, 0.310350};
    vector<vector<Scalar>> point_on_side1{p_hat_rim, p_band_top, p_band_bottom, p_neck, p_shirt_top_loop, p_shirt_bottom_loop, p_emblem_bottom, p_emblem_top, p_band_top_left1, p_band_top_left2, p_band_top_left3, p_band_top_left4, p_band_top_left5, p_band_top_left6, p_band_bottom_left1, p_band_bottom_left2, p_band_bottom_left3, p_band_bottom_left4, p_band_top_right1, p_band_top_right2, p_band_top_right3, p_band_top_right4, p_band_top_right5, p_band_top_right6, p_band_bottom_right1, p_band_bottom_right2, p_band_bottom_right3, p_band_bottom_right4, p_bottom/*, p_shirt0*/};

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Scalar dx = 0.003125;

    g_ibc_surface_specs.resize(ibc_files.size());
    for(size_t c = 0; c < ibc_files.size(); ++c)
    {
        g_ibc_surface_specs[c].SetSurface(ibc_files[c]);

        IBCMetaData meta;
        meta.is_oriented = true;
        meta.on_ibc_tolerance = 0.1 * dx * dx;
        meta.point_on_side1 = point_on_side1[c];
        meta.boundary_type = ibc_type[c];

        g_ibc_surface_specs[c].SetIBCMetaData(meta);
    }
    
    DiffusionCurves(ibc_colors, dx);

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
