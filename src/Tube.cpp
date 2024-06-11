#include "Tube.h"
#include "Defines.h"

#ifdef USE_SPARSE_GRID
#include "SPGrid_Initializer.h"
#include "CPM_Circle_Rasterizer.h"
#include "CPM_Shifted_Sphere_Rasterizer.h"
#include "CPM_Shifted_Hemisphere_Rasterizer.h"
#include "CPM_Mesh_Rasterizer.h"
#include "CPM_Rasterizer.h"

#include "Influence_Iterator.h"
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy_Iterator.h>
#include <nova/Tools/Log/Log.h>
#endif

namespace cpm {

Tube::Tube()
#ifdef USE_SPARSE_GRID
: hierarchy2(nullptr), hierarchy3(nullptr)
#endif
{
#ifdef TRACK_WHERE
    std::cout << "Tube::Tube(0)" << std::endl;
#endif
}

Tube::~Tube()
{
#ifdef TRACK_WHERE
    std::cout << "Tube::~Tube()" << std::endl;
#endif
#ifdef USE_SPARSE_GRID 
    if(hierarchy2 != nullptr) {
        delete hierarchy2;
        hierarchy2 = nullptr;
    }
    if(hierarchy3 != nullptr) {
        delete hierarchy3;
        hierarchy3 = nullptr;
    }
#endif
}


Tube::Tube(SurfaceSpecifier &surface_specs, size_t interp_deg, Scalar dx)
#ifdef USE_SPARSE_GRID
: hierarchy2(nullptr), hierarchy3(nullptr)
#endif
{
#ifdef TRACK_WHERE
    std::cout << "Tube::Tube(1)" << std::endl;
#endif
    ConstructTube(surface_specs, interp_deg, dx);
}


void Tube::Clear()
{
    m_surface.Clear();
    m_grid.Clear();

    m_base_points.clear();
    m_is_bpt_visited.clear();

    m_tube_coords.clear(); // coordinates of the tube nodes
    
    m_tube_nodes.clear(); // grid nodes in the computational tube

#ifdef USE_SPARSE_GRID
#ifdef COMPARE_GRIDS
    m_tube_indices_on_grid.clear(); // indicies with respect to tube nodes, but stored for each rectangular grid index
#endif
#endif
#ifndef USE_SPARSE_GRID
    m_tube_indices_on_grid.clear(); // indicies with respect to tube nodes, but stored for each rectangular grid index
#endif

    m_num_nodes = 0; // total number of grid points in computational tube

    m_cpx.clear(); // closest points to nodes in the computational tube
    m_dist.clear();
    m_bdy.clear(); // tag for which boundary the cp in tube belongs to, 0 if not a boundary cp
}


void Tube::ConstructTube(SurfaceSpecifier &surface_specs, size_t interp_deg, Scalar dx)
{
#ifdef TRACK_WHERE
    std::cout << "Tube::ConstructTube(1)" << std::endl;
#endif
    m_surface.SetSurfaceSpecs(surface_specs);
    ConstructTube(interp_deg, dx);
}


void Tube::replaceClosestPoint(size_t index, vector<Scalar> &cpx)
{ 
#ifdef TRACK_WHERE
    std::cout << "Tube::replaceClosestPoint(size_t index, vector<Scalar> &cpx)" << std::endl;
#endif
    for(size_t d = 0; d < cpx.size(); ++d)
    {
        m_cpx[index][d] = cpx[d];
    } 
}


////////////////////////////////////////////////////////////////////////////////////
// Private Functions
////////////////////////////////////////////////////////////////////////////////////


void Tube::ConstructTube(size_t interp_deg, Scalar dx)
{
#ifdef TRACK_WHERE
    std::cout << "Tube::ConstructTube(3)" << std::endl;
#endif
    m_interp_deg = interp_deg;
    m_grid.SetDim(m_surface.Dim());
    m_grid.SetDx(dx);
    m_order = 2;
    
    ComputeTubeRadius();

    ConstructTube();
}


void Tube::ComputeTubeRadius()
{
#ifdef TRACK_WHERE
    std::cout << "Tube::ComputeTubeRadius()" << std::endl;
#endif
    m_tube_radius_factor = 1.0001 * sqrt( (m_grid.dim()-1) * pow(0.5 * (m_interp_deg+1), 2) + pow(0.5 * (m_order + m_interp_deg + 1), 2) );
    m_tube_radius = m_tube_radius_factor * m_grid.dx();
}


// Set up rectangular embedding space grid within bounding box
void Tube::ConstructTube()
{
#ifdef TRACK_WHERE
    std::cout << "Tube::ConstructTube(6)" << std::endl;
#endif

    size_t pad = 5;
    size_t spacings_covering_bbox = ceil( (m_surface.boundingBox()[1] - m_surface.boundingBox()[0]) / m_grid.dx() ); // see how many dx spacings we need starting from bBox[0] to get past bBox[1]
    size_t n = ceil( log2(spacings_covering_bbox + 2 * pad) ); // see what power of 2^n is enough spacings to cover the bBox with padding on each side

    size_t spacing_diff = ceil( 0.5 * (pow(2, n) - spacings_covering_bbox - 2 * pad) ); // how many extra spacings you must add on one side of the bBox plus padding

    vector<Scalar> xstart(m_grid.dim(), m_surface.boundingBox()[0] - (pad + spacing_diff) * m_grid.dx()); 
    m_grid.SetXstart(xstart);

    m_grid.SetNumberGridPoints(n);
#ifdef USE_SPARSE_GRID
    if(m_grid.dim() == 2)
    {
        T_INDEX2 counts = T_INDEX2(int(pow(2, n)));
        TV2 min_corner, max_corner;
        min_corner(0) = xstart[0];
        min_corner(1) = xstart[1];
        max_corner = min_corner + float(m_grid.dx()) * TV2(pow(2, n));
        T_Range2 domain = T_Range2(min_corner, max_corner);
        T_Grid2 grid = T_Grid2(counts, domain);
        int base_levels = 1;
        
        if(hierarchy2 != nullptr) {
            delete hierarchy2;
            hierarchy2 = nullptr;
        }
        hierarchy2 = new Hierarchy2(grid, base_levels);

        Nova::Vector<double, 2> min_cornerd;
        min_cornerd(0) = (double)xstart[0];
        min_cornerd(1) = (double)xstart[1];

        using CPM_Rasterizer = CPM_Rasterizer<Struct_type, data_type, 2>;
        CPM_Rasterizer cpm_rasterizer(*hierarchy2, m_surface, m_tube_radius, m_grid.dx(), min_cornerd);
        for(Grid_Hierarchy_Iterator<2, CPM_Rasterizer> iterator(hierarchy2->Lattice(base_levels-1).Node_Indices(), base_levels - 1, cpm_rasterizer);iterator.Valid();iterator.Next());
        hierarchy2->Update_Block_Offsets();

        // if(m_surface.SurfaceType().compare("Circle") == 0)
        // {
        //     using CPM_Rasterizer = CPM_Circle_Rasterizer<Struct_type, data_type, 2>;
        //     CPM_Rasterizer cpm_rasterizer(*hierarchy2, m_surface.surfaceParams()[0], m_tube_radius, m_grid.dx(), min_cornerd);
        //     for(Grid_Hierarchy_Iterator<2, CPM_Rasterizer> iterator(hierarchy2->Lattice(base_levels - 1).Node_Indices(), base_levels - 1, cpm_rasterizer); iterator.Valid(); iterator.Next());
        //     hierarchy2->Update_Block_Offsets();
        // }
    }
    else if(m_grid.dim() == 3)
    {
        T_INDEX3 counts = T_INDEX3(int(pow(2, n)));
        TV3 min_corner, max_corner;
        min_corner(0) = xstart[0];
        min_corner(1) = xstart[1];
        min_corner(2) = xstart[2];
        max_corner = min_corner + float(m_grid.dx()) * TV3(pow(2, n));
        T_Range3 domain = T_Range3(min_corner, max_corner);
        T_Grid3 grid = T_Grid3(counts, domain);
        int base_levels = 1;

        if(hierarchy3 != nullptr) {
            delete hierarchy3;
            hierarchy3 = nullptr;
        }
        hierarchy3= new Hierarchy3(grid, base_levels);

        Nova::Vector<double, 3> min_cornerd;
        min_cornerd(0) = (double)xstart[0];
        min_cornerd(1) = (double)xstart[1];
        min_cornerd(2) = (double)xstart[2];

        if(m_surface.SurfaceType().compare("Sphere") == 0)
        {
            using CPM_Rasterizer = CPM_Shifted_Sphere_Rasterizer<Struct_type, data_type, 3>;
            CPM_Rasterizer cpm_rasterizer(*hierarchy3, m_surface.surfaceParams()[0], m_tube_radius, m_grid.dx(), min_cornerd);
            for(Grid_Hierarchy_Iterator<3, CPM_Rasterizer> iterator(hierarchy3->Lattice(base_levels-1).Node_Indices(), base_levels - 1, cpm_rasterizer);iterator.Valid();iterator.Next());
            hierarchy3->Update_Block_Offsets();
        }
        else if(m_surface.SurfaceType().compare("Hemisphere") == 0)
        {
            using CPM_Rasterizer = CPM_Shifted_Hemisphere_Rasterizer<Struct_type, data_type, 3>;
            CPM_Rasterizer cpm_rasterizer(*hierarchy3, m_surface.surfaceParams()[0], m_tube_radius, m_grid.dx(), min_cornerd);
            for(Grid_Hierarchy_Iterator<3, CPM_Rasterizer> iterator(hierarchy3->Lattice(base_levels-1).Node_Indices(), base_levels - 1, cpm_rasterizer);iterator.Valid();iterator.Next());
            hierarchy3->Update_Block_Offsets();
        }
        else if(m_surface.SurfaceType().compare("Triangulation") == 0)
        {
            if(m_surface.SurfaceSpecs().CoarseMesh().nVertices() == 0)
            {
                using CPM_Rasterizer = CPM_Mesh_Rasterizer<Struct_type, data_type, 3>;
                CPM_Rasterizer cpm_rasterizer(*hierarchy3, m_surface, m_tube_radius, m_grid.dx(), min_cornerd);
                for(Grid_Hierarchy_Iterator<3, CPM_Rasterizer> iterator(hierarchy3->Lattice(base_levels-1).Node_Indices(), base_levels - 1, cpm_rasterizer);iterator.Valid();iterator.Next());
                hierarchy3->Update_Block_Offsets();
            }
            else
            {
                std::cout << "# vertices: " << m_surface.SurfaceSpecs().CoarseMesh().nVertices() << std::endl;
                std::cout << "# faces: " << m_surface.SurfaceSpecs().CoarseMesh().nFaces() << std::endl;

                using T_Influence_Iterator = Influence_Iterator<T, 3, T_INDEX3>;
                T_INDEX3 bbox_min, bbox_max;
                int tr = int(m_tube_radius_factor) + 13 * n;
                // int tr = int(m_tube_radius_factor) + 30 * n;
                for(int f = 0; f < m_surface.SurfaceSpecs().CoarseMesh().nFaces(); ++f)
                {
                    bbox_min = T_INDEX3(1000000); 
                    bbox_max = T_INDEX3(-1000000);
                    for(int v = 0; v < m_surface.SurfaceSpecs().CoarseMesh().faces[f].size(); ++v)
                    {
                        TV3 vpos = TV3({data_type(m_surface.SurfaceSpecs().CoarseMesh().vertices[m_surface.SurfaceSpecs().CoarseMesh().faces[f][v]][0]), data_type(m_surface.SurfaceSpecs().CoarseMesh().vertices[m_surface.SurfaceSpecs().CoarseMesh().faces[f][v]][1]), data_type(m_surface.SurfaceSpecs().CoarseMesh().vertices[m_surface.SurfaceSpecs().CoarseMesh().faces[f][v]][2])});
                        T_INDEX3 closest_node = grid.Closest_Node(vpos);
                    #ifdef OUTPUT_INFO
                        Log::cout << "cn: " << closest_node << std::endl;
                    #endif
                        for(int axis = 0; axis < 3; ++axis)
                        {
                            bbox_min(axis) = std::min(bbox_min(axis), closest_node(axis));
                            bbox_max(axis) = std::max(bbox_max(axis), closest_node(axis));
                        }
                    }
                    bbox_min -= T_INDEX3(tr + 1);
                    bbox_max += T_INDEX3(tr + 1);
                    for(int axis = 0; axis < 3; ++ axis)
                    {
                        bbox_min(axis) = std::max(1, bbox_min(axis));
                        bbox_max(axis) = std::min(grid.counts(axis), bbox_max(axis));
                    }
                    for (T_Influence_Iterator iterator(bbox_min, bbox_max, T_INDEX3()); iterator.Valid(); iterator.Next())
                    {
                        hierarchy3->Activate_Cell(0, iterator.Current_Cell(), Node_Active);
                    }
                }
                hierarchy3->Update_Block_Offsets();
            }
        }
        else 
        {
            std::cout << "Otherwise" << std::endl;
            using CPM_Rasterizer = CPM_Rasterizer<Struct_type, data_type, 3>;
            CPM_Rasterizer cpm_rasterizer(*hierarchy3, m_surface, m_tube_radius, m_grid.dx(), min_cornerd);
            for(Grid_Hierarchy_Iterator<3, CPM_Rasterizer> iterator(hierarchy3->Lattice(base_levels-1).Node_Indices(), base_levels - 1, cpm_rasterizer);iterator.Valid();iterator.Next());
            hierarchy3->Update_Block_Offsets();
        }
    }
#endif
    BuildTube();
}


// build the computational tube based on implicit N1x, N1y, N1z uniform grid (but don't actually construct the full grid!)
void Tube::BuildTube()
{
#ifdef TRACK_WHERE
    std::cout << "Tube::BuildTube()" << std::endl;
#endif
    m_num_nodes = 0;

    // first find one grid node within the tube that is close to the surface
    vector<size_t> idx = FindStartingPoint(); 
    size_t index = m_grid.LinearIndexFromTriplet(idx);
    // add this first node to the tube
#ifdef USE_SPARSE_GRID
    if(dim() == 2)
    {
        SPGrid_Initializer<Struct_type, data_type, 2>(hierarchy2->Allocator(0), hierarchy2->Blocks(0), &Struct_type::ch_tube_id_lower, (data_type)-1.);
        SPGrid_Initializer<Struct_type, data_type, 2>(hierarchy2->Allocator(0), hierarchy2->Blocks(0), &Struct_type::ch_tube_id_higher, (data_type)0.);
    }
    else if(dim() == 3)
    {
        SPGrid_Initializer<Struct_type, data_type, 3>(hierarchy3->Allocator(0), hierarchy3->Blocks(0), &Struct_type::ch_tube_id_lower, (data_type)-1.);
        SPGrid_Initializer<Struct_type, data_type, 3>(hierarchy3->Allocator(0), hierarchy3->Blocks(0), &Struct_type::ch_tube_id_higher, (data_type)0.);
    }
#ifdef COMPARE_GRIDS
    m_tube_indices_on_grid.assign(m_grid.N(), -1); // vector that returns the index in the tube from the index in the implicit rectangular grid. Also used to tell you if you visited the node already for the breadth-first traversal below
#endif
#else
    m_tube_indices_on_grid.assign(m_grid.N(), -1); // vector that returns the index in the tube from the index in the implicit rectangular grid. Also used to tell you if you visited the node already for the breadth-first traversal below
#endif
    vector<Scalar> x(m_grid.dim());
    vector<Scalar> cpx(m_grid.dim());
    Scalar dist;
    size_t bdy;
    
    x = m_grid.CoordinateFromTriplet(idx);
    m_surface.ClosestPoint(x, cpx, dist, bdy);

    AddTubeNode(index, x, cpx, dist, bdy);

    // do a breadth-first traversal to construct the tube
    // grow outward by moving to neighbouring grid nodes within the tube
    list<size_t> queue;
    queue.push_back(m_tube_nodes[0].recIndex());
 
#ifdef USE_SPARSE_GRID
    if(dim() == 2)
    {
        auto tube_id_lower = hierarchy2->Channel(0, &Struct_type::ch_tube_id_lower);
        auto tube_id_higher = hierarchy2->Channel(0, &Struct_type::ch_tube_id_higher);
        while(!queue.empty())
        {
            for(size_t shift_dim = 0; shift_dim < m_grid.dim(); ++shift_dim)
            {
                for(size_t shift : {1, -1})
                {
                    size_t nbr = m_grid.NeighbourLinearIndex(queue.front(), shift, shift_dim);
                    vector<size_t> nidx = ID2Vector(nbr);
                    size_t idxx = nidx[0], idxy = nidx[1];
                    ssize_t cur_id = Merge_Digits(tube_id_higher(idxx, idxy), tube_id_lower(idxx, idxy));
                    if(cur_id == -1)
                    {   
                    #ifdef COMPARE_GRIDS
                        if(m_tube_indices_on_grid[nbr] != -1)
                        {
                            Log::cout << "DIFF 2" << std::endl;
                            getchar();
                        }
                    #endif
                        tube_id_lower(idxx, idxy) = data_type(-2);
                        tube_id_higher(idxx, idxy) = data_type(0);
                    #ifdef COMPARE_GRIDS
                        m_tube_indices_on_grid[nbr] = -2; // give -2 to nodes that are visited but not within the tube
                    #endif
                        idx = m_grid.TripletFromLinearIndex(nbr);
                        x = m_grid.CoordinateFromTriplet(idx);
                        m_surface.ClosestPoint(x, cpx, dist, bdy);

                        if(abs(dist) <= m_tube_radius) // is neighbour inside the tube radius
                        {
                            queue.push_back(nbr); // add this neighbour to the queue to process its neighbours after

                            AddTubeNode(nbr, x, cpx, dist, bdy);
                        }
                    }
                }
            }
            queue.pop_front();
        }
    }
    else if(dim() == 3)
    {
        auto tube_id_lower = hierarchy3->Channel(0, &Struct_type::ch_tube_id_lower);
        auto tube_id_higher = hierarchy3->Channel(0, &Struct_type::ch_tube_id_higher);
        while(!queue.empty())
        {
            for(size_t shift_dim = 0; shift_dim < m_grid.dim(); ++shift_dim)
            {
                for(size_t shift : {1, -1})
                {
                    size_t nbr = m_grid.NeighbourLinearIndex(queue.front(), shift, shift_dim);
                    vector<size_t> nidx = ID2Vector(nbr);
                    size_t idxx = nidx[0], idxy = nidx[1], idxz = nidx[2];
                    ssize_t cur_id = Merge_Digits(tube_id_higher(idxx, idxy, idxz), tube_id_lower(idxx, idxy, idxz));
                    if(cur_id == -1)
                    {   
                    #ifdef COMPARE_GRIDS
                        if(m_tube_indices_on_grid[nbr] != -1)
                        {
                            Log::cout << "DIFF 2" << std::endl;
                            getchar();
                        }
                    #endif
                        tube_id_lower(idxx, idxy, idxz) = data_type(-2);
                        tube_id_higher(idxx, idxy, idxz) = data_type(0);
                    #ifdef COMPARE_GRIDS
                        m_tube_indices_on_grid[nbr] = -2; // give -2 to nodes that are visited but not within the tube
                    #endif
                        idx = m_grid.TripletFromLinearIndex(nbr);
                        x = m_grid.CoordinateFromTriplet(idx);

                        m_surface.ClosestPoint(x, cpx, dist, bdy);

                        if(abs(dist) <= m_tube_radius) // is neighbour inside the tube radius
                        {
                            queue.push_back(nbr); // add this neighbour to the queue to process its neighbours after

                            AddTubeNode(nbr, x, cpx, dist, bdy);
                        }
                    }
                }
            }
            queue.pop_front();
        }
    }
#else
    while(!queue.empty())
    {
        for(size_t shift_dim = 0; shift_dim < m_grid.dim(); ++shift_dim)
        {
            for(size_t shift : {1, -1})
            {
                size_t nbr = m_grid.NeighbourLinearIndex(queue.front(), shift, shift_dim);
                if(m_tube_indices_on_grid[nbr] == -1)
                {
                    m_tube_indices_on_grid[nbr] = -2; // give -2 to nodes that are visited but not within the tube
                    idx = m_grid.TripletFromLinearIndex(nbr);
                    x = m_grid.CoordinateFromTriplet(idx);

                    m_surface.ClosestPoint(x, cpx, dist, bdy);
                    
                    if(abs(dist) <= m_tube_radius) // is neighbour inside the tube radius
                    {
                        queue.push_back(nbr); // add this neighbour to the queue to process its neighbours after

                        AddTubeNode(nbr, x, cpx, dist, bdy);
                    }
                }
            }
        }
        queue.pop_front();
    }
#endif
    
    // now set the neighbours of all the nodes in the tube
    SetNeighbours();
    m_base_points.clear();
    m_base_points = ComputeBasePoints(m_cpx);
}


vector<size_t> Tube::ComputeBasePoints(const vector<vector<Scalar>> &xq)
{
#ifdef TRACK_WHERE
    std::cout << "Tube::ComputeBasePoints()" << std::endl;
#endif
    vector<size_t> base_points(xq.size()); // vector of base point index for each query point xq

    for(size_t i = 0; i < xq.size(); ++i)
    {
        vector<size_t> bpt_triplet = InterpBasePointTriplet(xq[i]); // rectangular grid indices
        size_t bpt_idx = m_grid.LinearIndexFromTriplet(bpt_triplet); // rectangular grid linear index
        base_points[i] = TubeIndexFromRecIndex(bpt_idx); // tube index

        if(base_points[i] < 0 || base_points[i] >= m_num_nodes)
        {
            throw std::runtime_error("Tube::ComputeBasePoints - Interpolation base point is not within the tube. Maybe a query point is outside the tube? Or maybe a problem with the tube construction based on tube radius estimate.");
        }
    }
    return base_points;
}


// Asymmetric picture
// p=0:   B====
// p=1:   B====x
// p=2:   B    x====x
// p=3:   B    x====x    x
// p=4:   B    x    x====x    x
// p=5:   B    x    x====x    x    x
// etc
// Symmetric picture
// p=0: ==B==
// p=1:   B====x
// p=2:   B  ==x==  x
// p=3:   B    x====x    x
// p=4:   B    x  ==x==  x    x
// p=5:   B    x    x====x    x    x
// etc
vector<size_t> Tube::InterpBasePointTriplet(const vector<Scalar> &xq)
{    
    vector<size_t> triplet(m_grid.dim());

    if(p() % 2 == 0)
    {
        for(size_t d = 0; d < m_grid.dim(); ++d)
        {
            // change below round to floor will give asymmetric interpolation stencil, but this breaks the stencils for even p in 3D, had to increase tube radius for it to work
            triplet[d] = round((xq[d] - xstart()[d]) / m_grid.dx()); // index of nearest lower left (actually, it is lower left after the shifting below right?) grid point in the interpolation stencil
        }
    }
    else
    {
        for(size_t d = 0; d < m_grid.dim(); ++d)
        {
            triplet[d] = floor((xq[d] - xstart()[d]) / m_grid.dx()); // index of nearest lower left (actually, it is lower left after the shifting below right?) grid point in the interpolation stencil
        }
    }
    
    // move neighbours from the base point above to get to the correct lower left corner of the interpolation stencil, the above only gives a nearest grid point to the query point
    size_t num_shifts = (p() % 2 == 0) ? p()/2 : (p()-1)/2;
    for(size_t d = 0; d < m_grid.dim(); ++d)
    {
        triplet[d] -= num_shifts;  // move negative one neighbour, p/2 times, in each direction
        assert(triplet[d] >= 0 && triplet[d] < m_grid.N1d()); // grid point is outside the virtual rectangular grid
    }

    return triplet;
}


size_t Tube::InterpBasePointLinearIndex(vector<Scalar> &xq)
{  
    return m_grid.LinearIndexFromTriplet(InterpBasePointTriplet(xq));
}


vector<size_t> Tube::FindStartingPoint() // TO DO: speed this up by just computing a closest point to the surface from some starting x, then find the nearest grid point to the cp(x), e.g., using floor. Like is done for the ibc
{
    vector<Scalar> cpx(m_grid.dim());
    Scalar dist;
    size_t bdy;

    m_surface.ClosestPoint(m_grid.xstart(), cpx, dist, bdy);

    vector<size_t> closest_idx = m_grid.NearestTriplet(cpx);

    return closest_idx;
}


void Tube::AddTubeNode(size_t rec_index, vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy)
{
#ifdef USE_SPARSE_GRID
    if(m_grid.dim() == 2)
    {
        auto tube_id_lower = hierarchy2->Channel(0, &Struct_type::ch_tube_id_lower);
        auto tube_id_higher = hierarchy2->Channel(0, &Struct_type::ch_tube_id_higher);

        vector<size_t> index = ID2Vector(rec_index);
        size_t idxx = index[0], idxy = index[1];

        ssize_t cur_id = Merge_Digits(tube_id_higher(idxx, idxy), tube_id_lower(idxx, idxy));
        tube_id_higher(idxx, idxy) = data_type(ssize_t(m_num_nodes / 1000000));
        tube_id_lower(idxx, idxy) = data_type(ssize_t(m_num_nodes) - ssize_t(m_num_nodes/1000000) * ssize_t(1000000));

    #ifdef COMPARE_GRIDS
        m_tube_indices_on_grid[rec_index] = m_num_nodes;
    #endif
    }
    if(m_grid.dim() == 3)
    {
        auto tube_id_lower = hierarchy3->Channel(0, &Struct_type::ch_tube_id_lower);
        auto tube_id_higher = hierarchy3->Channel(0, &Struct_type::ch_tube_id_higher);

        vector<size_t> index = ID2Vector(rec_index);
        size_t idxx = index[0], idxy = index[1], idxz = index[2];

        ssize_t cur_id = Merge_Digits(tube_id_higher(idxx, idxy, idxz), tube_id_lower(idxx, idxy, idxz));
        tube_id_higher(idxx, idxy, idxz) = data_type(ssize_t(m_num_nodes / 1000000));
        tube_id_lower(idxx, idxy, idxz) = data_type(ssize_t(m_num_nodes) - ssize_t(m_num_nodes/1000000) * ssize_t(1000000));

    #ifdef COMPARE_GRIDS
        m_tube_indices_on_grid[rec_index] = m_num_nodes;
    #endif
    }
#else
    m_tube_indices_on_grid[rec_index] = m_num_nodes;
#endif

    Node tube_node;
    tube_node.setRecIndex(rec_index);
    tube_node.setTubeIndex(m_num_nodes);

    m_tube_nodes.push_back(tube_node);
    m_tube_coords.push_back(x);

    m_cpx.push_back(cpx);
    m_dist.push_back(dist);
    m_bdy.push_back(bdy);

    m_num_nodes++;
}

#ifdef USE_SPARSE_GRID
vector<size_t> Tube::ID2Vector(size_t id) const
{
    std::vector<size_t> idx = TripletFromLinearIndex(id);
    for(int axis=0;axis<dim();++axis)
    {
        idx[axis] += 1;
    }
    return idx;
}
#endif


void Tube::SetNeighbours()
{
#ifdef TRACK_WHERE
    std::cout << "Tube::SetNeighbours" << std::endl;
#endif
#ifdef USE_SPARSE_GRID
    if(dim() == 2)
    {
        auto tube_id_lower = hierarchy2->Channel(0, &Struct_type::ch_tube_id_lower);
        auto tube_id_higher = hierarchy2->Channel(0, &Struct_type::ch_tube_id_higher);
        auto flags = hierarchy2->Channel(0, &Struct_type::flags);
        vector<Scalar> shift{1.0, -1.0};

        for(size_t i = 0; i < m_num_nodes; ++i)
        {
            for(size_t d = 0; d < m_grid.dim(); ++d)
            {
                for(size_t shift_idx = 0; shift_idx < 2; ++shift_idx)
                {
                    size_t nbr = m_grid.NeighbourLinearIndex(m_tube_nodes[i].recIndex(), shift[shift_idx], d);
                    vector<size_t> nidx = ID2Vector(nbr);
                    size_t idxx = nidx[0], idxy = nidx[1];
                    ssize_t cur_id = Merge_Digits(tube_id_higher(idxx, idxy), tube_id_lower(idxx, idxy));
                    if((flags(idxx, idxy) & Node_Active) && cur_id >= 0 && cur_id < m_num_nodes)
                    {
                        m_tube_nodes[i].setNeighbour(2 * d + shift_idx, &m_tube_nodes[cur_id]);
                    }
                }
            }
        }
    }
    else if(dim() == 3)
    {
        auto tube_id_lower = hierarchy3->Channel(0, &Struct_type::ch_tube_id_lower);
        auto tube_id_higher = hierarchy3->Channel(0, &Struct_type::ch_tube_id_higher);
        auto flags = hierarchy3->Channel(0, &Struct_type::flags);
        vector<Scalar> shift{1.0, -1.0};

        for(size_t i = 0; i < m_num_nodes; ++i)
        {
            for(size_t d = 0; d < m_grid.dim(); ++d)
            {
                for(size_t shift_idx = 0; shift_idx < 2; ++shift_idx)
                {
                    size_t nbr = m_grid.NeighbourLinearIndex(m_tube_nodes[i].recIndex(), shift[shift_idx], d);
                    vector<size_t> nidx = ID2Vector(nbr);
                    size_t idxx = nidx[0], idxy = nidx[1], idxz = nidx[2];
                    ssize_t cur_id = Merge_Digits(tube_id_higher(idxx, idxy, idxz), tube_id_lower(idxx, idxy, idxz));
                    if((flags(idxx, idxy, idxz) & Node_Active) && cur_id >= 0 && cur_id < m_num_nodes)
                    {
                        m_tube_nodes[i].setNeighbour(2 * d + shift_idx, &m_tube_nodes[cur_id]);
                    }
                }
            }
        }
    }
#else
    vector<Scalar> shift{1.0, -1.0};

    for(size_t i = 0; i < m_num_nodes; ++i)
    {
        for(size_t d = 0; d < m_grid.dim(); ++d)
        {
            for(size_t shift_idx = 0; shift_idx < 2; ++shift_idx)
            {
                size_t nbr = m_grid.NeighbourLinearIndex(m_tube_nodes[i].recIndex(), shift[shift_idx], d);

                if(m_tube_indices_on_grid[nbr] >= 0 && m_tube_indices_on_grid[nbr] < m_num_nodes)
                {
                    m_tube_nodes[i].setNeighbour(2 * d + shift_idx, &m_tube_nodes[m_tube_indices_on_grid[nbr]]);
                }
            }
        }
    }
#endif
}

} // namespace cpm