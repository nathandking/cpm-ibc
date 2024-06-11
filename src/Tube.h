#pragma once

#include <vector>
#include <iostream>
#include <limits>
#include <list>

#include "Scalar.h"
#include "Node.h"
#include "Surface.h"
#include "Grid.h"

#ifdef USE_SPARSE_GRID
#include "CPM_Data.h"
#include <nova/Dynamics/Hierarchy/Grid_Hierarchy.h>
using namespace SPGrid;
using namespace Nova;
#endif

using namespace std;

namespace cpm {

typedef float data_type;

class Tube {

#ifdef USE_SPARSE_GRID
    using Struct_type       = CPM_Data<data_type>;
    using TV2               = Nova::Vector<data_type, 2>;
    using TV3               = Nova::Vector<data_type, 3>;
    using T_INDEX2          = Nova::Vector<int, 2>;
    using T_INDEX3          = Nova::Vector<int, 3>;
    using T_Range2          = Nova::Range<data_type, 2>;
    using T_Range3          = Nova::Range<data_type, 3>;
    using T_Grid2           = Nova::Grid<data_type, 2>;
    using T_Grid3           = Nova::Grid<data_type, 3>;
    using Hierarchy2        = Grid_Hierarchy<Struct_type, data_type, 2>;
    using Hierarchy3        = Grid_Hierarchy<Struct_type, data_type, 3>;
#endif

    public:
        Tube();
        ~Tube();
        Tube(SurfaceSpecifier &surface_specs, size_t interp_deg, Scalar dx);

        void Clear();

        void ConstructTube(SurfaceSpecifier &surface_specs, size_t interp_deg, Scalar dx);
        
        size_t dim() const { return m_grid.dim(); };
        size_t p() const { return m_interp_deg; };
        size_t nNodes() const { return m_num_nodes; };
        Scalar dx() const { return m_grid.dx(); };
        size_t N1d() const { return m_grid.N1d(); };
        size_t n() const { return m_grid.n(); }
        vector<Scalar> xstart() const { return m_grid.xstart(); };
        Scalar TubeRadiusFactor() const { return m_tube_radius_factor; };
        Scalar TubeRadius() const { return m_tube_radius; };

        Node TubeNode(size_t index) const { return m_tube_nodes[index]; };

        void ClosestPoint(const vector<Scalar> &x, vector<Scalar> &cpx){ m_surface.ClosestPoint(x, cpx); };
        
        const vector<vector<Scalar>>& x() const 
        { 
            return m_tube_coords; 
        };

        ssize_t TubeIndexFromRecIndex(size_t grid_index)
        {
        #ifdef USE_SPARSE_GRID
            if(dim() == 2)
            {
                auto tube_id_lower = hierarchy2->Channel(0, &Struct_type::ch_tube_id_lower);
                auto tube_id_higher = hierarchy2->Channel(0, &Struct_type::ch_tube_id_higher);
                vector<size_t> idx = ID2Vector(grid_index);
                size_t idxx = idx[0], idxy = idx[1];
                return ssize_t(tube_id_lower(idxx, idxy)) + ssize_t(1000000) * ssize_t(tube_id_higher(idxx, idxy));
            }
            else if(dim() == 3)
            {
                auto tube_id_lower = hierarchy3->Channel(0, &Struct_type::ch_tube_id_lower);
                auto tube_id_higher = hierarchy3->Channel(0, &Struct_type::ch_tube_id_higher);
                vector<size_t> idx = ID2Vector(grid_index);
                size_t idxx = idx[0], idxy = idx[1], idxz = idx[2];
                return ssize_t(tube_id_lower(idxx, idxy, idxz)) + ssize_t(1000000) * ssize_t(tube_id_higher(idxx, idxy, idxz));
            }
            
        #else
            return m_tube_indices_on_grid[grid_index];
        #endif
        }

        vector<size_t> TripletFromLinearIndex(size_t linear_rec_index) const
        {
            return m_grid.TripletFromLinearIndex(linear_rec_index);
        }

        size_t LinearIndexFromTriplet(vector<size_t> triplet)
        {
            return m_grid.LinearIndexFromTriplet(triplet);
        }

        size_t NearestTubeIndex(const vector<Scalar> &xq)
        {
            vector<size_t> triplet = m_grid.NearestTriplet(xq);
            size_t linear_idx = m_grid.LinearIndexFromTriplet(triplet);
            return TubeIndexFromRecIndex(linear_idx);
        }

        ssize_t Merge_Digits(ssize_t hd, ssize_t ld)
        {
            return ssize_t(1000000) * hd + ld;
        }

        const Surface& surface() const
        {
            return m_surface;
        };

        const vector<size_t>& bdy() const 
        { 
            return m_bdy; 
        };

        const vector<vector<Scalar>>& cpx() const
        { 
            return m_cpx; 
        };

        void replaceClosestPoint(size_t index, vector<Scalar> &cpx);

        const vector<Scalar>& dist() const
        { 
            return m_dist; 
        };

        const vector<size_t>& bpt() const 
        { 
            return m_base_points; 
        };

        vector<size_t> ComputeBasePoints(const vector<vector<Scalar>> &xq);

        void cpBar(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist){ m_surface.cpBar(x, cpx, dist); };

        void OverideParameterizationPoints(const vector<vector<Scalar>> xp){ m_surface.OverideParameterizationPoints(xp); };
        void OverideParameterizationFaces(const vector<vector<size_t>> faces){ m_surface.OverideParameterizationFaces(faces); };
        
        const Scalar IBCParameter(const vector<Scalar>& x){ return m_surface.IBCParameter(x); };

    private:

        Surface m_surface;
        Grid m_grid;

        void ConstructTube(size_t interp_deg, Scalar dx);

        void ComputeTubeRadius();
        void ConstructTube();
        void BuildTube();

#ifdef USE_SPARSE_GRID
        Hierarchy2* hierarchy2;
        Hierarchy3* hierarchy3;
        vector<size_t> ID2Vector(size_t id) const;
#endif

        vector<size_t> InterpBasePointTriplet(const vector<Scalar> &xq);
        size_t InterpBasePointLinearIndex(vector<Scalar> &xq);
        
        vector<size_t> m_base_points;
        vector<bool> m_is_bpt_visited;

        vector<size_t> FindStartingPoint();
        void AddTubeNode(size_t rec_index, vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy);
        void SetNeighbours();

        size_t m_interp_deg; // polynomial interpolation order
        size_t m_order; // Laplacian order, used for tube radius calculation. Might depend on gradient order or largest FD stencil

        vector<vector<Scalar>> m_tube_coords; // coordinates of the tube nodes
        vector<Node> m_tube_nodes; // grid nodes in the computational tube
#ifdef USE_SPARSE_GRID
    #ifdef COMPARE_GRIDS
        vector<ssize_t> m_tube_indices_on_grid; // indicies with respect to tube nodes, but stored for each rectangular grid index
    #endif
#endif
#ifndef USE_SPARSE_GRID
        vector<ssize_t> m_tube_indices_on_grid; // indicies with respect to tube nodes, but stored for each rectangular grid index
#endif

        // tubing
        Scalar m_tube_radius_factor; // how many dx spacings should the tube radius be
        Scalar m_tube_radius;
        size_t m_num_nodes; // total number of grid points in computational tube
        
        vector<Scalar> m_min_corner;
        vector<Scalar> m_max_corner;

        vector<vector<Scalar>> m_cpx; // closest points to nodes in the computational tube
        vector<Scalar> m_dist;
        vector<size_t> m_bdy; // tag for which boundary the cp in tube belongs to, 0 if not a boundary cp
};

} // namespace cpm