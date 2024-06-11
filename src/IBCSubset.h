#pragma once

#include <vector>
#include <iostream>

#include "VectorMath.h"
#include "Surface.h"
#include "Tube.h"
#include "Visualization.h"
#include "IBCMetaData.h"
#include "Interpolation.h"
#include "Geometry.h"

using namespace std;

namespace cpm {

class IBCSubset {
    
    public:
        IBCSubset();

        void ComputeSubset(Surface& ibc_surface, Tube &surface_tube);
        
        const vector<size_t>& SetIndexFromSubsetIndex() const
        {
            return m_set_idx_from_subset_idx;
        }

        const vector<ssize_t>& SubsetIndexFromSetIndex() const
        {
            return m_subset_idx_from_set_idx;
        }

        const size_t size() const
        { 
            return m_set_idx_from_subset_idx.size(); 
        };

        const vector<vector<Scalar>>& cpx() const
        { 
            return m_cpx_c; 
        };

        const vector<vector<Scalar>>& cp_diff() const
        { 
            return m_cp_diff; 
        };

        const vector<vector<Scalar>>& surface_normals_on_tube() const
        { 
            return m_surface_normals_on_tube; 
        };

        const vector<vector<Scalar>>& surface_normals() const
        { 
            return m_surface_normals; 
        };

        const vector<vector<Scalar>>& tangents() const
        { 
            return m_tangents; 
        };

        const vector<vector<Scalar>>& binormals() const
        { 
            return m_binormals; 
        };

        const vector<vector<Scalar>>& directions() const
        { 
            return m_directions; 
        };

        const vector<Scalar>& dist() const
        { 
            return m_dist_c; 
        };

        const vector<size_t>& bdy() const 
        { 
            return m_bdy; 
        };

        const IBCMetaData& meta() const { return m_ibc_meta; };

    private:
        Geometry m_geom;

        IBCMetaData m_ibc_meta;

        vector<vector<Scalar>> m_cpx_c;
        vector<vector<Scalar>> m_cp_diff;
        vector<Scalar> m_dist_c;
        vector<size_t> m_bdy;

        vector<vector<Scalar>> m_surface_normals_on_tube;
        vector<vector<Scalar>> m_surface_normals;
        vector<vector<Scalar>> m_tangents;
        vector<vector<Scalar>> m_binormals;
        vector<vector<Scalar>> m_directions;

        vector<ssize_t> m_subset_idx_from_set_idx;
        vector<size_t> m_set_idx_from_subset_idx;
        
        void ComputeSubsetNearIBC(Surface& ibc_surface, Tube& surface_tube);
        size_t FindStartingPoint(Surface& ibc_surface, Tube& surface_tube);
        
        void ComputeInverseMapping(size_t set_size, const vector<size_t> &set_from_subset, vector<ssize_t> &subset_from_set);

        vector<size_t> ComputeSubsetNearIBCBruteForce(Surface& ibc_surface, Tube &surface_tube);

        void ComputeCPdiff(Surface &ibc_surface, Tube &surface_tube);
        void ComputeFrameBasedDirections(Tube &surface_tube);
};

} // namespace cpm