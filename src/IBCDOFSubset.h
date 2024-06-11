#pragma once

#include <vector>
#include <iostream>

#include "VectorMath.h"
#include "Surface.h"
#include "Tube.h"
#include "IBCMetaData.h"
#include "IBCSubset.h"
#include "Interpolation.h"
#include <Eigen/Sparse>

using namespace std;

namespace cpm {

class IBCDOFSubset {
    
    public:
        IBCDOFSubset();

        void ComputeDOFSubset(Surface& ibc_surface, IBCSubset &ibc_subset, Tube &surface_tube);
        
        const vector<size_t>& SetIndexFromSubsetIndex() const
        {
            return m_set_idx_from_subset_idx;
        }

        const vector<ssize_t>& SubsetIndexFromSetIndex() const
        {
            return m_subset_idx_from_set_idx;
        }

        const vector<size_t>& DirichletSetIndexFromSubsetIndex() const
        {
            return m_dirichlet_set_idx_from_subset_idx;
        }

        const vector<ssize_t>& DirichletSubsetIndexFromSetIndex() const
        {
            return m_dirichlet_subset_idx_from_set_idx;
        }

        const size_t size() const
        { 
            return m_set_idx_from_subset_idx.size(); 
        };

        const vector<vector<Scalar>>& cpx() const
        { 
            return m_cpx_c; 
        };

        const vector<vector<Scalar>>& cpbar() const
        { 
            return m_cp_c_bar; 
        };

        const vector<vector<Scalar>>& cp_diff() const
        { 
            return m_cp_diff; 
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

        const vector<size_t>& which_side() const 
        { 
            return m_which_side; 
        };  

        const SpMat& cpExtensionMatrix() const 
        { 
            return m_E; 
        };

        const IBCMetaData& meta() const { return m_ibc_meta; };

    private:
        IBCMetaData m_ibc_meta;

        vector<vector<Scalar>> m_cp_c_bar;

        vector<vector<Scalar>> m_cpx_c;
        vector<vector<Scalar>> m_cp_diff;
        vector<Scalar> m_dist_c;
        vector<size_t> m_bdy;
        vector<size_t> m_which_side;

        vector<vector<Scalar>> m_surface_normals;
        vector<vector<Scalar>> m_tangents;
        vector<vector<Scalar>> m_binormals;
        vector<vector<Scalar>> m_directions;

        vector<ssize_t> m_subset_idx_from_set_idx;
        vector<size_t> m_set_idx_from_subset_idx;

        vector<ssize_t> m_dirichlet_subset_idx_from_set_idx;
        vector<size_t> m_dirichlet_set_idx_from_subset_idx;

        SpMat m_E;
        
        void ComputeInverseMapping(size_t set_size, const vector<size_t> &set_from_subset, vector<ssize_t> &subset_from_set);

        void ComputeExtensionMatrix(Tube& surface_tube);

        void OrientSubset(Surface& ibc_surface, Tube &surface_tube);
};

} // namespace cpm