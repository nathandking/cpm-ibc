#pragma once

#include <vector>
#include <iostream>

#include "VectorMath.h"
#include "Surface.h"
#include "Tube.h"
#include "IBCSubset.h"
#include "IBCDOFSubset.h"
#include "IBCMetaData.h"
#include "Geometry.h"

using namespace std;

namespace cpm {

class IBCSubsets {
    
    public:
        IBCSubsets();
        
        void Initialize(Tube &surface_tube, vector<SurfaceSpecifier> ibc_surface_specs);

        const IBCMetaData& meta(size_t ibc_index) const 
        { 
            return m_ibc_DOF_subset[ibc_index].meta();
        }; 

        size_t TotalIBCDOFs(){ return m_total_ibc_DOFs; };

        size_t NumIBCs(){ return m_num_ibcs; };

        size_t IBCColumnStartIndex(size_t ibc_index){ return m_ibc_starting_column[ibc_index]; };

        size_t DirichletIBCColumnStartIndex(size_t ibc_index){ return m_dirichlet_ibc_starting_column[ibc_index]; };

        const Surface& surface(size_t ibc_index) const
        {
            return m_ibc_surface[ibc_index];
        };

        Scalar IBCParameter(size_t ibc_index, const vector<Scalar>& x)
        {
            return m_ibc_surface[ibc_index].IBCParameter(x);
        }

        vector<Scalar> IBCTangent(size_t ibc_index, const vector<Scalar>& x)
        {
            return m_ibc_surface[ibc_index].IBCTangent(x);
        }

        const vector<vector<Scalar>>& xp(size_t ibc_index) const
        {
            return m_ibc_surface[ibc_index].xp();
        }

        void ClosestPoint(size_t ibc_index, const vector<Scalar> &x, vector<Scalar> &cpx){ m_ibc_surface[ibc_index].ClosestPoint(x, cpx); };
        void ClosestPoint(size_t ibc_index, const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy){ m_ibc_surface[ibc_index].ClosestPoint(x, cpx, dist, bdy); };

        const size_t& IdentityRowsStart() const 
        { 
            return m_identity_rows_start; 
        };

        const size_t& Num2ndOrderDirichletRows() const 
        { 
            return m_num_2nd_order_dirichlet_rows; 
        };

        const vector<IBCSubset>& TaggingSubset() const 
        { 
            return m_ibc_subset; 
        };

        const vector<IBCDOFSubset>& DOFSubset() const 
        { 
            return m_ibc_DOF_subset; 
        };

    private:

        size_t m_num_ibcs;
        vector<Surface> m_ibc_surface;
        vector<IBCSubset> m_ibc_subset;
        vector<IBCDOFSubset> m_ibc_DOF_subset;

        vector<size_t> m_ibc_starting_column;
        vector<size_t> m_dirichlet_ibc_starting_column;
        size_t m_total_ibc_DOFs;
        size_t m_identity_rows_start;
        size_t m_num_2nd_order_dirichlet_rows;

        void Initialize(Tube &surface_tube);

        void ComputeSubsets(Tube& surface_tube);
        void IBCStartingColumns(Tube& surface_tube);
        void TotalDOFs();

        void VisualizeSperp(Tube &surface_tube, size_t ibc_index);
        
        void VisualizeSubsets(Tube &surface_tube, size_t ibc_index);

        template <typename T>
        void VisualizeSubset(Tube &surface_tube, const T &subset, string subset_name, size_t ibc_index);
};

} // namespace cpm