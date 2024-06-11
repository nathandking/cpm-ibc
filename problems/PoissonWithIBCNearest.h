#pragma once

#include "Scalar.h"
#include "Tube.h"
#include "Helpers.h"
#include "IBCMatrixManipulations.h"
#include "IBCMetaData.h"

#include <Eigen/Sparse> 

using namespace std;

namespace cpm {

class PoissonWithIBCNearest{
    public:
        PoissonWithIBCNearest(SurfaceSpecifier &surface_specs, vector<SurfaceSpecifier> &ibc_specs, Scalar dx, bool freeze_nearest, size_t interp_deg);

        VectorX GetRHS(function<Scalar(const vector<Scalar>&)> f, function<Scalar(const vector<Scalar>&)> surface_exact);

        SpMat GetPlotInterpMatrix();

        size_t dim(){ return m_surface_tube.dim(); };

        const vector<vector<Scalar>>& cpx() const
        { 
            return m_surface_tube.cpx(); 
        };

        const vector<size_t>& bdy() const 
        { 
            return m_surface_tube.bdy(); 
        };

        const vector<Scalar>& thetap() const
        { 
            return m_surface_tube.surface().thetap(); 
        };

        const vector<Scalar>& phip() const
        { 
            return m_surface_tube.surface().phip(); 
        };

        const vector<vector<Scalar>>& xp() const
        { 
            return m_surface_tube.surface().xp(); 
        };

        const vector<vector<size_t>>& faces() const
        { 
            return m_surface_tube.surface().faces(); 
        };

        Scalar IBCParameter(size_t ibc_index, const vector<Scalar>& x)
        {
            return m_ibc.IBCParameter(ibc_index, x);
        }

        const vector<vector<Scalar>>& ibc_xp(size_t ibc_index) const
        { 
            return m_ibc.surface(ibc_index).xp(); 
        };

        const vector<vector<size_t>>& ibc_faces(size_t ibc_index) const
        { 
            return m_ibc.surface(ibc_index).faces(); 
        };

        const vector<Scalar>& surfaceParams() const
        {
            return m_surface_tube.surface().surfaceParams();
        };

        Node node(size_t index) const
        {
            return m_surface_tube.TubeNode(index);
        }

        const Tube& tube() const
        {
            return m_surface_tube;
        }

        const size_t n() const
        {
            return m_surface_tube.n();
        }

        const vector<vector<size_t>>& GetIBCIDX() const
        {
            return m_ibc_idx;
        }

        size_t DOFs(){ return m_surface_tube.nNodes() + m_ibc.TotalIBCDOFs(); };

        size_t NumIBCs(){ return m_ibc.NumIBCs(); };

        // integration
        const IBCSubsets& ibc()
        {
            return m_ibc;
        }

        const bool& FreezeNearest() const 
        {
            return m_freeze_nearest;
        }

        const SpMat GetLfull() const
        {
            return m_L;
        }

        const SpMat GetEfull() const
        {
            return m_E;
        }

        const vector<bool>& GetIdentityRows() const
        {
            return m_identity_rows;
        }
        
    private:
        IBCSubsets m_ibc;
        
        vector<bool> m_identity_rows;
        bool m_freeze_nearest;
        vector<vector<size_t>> m_ibc_idx;

        void Initialize();

        Helpers m_h;

        Tube m_surface_tube;

        SpMat m_L;
        SpMat m_E; 
};

} // namespace cpm
