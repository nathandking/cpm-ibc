#pragma once

#include "Scalar.h"
#include "Tube.h"
#include "Helpers.h"
#include "IBCMatrixManipulations.h"
#include "IBCMetaData.h"

#include <Eigen/Sparse>

using namespace std;

namespace cpm {

class PoissonWithIBC{
    public:
        PoissonWithIBC(SurfaceSpecifier &surface_specs, vector<SurfaceSpecifier> &ibc_specs, Scalar dx, size_t interp_deg);
        
        VectorX GetRHS(function<Scalar(const vector<Scalar>&)> f, function<Scalar(const vector<Scalar>&)> surface_exact);
        VectorX GetRHS(function<Scalar(const vector<Scalar>&)> f, function<Scalar(size_t &, const size_t &, const vector<Scalar>&)> ibc_exact);

        SpMat GetPlotInterpMatrix();
        SpMat GetPlotInterpMatrix(vector<vector<size_t>> &director_set_index_from_subset_index, vector<vector<size_t>> &director_which_side, vector<vector<bool>> &director_is_on_ibc);

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

        size_t DOFs(){ return m_surface_tube.nNodes() + m_ibc.TotalIBCDOFs(); };

        size_t NumIBCs(){ return m_ibc.NumIBCs(); };

        // integration
        const IBCSubsets& ibc()
        {
            return m_ibc;
        }

        const SpMat GetLfull() const
        {
            return m_Lfull;
        }

        const SpMat GetEfull() const
        {
            return m_Efull;
        }

        const vector<bool>& GetIdentityRows() const
        {
            return m_identity_rows;
        }
        
        void BuildExtensionMatrix();
        void OverideParameterizationPoints(const vector<vector<Scalar>> xp){ m_surface_tube.OverideParameterizationPoints(xp); };

    private:

        chrono::time_point<chrono::system_clock> m_start, m_end;
        chrono::duration<Scalar> m_elapsed_seconds;

        IBCSubsets m_ibc;
        IBCMatrixManipulations m_ibc_mat_manip;

        void Initialize();

        Helpers m_h;

        Tube m_surface_tube;

        SpMat m_L;
        SpMat m_E; 

        SpMat m_Lfull;
        SpMat m_Efull;

        vector<bool> m_identity_rows;
};

} // namespace cpm
