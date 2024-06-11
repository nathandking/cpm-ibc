#pragma once

#include "Scalar.h"
#include "Tube.h"
#include "Helpers.h"
#include "IBCMetaData.h"
#include "IBCMatrixManipulations.h"

#include "Solver.h"

#include <Eigen/Sparse>

using namespace std;

namespace cpm {

class ReactionDiffusionWithIBC{
    public:
        ReactionDiffusionWithIBC(SurfaceSpecifier &surface_specs, vector<SurfaceSpecifier> &ibc_specs, Scalar dx, Scalar dt, size_t interp_deg);

        void TimeStep(VectorX &ufull, VectorX &vfull, function<Scalar(size_t &, const size_t &, const vector<Scalar>&)> ibc_exact_u, function<Scalar(size_t &, const size_t &, const vector<Scalar>&)> ibc_exact_v);

        SpMat GetPlotInterpMatrix();
        SpMat GetInterpMatrix(vector<vector<Scalar>> &xp);
        
        void ClosestPoint(const vector<Scalar> &x, vector<Scalar> &cpx){ m_surface_tube.ClosestPoint(x, cpx); };

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

        vector<Scalar> IBCTangent(size_t ibc_index, const vector<Scalar>& x)
        {
            return m_ibc.IBCTangent(ibc_index, x);
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


        IBCSubsets& ibc()
        {
            return m_ibc;
        }

        Scalar dt(){ return m_dt; };

        size_t nNodes(){ return m_surface_tube.nNodes(); };

        size_t DOFs(){ return m_surface_tube.nNodes() + m_ibc.TotalIBCDOFs(); };

        size_t NumIBCs(){ return m_ibc.NumIBCs(); };

    private:

        chrono::time_point<chrono::system_clock> m_start, m_end;
        chrono::duration<Scalar> m_elapsed_seconds;

        IBCSubsets m_ibc;
        IBCMatrixManipulations m_ibc_mat_manip;

        void Initialize();

        Helpers m_h;

        Scalar m_dt;

        Tube m_surface_tube;

        SpMat m_I;
        SpMat m_L;
        SpMat m_E; 
        SpMat m_M;
        SpMat m_Lfull;
        SpMat m_Efull; 
        SpMat m_A;

        vector<bool> m_identity_rows;

        void SetDirichletValue(VectorX &ufull, function<Scalar(size_t &, const size_t &, const vector<Scalar>&)> ibc_exact);

        Solver sv;
};

} // namespace cpm
