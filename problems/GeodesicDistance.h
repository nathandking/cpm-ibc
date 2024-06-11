#pragma once

#include "Scalar.h"
#include "Tube.h"
#include "IBCMatrixManipulations.h"
#include "IBCSubsets.h"
#include "Helpers.h"
#include "FDMatrices.h"
#include "Defines.h"
#include "Solver.h"

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/SparseCore>

using namespace std;

namespace cpm {

class GeodesicDistance{
    public:
        GeodesicDistance(bool visualize = false);

        VectorX ComputeDistance(Tube &surface_tube, IBCSubsets &ibc, Scalar dt);

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

        Helpers m_h;
        bool m_visualize;
        
        IBCMatrixManipulations m_ibc_mat_manip;

        SpMat m_L;
        SpMat m_E;
        SpMat m_Lfull;
        SpMat m_Efull;
        vector<SpMat> m_Dc;
        vector<SpMat> m_Dcfull;
        vector<bool> m_identity_rows;

        VectorX Step1(Tube &surface_tube, IBCSubsets &ibc, Scalar dt);
        vector<VectorX> Step2(VectorX &ufull, Tube &surface_tube, IBCSubsets &ibc);
        VectorX ComputeDivergence(vector<VectorX> &gradu, Tube &surface_tube, IBCSubsets &ibc);
        VectorX Step3(VectorX &div, Tube &surface_tube, IBCSubsets &ibc);

        VectorX ExactMatchStep1ScalarField(Tube &surface_tube, IBCSubsets &ibc);
        vector<VectorX> ExactGradientSpherePoint(Tube &surface_tube, IBCSubsets &ibc);
        VectorX ExactDivergenceSpherePoint(Tube &surface_tube, IBCSubsets &ibc);
        VectorX ExactDivergenceCirclePoint(Tube &surface_tube, IBCSubsets &ibc);

        vector<vector<Scalar>> m_surface_normals;
};

} // namespace cpm
