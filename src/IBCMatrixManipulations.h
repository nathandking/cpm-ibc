#pragma once

#include <vector>
#include <iostream>

#include "Scalar.h"
#include "VectorMath.h"
#include "Tube.h"
#include "IBCSubsets.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace std;


namespace cpm {

class IBCMatrixManipulations {
    public:
        IBCMatrixManipulations();

        SpMat AddIBCExtension(SpMat &A, Tube &surface_tube, IBCSubsets &ibc, const vector<vector<Scalar>> &director_x);
        SpMat AddIBCFiniteDifference(SpMat &A, Tube &surface_tube, IBCSubsets &ibc, const vector<vector<Scalar>> &director_x);
        SpMat AddIBCPlotInterpolation(SpMat &A, Tube &surface_tube, IBCSubsets &ibc, const vector<vector<Scalar>> &director_x);
        SpMat AddIBCPlotInterpolation(SpMat &A, Tube &surface_tube, IBCSubsets &ibc, const vector<vector<Scalar>> &director_x, vector<vector<size_t>> &director_set_index_from_subset_index, vector<vector<size_t>> &director_which_side, vector<vector<bool>> &director_is_on_ibc);

        SpMat TagStencils(SpMat &A, Tube &surface_tube, IBCSubsets &ibc, const vector<vector<Scalar>> &director_cp_diff, const vector<ssize_t>& director_subset_index_from_set_index, size_t ibc_index, Scalar flip_scale = 1.0);
        
    private:

        // Functions for matrix manipulation based on normal crossings of closest points compared to stencil director
        void MoveToIBCColumns(SpMat &A, vector<SpMat> &A_tags, IBCSubsets &ibc, vector<T> &coeffs);
        void AddExtensionIBCRows(Tube &surface_tube, IBCSubsets &ibc, vector<T> &coeffs);
        void AddFiniteDifferenceIBCRows(SpMat &A, vector<SpMat> &A_tags, vector<Eigen::Triplet<Scalar>> &coeffs, IBCSubsets &ibc);

        void PlotPointsCPDiff(Tube &surface_tube, IBCSubsets &ibc, vector<vector<Scalar>> &subset_cp_diff, vector<size_t> &subset_bdy, vector<ssize_t> &subset_index_from_set_index, vector<size_t> &set_index_from_subset_index, vector<size_t> &which_side, vector<bool> &is_on_ibc, size_t ibc_index);
};

} // namespace cpm