#pragma once

#include <iostream>
#include <vector>
#include "math.h"

#include "Scalar.h"
#include "Tube.h"
#include "Interpolation.h"
#include "FDMatrices.h"

#include <Eigen/Sparse>

using namespace std;

namespace cpm {

class Geometry{
    public:
        
        Geometry();
        
        template <typename M, typename V>
        void EigenDecompositionOfJcp(Tube &t, vector<M> &eigenvectors, vector<V> &eigenvalues);

        template <typename M, typename V>
        void EigenDecompositionOfJcp(Tube &t, const vector<vector<Scalar>> &subset_cpx, const vector<size_t> &set_idx_from_subset_idx, const vector<ssize_t> &subset_idx_from_set_idx, vector<M> &eigenvectors, vector<V> &eigenvalues);

        template <typename M, typename V>
        void EigenDecomposition(Tube &t, vector<vector<VectorX>> &J, vector<M> &eigenvectors, vector<V> &eigenvalues);

        vector<vector<Scalar>> UnorientedSurfaceNormals(Tube &t, bool is_normalized = true, Scalar zero_tol = 1e-15);

        void InterpolateUnorientedVectorFromNeighbours(Tube &t, vector<vector<Scalar>> &vectors, size_t interp_index, Scalar zero_tol = 1e-15);
        void InterpolateUnorientedVectorFromNeighbours(Tube &t, vector<vector<Scalar>> &subset_vectors, const vector<size_t> &set_idx_from_subset_idx, const vector<ssize_t> &subset_idx_from_set_idx, size_t subset_interp_index, Scalar zero_tol = 1e-15);

        vector<vector<Scalar>> InterpolateUnorientedVectors(Tube &t, const vector<vector<Scalar>> &vectors, const vector<vector<Scalar>> &xq, bool is_normalized = true);
        vector<vector<Scalar>> InterpolateUnorientedVectors(Tube &t, const vector<size_t> &set_idx_from_subset_idx, const vector<ssize_t> &subset_idx_from_set_idx, const vector<vector<Scalar>> &vectors, const vector<vector<Scalar>> &xq, bool is_normalized = true);

        void ApplyInterpolationToUnorientedVectors(const SpMat &E, const vector<vector<Scalar>> &v, vector<vector<Scalar>> &vq);

    private:

};

} // namespace cpm