#pragma once

#include <iostream>
#include <vector>
#include "math.h"

#include "Scalar.h"
#include "Tube.h"

#include <Eigen/Sparse>

using namespace std;

namespace cpm {

class Interpolation{
    public:
        
        Interpolation(Tube &t); // for interpolation onto closest points of tube
        Interpolation(Tube &t, const vector<vector<Scalar>> &xquery); // for interpolation of data at any point within the tube from data on the tube

        void BuildInterpolationMatrix(Tube &t, SpMat &E);
        void BuildInterpolationMatrix(Tube &t, const vector<ssize_t> &subset_idx_from_set_idx, SpMat &E);

    private:

        int m_stencil_size;
        size_t m_Nq;
        vector<vector<Scalar>> m_xq;
        vector<vector<vector<Scalar>>> m_w; // interpolation weights for each dimension
        vector<Scalar> m_uniform_grid_w;
        vector<size_t> m_Ibpt;

        // computing weights
        void UniformGridWeights1D();

        void BuildInterpolationWeights1D(vector<Scalar> &x, Scalar &xq_subset, vector<Scalar> &w);

        void BuildInterpolationWeights(Tube &t);

        // testing
        void TestRowSumOne(SpMat &E);
};

} // namespace cpm