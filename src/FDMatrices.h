#pragma once

#include "Scalar.h"
#include "Tube.h"

#include <Eigen/Sparse>

using namespace std;

namespace cpm {

class FDMatrices{
    public:
        FDMatrices();

        void BuildLaplacianMatrix(Tube &t, SpMat &L);
        void BuildLaplacianMatrix(Tube &t, const vector<size_t> &set_idx_from_subset_idx, const vector<ssize_t> &subset_idx_from_set_idx, SpMat &L);

        void BuildGradientMatrices(Tube &t, vector<SpMat> &Dc);
        void BuildGradientMatrices(Tube &t, const vector<size_t> &set_idx_from_subset_idx, const vector<ssize_t> &subset_idx_from_set_idx, vector<SpMat> &Dc);
};

} // namespace cpm
