#pragma once

#include "Scalar.h"
#include "Tube.h"
#include "Helpers.h"
#include "FDMatrices.h"
#include "Geometry.h"
#include "Solver.h"

#ifdef POLYSCOPE
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#endif

#include <Eigen/Sparse>

using namespace std;

namespace cpm {

class GeodesicDistanceNearest{
    public:
        GeodesicDistanceNearest(bool visualize = false);

        VectorX ComputeDistance(Tube &surface_tube, vector<vector<size_t>>& ibc_idx, vector<vector<size_t>> &inside_ibc_idx, vector<vector<Scalar>> &inside_ibc_idx_dist, Scalar dt);

    private:

        Geometry m_geom;

        Helpers m_h;
        bool m_visualize;

        SpMat m_L;
        SpMat m_E;
        vector<SpMat> m_Dc;
        vector<bool> m_identity_rows;

        VectorX Step1(Tube &surface_tube, vector<vector<size_t>> &ibc_idx, vector<vector<size_t>> &inside_ibc_idx, vector<vector<Scalar>> &inside_ibc_idx_dist, Scalar dt);
        vector<VectorX> Step2(VectorX &ufull, Tube &surface_tube);
        VectorX Step3(VectorX &div, Tube &surface_tube, vector<vector<size_t>> &ibc_idx);

        VectorX ComputeDivergence(vector<VectorX> &gradu, Tube &surface_tube);

        vector<vector<Scalar>> m_surface_normals;
        void ComputeSurfaceNormals(Tube &surface_tube);
        vector<vector<Scalar>> InterpolateUnorientedVectors(Tube &surface_tube, const vector<vector<Scalar>> &tube_vectors, const vector<vector<Scalar>> &xq, bool is_normalized);
};

} // namespace cpm
