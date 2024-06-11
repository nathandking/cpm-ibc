#pragma once

#ifdef POLYSCOPE
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/volume_mesh.h"
#endif

#include "Scalar.h"
#include "VectorMath.h"
#include "Tube.h"

#include <Eigen/Sparse>

namespace cpm {

class Visualization {
    public:
        Visualization();
        
    private:

};

std::vector<Scalar> PlaneVertex(const std::vector<Scalar> &n, const std::vector<Scalar> &center, Scalar radius, Scalar corner_scale1, Scalar corner_scale2);
void PlotPlane(const std::vector<Scalar> &n, const std::vector<Scalar> &center, Scalar radius, std::string planeName);
void PlotAABB(const std::vector<std::vector<Scalar>> &bounds, std::string AABBName);
void VisualizeTubeNeighbours(Tube &t, std::string plot_name);
void VisualizeTubeVolume(Tube &t, std::string plot_name);
void Visualize(const vector<vector<Scalar>> &points, std::string plot_name);
void VisualizeCurve(const vector<vector<Scalar>> &points, std::string plot_name);
void VisualizeTube(Tube &t, std::string plot_name);
vector<vector<Scalar>> VisualizeTubeSubset(Tube &t, const vector<size_t> &set_idx_from_subset_idx, std::string plot_name);
void VisualizeTubeSubset(Tube &t, const vector<size_t> &set_idx_from_subset_idx, const vector<vector<Scalar>> &subset_cpc, std::string plot_name);

void GetStencilIndicesAndEntriesPerRow(const SpMat &A, vector<vector<size_t>> &stencil_indices, vector<vector<Scalar>> &stencil_entries);
vector<size_t> GetStencilsInvolvingIBCs(const vector<vector<size_t>> &stencil_indices, size_t num_surface_tube_nodes);
void VisualizeStencil(string stencil_name, size_t row_index, const vector<vector<size_t>> &stencil_indices, const vector<vector<Scalar>> &stencil_entries, const vector<vector<Scalar>> &row_x, const vector<vector<Scalar>> &col_x, size_t num_surface_tube_nodes);

template <typename T>
void VisualizeTubeSubsetAddScalarQuantity(const vector<T> &scalar_quantity, std::string plot_name, std::string quantity_name);

template <typename T>
void VisualizeTubeSubsetAddScalarQuantity(const vector<T> &scalar_quantity, const vector<size_t> &set_idx_from_subset_idx, std::string plot_name, std::string quantity_name);

template <typename T>
void VisualizeAddScalarQuantity(const vector<T> &scalar_quantity, std::string plot_name, std::string quantity_name);

template <typename T>
void VisualizeTubeSubsetAddVectorQuantity(const vector<vector<T>> &vector_quantity, std::string plot_name, std::string quantity_name);

template <typename T>
void VisualizeTubeSubsetAddVectorQuantity(const vector<vector<T>> &vector_quantity, const vector<size_t> &set_idx_from_subset_idx, std::string plot_name, std::string quantity_name);

template <typename T>
void VisualizeAddVectorQuantity(const vector<vector<T>> &vector_quantity, std::string plot_name, std::string quantity_name);

template <typename T>
void VisualizeTubeAddScalarQuantity(const vector<T> &scalar_quantity, std::string plot_name, std::string quantity_name);

template <typename T>
void VisualizeTubeAddVectorQuantity(const vector<vector<T>> &vector_quantity, std::string plot_name, std::string quantity_name);

} // namespace cpm