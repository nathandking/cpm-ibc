#include "Visualization.h"

namespace cpm {

Visualization::Visualization()
{

}


std::vector<Scalar> PlaneVertex(const std::vector<Scalar> &n, const std::vector<Scalar> &center, Scalar radius, Scalar corner_scale1, Scalar corner_scale2)
{
    Scalar d = -DotProduct(n, center);
    std::vector<Scalar> vert(3, 0.0);
    if(abs(n[0]) > 1e-1)
    {
        vert[1] = center[1] + corner_scale1 * radius;
        vert[2] = center[2] + corner_scale2 * radius;
        vert[0] = (-n[1] * vert[1] - n[2] * vert[2] - d) / n[0]; 
    }
    else if(abs(n[1]) > 1e-1)
    {
        vert[0] = center[0] + corner_scale1 * radius;
        vert[2] = center[2] + corner_scale2 * radius;
        vert[1] = (-n[0] * vert[0] - n[2] * vert[2] - d) / n[1]; 
    }
    else if(abs(n[2]) > 1e-1)
    {
        vert[0] = center[0] + corner_scale1 * radius;
        vert[1] = center[1] + corner_scale2 * radius;
        vert[2] = (-n[0] * vert[0] - n[1] * vert[1] - d) / n[2]; 
    }
    else
    {
        std::cout << "normal to plane is zero std::vector" << std::endl;
    }

    return vert;
}


void PlotPlane(const std::vector<Scalar> &n, const std::vector<Scalar> &center, Scalar radius, std::string planeName)
{           
#ifdef POLYSCOPE
    std::vector<std::vector<Scalar>> vert(4, std::vector<Scalar>(3));
    vert[0] =  PlaneVertex(n, center, radius, 1.0, 1.0);
    vert[1] =  PlaneVertex(n, center, radius, -1.0, 1.0);
    vert[2] =  PlaneVertex(n, center, radius, -1.0, -1.0);
    vert[3] =  PlaneVertex(n, center, radius, 1.0, -1.0);

    std::vector<std::vector<size_t>> face(1, std::vector<size_t>(4));
    face[0][0] = 0; face[0][1] = 1; face[0][2] = 2; face[0][3] = 3; 

    polyscope::registerSurfaceMesh(planeName, vert, face);
#endif
}


void PlotAABB(const std::vector<std::vector<Scalar>> &bounds, std::string AABBName)
{           
#ifdef POLYSCOPE
    if(bounds[0].size() == 2)
    {
        std::vector<std::vector<Scalar>> vert(4, std::vector<Scalar>(3));
        vert[0] = bounds[0];
        vert[2] = bounds[1];
        vert[1][0] = bounds[1][0];
        vert[1][1] = bounds[0][1];
        vert[3][0] = bounds[0][0];
        vert[3][1] = bounds[1][1];  

        std::vector<std::vector<size_t>> edges(4, std::vector<size_t>(4));
        edges[0][0] = 0; edges[0][1] = 1; 
        edges[1][0] = 1; edges[1][1] = 2;
        edges[2][0] = 2; edges[2][1] = 3;
        edges[3][0] = 3; edges[3][1] = 0;

        glm::vec3 black{0,0,0};
        polyscope::registerCurveNetwork2D(AABBName, vert, edges)->setColor(black);
    }
    else if(bounds[0].size() == 3)
    {
        std::vector<std::vector<Scalar>> vert(8, std::vector<Scalar>(3));
        vert[0] = bounds[0];
        
        vert[1][0] = bounds[1][0];
        vert[1][1] = bounds[0][1];
        vert[1][2] = bounds[0][2];
        
        vert[2][0] = bounds[1][0];
        vert[2][1] = bounds[1][1];
        vert[2][2] = bounds[0][2];  

        vert[3][0] = bounds[0][0];
        vert[3][1] = bounds[1][1];
        vert[3][2] = bounds[0][2];  

        vert[4][0] = bounds[0][0];
        vert[4][1] = bounds[0][1];
        vert[4][2] = bounds[1][2]; 

        vert[5][0] = bounds[1][0];
        vert[5][1] = bounds[0][1];
        vert[5][2] = bounds[1][2];

        vert[6] = bounds[1];

        vert[7][0] = bounds[0][0];
        vert[7][1] = bounds[1][1];
        vert[7][2] = bounds[1][2]; 


        std::vector<std::vector<size_t>> edges(12, std::vector<size_t>(4));
        edges[0][0] = 0; edges[0][1] = 1; 
        edges[1][0] = 1; edges[1][1] = 2;
        edges[2][0] = 2; edges[2][1] = 3;
        edges[3][0] = 3; edges[3][1] = 0;

        edges[4][0] = 4; edges[4][1] = 5; 
        edges[5][0] = 5; edges[5][1] = 6;
        edges[6][0] = 6; edges[6][1] = 7;
        edges[7][0] = 7; edges[7][1] = 4;

        edges[8][0] = 0; edges[8][1] = 4; 
        edges[9][0] = 1; edges[9][1] = 5;
        edges[10][0] = 2; edges[10][1] = 6;
        edges[11][0] = 3; edges[11][1] = 7;

        glm::vec3 black{0,0,0};
        polyscope::registerCurveNetwork(AABBName, vert, edges)->setColor(black);
    }
#endif
}


void VisualizeTubeNeighbours(Tube &t, std::string plot_name)
{
#ifdef POLYSCOPE
    vector<size_t> edge(2);
    vector<vector<size_t>> edges;
    for(size_t i = 0; i < t.nNodes(); ++i)
    {
        for(size_t nbr = 0; nbr < 2 * t.dim(); ++nbr)
        {
            if(t.TubeNode(i).neighbour(nbr) != nullptr)
            {
                edge[0] = i;
                edge[1] = t.TubeNode(i).neighbour(nbr)->TubeIndex();
                edges.push_back(edge);
            }
        }
    }
    
    if(t.dim() == 2)
    {
        polyscope::registerCurveNetwork2D(plot_name + " neighbours", t.x(), edges);
    }
    else if(t.dim() == 3)
    {
        polyscope::registerCurveNetwork(plot_name + " neighbours", t.x(), edges);
    }
    polyscope::getCurveNetwork(plot_name + " neighbours")->addNodeScalarQuantity("bdy", t.bdy());
#endif
}


void VisualizeTubeVolume(Tube &t, std::string plot_name)
{
#ifdef POLYSCOPE
    vector<size_t> hex(8);
    vector<vector<size_t>> hexes;
    for(size_t i = 0; i < t.nNodes(); ++i)
    {
        vector<size_t> plane_corner{i};
        if(t.TubeNode(i).neighbour(2/*+y*/) != nullptr)
        {
            plane_corner.push_back(t.TubeNode(i).neighbour(2/*+y*/)->TubeIndex());
            
            bool all_nbrs_exist = true;
            vector<size_t> nbr_plane_order{0/*+x*/, 5/*-z*/, 1/*-x*/};  
            for(size_t c = 0; c < 2; ++c)
            {
                hex[4 * c] = plane_corner[c];

                size_t nbr_tube_idx = plane_corner[c];
                for(size_t v = 0; v < nbr_plane_order.size(); ++v)
                {
                    if(t.TubeNode(nbr_tube_idx).neighbour(nbr_plane_order[v]) != nullptr)
                    {
                        nbr_tube_idx = t.TubeNode(nbr_tube_idx).neighbour(nbr_plane_order[v])->TubeIndex();
                        hex[4 * c + v + 1] = nbr_tube_idx;
                    }
                    else
                    {
                        all_nbrs_exist = false;
                    }
                }
            }

            if(all_nbrs_exist)
            {
                hexes.push_back(hex);
            }
        }
    }
    
    polyscope::registerHexMesh(plot_name + " Volume", t.x(), hexes)->addVertexScalarQuantity("Distance", t.dist());
#endif
}


void Visualize(const vector<vector<Scalar>> &points, std::string plot_name)
{
#ifdef POLYSCOPE
    if(points[0].size() == 2)
    {
        polyscope::registerPointCloud2D(plot_name, points);
    }
    else if(points[0].size() == 3)
    {
        polyscope::registerPointCloud(plot_name, points);
    }
#endif
}


void VisualizeCurve(const vector<vector<Scalar>> &points, std::string plot_name)
{
#ifdef POLYSCOPE
    if(points[0].size() == 2)
    {
        polyscope::registerCurveNetworkLine2D(plot_name, points);
    }
    else if(points[0].size() == 3)
    {
        polyscope::registerCurveNetworkLine(plot_name, points);
    }
#endif
}


void VisualizeTube(Tube &t, std::string plot_name)
{
    Visualize(t.x(), plot_name + " x");
    Visualize(t.cpx(), plot_name + " cpx");

    // visualize vectors from cp(x) to x
    vector<vector<Scalar>> cpx_to_x(t.nNodes());
    for(size_t i = 0; i < t.nNodes(); ++i)
    {
        cpx_to_x[i] = t.x()[i] - t.cpx()[i];
    }
    VisualizeAddVectorQuantity(cpx_to_x, plot_name + " cpx", "x - cp(x)");
}


vector<vector<Scalar>> VisualizeTubeSubset(Tube &t, const vector<size_t> &set_idx_from_subset_idx, std::string plot_name)
{
    size_t subset_size = set_idx_from_subset_idx.size();
    vector<vector<Scalar>> subset_x;
    if(subset_size > 0)
    {
        subset_x.resize(subset_size, vector<Scalar>(t.dim()));
        vector<vector<Scalar>> subset_cps(subset_size, vector<Scalar>(t.dim()));
        for(size_t i = 0; i < subset_size; ++i)
        {
            subset_x[i] = t.x()[set_idx_from_subset_idx[i]];
            subset_cps[i] = t.cpx()[set_idx_from_subset_idx[i]];
        }

        Visualize(subset_x, plot_name + " x");
        Visualize(subset_cps, plot_name + " cps");

        // visualize vectors from cp(x) to x
        vector<vector<Scalar>> cps_to_x(subset_size);
        for(size_t i = 0; i < subset_size; ++i)
        {
            cps_to_x[i] = subset_x[i] - subset_cps[i];
        }
        VisualizeAddVectorQuantity(cps_to_x, plot_name + " cps", "x - cps(x)");
    }

    return subset_x;
}


void VisualizeTubeSubset(Tube &t, const vector<size_t> &set_idx_from_subset_idx, const vector<vector<Scalar>> &subset_cpc, std::string plot_name)
{
    vector<vector<Scalar>> subset_x = VisualizeTubeSubset(t, set_idx_from_subset_idx, plot_name);

    size_t subset_size = set_idx_from_subset_idx.size();
    Visualize(subset_cpc, plot_name + " cpc");

    // visualize vectors from cp(x) to x
    vector<vector<Scalar>> cpc_to_x(subset_size);
    for(size_t i = 0; i < subset_size; ++i)
    {
        cpc_to_x[i] = subset_x[i] - subset_cpc[i];
    }
    VisualizeAddVectorQuantity(cpc_to_x, plot_name + " cpc", "x - cpc(x)");
}


void GetStencilIndicesAndEntriesPerRow(const SpMat &A, vector<vector<size_t>> &stencil_indices, vector<vector<Scalar>> &stencil_entries)
{
    stencil_indices.resize(A.rows());
    stencil_entries.resize(A.rows()); 

    // put the nonzeros for each row of the matrix in vectors for easy access
    for(size_t k=0; k < A.outerSize(); ++k)
    {
        for(SpMat::InnerIterator it(A,k); it; ++it)
        {
            stencil_indices[it.row()].push_back(it.col());
            stencil_entries[it.row()].push_back(it.value());
        }
    }
}


vector<size_t> GetStencilsInvolvingIBCs(const vector<vector<size_t>> &stencil_indices, size_t num_surface_tube_nodes)
{
    vector<size_t> involves_ibc_subset;
    for(size_t row = 0; row < stencil_indices.size(); ++row)
    {
        bool involves_ibc = false;
        for(size_t col = 0; col < stencil_indices[row].size(); ++col)
        {
            if(stencil_indices[row][col] >= num_surface_tube_nodes)
            {
                involves_ibc = true;
            }
        }
        
        if(involves_ibc)
        {
            involves_ibc_subset.push_back(row);
        }
    }

    return involves_ibc_subset;
}


// plot the stencil grid points as well as the corresponding stencil director for that row (the FD center or query point for interpolation)
void VisualizeStencil(string stencil_name, size_t row_index, const vector<vector<size_t>> &stencil_indices, const vector<vector<Scalar>> &stencil_entries, const vector<vector<Scalar>> &row_x, const vector<vector<Scalar>> &col_x, size_t num_surface_tube_nodes)
{
#ifdef POLYSCOPE
    vector<vector<Scalar>> stencil_x(stencil_indices[row_index].size());
    vector<Scalar> stencil_x_entry(stencil_indices[row_index].size());
    vector<vector<Scalar>> director_x(1); // there is only one stencil director
    director_x[0] = row_x[row_index];

    vector<Scalar> is_ibc_index(stencil_indices[row_index].size(), 0.0);

    for(size_t col = 0; col < stencil_indices[row_index].size(); ++col)
    {
        stencil_x[col] = col_x[stencil_indices[row_index][col]];
        stencil_x_entry[col] = stencil_entries[row_index][col];
        if(stencil_indices[row_index][col] >= num_surface_tube_nodes)
        {
            is_ibc_index[col] = 1.0;
        }
    }

    if(stencil_x[0].size() == 2)
    {
        polyscope::registerPointCloud2D(stencil_name, stencil_x)->addScalarQuantity("Entry", stencil_x_entry);
        polyscope::registerPointCloud2D(stencil_name + " Director", director_x);
    }
    else if(stencil_x[0].size() == 3)
    {
        polyscope::registerPointCloud(stencil_name, stencil_x)->addScalarQuantity("Entry", stencil_x_entry);
        polyscope::registerPointCloud(stencil_name + " Director", director_x);
    }
    polyscope::getPointCloud(stencil_name)->addScalarQuantity("DOF type", is_ibc_index)->setEnabled(true);
#endif
}


template <typename T>
void VisualizeTubeSubsetAddScalarQuantity(const vector<T> &scalar_quantity, std::string plot_name, std::string quantity_name)
{
#ifdef POLYSCOPE
    VisualizeAddScalarQuantity(scalar_quantity, plot_name + " x", quantity_name);
    VisualizeAddScalarQuantity(scalar_quantity, plot_name + " cps", quantity_name);
    if(polyscope::getPointCloud(plot_name + " cpc") != nullptr)
    {
        VisualizeAddScalarQuantity(scalar_quantity, plot_name + " cpc", quantity_name);
    }
#endif
}


template <typename T>
void VisualizeTubeSubsetAddScalarQuantity(const vector<T> &scalar_quantity, const vector<size_t> &set_idx_from_subset_idx, std::string plot_name, std::string quantity_name)
{
#ifdef POLYSCOPE
    size_t subset_size = set_idx_from_subset_idx.size();
    vector<T> subset_scalar_quantity(subset_size);
    for(size_t i = 0; i < subset_size; ++i)
    {
        subset_scalar_quantity[i] = scalar_quantity[set_idx_from_subset_idx[i]];
    }

    VisualizeAddScalarQuantity(subset_scalar_quantity, plot_name + " x", quantity_name);
    VisualizeAddScalarQuantity(subset_scalar_quantity, plot_name + " cps", quantity_name);
    if(polyscope::getPointCloud(plot_name + " cpc") != nullptr)
    {
        VisualizeAddScalarQuantity(subset_scalar_quantity, plot_name + " cpc", quantity_name);
    }
#endif
}


template <typename T>
void VisualizeAddScalarQuantity(const vector<T> &scalar_quantity, std::string plot_name, std::string quantity_name)
{
#ifdef POLYSCOPE
    polyscope::getPointCloud(plot_name)->addScalarQuantity(quantity_name, scalar_quantity);
#endif
}


template <typename T>
void VisualizeTubeSubsetAddVectorQuantity(const vector<vector<T>> &vector_quantity, std::string plot_name, std::string quantity_name)
{
#ifdef POLYSCOPE
    VisualizeAddVectorQuantity(vector_quantity, plot_name + " x", quantity_name);
    VisualizeAddVectorQuantity(vector_quantity, plot_name + " cps", quantity_name);
    if(polyscope::getPointCloud(plot_name + " cpc") != nullptr)
    {
        VisualizeAddVectorQuantity(vector_quantity, plot_name + " cpc", quantity_name);
    }
#endif
}


template <typename T>
void VisualizeTubeSubsetAddVectorQuantity(const vector<vector<T>> &vector_quantity, const vector<size_t> &set_idx_from_subset_idx, std::string plot_name, std::string quantity_name)
{
#ifdef POLYSCOPE
    size_t subset_size = set_idx_from_subset_idx.size();
    vector<vector<T>> subset_vector_quantity(subset_size);
    for(size_t i = 0; i < subset_size; ++i)
    {
        subset_vector_quantity[i] = vector_quantity[set_idx_from_subset_idx[i]];
    }

    VisualizeAddVectorQuantity(subset_vector_quantity, plot_name + " x", quantity_name);
    VisualizeAddVectorQuantity(subset_vector_quantity, plot_name + " cps", quantity_name);
    if(polyscope::getPointCloud(plot_name + " cpc") != nullptr)
    {
        VisualizeAddVectorQuantity(subset_vector_quantity, plot_name + " cpc", quantity_name);
    }
#endif
}


template <typename T>
void VisualizeAddVectorQuantity(const vector<vector<T>> &vector_quantity, std::string plot_name, std::string quantity_name)
{
#ifdef POLYSCOPE
    Scalar max_vec_norm = 0;
    for(size_t i = 0; i < vector_quantity.size(); ++i)
    {
        if(Norm(vector_quantity[i]) > max_vec_norm)
        {
            max_vec_norm = Norm(vector_quantity[i]);
        }
    }
    
    if(vector_quantity[0].size() == 2)
    {
        polyscope::getPointCloud(plot_name)->addVectorQuantity2D(quantity_name, vector_quantity)->setVectorLengthScale(max_vec_norm, false);
    }
    else if(vector_quantity[0].size() == 3)
    {
        polyscope::getPointCloud(plot_name)->addVectorQuantity(quantity_name, vector_quantity)->setVectorLengthScale(max_vec_norm, false);
    }
#endif
}


template <typename T>
void VisualizeTubeAddScalarQuantity(const vector<T> &scalar_quantity, std::string plot_name, std::string quantity_name)
{
    VisualizeAddScalarQuantity(scalar_quantity, plot_name + " x", quantity_name);
    VisualizeAddScalarQuantity(scalar_quantity, plot_name + " cpx", quantity_name);
}


template <typename T>
void VisualizeTubeAddVectorQuantity(const vector<vector<T>> &vector_quantity, std::string plot_name, std::string quantity_name)
{
    VisualizeAddVectorQuantity(vector_quantity, plot_name + " x", quantity_name);
    VisualizeAddVectorQuantity(vector_quantity, plot_name + " cpx", quantity_name);
}


template void VisualizeTubeSubsetAddScalarQuantity<ssize_t>(const vector<ssize_t> &scalar_quantity, std::string plot_name, std::string quantity_name);
template void VisualizeTubeSubsetAddScalarQuantity<size_t>(const vector<size_t> &scalar_quantity, std::string plot_name, std::string quantity_name);
template void VisualizeTubeSubsetAddScalarQuantity<Scalar>(const vector<Scalar> &scalar_quantity, std::string plot_name, std::string quantity_name);

template void VisualizeTubeSubsetAddScalarQuantity<ssize_t>(const vector<ssize_t> &scalar_quantity, const vector<size_t> &set_idx_from_subset_idx, std::string plot_name, std::string quantity_name);
template void VisualizeTubeSubsetAddScalarQuantity<size_t>(const vector<size_t> &scalar_quantity, const vector<size_t> &set_idx_from_subset_idx, std::string plot_name, std::string quantity_name);
template void VisualizeTubeSubsetAddScalarQuantity<Scalar>(const vector<Scalar> &scalar_quantity, const vector<size_t> &set_idx_from_subset_idx, std::string plot_name, std::string quantity_name);

template void VisualizeAddScalarQuantity<ssize_t>(const vector<ssize_t> &scalar_quantity, std::string plot_name, std::string quantity_name);
template void VisualizeAddScalarQuantity<size_t>(const vector<size_t> &scalar_quantity, std::string plot_name, std::string quantity_name);
template void VisualizeAddScalarQuantity<Scalar>(const vector<Scalar> &scalar_quantity, std::string plot_name, std::string quantity_name);

template void VisualizeTubeSubsetAddVectorQuantity<ssize_t>(const vector<vector<ssize_t>> &vector_quantity, std::string plot_name, std::string quantity_name);
template void VisualizeTubeSubsetAddVectorQuantity<size_t>(const vector<vector<size_t>> &vector_quantity, std::string plot_name, std::string quantity_name);
template void VisualizeTubeSubsetAddVectorQuantity<Scalar>(const vector<vector<Scalar>> &vector_quantity, std::string plot_name, std::string quantity_name);

template void VisualizeTubeSubsetAddVectorQuantity<ssize_t>(const vector<vector<ssize_t>> &vector_quantity, const vector<size_t> &set_idx_from_subset_idx, std::string plot_name, std::string quantity_name);
template void VisualizeTubeSubsetAddVectorQuantity<size_t>(const vector<vector<size_t>> &vector_quantity, const vector<size_t> &set_idx_from_subset_idx, std::string plot_name, std::string quantity_name);
template void VisualizeTubeSubsetAddVectorQuantity<Scalar>(const vector<vector<Scalar>> &vector_quantity, const vector<size_t> &set_idx_from_subset_idx, std::string plot_name, std::string quantity_name);

template void VisualizeAddVectorQuantity<ssize_t>(const vector<vector<ssize_t>> &vector_quantity, std::string plot_name, std::string quantity_name);
template void VisualizeAddVectorQuantity<size_t>(const vector<vector<size_t>> &vector_quantity, std::string plot_name, std::string quantity_name);
template void VisualizeAddVectorQuantity<Scalar>(const vector<vector<Scalar>> &vector_quantity, std::string plot_name, std::string quantity_name);

template void VisualizeTubeAddScalarQuantity<ssize_t>(const vector<ssize_t> &scalar_quantity, std::string plot_name, std::string quantity_name);
template void VisualizeTubeAddScalarQuantity<size_t>(const vector<size_t> &scalar_quantity, std::string plot_name, std::string quantity_name);
template void VisualizeTubeAddScalarQuantity<Scalar>(const vector<Scalar> &scalar_quantity, std::string plot_name, std::string quantity_name);

template void VisualizeTubeAddVectorQuantity<ssize_t>(const vector<vector<ssize_t>> &vector_quantity, std::string plot_name, std::string quantity_name);
template void VisualizeTubeAddVectorQuantity<size_t>(const vector<vector<size_t>> &vector_quantity, std::string plot_name, std::string quantity_name);
template void VisualizeTubeAddVectorQuantity<Scalar>(const vector<vector<Scalar>> &vector_quantity, std::string plot_name, std::string quantity_name);


} // namespace cpm