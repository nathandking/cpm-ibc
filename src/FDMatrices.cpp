#include "FDMatrices.h"

namespace cpm {

FDMatrices::FDMatrices()
{

}


void FDMatrices::BuildLaplacianMatrix(Tube &t, SpMat &L)
{
    vector<ssize_t> subset_idx_from_set_idx(t.nNodes());
    vector<size_t> set_idx_from_subset_idx(t.nNodes());
    for(size_t i = 0; i < t.nNodes(); ++i)
    {
        subset_idx_from_set_idx[i] = i;
        set_idx_from_subset_idx[i] = i;
    }

    BuildLaplacianMatrix(t, set_idx_from_subset_idx, subset_idx_from_set_idx, L);
}


void FDMatrices::BuildLaplacianMatrix(Tube &t, const vector<size_t> &set_idx_from_subset_idx, const vector<ssize_t> &subset_idx_from_set_idx, SpMat &L)
{
    vector<T> coeffs; // list of non-zeros coefficients, (row, col, value) triplets to be used to construct sparse matrix afterward
    for(size_t row = 0; row < set_idx_from_subset_idx.size(); ++row)
    {
        coeffs.push_back(T(row, row, -2.0 * (Scalar)t.dim())); // diagonal entry

        size_t set_idx = set_idx_from_subset_idx[row];
        for(size_t nbr = 0; nbr < 2 * t.dim(); ++nbr)
        {
            if(t.TubeNode(set_idx).neighbour(nbr) != nullptr)
            {
                size_t nbr_idx = t.TubeNode(set_idx).neighbour(nbr)->TubeIndex();
                ssize_t nbr_subset_idx = subset_idx_from_set_idx[nbr_idx];
                if(nbr_subset_idx >= 0)
                {
                    coeffs.push_back(T(row, nbr_subset_idx, 1.0)); // neighour entries
                }
            }
        }
    }
    
    L.setFromTriplets(coeffs.begin(), coeffs.end());
    L /= (t.dx() * t.dx());
}


void FDMatrices::BuildGradientMatrices(Tube &t, vector<SpMat> &Dc)
{
    vector<ssize_t> subset_idx_from_set_idx(t.nNodes());
    vector<size_t> set_idx_from_subset_idx(t.nNodes());
    for(size_t i = 0; i < t.nNodes(); ++i)
    {
        subset_idx_from_set_idx[i] = i;
        set_idx_from_subset_idx[i] = i;
    }

    BuildGradientMatrices(t, set_idx_from_subset_idx, subset_idx_from_set_idx, Dc);
}


void FDMatrices::BuildGradientMatrices(Tube &t, const vector<size_t> &set_idx_from_subset_idx, const vector<ssize_t> &subset_idx_from_set_idx, vector<SpMat> &Dc)
{
    vector<vector<T>> coeffs(t.dim()); // list of non-zeros coefficients for each dimension in a vector, (row, col, value) triplets to be used to construct sparse matrix afterward

    for(size_t row = 0; row < set_idx_from_subset_idx.size(); ++row)
    {
        size_t set_idx = set_idx_from_subset_idx[row];
        for(size_t nbr = 0; nbr < 2 * t.dim(); ++nbr)
        {
            if(t.TubeNode(set_idx).neighbour(nbr) != nullptr)
            {
                size_t nbr_idx = t.TubeNode(set_idx).neighbour(nbr)->TubeIndex();
                ssize_t nbr_subset_idx = subset_idx_from_set_idx[nbr_idx];
                if(nbr_subset_idx >= 0)
                {
                    coeffs[nbr/2].push_back(T(row, nbr_subset_idx, pow(-1, nbr) * 0.5)); // neighour entry
                }
            }
        }
    }
    
    for(int d = 0; d < t.dim(); ++d)
    {
        Dc[d].setFromTriplets(coeffs[d].begin(), coeffs[d].end());
        Dc[d] /= t.dx();
    }
}


} // namespace cpm