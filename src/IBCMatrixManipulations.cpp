#include "IBCMatrixManipulations.h"

namespace cpm {


IBCMatrixManipulations::IBCMatrixManipulations()
{

}


SpMat IBCMatrixManipulations::AddIBCExtension(SpMat &A, Tube &surface_tube, IBCSubsets &ibc, const vector<vector<Scalar>> &director_x)
{
    // A.makeCompressed();

    vector<SpMat> A_tags(ibc.NumIBCs());
    for(size_t c = 0; c < ibc.NumIBCs(); ++c)
    {
        A_tags[c] = TagStencils(A, surface_tube, ibc, ibc.TaggingSubset()[c].directions(), ibc.TaggingSubset()[c].SubsetIndexFromSetIndex(), c);
    }

    vector<T> coeffs;
    MoveToIBCColumns(A, A_tags, ibc, coeffs);

    AddExtensionIBCRows(surface_tube, ibc, coeffs);

    SpMat Afull(A.rows() + ibc.TotalIBCDOFs() + ibc.Num2ndOrderDirichletRows(), A.cols() + ibc.TotalIBCDOFs() + ibc.Num2ndOrderDirichletRows());
    Afull.setFromTriplets(coeffs.begin(), coeffs.end());

    // clear memory of original A after done with it
    // A.resize(0,0);
    // A.data().squeeze();

    return Afull;
}


SpMat IBCMatrixManipulations::AddIBCFiniteDifference(SpMat &A, Tube &surface_tube, IBCSubsets &ibc, const vector<vector<Scalar>> &director_x)
{
    vector<SpMat> A_tags(ibc.NumIBCs());
    for(size_t c = 0; c < ibc.NumIBCs(); ++c)
    {
        A_tags[c] = TagStencils(A, surface_tube, ibc, ibc.TaggingSubset()[c].directions(), ibc.TaggingSubset()[c].SubsetIndexFromSetIndex(), c);
    }

    vector<Eigen::Triplet<Scalar>> coeffs;
    MoveToIBCColumns(A, A_tags, ibc, coeffs);

    AddFiniteDifferenceIBCRows(A, A_tags, coeffs, ibc);

    for(size_t c = 0; c < ibc.NumIBCs(); ++c)
    {
        for(size_t i = 0; i < ibc.DOFSubset()[c].size(); ++i)
        {
            const size_t bc_type_idx = ibc.meta(c).is_oriented ? ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
            if(ibc.meta(c).boundary_type[bc_type_idx] == 0 && ibc.meta(c).boundary_order == 2) // add in identity rows
            {
                coeffs.emplace_back(ibc.IdentityRowsStart() + ibc.DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + ibc.DirichletIBCColumnStartIndex(c), ibc.IdentityRowsStart() + ibc.DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + ibc.DirichletIBCColumnStartIndex(c), 1.0);
            }
        }
    }

    SpMat Afull(A.rows() + ibc.TotalIBCDOFs() + ibc.Num2ndOrderDirichletRows(), A.cols() + ibc.TotalIBCDOFs() + ibc.Num2ndOrderDirichletRows());
    Afull.setFromTriplets(coeffs.begin(), coeffs.end());

    // clear memory of original A after done with it
    // A.resize(0,0);
    // A.data().squeeze();

    return Afull;
}


SpMat IBCMatrixManipulations::AddIBCPlotInterpolation(SpMat &A, Tube &surface_tube, IBCSubsets &ibc, const vector<vector<Scalar>> &director_x)
{
    vector<vector<size_t>> director_set_index_from_subset_index(ibc.NumIBCs());
    vector<vector<size_t>> director_which_side(ibc.NumIBCs());
    vector<vector<bool>> director_is_on_ibc(ibc.NumIBCs());

    return AddIBCPlotInterpolation(A, surface_tube, ibc, director_x, director_set_index_from_subset_index, director_which_side, director_is_on_ibc);
}


SpMat IBCMatrixManipulations::AddIBCPlotInterpolation(SpMat &A, Tube &surface_tube, IBCSubsets &ibc, const vector<vector<Scalar>> &director_x, vector<vector<size_t>> &director_set_index_from_subset_index, vector<vector<size_t>> &director_which_side, vector<vector<bool>> &director_is_on_ibc)
{
    vector<vector<vector<Scalar>>> director_directions(ibc.NumIBCs());
    vector<vector<size_t>> director_bdy(ibc.NumIBCs());
    vector<vector<ssize_t>> director_subset_index_from_set_index(ibc.NumIBCs());
    director_set_index_from_subset_index.resize(ibc.NumIBCs());
    director_which_side.resize(ibc.NumIBCs());
    director_is_on_ibc.resize(ibc.NumIBCs());
    vector<SpMat> A_tags(ibc.NumIBCs());
    for(size_t c = 0; c < ibc.NumIBCs(); ++c)
    {
        PlotPointsCPDiff(surface_tube, ibc, director_directions[c], director_bdy[c], director_subset_index_from_set_index[c], director_set_index_from_subset_index[c], director_which_side[c], director_is_on_ibc[c], c);
        A_tags[c] = TagStencils(A, surface_tube, ibc, director_directions[c], director_subset_index_from_set_index[c], c);
    }

    vector<Eigen::Triplet<Scalar>> coeffs;
    MoveToIBCColumns(A, A_tags, ibc, coeffs);

    SpMat Afull(A.rows(), A.cols() + ibc.TotalIBCDOFs() + ibc.Num2ndOrderDirichletRows());
    Afull.setFromTriplets(coeffs.begin(), coeffs.end());

    return Afull;
}


////////////////////////////////////////////////////////////////////////////////////
// Private Functions
////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for matrix manipulation based on normal crossings of grid points compared to stencil director
///////////////////////////////////////////////////////////////////////////////////////////////////////////


void IBCMatrixManipulations::AddExtensionIBCRows(Tube &surface_tube, IBCSubsets &ibc, vector<T> &coeffs)
{   
    for(size_t c = 0; c < ibc.NumIBCs(); ++c)
    {
        SpMat E_ibc = ibc.DOFSubset()[c].cpExtensionMatrix();
        vector<T> ibc_coeffs;

        vector<ssize_t> subset_indices(ibc.DOFSubset()[c].size()); // the rows of E_ibc is just one per cp_c, so you do not use SubsetIndexFromSetIndex here
        for(size_t i = 0; i < ibc.DOFSubset()[c].size(); ++i)
        {
            subset_indices[i] = i;
        }
        SpMat E_tags = TagStencils(E_ibc, surface_tube, ibc, ibc.DOFSubset()[c].directions(), subset_indices, c, -1.0);
        
        // move to ibc columns
        for (size_t k=0; k < E_ibc.outerSize(); ++k)
        {
            for (SpMat::InnerIterator it(E_ibc,k); it; ++it)
            {
                const size_t bc_type_idx = ibc.meta(c).is_oriented ? ibc.DOFSubset()[c].which_side()[it.row()] - 1 : 0;
                const Scalar boundary_type_scale = (ibc.meta(c).boundary_type[bc_type_idx] == 0) ? -1.0 : 1.0;
                if(ibc.meta(c).boundary_type[bc_type_idx] == 1 /* Neumann */ || (ibc.meta(c).boundary_type[bc_type_idx] == 0 /* Dirichlet */ && ibc.meta(c).boundary_order == 2))
                {
                    if(E_tags.coeffRef(it.row(), it.col()) == 2.0 && ibc.DOFSubset()[c].SubsetIndexFromSetIndex()[it.col()] >= 0)
                    {
                        coeffs.emplace_back(it.row() + ibc.IBCColumnStartIndex(c), ibc.DOFSubset()[c].SubsetIndexFromSetIndex()[it.col()] + ibc.IBCColumnStartIndex(c), boundary_type_scale * it.value()); 
                    }
                    else
                    {
                        coeffs.emplace_back(it.row() + ibc.IBCColumnStartIndex(c), it.col(), boundary_type_scale * it.value());
                    }
                }
            }
        }

        if(ibc.meta(c).boundary_order == 2) // only for 2nd order Dirichlet
        {
            for(size_t i = 0; i < ibc.DOFSubset()[c].size(); ++i)
            {
                const size_t bc_type_idx = ibc.meta(c).is_oriented ? ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
                if(ibc.meta(c).boundary_type[bc_type_idx] == 0 /* Dirichlet */)
                {
                    coeffs.emplace_back(i + ibc.IBCColumnStartIndex(c), ibc.IdentityRowsStart() + ibc.DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + ibc.DirichletIBCColumnStartIndex(c), 2.0);
                }
            }

            // add in identity rows for Dirichlet interior BCs
            for(size_t i = 0; i < ibc.DOFSubset()[c].size(); ++i)
            {
                const size_t bc_type_idx = ibc.meta(c).is_oriented ? ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
                if(ibc.meta(c).boundary_type[bc_type_idx] == 0)
                {
                    coeffs.emplace_back(ibc.IdentityRowsStart() + ibc.DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + ibc.DirichletIBCColumnStartIndex(c), ibc.IdentityRowsStart() + ibc.DOFSubset()[c].DirichletSubsetIndexFromSetIndex()[i] + ibc.DirichletIBCColumnStartIndex(c), 1.0);
                }
            }
        }
        else
        {
            for(size_t i = 0; i < ibc.DOFSubset()[c].size(); ++i)
            {
                const size_t bc_type_idx = ibc.meta(c).is_oriented ? ibc.DOFSubset()[c].which_side()[i] - 1 : 0;
                if(ibc.meta(c).boundary_type[bc_type_idx] == 0)
                {
                    coeffs.emplace_back(ibc.IdentityRowsStart() + i + ibc.IBCColumnStartIndex(c), ibc.IdentityRowsStart() + i + ibc.IBCColumnStartIndex(c), 1.0);
                }
            }
        }
    }
}


void IBCMatrixManipulations::MoveToIBCColumns(SpMat &A, vector<SpMat> &A_tags, IBCSubsets &ibc, vector<T> &coeffs)
{   
    for (size_t k=0; k < A.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(A,k); it; ++it)
        {        
            size_t is_ibc_point = 0;
            for(size_t c = 0; c < ibc.NumIBCs(); ++c)
            {
                if(A_tags[c].coeffRef(it.row(), it.col()) == 2.0 && ibc.DOFSubset()[c].SubsetIndexFromSetIndex()[it.col()] >= 0)
                {
                    ++is_ibc_point;
                    coeffs.emplace_back(it.row(), ibc.DOFSubset()[c].SubsetIndexFromSetIndex()[it.col()] + ibc.IBCColumnStartIndex(c), it.value()); 
                }
            }

            if(is_ibc_point > 1)
            {
                // cout << "IBCMatrixManipulations::MoveToIBCColumns - The ibcs are too close together" << endl;
            }
            
            if(is_ibc_point == 0)
            {
                coeffs.emplace_back(it.row(), it.col(), it.value());
            }
        }
    }
}


void IBCMatrixManipulations::AddFiniteDifferenceIBCRows(SpMat &A, vector<SpMat> &A_tags, vector<Eigen::Triplet<Scalar>> &coeffs, IBCSubsets &ibc)
{    
    for (size_t k=0; k < A.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(A,k); it; ++it)
        {
            size_t is_ibc_point = 0;
            for(size_t c = 0; c < ibc.NumIBCs(); ++c)
            {
                if(A_tags[c].coeffRef(it.row(), it.col()) == 2.0 && ibc.DOFSubset()[c].SubsetIndexFromSetIndex()[it.row()] >= 0)
                {
                    ++is_ibc_point;
                        coeffs.emplace_back(ibc.DOFSubset()[c].SubsetIndexFromSetIndex()[it.row()] + ibc.IBCColumnStartIndex(c), it.col(), it.value()); 
                }
                else
                {
                    if(ibc.DOFSubset()[c].SubsetIndexFromSetIndex()[it.row()] >= 0) // if your row has a corresponding ibc grid point
                    {
                        if(ibc.DOFSubset()[c].SubsetIndexFromSetIndex()[it.col()] >= 0) // if your column has a corresponding ibc grid point
                        {
                            coeffs.emplace_back(ibc.DOFSubset()[c].SubsetIndexFromSetIndex()[it.row()] + ibc.IBCColumnStartIndex(c), ibc.DOFSubset()[c].SubsetIndexFromSetIndex()[it.col()] + ibc.IBCColumnStartIndex(c), it.value());
                        }
                    }
                }
            }

            if(is_ibc_point > 1)
            {
                // cout << "IBCMatrixManipulations::AddFiniteDifferenceIBCRows - The ibcs are too close together" << endl;
            }
        }
    }
}


void IBCMatrixManipulations::PlotPointsCPDiff(Tube &surface_tube, IBCSubsets &ibc, vector<vector<Scalar>> &subset_directions, vector<size_t> &subset_bdy, vector<ssize_t> &subset_index_from_set_index, vector<size_t> &set_index_from_subset_index, vector<size_t> &which_side, vector<bool> &is_on_ibc, size_t ibc_index)
{
    vector<Scalar> cp_c(surface_tube.dim());
    vector<vector<Scalar>> subset_cp_c;
    vector<vector<Scalar>> subset_cp_diff;
    vector<Scalar> cp_diff(surface_tube.dim());
    Scalar dummy;
    size_t bdy;
    subset_index_from_set_index.assign(surface_tube.surface().xp().size(), -1);

    size_t subset_count = 0;
    for(size_t i = 0; i < surface_tube.surface().xp().size(); ++i)
    {
        ibc.ClosestPoint(ibc_index, surface_tube.surface().xp()[i], cp_c, dummy, bdy);

        cp_diff = surface_tube.surface().xp()[i] - cp_c;
        
        if(Norm(cp_diff) < 2.0 * surface_tube.TubeRadius())
        {
            subset_cp_c.push_back(cp_c);
            subset_cp_diff.push_back(cp_diff);
            subset_bdy.push_back(bdy);
            subset_index_from_set_index[i] = subset_count;
            set_index_from_subset_index.push_back(i);
            ++subset_count;
        }
    }

    // check if any cp_diff have length less than tolerance
    is_on_ibc.assign(subset_count, false);
    for(size_t i = 0; i < subset_count; ++i)
    {
        if(Norm(subset_cp_diff[i]) < ibc.meta(ibc_index).on_ibc_tolerance)
        {
            is_on_ibc[i] = true;
            size_t nearest_idx = surface_tube.NearestTubeIndex(surface_tube.surface().xp()[set_index_from_subset_index[i]]);
    
            if(ibc.TaggingSubset()[ibc_index].SubsetIndexFromSetIndex()[nearest_idx] >= 0) // should be included in the subset because it is on Sperp, except if bdy > 0
            {
                subset_cp_diff[i] = ibc.TaggingSubset()[ibc_index].cp_diff()[ibc.TaggingSubset()[ibc_index].SubsetIndexFromSetIndex()[nearest_idx]];
            }
        }
    }

    subset_directions.resize(subset_count, vector<Scalar>(surface_tube.dim()));
    if(ibc.meta(ibc_index).use_cp_diff_directions)
    {
        for(size_t i = 0; i < subset_count; ++i)
        {
            if(subset_bdy[i] > 0)
            {
                vector<Scalar> projection = 2.0 * subset_cp_c[i] - surface_tube.surface().xp()[set_index_from_subset_index[i]];

                vector<Scalar> cpBar_s(surface_tube.dim()); // not quite cpBar_s, it is actually cpBar_c below, but here it is cp_s of the projection point onto the end of the ibc
                vector<Scalar> cpBar_c(surface_tube.dim());
                
                surface_tube.ClosestPoint(projection, cpBar_s);
                ibc.ClosestPoint(ibc_index, projection, cpBar_c);
                subset_cp_diff[i] = cpBar_c - cpBar_s; // -ve cp_diff of projection point
            }
            subset_directions[i] = subset_cp_diff[i];
        }
    }
    else
    {
        // interpolate the binormals onto cpc
        Geometry geom;
        vector<vector<Scalar>> binormals = geom.InterpolateUnorientedVectors(surface_tube, ibc.TaggingSubset()[ibc_index].SetIndexFromSubsetIndex(), ibc.TaggingSubset()[ibc_index].SubsetIndexFromSetIndex(), ibc.TaggingSubset()[ibc_index].binormals(), subset_cp_c);

        for(size_t i = 0; i < subset_count; ++i)
        {
            subset_directions[i] = DotProduct(binormals[i], subset_cp_diff[i]) * binormals[i];
        }
    }

    vector<Scalar> ibc_parameter(subset_count);
    vector<size_t> sorted_idx;
    if(ibc.meta(ibc_index).is_oriented && subset_count > 0)
    {
        for(size_t i = 0; i < subset_count; ++i)
        {
            ibc_parameter[i] = ibc.IBCParameter(ibc_index, subset_cp_c[i]);
        }

        sorted_idx = SortedIndices(ibc_parameter);
    }

    which_side.assign(subset_count, 0);
    if(ibc.meta(ibc_index).is_oriented && subset_count > 0)
    {
        which_side[sorted_idx[0]] = 1;
        for(size_t i = 1; i < subset_count; ++i)
        {
            if(DotProduct(subset_directions[sorted_idx[i]], subset_directions[sorted_idx[i-1]]) < 0)
            {
                assert(which_side[sorted_idx[i-1]] != 0);
                if(which_side[sorted_idx[i-1]] == 2)
                {
                    which_side[sorted_idx[i]] = 1;
                }
                else
                {
                    which_side[sorted_idx[i]] = 2;
                }
            }
            else
            {
                which_side[sorted_idx[i]] = which_side[sorted_idx[i-1]];
            }
        }

        if(ibc.meta(ibc_index).point_on_side1.size() == surface_tube.dim())
        {
            // make sure the which_side tag is such that points near point_on_side1 have tag of which_side = 1
            size_t nearest_subset_idx = 0;
            Scalar nearest_dist = numeric_limits<Scalar>::max(); // find nearest subset index
            for(size_t i = 0; i < subset_count; ++i)
            {
                Scalar dist = Distance(ibc.meta(ibc_index).point_on_side1, surface_tube.surface().xp()[set_index_from_subset_index[i]]); 
                
                vector<Scalar> cpx_c(surface_tube.dim());
                
                ibc.ClosestPoint(ibc_index, ibc.meta(ibc_index).point_on_side1, cpx_c);
                vector<Scalar> cp_diff_side1 = ibc.meta(ibc_index).point_on_side1 - cpx_c;

                if(subset_count == 1 && DotProduct(cp_diff_side1, subset_cp_diff[i]) < 0)
                {
                    which_side[0] = 2;
                }

                if(dist < nearest_dist && Norm(subset_cp_diff[i]) >= 2 * surface_tube.dx() && DotProduct(cp_diff_side1, subset_cp_diff[i]) > 0)
                {
                    nearest_dist = dist;
                    nearest_subset_idx = i; 
                }
            }

            if(which_side[nearest_subset_idx] == 2 && subset_count != 1)
            {
                for(size_t i = 0; i < subset_count; ++i)
                {
                    if(which_side[i] == 1)
                    {
                        which_side[i] = 2;
                    }
                    else if(which_side[i] == 2)
                    {
                        which_side[i] = 1;
                    }
                }
            }
        }
    }

//     vector<vector<Scalar>> xp(subset_count, vector<Scalar>(surface_tube.dim()));
//     for(size_t i = 0; i < subset_count; ++i)
//     {
//         xp[i] = surface_tube.surface().xp()[set_index_from_subset_index[i]];
//     }
// #ifdef POLYSCOPE
//     polyscope::registerPointCloud2D("Plotting IBC Subset", xp)->addScalarQuantity("which side", which_side);
// #endif
}


// A - matrix with stencils defined for each row based on non-zero entries of the columns in each row
// stencil_director - a vector of vectors defining the direction of the stencil for testing which side a grid point is on of the ibc normal.
//                     Each row corresponds to the stencil direction for rows in A whose stencil center (finite-difference) or query point (interpolation) are within the ibc tube radius
// director_ibc_idx - a mapping from row index in A to the index in stencil_director for stencil centers or query points that are within the ibc tube radius
SpMat IBCMatrixManipulations::TagStencils(SpMat &A, Tube& surface_tube, IBCSubsets &ibc, const vector<vector<Scalar>> &director_directions, const vector<ssize_t>& director_subset_index_from_set_index, size_t ibc_index, Scalar flip_scale)
{
    vector<Eigen::Triplet<Scalar>> coeffs; 
    for (size_t k=0; k < A.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(A,k); it; ++it)
        {
            ssize_t director_idx = director_subset_index_from_set_index[it.row()];
            ssize_t stencil_point_idx = ibc.TaggingSubset()[ibc_index].SubsetIndexFromSetIndex()[it.col()];
           
            if(director_idx >= 0 && stencil_point_idx >= 0)
            {
                if(ibc.TaggingSubset()[ibc_index].bdy()[stencil_point_idx] == 0) // grid point not in the boundary set at the ends of a open curve (half ball)
                {
                    Scalar dot_product = flip_scale * DotProduct(director_directions[director_idx], ibc.TaggingSubset()[ibc_index].directions()[stencil_point_idx]);
                    if(dot_product < 0.0)
                    {
                        coeffs.emplace_back(it.row(), it.col(), 2.0); // move column to ibc column
                    }
                    else
                    {
                        coeffs.emplace_back(it.row(), it.col(), 1.0);
                    }
                }
                else
                {
                    coeffs.emplace_back(it.row(), it.col(), 1.0);
                }
            }
            else
            {
                coeffs.emplace_back(it.row(), it.col(), 1.0);
            }
        }
    }

    SpMat A_tags(A.rows(), A.cols());
    A_tags.setFromTriplets(coeffs.begin(), coeffs.end());

    return A_tags;
}


} // namespace cpm