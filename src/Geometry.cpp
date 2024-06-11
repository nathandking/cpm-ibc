#include "Geometry.h"

namespace cpm {


Geometry::Geometry()
{

}


template <typename M, typename V>
void Geometry::EigenDecompositionOfJcp(Tube &t, vector<M> &eigenvectors, vector<V> &eigenvalues)
{
    vector<ssize_t> subset_idx_from_set_idx(t.nNodes());
    vector<size_t> set_idx_from_subset_idx(t.nNodes());
    for(size_t i = 0; i < t.nNodes(); ++i)
    {
        subset_idx_from_set_idx[i] = i;
        set_idx_from_subset_idx[i] = i;
    }

    EigenDecompositionOfJcp(t, t.cpx(), set_idx_from_subset_idx, subset_idx_from_set_idx, eigenvectors, eigenvalues);
}


template <typename M, typename V>
void Geometry::EigenDecompositionOfJcp(Tube &t, const vector<vector<Scalar>> &subset_cpx, const vector<size_t> &set_idx_from_subset_idx, const vector<ssize_t> &subset_idx_from_set_idx, vector<M> &eigenvectors, vector<V> &eigenvalues)
{
    size_t subset_size = set_idx_from_subset_idx.size();

    // build gradient matrices
    vector<SpMat> Dc(t.dim());
    for(size_t d = 0; d < t.dim(); ++d)
    {
        Dc[d].resize(subset_size, subset_size);
    }

    FDMatrices fdmat;
    fdmat.BuildGradientMatrices(t, set_idx_from_subset_idx, subset_idx_from_set_idx, Dc); 

    // build closest point extension matrix
    SpMat E(subset_size, subset_size);
    Interpolation interp(t, subset_cpx); 
    interp.BuildInterpolationMatrix(t, subset_idx_from_set_idx, E);

    // store vectors for each dimension of the closest points
    vector<VectorX> cp(t.dim(), VectorX(subset_size)); 
    for(size_t i = 0; i < subset_size; ++i)
    {
        for(size_t d = 0; d < t.dim(); ++d)
        {
            cp[d][i] = subset_cpx[i][d];
        }
    }

    // compute the Jacobian of the closest point function
    vector<vector<VectorX>> Jcp(t.dim(), vector<VectorX>(t.dim(), VectorX(subset_size)));
    for(size_t row = 0; row < t.dim(); ++row)
    {
        for(size_t col = 0; col < t.dim(); ++col)
        {
            Jcp[row][col] = Dc[col] * cp[row];
            Jcp[row][col] = E * Jcp[row][col]; // perform closest point extension of Jacobian
        }
    }

    EigenDecomposition(t, Jcp, eigenvectors, eigenvalues);
}


template <typename M, typename V>
void Geometry::EigenDecomposition(Tube &t, vector<vector<VectorX>> &J, vector<M> &eigenvectors, vector<V> &eigenvalues)
{
    size_t size = J[0][0].size();

    Eigen::SelfAdjointEigenSolver<M> es;
    eigenvectors.resize(size);
    eigenvalues.resize(size);
    for(size_t i = 0; i < size; ++i)
    {
        M A;
        for(size_t row = 0; row < t.dim(); ++row)
        {
            for(size_t col = 0; col < t.dim(); ++col)
            {    
                A(row, col) = J[row][col][i];
            }
        }

        A = 0.5 * (A + A.transpose()); // matrix shoud have Jcp = Jcp^T, average to make symmetric if numerical errors

        es.compute(A);

        eigenvectors[i] = es.eigenvectors();
        eigenvalues[i] = es.eigenvalues();
    }
}


vector<vector<Scalar>> Geometry::UnorientedSurfaceNormals(Tube &t, bool is_normalized, Scalar zero_tol)
{
    vector<vector<Scalar>> normals(t.nNodes(), vector<Scalar>(t.dim()));
    for(size_t i = 0; i < t.nNodes(); ++i)
    {
        normals[i] = t.x()[i] - t.cpx()[i];
    }
    
    // if x is too close to the surface, interpolate the normal from neighbouring grid points
    for(size_t i = 0; i < t.nNodes(); ++i)
    {
        if(Norm(normals[i]) < zero_tol)
        {
            InterpolateUnorientedVectorFromNeighbours(t, normals, i, zero_tol);
        }
    }

    if(is_normalized)
    {
        for(size_t i = 0; i < t.nNodes(); ++i)
        {
            Normalize(normals[i], zero_tol);
        }
    }

    return normals;
}


void Geometry::InterpolateUnorientedVectorFromNeighbours(Tube &t, vector<vector<Scalar>> &vectors, size_t interp_index, Scalar zero_tol)
{
    vector<ssize_t> subset_idx_from_set_idx(t.nNodes());
    vector<size_t> set_idx_from_subset_idx(t.nNodes());
    for(size_t i = 0; i < t.nNodes(); ++i)
    {
        subset_idx_from_set_idx[i] = i;
        set_idx_from_subset_idx[i] = i;
    }

    InterpolateUnorientedVectorFromNeighbours(t, vectors, set_idx_from_subset_idx, subset_idx_from_set_idx, interp_index, zero_tol);
}


void Geometry::InterpolateUnorientedVectorFromNeighbours(Tube &t, vector<vector<Scalar>> &subset_vectors, const vector<size_t> &set_idx_from_subset_idx, const vector<ssize_t> &subset_idx_from_set_idx, size_t subset_interp_index, Scalar zero_tol)
{
    vector<Scalar> accumulated_vector(t.dim(), 0.0);
    vector<Scalar> first_direction(t.dim());
    size_t num_vec_accumulated = 0;
    for(size_t nbr = 0; nbr < 2 * t.dim(); ++nbr)
    {
        size_t index = set_idx_from_subset_idx[subset_interp_index];
        if(t.TubeNode(index).neighbour(nbr) != nullptr)
        {
            size_t nbr_idx = t.TubeNode(index).neighbour(nbr)->TubeIndex();
            ssize_t subset_nbr_idx = subset_idx_from_set_idx[nbr_idx];
            if(subset_nbr_idx >= 0)
            {
                if(Norm(subset_vectors[subset_nbr_idx]) >= zero_tol)
                {
                    if(num_vec_accumulated == 0)
                    {
                        first_direction = subset_vectors[subset_nbr_idx];
                    }

                    if(DotProduct(first_direction, subset_vectors[subset_nbr_idx]) > 0)
                    {
                        accumulated_vector += subset_vectors[subset_nbr_idx];
                    }
                    else
                    {
                        accumulated_vector -= subset_vectors[subset_nbr_idx];
                    }

                    ++num_vec_accumulated;
                }
            }
        }
    }

    if(num_vec_accumulated > 0)
    {
        subset_vectors[subset_interp_index] = accumulated_vector / num_vec_accumulated;
    } 
    else // else, I should probably use the next neighbours, but for now it just doesn't change the vector
    {
        cout << "We have a problem with the cp_diff on S_perp in IBCSubset" << endl;
    }
}


vector<vector<Scalar>> Geometry::InterpolateUnorientedVectors(Tube &t, const vector<vector<Scalar>> &vectors, const vector<vector<Scalar>> &xq, bool is_normalized)
{
    vector<ssize_t> subset_idx_from_set_idx(t.nNodes());
    vector<size_t> set_idx_from_subset_idx(t.nNodes());
    for(size_t i = 0; i < t.nNodes(); ++i)
    {
        subset_idx_from_set_idx[i] = i;
        set_idx_from_subset_idx[i] = i;
    }

    return InterpolateUnorientedVectors(t, set_idx_from_subset_idx, subset_idx_from_set_idx, vectors, xq, is_normalized);
}


vector<vector<Scalar>> Geometry::InterpolateUnorientedVectors(Tube &t, const vector<size_t> &set_idx_from_subset_idx, const vector<ssize_t> &subset_idx_from_set_idx, const vector<vector<Scalar>> &subset_vectors,  const vector<vector<Scalar>> &xq, bool is_normalized)
{
    // build interpolation matrix for interpolating from the tube points onto query points xq
    size_t subset_size = set_idx_from_subset_idx.size();
    SpMat E(xq.size(), subset_size);
    Interpolation interp(t, xq);
    interp.BuildInterpolationMatrix(t, subset_idx_from_set_idx, E);

    // now apply interpolation one weight at a time, checking if there is a flip in the sign of the vector, since they are unoriented
    vector<vector<Scalar>> vq(xq.size(), vector<Scalar>(subset_vectors[0].size(), 0.0));
    ApplyInterpolationToUnorientedVectors(E, subset_vectors, vq);

    if(is_normalized)
    {
        // normalize the accumulated vector
        for(size_t i = 0; i < vq.size(); ++i)
        {
            Normalize(vq[i], 0.0);
        }
    }
    else
    {
        // divide by the number of accumlated vectors, i.e., number of points in the interpolation stencil
        size_t num_interp_points = pow(t.p() + 1, t.dim());
        for(size_t i = 0; i < vq.size(); ++i)
        {
            vq[i] /= num_interp_points;
        }
    }

    return vq;
}


void Geometry::ApplyInterpolationToUnorientedVectors(const SpMat &E, const vector<vector<Scalar>> &v, vector<vector<Scalar>> &vq)
{
    vector<bool> first_row_visit(E.rows(), true);
    vector<vector<Scalar>> first_direction(E.rows(), vector<Scalar>(v[0].size(), 0.0));
    for (size_t k=0; k < E.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(E,k); it; ++it)
        {   
            if(first_row_visit[it.row()])
            {
                first_direction[it.row()] = v[it.col()];
                first_row_visit[it.row()] = false;
            }

            Scalar dot_product = DotProduct(first_direction[it.row()], v[it.col()]);
            assert(dot_product != 0.0);

            if(dot_product > 0.0)
            {
                vq[it.row()] += it.valueRef() * v[it.col()];
            }
            else // flip vector direction before using it in the interpolation
            {
                vq[it.row()] -= it.valueRef() * v[it.col()];
            }
        }
    }
}


template void Geometry::EigenDecompositionOfJcp<Matrix2, Vector2>(Tube &t, vector<Matrix2> &eigenvectors, vector<Vector2> &eigenvalues);
template void Geometry::EigenDecompositionOfJcp<Matrix3, Vector3>(Tube &t, vector<Matrix3> &eigenvectors, vector<Vector3> &eigenvalues);

template void Geometry::EigenDecompositionOfJcp<Matrix2, Vector2>(Tube &t, const vector<vector<Scalar>> &subset_cpx, const vector<size_t> &set_idx_from_subset_idx, const vector<ssize_t> &subset_idx_from_set_idx, vector<Matrix2> &eigenvectors, vector<Vector2> &eigenvalues);
template void Geometry::EigenDecompositionOfJcp<Matrix3, Vector3>(Tube &t, const vector<vector<Scalar>> &subset_cpx, const vector<size_t> &set_idx_from_subset_idx, const vector<ssize_t> &subset_idx_from_set_idx, vector<Matrix3> &eigenvectors, vector<Vector3> &eigenvalues);

template void Geometry::EigenDecomposition<Matrix2, Vector2>(Tube &t, vector<vector<VectorX>> &J, vector<Matrix2> &eigenvectors, vector<Vector2> &eigenvalues);
template void Geometry::EigenDecomposition<Matrix3, Vector3>(Tube &t, vector<vector<VectorX>> &J, vector<Matrix3> &eigenvectors, vector<Vector3> &eigenvalues);


} // namespace cpm