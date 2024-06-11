#include "Interpolation.h"

namespace cpm {

// #define INTERPOLATION_TESTS


Interpolation::Interpolation(Tube &t)
{
    m_xq = t.cpx();
    m_Nq = t.nNodes(); 

    // find base point for each query point (lower left corner of interpolation stencil (hypercube))
    m_Ibpt = t.bpt();
    m_stencil_size = t.p() + 1;

    UniformGridWeights1D();

    BuildInterpolationWeights(t);
}


Interpolation::Interpolation(Tube &t, const vector<vector<Scalar>> &xquery)
{
    m_xq = xquery;
    m_Nq = m_xq.size(); 

    // find base point for each query point (lower left corner of interpolation stencil (hypercube))
    m_Ibpt = t.ComputeBasePoints(m_xq);
    m_stencil_size = t.p() + 1;

    UniformGridWeights1D();

    BuildInterpolationWeights(t);
}


// interpolation using all grid points in the tube
void Interpolation::BuildInterpolationMatrix(Tube &t, SpMat &E)
{
    vector<ssize_t> subset_idx_from_set_idx(t.nNodes());
    for(size_t i = 0; i < t.nNodes(); ++i)
    {
        subset_idx_from_set_idx[i] = i;
    }

    BuildInterpolationMatrix(t, subset_idx_from_set_idx, E);
}


// interpolation using a subset of grid points in the tube
void Interpolation::BuildInterpolationMatrix(Tube &t, const vector<ssize_t> &subset_idx_from_set_idx, SpMat &E)
{
    vector<T> coeffs; // list of non-zeros coefficients, (row, col, value) triplets to be used to construct sparse matrix afterward

    if(t.dim() == 2)
    {
        int xindex;
        int findex;
        for(int q = 0; q < m_Nq; ++q)
        {
            xindex = m_Ibpt[q];
            findex = m_Ibpt[q];
            for(int i = 0; i < m_stencil_size; ++i)
            {
                for(int j = 0; j < m_stencil_size; ++j)
                {
                    ssize_t subset_idx = subset_idx_from_set_idx[t.TubeNode(findex).TubeIndex()];
                    if(subset_idx >= 0)
                    {
                        coeffs.push_back(T(q, subset_idx, m_w[1][q][j] * m_w[0][q][i]));
                    }
                    findex = t.TubeNode(findex).neighbour(2)->TubeIndex(); // move one neighbour up
                }
                findex = t.TubeNode(xindex).neighbour(0)->TubeIndex(); // move one neighbour right, from the current x bottom base point
                xindex = findex;
            }
        }
    }
    else if(t.dim() == 3)
    {
        int xindex;
        int yindex;
        int zindex;
        for(int q = 0; q < m_Nq; ++q)
        {
            xindex = m_Ibpt[q];
            yindex = m_Ibpt[q];
            zindex = m_Ibpt[q];
            for(int i = 0; i < m_stencil_size; ++i)
            {
                for(int j = 0; j < m_stencil_size; ++j)
                {
                    for(int k = 0; k < m_stencil_size; ++k)
                    {
                        ssize_t subset_idx = subset_idx_from_set_idx[t.TubeNode(zindex).TubeIndex()];
                        if(subset_idx >= 0)
                        {
                            coeffs.push_back(T(q, subset_idx, m_w[2][q][k] * m_w[1][q][j] * m_w[0][q][i]));
                        }
                        zindex = t.TubeNode(zindex).neighbour(4)->TubeIndex(); // move +z one neighbour
                    }
                    zindex = t.TubeNode(yindex).neighbour(2)->TubeIndex(); // move one neighbour up
                    yindex = zindex;
                }
                zindex = t.TubeNode(xindex).neighbour(0)->TubeIndex(); // move one neighbour right, from the current x bottom base point
                yindex = zindex;
                xindex = zindex;
            }
        }
    }
    else
    {
        cout << "Dimension of embedding space not implemented in Interpolation::BuildInterpolationMatrix" << endl;
    }
    E.setFromTriplets(coeffs.begin(), coeffs.end());

#ifdef INTERPOLATION_TESTS
    TestRowSumOne(E);
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////
//    COMPUTING WEIGHTS AND INTERPOLATION MATRIX
/////////////////////////////////////////////////////////////////////////////////////////

void Interpolation::UniformGridWeights1D()
{
    m_uniform_grid_w.resize(m_stencil_size);
    switch(m_stencil_size)
    {
        case 1:
            m_uniform_grid_w[0] = 1;
            break;
        case 2:
            m_uniform_grid_w[0] = 1;
            m_uniform_grid_w[1] = -1;
            break;
        case 3:
            m_uniform_grid_w[0] = 1;
            m_uniform_grid_w[1] = -2;
            m_uniform_grid_w[2] = 1;
            break;
        case 4:
            m_uniform_grid_w[0] = 1;
            m_uniform_grid_w[1] = -3;
            m_uniform_grid_w[2] = 3;
            m_uniform_grid_w[3] = -1;
            break;
        case 5:
            m_uniform_grid_w[0] = 1;
            m_uniform_grid_w[1] = -4;
            m_uniform_grid_w[2] = 6;
            m_uniform_grid_w[3] = -4;
            m_uniform_grid_w[4] = 1;
            break;
    }
}


void Interpolation::BuildInterpolationWeights1D(vector<Scalar> &x, Scalar &xq_subset, vector<Scalar> &w)
{
    // check to see if any x[j] = m_xq
    bool isInterpPoint = false;
    int interpPoint;
    for(int j = 0; j < m_stencil_size; ++j)
    {
        if(abs(xq_subset - x[j]) < 1e-17)
        {
            isInterpPoint = true;
            interpPoint = j;
        }
    }

    if(isInterpPoint)
    {
        // set all weights to zero except the one for m_xq = x[j]
        for(int j = 0; j < m_stencil_size; ++j)
        {
            w[j] = 0.0;
        }
        w[interpPoint] = 1.0;
    }
    else
    {
        // add dependence on the query point
        Scalar wqSum = 0.0;
        for(int j = 0; j < m_stencil_size; ++j)
        {
            w[j] = m_uniform_grid_w[j] / (xq_subset - x[j]);
            wqSum += w[j];
        }

        for(int j = 0; j < m_stencil_size; ++j)
        {
            w[j] /= wqSum;
        }
    }
}


// template <typename B>
void Interpolation::BuildInterpolationWeights(Tube &t)
{
    m_w.resize(t.dim(), vector<vector<Scalar>>(m_Nq, vector<Scalar>(m_stencil_size)));

    // compute interpolation weights
    int index;
    vector<Scalar> x(m_stencil_size);
    for(int d = 0; d < t.dim(); ++d) // for each direction x,y(,z)
    {
        for(int q = 0; q < m_Nq; ++q)
        {
            index = m_Ibpt[q];
            x[0] = t.x()[m_Ibpt[q]][d];
            for(int i = 1; i < m_stencil_size; ++i)
            {
                x[i] = t.x()[t.TubeNode(index).neighbour(2 * d)->TubeIndex()][d];
                index = t.TubeNode(index).neighbour(2 * d)->TubeIndex(); // set index to your neighbour, to move one neighbour at a time
            }
            BuildInterpolationWeights1D(x, m_xq[q][d], m_w[d][q]);
        }
    }
}


//////////////////////////////////////////////////////////////////////////////
//    TESTING
//////////////////////////////////////////////////////////////////////////////

void Interpolation::TestRowSumOne(SpMat &E)
{
    VectorX sumVec = E * VectorX::Ones(E.cols());
    if(!sumVec.isOnes())
    {
        cout << "Interpolation matrix rows do not all sum to one" << endl;
    }
}

} // namespace cpm