#include "Grid.h"

namespace cpm {

Grid::Grid()
{

}


Grid::Grid(size_t embed_dim, Scalar dx, std::vector<Scalar> xstart, size_t N1d)
{
        m_embed_dim = embed_dim;
        m_dx = dx;
        m_xstart = xstart;
        
        m_N1d = N1d;
        m_N = pow(m_N1d, m_embed_dim);
}


void Grid::Clear()
{
    m_xstart.clear(); // (x, y, z,...) of negative (bottom-left in 2D) corner of bounding box for rectangular grid
}


std::vector<Scalar> Grid::CoordinateFromTriplet(std::vector<size_t> &idx)
{
    assert(idx.size() == m_embed_dim);
    
    std::vector<Scalar> x(m_embed_dim);
    for(size_t d = 0; d < m_embed_dim; ++d)
    {
        x[d] = m_xstart[d] + idx[d] * m_dx;
    }

    return x;
}


std::vector<Scalar> Grid::CoordinateFromLinearIndex(size_t &index)
{
    std::vector<size_t> triplet = TripletFromLinearIndex(index);
    return CoordinateFromTriplet(triplet);
}


std::vector<size_t> Grid::TripletFromCoordinate(std::vector<Scalar> &x)
{
    assert(x.size() == m_embed_dim);

    std::vector<size_t> idx(m_embed_dim);
    for(size_t d = 0; d < m_embed_dim; ++d)
    {
        assert(x[d] >= m_xstart[d]);
        idx[d] =  (size_t) ( (x[d] - m_xstart[d]) / m_dx );
        assert(idx[d] < m_N1d);
    }

    return idx;
}


size_t Grid::LinearIndexFromCoordinate(std::vector<Scalar> &x)
{
    return LinearIndexFromTriplet(TripletFromCoordinate(x));
}


size_t Grid::NeighbourLinearIndex(size_t current_index, ssize_t shift, size_t shift_dim)
{
    size_t nbr_index;

    std::vector<size_t> triplet = TripletFromLinearIndex(current_index);
    triplet[shift_dim] += shift;

    assert(triplet[shift_dim] >= 0 && triplet[shift_dim] < m_N1d); // if false the padding around the surface needs to be increased

    nbr_index = LinearIndexFromTriplet(triplet);

    return nbr_index;
}


size_t Grid::NeighbourLinearIndex(std::vector<size_t> triplet, ssize_t shift, size_t shift_dim)
{
    triplet[shift_dim] += shift;

    assert(triplet[shift_dim] >= 0 && triplet[shift_dim] < m_N1d); // if false the padding around the surface needs to be increased

    return LinearIndexFromTriplet(triplet);
}


std::vector<size_t> Grid::TripletFromLinearIndex(size_t index) const
{
    assert(index >= 0 && index < m_N);

    std::vector<size_t> idx(m_embed_dim);

    if(m_embed_dim == 2)
    {
        idx[1] = index % m_N1d;
        idx[0] = (index - idx[1]) / m_N1d;
    }
    else if(m_embed_dim == 3)
    {
        idx[2] = (index % (m_N1d * m_N1d)) % m_N1d;
        idx[1] = ((index - idx[2]) / m_N1d) % m_N1d;
        idx[0] = (index - idx[2] - idx[1] * m_N1d) / (m_N1d * m_N1d);
    }

    for(size_t d = 0; d < m_embed_dim; ++d)
    {
        assert(idx[d] >= 0 && idx[d] < m_N1d);
    }

    return idx;
}


size_t Grid::LinearIndexFromTriplet(std::vector<size_t> idx)
{
    for(size_t d = 0; d < m_embed_dim; ++d)
    {
        assert(idx[d] >= 0 && idx[d] < m_N1d);
    }

    size_t index;

    if(m_embed_dim == 2)
    {
        index = idx[0] * m_N1d + idx[1]; // indicies run up the y direction and then move one x and up again
    }
    else if(m_embed_dim == 3)
    {
        index = idx[0] * m_N1d * m_N1d + idx[1] * m_N1d + idx[2];
    }

    assert(index >= 0 && index < m_N);
    
    return index;
}


std::vector<size_t> Grid::NearestTriplet(const std::vector<Scalar> &xq)
{
    std::vector<size_t> triplet(m_embed_dim);
    for(size_t d = 0; d < dim(); ++d)
    {
        triplet[d] = round((xq[d] - xstart()[d]) / dx());
    }

    return triplet;
}

} // namespace cpm