#pragma once

#include <vector>
#include "Scalar.h"

namespace cpm {

class Grid {
    public:
        Grid();
        Grid(size_t embed_dim, Scalar dx, std::vector<Scalar> xstart, size_t N1d);

        void Clear();

        std::vector<Scalar> CoordinateFromTriplet(std::vector<size_t> &idx);
        std::vector<Scalar> CoordinateFromLinearIndex(size_t &index);
        std::vector<size_t> TripletFromCoordinate(std::vector<Scalar> &x);
        size_t NeighbourLinearIndex(std::vector<size_t> triplet, ssize_t shift, size_t shift_dim);
        size_t NeighbourLinearIndex(size_t current_index, ssize_t shift, size_t shift_dim);
        
        std::vector<size_t> TripletFromLinearIndex(size_t index) const;
        size_t LinearIndexFromTriplet(std::vector<size_t> idx);
        size_t LinearIndexFromCoordinate(std::vector<Scalar> &x);

        std::vector<size_t> NearestTriplet(const std::vector<Scalar> &xq);

        size_t dim() const { return m_embed_dim; };
        Scalar dx() const { return m_dx; };
        size_t N1d() const { return m_N1d; };
        size_t N() { return m_N; };
        size_t n() const { return m_n;}
        std::vector<Scalar> xstart() const { return m_xstart; };

        void SetDim(size_t dim){ m_embed_dim = dim; };
        void SetDx(Scalar dx){ m_dx = dx; };
        void SetXstart(std::vector<Scalar> xstart){ m_xstart = xstart; };
        
        // find number of grid points in the implicit rectangular grid
        void SetNumberGridPoints(size_t n){ m_n = n;
                                            m_N1d = pow(2.0, n) + 1;
                                            m_N = pow(m_N1d, m_embed_dim); };

    private:

        // grid spacing
        Scalar m_dx; // grid spacing, NOTE: code uses same spacing in x,y,z,...

        // specify the bounding box
        std::vector<Scalar> m_xstart; // (x, y, z,...) of negative (bottom-left in 2D) corner of bounding box for rectangular grid

        size_t m_embed_dim; // dimension of the embedding space

        size_t m_n;
        size_t m_N1d; // number of grid points in each dimension
        size_t m_N; // total number of grid points
};


} // namespace cpm