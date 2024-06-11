#pragma once

#include "stddef.h"

namespace cpm {

class Node {
    public:
        Node();

        size_t recIndex(){ return m_rec_index; }; 
        size_t TubeIndex(){ return m_tube_index; }; 
        Node* neighbour(size_t nbr_index){ return m_neighbours[nbr_index]; };

        void setRecIndex(size_t rec_index){ m_rec_index = rec_index; }; 
        void setTubeIndex(size_t tube_index){ m_tube_index = tube_index; }; 
        void setNeighbour(size_t nbr_index, Node* neighbour){ m_neighbours[nbr_index] = neighbour; }; 

    private:
        size_t m_rec_index; // index in original rectangular grid
        size_t m_tube_index; // index in the computational tube, instead of index in original rectangular grid. Set to -1 if grid point is not in tube

        Node *m_neighbours[6]; // neighbours of current node, order is right(+x), left(-x), up(+y), down(-y), front(+z), back(-z).
};

} // namespace cpm