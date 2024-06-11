#include "Node.h"

namespace cpm {

Node::Node()
{
    for(int i = 0; i < 6; ++i) // TO DO: Make this work for nD, right now I only have an array enough for 3D. Also, this is wasteful in 2D since there are always 2 nullptr neighbours
    {
        m_neighbours[i] = nullptr;
    }
}

} // namespace cpm