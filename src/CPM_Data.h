#pragma once 
#include <stdint.h>

#include "Scalar.h"

namespace cpm{
template<class Scalar,class T_FLAGS=uint32_t>
struct CPM_Data
{
    typedef T_FLAGS Flags_type;
    T_FLAGS flags;          
    Scalar ch_tube_id_lower;
    Scalar ch_tube_id_higher;
};
}
