#pragma once

#include <nova/Tools/Log/Log.h>
#include <nova/SPGrid/Core/SPGrid_Allocator.h>
#include <nova/SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <nova/Dynamics/Utilities/SPGrid_Flags.h>
using namespace Nova;
namespace SPGrid{
template<class Struct_type,class T,int d>
class SPGrid_Initializer
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Allocator_type        = SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    SPGrid_Initializer(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel,const T& val)
    {Run(allocator,blocks,channel,val);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel,const T& val) const
    {
        auto data=allocator.template Get_Array<Struct_type,T>(channel);
        auto flags=allocator.template Get_Const_Array<Struct_type, unsigned>(&Struct_type::flags);
        auto channel_cleaner=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(Flags_type)) {
                if(flags(offset) & Node_Active)
                {
                    data(offset) = val;
                }
            }
        };
        SPGrid_Computations::Run_Parallel_Blocks(blocks,channel_cleaner);
    }
};

template<class Struct_type,class T,int d>
class SPGrid_Masked_Initializer
{
    using Flags_type            = typename Struct_type::Flags_type;
    using Allocator_type        = SPGrid_Allocator<Struct_type,d>;
    using Flag_array_mask       = typename Allocator_type::template Array_mask<unsigned>;

  public:
    SPGrid_Masked_Initializer(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel,const unsigned mask,const T& val)
    {Run(allocator,blocks,channel,mask,val);}

    void Run(Allocator_type& allocator,const std::pair<const uint64_t*,unsigned>& blocks,T Struct_type::* channel,const unsigned mask,const T& val) const
    {
        auto flags=allocator.template Get_Const_Array<Struct_type,unsigned>(&Struct_type::flags);
        auto data=allocator.template Get_Array<Struct_type,T>(channel);

        auto masked_clear=[&](uint64_t offset)
        {
            for(int e=0;e<Flag_array_mask::elements_per_block;++e,offset+=sizeof(T)) {
                if(flags(offset)&mask) {
                    data(offset)=val;
                }
            }
        };

        SPGrid_Computations::Run_Parallel_Blocks(blocks,masked_clear);
    }
};
}
