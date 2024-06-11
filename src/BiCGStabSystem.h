#pragma once

#include <chrono>

#include "Scalar.h"
#include "BiCGStabVector.h"
#include <Eigen/Sparse>
// #include <eigen3/Eigen/Sparse>

using namespace std;
using namespace std::chrono;
using namespace Eigen;

namespace cpm{

class BiCGStabSystem
{
    using VectorBase                    = BiCGStabVector;

  public:
    BiCGStabSystem() {}
    ~BiCGStabSystem() {}

    virtual void Multiply(const VectorBase& x,VectorBase& result) const = 0;   
    virtual void Precondition(const VectorBase& x,VectorBase& result) const = 0;    
    
   

    Scalar InnerProduct(const VectorBase& v1,const VectorBase& v2) const
    {
        return v1.data.dot(v2.data);
    }

    Scalar SquaredNorm(const VectorBase& v) const
    {
        return v.data.squaredNorm();
    }
};

}
