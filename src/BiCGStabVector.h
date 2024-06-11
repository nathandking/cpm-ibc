#pragma once

#include <cassert>

#include "Scalar.h"

// #include <eigen3/Eigen/Dense>
#include <Eigen/Dense>

using namespace Eigen;
namespace cpm{
class BiCGStabVector
{
    using T_Data = Eigen::Matrix<Scalar,Dynamic,1>;
  public:
    T_Data& data;

    BiCGStabVector(T_Data& data_input): data(data_input) 
    {
    }
    ~BiCGStabVector() {}

    const BiCGStabVector& operator=(const BiCGStabVector& bv)
    {Copy((Scalar)1,bv);return *this;}
    
    BiCGStabVector& operator+=(const BiCGStabVector& bv)
    {
        assert(bv.data.size()==data.size());
    #pragma omp parallel for
        for(int i=0;i<data.size();++i) data(i)+=bv.data(i);
        return *this;
    }

    BiCGStabVector& operator-=(const BiCGStabVector& bv)
    {
        assert(bv.data.size()==data.size());
    #pragma omp parallel for
        for(int i=0;i<data.size();++i) data(i)-=bv.data(i);
        return *this;
    }

    BiCGStabVector& operator*=(const Scalar a)
    {
    #pragma omp parallel for
        for(int i=0;i<data.size();++i) data(i)*=a;
        return *this;                                                
    }

    void Copy(const Scalar c,const BiCGStabVector& bv)
    {
        assert(bv.data.size()==data.size()); 
    #pragma omp parallel for
        for(int i=0;i<data.size();++i) data(i)=c*bv.data(i);
    }

    void Copy(const Scalar c1,const BiCGStabVector& bv1,const BiCGStabVector& bv2)
    {
        assert(bv1.data.size()==data.size()); assert(bv1.data.size()==data.size());
    #pragma omp parallel for
        for(int i=0;i<data.size();++i) data(i)=c1*bv1.data(i)+bv2.data(i);
    }

    Scalar Squared_Norm() const
    {
        return data.squaredNorm();
    }

    void Clear()
    {
    #pragma omp parallel for
        for(int i=0;i<data.size();++i) data(i)=(Scalar)0.;
    }
};
}
