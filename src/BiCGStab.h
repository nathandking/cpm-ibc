#pragma once

#include <iostream>

#include "Scalar.h"
#include "Defines.h"
#include "BiCGStabVector.h"
#include "BiCGStabSystem.h"


using namespace std;

namespace cpm{

class BiCGStab
{
    using VectorBase           = BiCGStabVector;
    using SystemBase           = BiCGStabSystem;

  public:
    BiCGStab() 
    {
    #ifdef TRACK_WHERE
        std::cout << "BiCGStab" << std::endl;
    #endif
    }

    ~BiCGStab() {}
//######################################################################
    bool Solve(const SystemBase& system,
            VectorBase& x,const VectorBase& b,
            VectorBase& r,VectorBase& r0,
            VectorBase& v,VectorBase& p,
            VectorBase& s,VectorBase& t,
            VectorBase& y,VectorBase& z,
            const Scalar tolerance,const int min_iterations,const int max_iterations, int& iters)
    {
    #ifdef TRACK_WHERE
        cout << "Solve" << endl;
    #endif
        r0.Clear();
        r=b;

        system.Multiply(x,r0); 

        r-=r0;
        r0=r;

        Scalar r0_sqnorm=system.SquaredNorm(r0);
        Scalar rhs_sqnorm=system.SquaredNorm(b);
        if(rhs_sqnorm==(Scalar)0.) {x.Clear(); return true;}
        Scalar rho=(Scalar)1.; Scalar alpha=(Scalar)1.; Scalar w=(Scalar)1.;
        v.Clear(); p.Clear();
        
        Scalar tol2=tolerance*tolerance*rhs_sqnorm;
        Scalar eps2=Eigen::NumTraits<Scalar>::epsilon()*Eigen::NumTraits<Scalar>::epsilon();
        iters=0; int restarts=0;
  
        while(system.SquaredNorm(r)>tol2 && iters<max_iterations){
            Scalar rho_old=rho;
            rho=system.InnerProduct(r0,r);
            if(std::abs(rho)<eps2*r0_sqnorm){
                // std::cout<<"restarting"<<std::endl;
                r0.Clear();
                r=b;
                system.Multiply(x, r0); 
                r-=r0;
                r0=r;
                rho=r0_sqnorm=system.SquaredNorm(r);
                if(restarts++==0) iters=0;
            }
            
            Scalar beta=(rho/rho_old)*(alpha/w);
            p.Copy(-w,v,p); p*=beta; 
            p+=r;

            system.Precondition(p, y);
            system.Multiply(y,v);
            alpha=rho/system.InnerProduct(r0,v);
            
            s.Copy(-alpha,v,r);

            system.Precondition(s, z);
            system.Multiply(z,t);
  
            Scalar val=system.SquaredNorm(t);            
            if(val>(Scalar)0.)
            {
                w=system.InnerProduct(t,s)/val;
            }
            else 
            {
                w=(Scalar)0.;
            }
        
            x.Copy(alpha,y,x);
            x.Copy(w,z,x);
            r.Copy(-w,t,s);
#ifdef SHOW_PROGRESS
            cout<<iters<<": "<<system.SquaredNorm(r)<<" vs "<<tol2<< endl;
#endif
            ++iters;
        }
#ifdef SHOW_ITERATION
        std::cout << "#iterations " << iters <<std::endl;
#endif
        return false;
    }
};
}
