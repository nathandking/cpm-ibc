#pragma once

#include <chrono>

#include "Defines.h"
#include "Scalar.h"
#include "BiCGStabVector.h"
#include "BiCGStabSystem.h"

#include <Eigen/Sparse>
// #include <eigen3/Eigen/Sparse>

using namespace std;
using namespace std::chrono;
using namespace Eigen;

namespace cpm{

class GeneralizedBiCGStabSystem: BiCGStabSystem
{
    using VectorBase                    = BiCGStabVector;
  public:
    Scalar m, n;
    Preconditioner precond;
    SpMat L;
    SpMat E;
    int precond_iters;
    int nD;
    size_t dofs;
    vector<bool> identity_rows;
    vector<Scalar> diag_entries;
    bool computed;
    
    GeneralizedBiCGStabSystem() {}
    ~GeneralizedBiCGStabSystem() {}

    void SetComponents(Scalar m_input, Scalar n_input, SpMat L_input, SpMat E_input, Preconditioner precond_input, 
                        const vector<bool>& identity_rows_input, int nD_input=0)
    {
#ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::SetComponents" << std::endl;
#endif
        m = m_input;
        n = n_input;
        L = L_input;
        E = E_input;
        precond = precond_input;
        identity_rows = identity_rows_input;
        nD = nD_input;
        diag_entries.resize(L.rows(), (Scalar)0.);
        computed = false;
        dofs = L.rows() - nD;
    }

    void SetPrecond(int precond_iters_input)
    {
    #ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::SetPrecond" << std::endl;
    #endif
        precond_iters = precond_iters_input;
    }

    void SetBoundaryConditions(VectorBase& v) const {}
    void ProjectNullspace(VectorBase& v) const {}
    void LdiagTimesX(const VectorX& x_array, VectorX& result_array, Scalar m_input, Scalar n_input) const 
    {
    #ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::LdiagTimesX" << std::endl;
    #endif
        result_array.setZero();
#ifdef OPENMP
    #pragma omp parallel for
#endif
        for (size_t k=0; k < dofs; ++k)
        {
            for (SpMat::InnerIterator it(L, k); it; ++it)
            {
                size_t row = it.row();
                size_t col = it.col();
                if(row == col)
                {
                    Scalar val = it.value();
                    result_array(row) += val * x_array(col);
                }
            }
        }
        for(int i=0;i<nD;++i)
        {
            result_array(dofs + i) = (m_input - n_input) * x_array(dofs + i);
        }
    }

    void LnodiagTimesX(const VectorX& x_array,VectorX& result_array) const
    {
    #ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::LnodiagTimesX" << std::endl;
    #endif
        result_array.setZero();
#ifdef OPENMP
    #pragma omp parallel for
#endif
        for (size_t k=0; k < L.outerSize(); ++k)
        {
            for (SpMat::InnerIterator it(L, k); it; ++it)
            {
                size_t row = it.row();
                size_t col = it.col();
                if(row != col)
                {
                    Scalar val = it.value();
                    result_array(row) += val * x_array(col);
                }
            }
        }
    }

    void Compute()
    {
    #ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::Compute" << std::endl;
    #endif
        computed = true;
        Scalar diag = (Scalar)0.;
        Scalar Ldiag = (Scalar)0.;
#ifdef OPENMP
    #pragma omp parallel for
#endif
        for(int i = 0; i < dofs; ++i) 
        {
            if(identity_rows[i])
            {
                diag_entries[i] = (Scalar)1.;
            }
            else 
            {
                diag = (Scalar)0.;
                Ldiag = (Scalar)0.;
                for(SpMat::InnerIterator itL(L, i); itL; ++itL) 
                {
                    int k = itL.col();

                    if(i == k) 
                    {
                        Ldiag = itL.value();
                    }
                    else 
                    {
                        Scalar L_val = itL.value();
                        for (SpMat::InnerIterator itE(E, k); itE; ++itE) 
                        {
                            int j = itE.col();
                            Scalar E_val = itE.value();

                            if (j == i) 
                            {
                                diag += L_val * E_val;
                            }
                        }
                    }
                }
                diag_entries[i] = m + n * (Ldiag + diag);
            }
        }
    }

    void FinalizeMultiplication(const VectorX& x_array,const VectorX& tmp_array,VectorX& result_array, Scalar m_input, Scalar n_input) const 
    {
    #ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::FinalizeMultiplication" << std::endl;
    #endif
#ifdef OPENMP
    #pragma omp parallel for
#endif
        for(int q = 0; q < dofs; ++q) 
        {
            if(identity_rows[q])
            {
                result_array(q) = x_array(q);
            }
            else 
            {
                result_array(q) = m_input * x_array(q) + n_input * (tmp_array(q) + result_array(q));
            }
        }
        for(int i=0;i<nD;++i)
        {
            result_array(dofs + i) = (m_input - n_input) * x_array(dofs + i);
        }
    }

    void Multiply(const VectorBase& x,VectorBase& result) const
    {
    #ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::Multiply" << std::endl;
    #endif
        VectorX& x_array                            = x.data;
        VectorX& result_array                       = result.data;
        VectorX tmp_array                           = VectorX::Zero(dofs + nD);

        Interpolation(x_array, tmp_array);

        result.Clear();
        LnodiagTimesX(tmp_array, result_array);
        
        tmp_array.setZero();
        LdiagTimesX(x_array, tmp_array, m, n);

        // Finalize
        // r = -Lnodiag * E * x + -Ldiag * x
        FinalizeMultiplication(x_array, tmp_array, result_array, m, n);
    }

    void Multiply(const VectorBase& x, VectorBase& result, Scalar m_input, Scalar n_input) const
    {
    #ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::Multiply" << std::endl;
    #endif
        VectorX& x_array                            = x.data;
        VectorX& result_array                       = result.data;
        VectorX tmp_array                           = VectorX::Zero(dofs + nD);

        Interpolation(x_array, tmp_array);

        result.Clear();
        LnodiagTimesX(tmp_array, result_array);
        
        tmp_array.setZero();
        LdiagTimesX(x_array, tmp_array, m_input, n_input);

        // Finalize
        // r = -Lnodiag * E * x + -Ldiag * x
        FinalizeMultiplication(x_array, tmp_array, result_array, m_input, n_input);
    }

    void Interpolation(const VectorX& x_array,VectorX& result_array) const 
    {
    #ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::Interpolation" << std::endl;
    #endif
        result_array.setZero();
#ifdef OPENMP
    #pragma omp parallel for
#endif
        for (size_t k=0; k < E.outerSize(); ++k)
        {
            if(identity_rows[k])
            {
                result_array(k) = x_array(k);
            }
            else 
            {
                for (SpMat::InnerIterator it(E, k); it; ++it)
                {
                    size_t row = it.row();
                    size_t col = it.col();
                    Scalar val = it.value();
                    result_array(row) += val * x_array(col);
                }
            }
        }
    }

    void Project(VectorBase& v) const {}

    Scalar InnerProduct(const VectorBase& v1,const VectorBase& v2) const
    {
    #ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::InnerProduct" << std::endl;
    #endif
        return v1.data.dot(v2.data);
    }

    Scalar SquaredNorm(const VectorBase& v) const
    {
    #ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::SquaredNorm" << std::endl;
    #endif
        return v.data.squaredNorm();
    }

    void Diagonal_Inverse(const VectorX& input, VectorX& output) const
    {
#ifdef TRACK_WHERE
        std::cout << "BiCGStabSystem_IEM::Diagonal_Inverse" << std::endl;
#endif
        if(computed)
        {
    #ifdef OPENMP
        #pragma omp parallel for
    #endif
            for (size_t i = 0; i < dofs; ++i)
            {
                if(identity_rows[i])
                {
                    output(i) = input(i);
                }
                else
                {
                    output(i) = input(i) / diag_entries[i];
                }
            }
        }
        else
        {
    #ifdef OPENMP
        #pragma omp parallel for
    #endif
            for (size_t k = 0; k < L.outerSize(); ++k)
            {
                if(identity_rows[k])
                {
                    output(k) = input(k);
                }
                else
                {
                    bool found = false;
                    for (SpMat::InnerIterator it(L, k); it; ++it)
                    {
                        size_t row = it.row();
                        size_t col = it.col();
                        if(row == col)
                        {
                            Scalar val = m + it.value() * n;
                            output(row) = input(col) / val;
                            found = true;
                            break;
                        }
                    }  
                    if(!found)
                    {
                        output(k) = input(k);
                    }
                }
            }
            for(int i=0;i<nD;++i)
            {
                output(dofs + i) = input(dofs + i) / (m - n);
            }
        }
    }

    void Identity(const VectorX& input, VectorX& output) const
    {
    #ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::Identity" << std::endl;
    #endif
        output = input;
    }

    void Precondition(const VectorBase& x, VectorBase& result) const
    {
    #ifdef TRACK_WHERE
        std::cout << "GeneralizedBiCGStabSystem::Precondition" << std::endl;
    #endif
        if(precond == Preconditioner::diagonal)
        {
            Diagonal_Inverse(x.data, result.data);
        }
        // A * result = x
        else if(precond == Preconditioner::jacobi)
        {
            VectorX tmp_array = VectorX::Zero(dofs + 2);
            VectorX res_array = VectorX::Zero(dofs + 2);
            result.data.setZero();
            const Scalar omega = (Scalar)2./(Scalar)3.;
            for(int i = 0; i < precond_iters; ++i)
            {
                // 1. compute r = x - A * result
                VectorBase tmp(tmp_array);
                Multiply(result, tmp);
                res_array = x.data - tmp.data;
                // 2. compute r = r / D
                Diagonal_Inverse(res_array, res_array);
                // 3. compute result = result + omega * r
                for(int q = 0; q < dofs; ++q) 
                {
                    result.data(q) += omega * res_array(q);
                }
            }
        }
        else if(precond == Preconditioner::identity)
        {
            Identity(x.data, result.data);
        }
    }
};
}
