#include "Solver.h"
#include "BiCGStab.h"
#include "GeneralizedBiCGStabSystem.h"
#include "BiCGStabVector.h"
namespace cpm {

Solver::Solver() {}
void Solver::Init(cpm::Scalar m, cpm::Scalar n, cpm::SpMat L, cpm::SpMat E, Preconditioner precond_input, 
                const vector<bool>& identity_rows_input, int nD)
{
#ifdef TRACK_WHERE
    std::cout << "GeneralizedBiCGStabSystem::Init" << std::endl;
#endif
    precond = precond_input;
    identity_rows = identity_rows_input;
#ifndef CUSTOM_SOLVER
    Helpers m_h;
    cpm::SpMat M = -m_h.LaplacianSharp(L, E);
    for(int i = M.rows() - 1; i >= M.rows() - nD; --i)
    {
        M.coeffRef(i, i) = (Scalar)1.0;
    }
    cpm::SpMat I(M.rows(), M.cols());
    I.setIdentity();
    A = m * I - n * M;
    m_h.SetIdentityRows(A, identity_rows);
    A.makeCompressed();

#ifdef ESL
    slu_solver.compute(A);
    if(slu_solver.info()!=Eigen::Success) {
        std::cout << "SparseLU.compute() failed" << std::endl;
    }
#endif
#ifdef EBCG
    precond = precond_input;
    if(precond == Preconditioner::diagonal)
    {
        ebcg_diag_solver.compute(A);
        if(ebcg_diag_solver.info()!=Eigen::Success) {
            std::cout << "BiCGStab.compute() failed" << std::endl;
        }    
        ebcg_diag_solver.setTolerance(1.e-10);
    }
    else if(precond == Preconditioner::jacobi)
    {
        ebcg_jac_solver.compute(A);
        if(ebcg_jac_solver.info()!=Eigen::Success) {
            std::cout << "BiCGStab.compute() failed" << std::endl;
        }    
        ebcg_jac_solver.setTolerance(1.e-10);
    }
#endif

#else
    gs.SetComponents(m, n, L, E, precond, identity_rows, nD);
    size_t DOFs = L.rows();
    r_vec=cpm::VectorX::Zero(DOFs);
    v_vec=cpm::VectorX::Zero(DOFs);
    p_vec=cpm::VectorX::Zero(DOFs);
    s_vec=cpm::VectorX::Zero(DOFs);
    t_vec=cpm::VectorX::Zero(DOFs);
    y_vec=cpm::VectorX::Zero(DOFs);
    z_vec=cpm::VectorX::Zero(DOFs);
    r0_vec=cpm::VectorX::Zero(DOFs);
    bicgstab_system = (BiCGStabSystem*)&gs;
#endif
}
cpm::VectorX Solver::Solve(cpm::VectorX& rhs)
{
    cpm::VectorX u = cpm::VectorX::Zero(rhs.size());
#ifndef CUSTOM_SOLVER
#ifdef ESL
    u = slu_solver.solve(rhs);
    if(slu_solver.info()!=Eigen::Success) {
        std::cout << "SparseLU.solve() failed" << std::endl;
    }
#endif
#ifdef EBCG
    if(precond == Preconditioner::diagonal)
    {
        u = ebcg_diag_solver.solve(rhs);
        if(ebcg_diag_solver.info()!=Eigen::Success) {
            std::cout << "BiCGSTAB.solve() failed" << std::endl;
        }
    }
    else if(precond == Preconditioner::jacobi)
    {
        u = ebcg_jac_solver.solve(rhs);
        if(ebcg_jac_solver.info()!=Eigen::Success) {
            std::cout << "BiCGSTAB.solve() failed" << std::endl;
        }
    }
#endif
#else
    size_t DOFs = rhs.size();
    BiCGStab bicgstab_solver;
    u=cpm::VectorX::Zero(DOFs);
    BiCGStabVector u_V(u);
    BiCGStabVector b_V(rhs);
    BiCGStabVector r_V(r_vec);
    BiCGStabVector r0_V(r0_vec);
    BiCGStabVector v_V(v_vec);
    BiCGStabVector p_V(p_vec);
    BiCGStabVector s_V(s_vec);
    BiCGStabVector t_V(t_vec);
    BiCGStabVector y_V(y_vec);
    BiCGStabVector z_V(z_vec);
    int iters = 0;
    bicgstab_solver.Solve((*bicgstab_system),u_V,b_V,r_V,r0_V,v_V,p_V,s_V,t_V,y_V,z_V,(Scalar)1.e-10,1,100000,iters);
#endif

    return u;
}
void Solver::SetPrecond(int precond_iters)
{
#ifdef CUSTOM_SOLVER
    gs.SetPrecond(precond_iters);
#endif
}

void Solver::Multiply(cpm::VectorX& x, cpm::VectorX& b, cpm::Scalar m, cpm::Scalar n,
                        cpm::SpMat L, cpm::SpMat E)
{
#ifdef CUSTOM_SOLVER
    BiCGStabVector u_V(x);
    BiCGStabVector b_V(b);
    gs.Multiply(u_V, b_V, m, n);
#else
    Helpers m_h;
    cpm::SpMat M = m_h.LaplacianSharp(L, E);
    cpm::SpMat I(M.rows(), M.cols());
    I.setIdentity();
    SpMat Acur = m * I + n * M;
    m_h.SetIdentityRows(Acur, identity_rows);
    b = Acur * x;
#endif
}

} // namespace cpm
