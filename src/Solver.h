#pragma once

#include <vector>
#include <iostream>

#include "Defines.h"
#include "Scalar.h"
#include "Tube.h"
#include "Helpers.h"

#include "BiCGStab.h"
#include "BiCGStabVector.h"
#include "GeneralizedBiCGStabSystem.h"

using namespace std;

namespace cpm {

class Solver {
    public:
        Solver();
        void Init(cpm::Scalar m, cpm::Scalar n, cpm::SpMat L_input, cpm::SpMat E_input, Preconditioner precond_input, 
            const vector<bool>& identity_rows_input, int nD=0);
        cpm::VectorX Solve(cpm::VectorX& rhs);
        void SetPrecond(int precond_iters);
        void Multiply(cpm::VectorX& x, cpm::VectorX& b, cpm::Scalar m, cpm::Scalar n, cpm::SpMat L, cpm::SpMat E);
    private:
        Preconditioner precond;
        int precond_iters;
        vector<bool> identity_rows;
    #ifdef CUSTOM_SOLVER
        BiCGStab bicgstab_solver;
        GeneralizedBiCGStabSystem gs;
        BiCGStabSystem* bicgstab_system;
        cpm::VectorX r_vec;
        cpm::VectorX v_vec;
        cpm::VectorX p_vec;
        cpm::VectorX s_vec;
        cpm::VectorX t_vec;
        cpm::VectorX y_vec;
        cpm::VectorX z_vec;
        cpm::VectorX r0_vec;
    #else
        SpMat A;
    #ifdef ESL
        Eigen::SparseLU<SpMat> slu_solver;
    #endif
    #ifdef EBCG
        Eigen::BiCGSTAB<SpMat> ebcg_diag_solver;
        Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<cpm::Scalar>> ebcg_jac_solver;
    #endif
    #endif
};

} // namespace cpm