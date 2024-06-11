#pragma once

#include "Defines.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#define USE_DOUBLES

namespace cpm {

#ifdef USE_DOUBLES
typedef double Scalar;
#else
typedef float Scalar;
#endif

typedef Eigen::Triplet<cpm::Scalar> T;

typedef Eigen::Matrix<cpm::Scalar, -1, 1> VectorX;
typedef Eigen::Matrix<cpm::Scalar, 2, 1> Vector2;
typedef Eigen::Matrix<cpm::Scalar, 3, 1> Vector3;

typedef Eigen::Matrix<cpm::Scalar, 2, 2> Matrix2;
typedef Eigen::Matrix<cpm::Scalar, 3, 3> Matrix3;
typedef Eigen::Matrix<cpm::Scalar, 4, 4> Matrix4;

#ifdef CUSTOM_SOLVER
typedef Eigen::SparseMatrix<cpm::Scalar, Eigen::RowMajor> SpMat;
#else
typedef Eigen::SparseMatrix<cpm::Scalar, Eigen::ColMajor> SpMat;
#endif

} // namespace cpm