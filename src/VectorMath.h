#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include "Scalar.h"

namespace cpm {

class VectorMath {
    public:
        VectorMath();
        
    private:

};

template <typename T>
Scalar DotProduct(const std::vector<T> &u, const std::vector<T> &v);

template <typename T>
Scalar Norm(const std::vector<T> &u);

template <typename T>
Scalar NormSquared(const std::vector<T> &u);

void Normalize(std::vector<Scalar> &u, Scalar zero_tol = 0.0);
std::vector<Scalar> CrossProduct(const std::vector<Scalar> &u, const std::vector<Scalar> &v);  
Scalar Distance(const std::vector<Scalar> &u, const std::vector<Scalar> &v);

std::vector<size_t> SortedIndices(const std::vector<Scalar> &u);

template<typename T>
bool isSortedVectorsEqual(std::vector<T> &first, std::vector<T> &second);

std::vector<Scalar> operator+(const std::vector<Scalar>& u, const std::vector<Scalar>& v);
std::vector<Scalar> operator-(const std::vector<Scalar>& u, const std::vector<Scalar>& v);
std::vector<Scalar>& operator+=(std::vector<Scalar>& u, const std::vector<Scalar>& v);
std::vector<Scalar>& operator-=(std::vector<Scalar>& u, const std::vector<Scalar>& v);
std::vector<Scalar>& operator*=(std::vector<Scalar>& u, const Scalar& scale);
std::vector<Scalar>& operator/=(std::vector<Scalar>& u, const Scalar& scale);
std::vector<Scalar>& operator+=(std::vector<Scalar>& u, const Scalar& constant);
std::vector<Scalar>& operator-=(std::vector<Scalar>& u, const Scalar& constant);
std::vector<Scalar> operator*(const Scalar& scale, const std::vector<Scalar>& u);
std::vector<Scalar> operator/(const std::vector<Scalar>& u, const Scalar& scale);

} // namespace cpm