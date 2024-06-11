#include "VectorMath.h"

namespace cpm {

VectorMath::VectorMath()
{

}


template <typename T>
Scalar DotProduct(const std::vector<T> &u, const std::vector<T> &v)
{
    assert(u.size() == v.size());

    Scalar dot_product = 0.0;
    for(size_t d = 0; d < u.size(); ++d)
    {
        dot_product += u[d] * v[d];
    }

    return dot_product;
}


template <typename T>
Scalar Norm(const std::vector<T> &u)
{
    return sqrt(NormSquared(u));
}


template <typename T>
Scalar NormSquared(const std::vector<T> &u)
{
    return DotProduct(u, u);
}


void Normalize(std::vector<Scalar> &u, Scalar zero_tol)
{
    Scalar norm_u = Norm(u);

    if(norm_u > zero_tol)
    {
        for(size_t d = 0; d < u.size(); ++d)
        {
            u[d] /= norm_u;
        }
    }
    else
    {
        for(size_t d = 0; d < u.size(); ++d)
        {
            u[d] = 0.0;
        }
    }
}


std::vector<Scalar> CrossProduct(const std::vector<Scalar> &u, const std::vector<Scalar> &v)
{
    assert(u.size() == v.size());
    assert(u.size() == 3);

    std::vector<Scalar> cross(u.size());

    cross[0] = u[1] * v[2] - u[2] * v[1];
    cross[1] = u[2] * v[0] - u[0] * v[2];
    cross[2] = u[0] * v[1] - u[1] * v[0];

    return cross;
}


// distance between two points u and v, i.e., Euclidean distance ||u - v||_2
Scalar Distance(const std::vector<Scalar> &u, const std::vector<Scalar> &v)
{
    assert(u.size() == v.size());

    Scalar distance = 0.0;
    for(size_t d = 0; d < u.size(); ++d) 
    {
        distance += pow(u[d] - v[d], 2);
    }
    distance = sqrt(distance);

    return distance;
}


// get the indices of u in order of sorted vector from smallest to largest
std::vector<size_t> SortedIndices(const std::vector<Scalar> &u)
{
    std::vector<size_t> sorted_indices(u.size());

    iota(sorted_indices.begin(), sorted_indices.end(), 0);

    sort(sorted_indices.begin(), sorted_indices.end(), [&u](size_t i1, size_t i2) {return u[i1] < u[i2];});

    return sorted_indices;
}


// elementwise addition of vectors
std::vector<Scalar> operator+(const std::vector<Scalar>& u, const std::vector<Scalar>& v)
{
    assert(u.size() == v.size());
    std::vector<Scalar> result;

    result.reserve(u.size());

    std::transform(u.begin(), u.end(), v.begin(), std::back_inserter(result), std::plus<Scalar>());

    return result;
}


// elementwise subtraction of vectors
std::vector<Scalar> operator-(const std::vector<Scalar>& u, const std::vector<Scalar>& v)
{
    assert(u.size() == v.size());
    std::vector<Scalar> result;

    result.reserve(u.size());

    std::transform(u.begin(), u.end(), v.begin(), std::back_inserter(result), std::minus<Scalar>());

    return result;
}


// elementwise addition of vectors
std::vector<Scalar>& operator+=(std::vector<Scalar>& u, const std::vector<Scalar>& v)
{
    assert(u.size() == v.size());

    std::transform(u.begin(), u.end(), v.begin(), u.begin(), std::plus<Scalar>());

    return u;
}


// elementwise subtraction of vectors
std::vector<Scalar>& operator-=(std::vector<Scalar>& u, const std::vector<Scalar>& v)
{
    assert(u.size() == v.size());

    std::transform(u.begin(), u.end(), v.begin(), u.begin(), std::minus<Scalar>());

    return u;
}


// elementwise multiplication of vector by scalar
std::vector<Scalar>& operator*=(std::vector<Scalar>& u, const Scalar& scale)
{
    for(size_t d = 0; d < u.size(); ++d) 
    {
        u[d] *= scale;
    }

    return u;
}


// elementwise division of vector by scalar
std::vector<Scalar>& operator/=(std::vector<Scalar>& u, const Scalar& scale)
{
    for(size_t d = 0; d < u.size(); ++d) 
    {
        u[d] /= scale;
    }

    return u;
}


// elementwise addition of vector by scalar
std::vector<Scalar>& operator+=(std::vector<Scalar>& u, const Scalar& constant)
{
    for(size_t d = 0; d < u.size(); ++d) 
    {
        u[d] += constant;
    }

    return u;
}


// elementwise subtraction of vector by scalar
std::vector<Scalar>& operator-=(std::vector<Scalar>& u, const Scalar& constant)
{
    for(size_t d = 0; d < u.size(); ++d) 
    {
        u[d] -= constant;
    }

    return u;
}


std::vector<Scalar> operator*(const Scalar& scale, const std::vector<Scalar>& u)
{
    std::vector<Scalar> result(u.size());
    for(size_t d = 0; d < u.size(); ++d) 
    {
        result[d] = scale * u[d];
    }

    return result;
}


std::vector<Scalar> operator/(const std::vector<Scalar>& u, const Scalar& scale)
{
    std::vector<Scalar> result(u.size());
    for(size_t d = 0; d < u.size(); ++d) 
    {
        result[d] = u[d] / scale;
    }

    return result;
}


template<typename T>
bool isSortedVectorsEqual(std::vector<T> &first, std::vector<T> &second)
{
    if (first.size() != second.size()) {
        return false;
    }
 
    std::sort(first.begin(), first.end());
    std::sort(second.begin(), second.end());
 
    return first == second;
}

template bool isSortedVectorsEqual<size_t>(std::vector<size_t> &first, std::vector<size_t> &second);

template Scalar DotProduct<Scalar>(const std::vector<Scalar> &u, const std::vector<Scalar> &v);
template Scalar DotProduct<ssize_t>(const std::vector<ssize_t> &u, const std::vector<ssize_t> &v);
template Scalar DotProduct<size_t>(const std::vector<size_t> &u, const std::vector<size_t> &v);

template Scalar Norm<Scalar>(const std::vector<Scalar> &u);
template Scalar Norm<ssize_t>(const std::vector<ssize_t> &u);
template Scalar Norm<size_t>(const std::vector<size_t> &u);

template Scalar NormSquared<Scalar>(const std::vector<Scalar> &u);
template Scalar NormSquared<ssize_t>(const std::vector<ssize_t> &u);
template Scalar NormSquared<size_t>(const std::vector<size_t> &u);

} // namespace cpm