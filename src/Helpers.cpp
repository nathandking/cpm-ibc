#include "Defines.h"
#include "Helpers.h"
namespace cpm {

Helpers::Helpers()
{

}


template <typename T>
void Helpers::WriteCSV(T &data, string filename)
{
    std::ofstream myfile;
    myfile.open (filename);
    for(size_t i = 0; i < data.size()-1; ++i)
    {
        myfile << data[i] << ",";
    }
    myfile << data[data.size() - 1] << "\n";
    myfile.close();
}


void Helpers::WriteCSV(vector<vector<Scalar>> &data, string filename)
{
    std::ofstream myfile;
    myfile.open (filename);
    for(size_t row = 0; row < data.size(); ++row)
    {
        for(size_t col = 0; col < data[row].size()-1; ++col)
        {
            myfile << data[row][col] << ",";
        }
        myfile << data[row][data[row].size() - 1] << "\n";
    }
    myfile.close();
}


template <typename T>
void Helpers::ReadCSV(string filename, T &data)
{
    vector<vector<Scalar>> temp;

    ReadCSV(filename, temp);

    if(temp.size() != 1)
    {
        throw std::runtime_error("In Helpers::ReadCSV(), " + filename + " does not contain one line of comma separated values");
    }

    data.resize(temp[0].size());
    for(size_t i = 0; i < temp[0].size(); ++i)
    {
        data[i] = temp[0][i];
    }
}


void Helpers::ReadCSV(string filename, vector<vector<Scalar>> &data)
{
    ifstream is(filename);
    if(is.is_open())
    {
        std::string line;
        while(getline(is, line))
        {
            std::stringstream ss(line);
            std::string token;
            
            vector<Scalar> temp;
            while(getline(ss, token, ','))
            {
                temp.push_back(stod(token));
            }

            data.push_back(temp);
        }
        is.close();
    }
}


// construct grid, compute closest points, build discrete Laplacian L, and build closest point extension E
void Helpers::SetupMatrices(Tube &t, SpMat &E, SpMat &L)
{
    L.resize(t.nNodes(), t.nNodes());
    GetLaplacianMatrix(t, L);

    E.resize(t.cpx().size(), t.nNodes());
    GetExtensionMatrix(t, E);
}


// construct plotting interpolation matrix
void Helpers::GetPlotInterpMatrix(Tube &t, SpMat &plotE)
{
    plotE.resize(t.surface().xp().size(), t.nNodes());

    Interpolation plotInterp(t, t.surface().xp());
    plotInterp.BuildInterpolationMatrix(t, plotE);
}


// Numerically stable discrete Laplace-Beltrami operator, see Macdonald & Ruuth 2009 (section 2.2.3 of https://steveruuth.org/wp-content/uploads/2020/10/icpm.pdf)
SpMat Helpers::LaplacianSharp(SpMat &L, SpMat &E)
{
    SpMat diagL;
    diagL = L.diagonal().asDiagonal();

    SpMat tmp = L - diagL;
    SpMat tmp1 = tmp * E;
    return diagL + tmp1;
}


SpMat Helpers::LaplacianSharpWithDirichletConditions2ndOrder(const SpMat &L, const SpMat &E, const vector<size_t> &bdy, Scalar scale)
{
    size_t rows = bdy.size();
    
    // get coefficients of E without boundary contribution
    vector<Eigen::Triplet<Scalar>> coeffs; 
    for (size_t k=0; k < E.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(E,k); it; ++it)
        {
            coeffs.emplace_back(it.row(), it.col(), it.value());
        }
    }

    // add in contribution from prescribed Dirichlet value
    for(size_t i = 0; i < rows; ++i)
    {
        if(bdy[i] > 0)
        {
            coeffs.emplace_back(i, rows + bdy[i] - 1, 2.0); // add in factor 2 for multiplying by the prescribe Dirichlet value to all the rows that are a boundary node
        }
    }

    SpMat Efull(rows + 2, rows + 2);
    Efull.setFromTriplets(coeffs.begin(), coeffs.end());

    // make Laplacian with two extra zero rows
    vector<Eigen::Triplet<Scalar>> Lcoeffs; // get coefficients of E without boundary contribution
    for (size_t k=0; k < L.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(L,k); it; ++it)
        {
            Lcoeffs.emplace_back(it.row(), it.col(), it.value());
        }
    }

    SpMat Lfull(rows + 2, rows + 2);
    Lfull.setFromTriplets(Lcoeffs.begin(), Lcoeffs.end());

    SpMat A = scale * LaplacianSharp(Lfull, Efull);
    A.coeffRef(rows, rows) = 1.0;
    A.coeffRef(rows + 1, rows + 1) = 1.0;

    return A;
}

void Helpers::LaplacianSharpWithDirichletConditions2ndOrder(const SpMat &L, const SpMat &E, SpMat &Lfull, SpMat &Efull, const vector<size_t> &bdy)
{
    size_t rows = bdy.size();
    
    // get coefficients of E without boundary contribution
    vector<Eigen::Triplet<Scalar>> coeffs; 
    for (size_t k=0; k < E.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(E,k); it; ++it)
        {
            coeffs.emplace_back(it.row(), it.col(), it.value());
        }
    }

    // add in contribution from prescribed Dirichlet value
    for(size_t i = 0; i < rows; ++i)
    {
        if(bdy[i] > 0)
        {
            coeffs.emplace_back(i, rows + bdy[i] - 1, 2.0); // add in factor 2 for multiplying by the prescribe Dirichlet value to all the rows that are a boundary node
        }
    }
    
    // Efull.resize(rows + 2, rows + 2);
    // Efull.setFromTriplets(coeffs.begin(), coeffs.end());
    
    SpMat E_tmp;
    E_tmp.resize(rows + 2, rows + 2);
    E_tmp.setFromTriplets(coeffs.begin(), coeffs.end());

    // make Laplacian with two extra zero rows
    vector<Eigen::Triplet<Scalar>> Lcoeffs; // get coefficients of E without boundary contribution
    for (size_t k=0; k < L.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(L,k); it; ++it)
        {
            Lcoeffs.emplace_back(it.row(), it.col(), it.value());
        }
    }
    // Lfull.resize(rows + 2, rows + 2);
    // Lfull.setFromTriplets(Lcoeffs.begin(), Lcoeffs.end());

    SpMat L_tmp;
    L_tmp.resize(rows + 2, rows + 2);
    L_tmp.setFromTriplets(Lcoeffs.begin(), Lcoeffs.end());

    Lfull = L_tmp;
    Efull = E_tmp;
}


SpMat Helpers::AddBdyColumnsAndZeroRows(const SpMat &A, const vector<vector<Scalar>> &bdy_values)
{
    size_t rows = A.rows();
    size_t num_bdries = bdy_values[0].size();

    // get coefficients of A
    vector<Eigen::Triplet<Scalar>> coeffs; 
    for (size_t k=0; k < A.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(A,k); it; ++it)
        {
            coeffs.emplace_back(it.row(), it.col(), it.value());
        }
    }

    // add in contribution from prescribed boundary value in the new column
    for(size_t i = 0; i < rows; ++i)
    {
        for(size_t nbdy = 0; nbdy < num_bdries; ++nbdy)
        {
            if(bdy_values[i][nbdy] != 0)
            {
                coeffs.emplace_back(i, rows + nbdy, bdy_values[i][nbdy]);
            }
        }
    }

    SpMat Afull(rows + num_bdries, rows + num_bdries);
    Afull.setFromTriplets(coeffs.begin(), coeffs.end());

    return Afull;
}


SpMat Helpers::AddZeroColumnsAndRows(const SpMat &A, size_t full_size)
{
    vector<Eigen::Triplet<Scalar>> coeffs; 
    for (size_t k=0; k < A.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(A,k); it; ++it)
        {
            coeffs.emplace_back(it.row(), it.col(), it.value());
        }
    }

    SpMat Afull(full_size, full_size);
    Afull.setFromTriplets(coeffs.begin(), coeffs.end());

    return Afull;
}


void Helpers::ReplaceClosestPointWithcpBar(Tube &t, const vector<size_t> &bdy)
{
    vector<Scalar> x(t.dim());
    vector<Scalar> cpx(t.dim());
    Scalar dist;
    for(size_t i = 0; i < t.cpx().size(); ++i)
    {
        if(bdy[i] > 0) // cpBar(x) = cp(x) for interior closest points, so only do computation for boundary points
        {
            for(size_t d = 0; d < t.dim(); ++d)
            {
                x[d] = t.x()[i][d];
                cpx[d] = t.cpx()[i][d];
            }

            Scalar dist;
            t.cpBar(x, cpx, dist);

            t.replaceClosestPoint(i, cpx);
        }
    }
}


// Implicit Euler time-stepping matrix
SpMat Helpers::ImplicitEulerMatrix(SpMat &L, SpMat &E, Scalar &dt)
{
    SpMat M = LaplacianSharp(L, E);
    SpMat I(M.rows(), M.cols());
    I.setIdentity();

    return I - dt * M;
}


// run a convergence study for a function f that returns an error value for all dx specified
// nStart: coarsest grid resolution power
// numLevels: number different grid resolutions
// f: any function that returns an error and a given grid resolution
Scalar Helpers::ConvergenceStudy(size_t nStart, size_t numLevels, function<vector<Scalar>(size_t&)> f)
{
    // compute all grid division powers for numLevels
    vector<size_t> gridDivisionPower(numLevels);
    for(size_t i = 0; i < numLevels; ++i)
    {
        gridDivisionPower[i] = i + nStart;
    }

    // compute error for each dx value for the function f
    vector<vector<Scalar>> error(numLevels, vector<Scalar>(2));
    for(size_t i = 0; i < numLevels; ++i)
    {
        error[i] = f(gridDivisionPower[i]); // must return error and dx 
    }
    
    vector<Scalar> order = ConvergenceOrder(error);
#ifdef DEMO
    cout << "dx " << error[0][1] << endl;
    cout << "error " << error[0][0] << endl;
    return -1.;
#else
    // output table of results
    cout << setw(8) << "dx" << '\t' << setw(8) << "error" << '\t' << setw(8) << "order" << endl;
    for(size_t i = 0; i < numLevels; ++i)
    {
        cout << setw(8) << error[i][1] << '\t' << setw(8) << error[i][0] << '\t' << setw(8) << order[i] << endl;
    }

    Scalar avg_conv_order = AverageConvergenceOrder(error);
    cout << "Average Convergence Order = " << avg_conv_order << endl;

    return avg_conv_order;
#endif
}


Scalar Helpers::ConvergenceStudy(Scalar dxStart, size_t numLevels, function<Scalar(Scalar&)> f)
{
    // compute all grid division powers for numLevels
    vector<Scalar> dx(numLevels);
    for(size_t i = 0; i < numLevels; ++i)
    {
        dx[i] = pow(0.5, i) * dxStart;
    }

    // compute error for each dx value for the function f
    vector<vector<Scalar>> error(numLevels, vector<Scalar>(2));
    for(size_t i = 0; i < numLevels; ++i)
    {
        error[i][0] = f(dx[i]); // must return error and dx 
        error[i][1] = dx[i];
    }
    
    vector<Scalar> order = ConvergenceOrder(error);

#ifdef DEMO
    cout << "dx " << error[0][1] << endl;
    cout << "error " << error[0][0] << endl;
    return -1.;
#else
    // output table of results
    cout << setw(8) << "dx" << '\t' << setw(8) << "error" << '\t' << setw(8) << "order" << endl;
    for(size_t i = 0; i < numLevels; ++i)
    {
        cout << setw(8) << error[i][1] << '\t' << setw(8) << error[i][0] << '\t' << setw(8) << order[i] << endl;
    }

    Scalar avg_conv_order = AverageConvergenceOrder(error);
    cout << "Average Convergence Order = " << avg_conv_order << endl;

    return avg_conv_order;
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//          Private Functions
/////////////////////////////////////////////////////////////////////////////////////////////////

// build Laplacian matrix for grid points within tube
void Helpers::GetLaplacianMatrix(Tube &t, SpMat &L)
{
    FDMatrices fdmat;
    fdmat.BuildLaplacianMatrix(t, L);
}


// interpolation on the tubed grid
void Helpers::GetExtensionMatrix(Tube &t, SpMat &E)
{
    Interpolation interp(t); 
    interp.BuildInterpolationMatrix(t, E);
}


// Make rows in matrix rows of the identity matrix when identity_rows[i] == true
void Helpers::SetIdentityRows(SpMat &A, const vector<bool> &identity_rows)
{
    // Note: this iterates over columns, and care must be taken for columns that have no diagonal entry
    for (size_t k=0; k < A.outerSize(); ++k)
    {
        bool diagonal_entry_exists = false;
        for (SpMat::InnerIterator it(A,k); it; ++it)
        {
            if(identity_rows[it.row()])
            {
                if(it.row() == it.col())
                {
                    diagonal_entry_exists = true;
                    it.valueRef() = 1.0;
                }
                else
                {
                    it.valueRef() = 0.0;
                }
            }
        }

        if(!diagonal_entry_exists && identity_rows[k]) // if the diagonal entry is zero we must insert a new element into the sparse matrix
        {
            A.insert(k,k) = 1.0;
        }
    }
    A.prune(0,0);  

    // this can be much much slower (6.5 hours vs 20 minutes for the above on the sphere with circle ibc) 
    // for(size_t i = 0; i < identity_rows.size(); ++i)
    // {
    //     if(identity_rows[i])
    //     {
    //         A.row(i) *= 0.0;
    //         A.coeffRef(i,i) = 1.0;
    //     }
    // }
    // A.prune(0,0);
}


// Compute the order of convergence
vector<Scalar> Helpers::ConvergenceOrder(vector<vector<Scalar>> &error)
{
    vector<Scalar> order(error.size());
    order[0] = NAN;
    for(size_t i = 1; i < error.size(); ++i)
    {
        order[i] = log(error[i-1][0]/error[i][0]) / log(error[i-1][1]/error[i][1]);
    }

    return order;
}


Scalar Helpers::AverageConvergenceOrder(vector<vector<Scalar>> &error)
{
    vector<Scalar> order = ConvergenceOrder(error);
    return accumulate(order.begin() + 1, order.end(), 0.0) / (order.size() - 1);
}

template void Helpers::WriteCSV<cpm::VectorX>(cpm::VectorX &data, string filename);
template void Helpers::WriteCSV<vector<size_t>>(vector<size_t> &data, string filename);

template void Helpers::ReadCSV<cpm::VectorX>(string filename, cpm::VectorX &data);
template void Helpers::ReadCSV<vector<size_t>>(string filename, vector<size_t> &data);

} // namespace cpm