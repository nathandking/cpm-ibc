#pragma once

#include <vector>
#include <iostream>
#include <iomanip> // for setw
#include <numeric> // for accumulate

#include "Scalar.h"
#include "Tube.h"
#include "FDMatrices.h"
#include "Interpolation.h"

#include "SimplePolygonMesh.h"

using namespace std;


namespace cpm {

class Helpers {
    public:
        Helpers();

        template <typename T>
        void WriteCSV(T &data, string filename);
        
        void WriteCSV(vector<vector<Scalar>> &data, string filename);

        template <typename T>
        void ReadCSV(string filename, T &data);

        void ReadCSV(string filename, vector<vector<Scalar>> &data);

        void SetupMatrices(Tube &t, SpMat &E, SpMat &L);

        void GetPlotInterpMatrix(Tube &t, SpMat &plotE);

        SpMat LaplacianSharp(SpMat &L, SpMat &E);
        SpMat LaplacianSharpWithDirichletConditions2ndOrder(const SpMat &L, const SpMat &E, const vector<size_t> &bdy, Scalar scale);
        void LaplacianSharpWithDirichletConditions2ndOrder(const SpMat &L, const SpMat &E, SpMat &Lfull, SpMat &Efull, const vector<size_t> &bdy);

        void ReplaceClosestPointWithcpBar(Tube &t, const vector<size_t> &bdy);

        SpMat ImplicitEulerMatrix(SpMat &L, SpMat &E, Scalar &dt);

        void SetIdentityRows(SpMat &A, const vector<bool> &identityRows);

        Scalar ConvergenceStudy(size_t nStart, size_t numLevels, function<vector<Scalar>(size_t&)> f);
        Scalar ConvergenceStudy(Scalar dxStart, size_t numLevels, function<Scalar(Scalar&)> f);

        void GetLaplacianMatrix(Tube &t, SpMat &L);
        void GetExtensionMatrix(Tube &t, SpMat &E);

    private:        

        SpMat AddBdyColumnsAndZeroRows(const SpMat &A, const vector<vector<Scalar>> &bdy_values);
        SpMat AddZeroColumnsAndRows(const SpMat &A, size_t full_size);

        vector<Scalar> ConvergenceOrder(vector<vector<Scalar>> &error);
        Scalar AverageConvergenceOrder(vector<vector<Scalar>> &error);
};

} // namespace cpm