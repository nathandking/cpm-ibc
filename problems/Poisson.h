#pragma once

#include "Scalar.h"
#include "Tube.h"
#include "Helpers.h"
#include "Visualization.h"

#include <Eigen/Sparse>

using namespace std;

namespace cpm {

class Poisson{
    public:
        Poisson(SurfaceSpecifier &surface_specs, Scalar dx, size_t bdy_order, vector<size_t> &which_bdy);
        Poisson(SurfaceSpecifier &surface_specs, Scalar dx);

        VectorX GetRHS(function<Scalar(const vector<Scalar>&)> f, function<Scalar(const vector<Scalar>&)> exact);
        VectorX GetRHS(function<Scalar(const Scalar&, const Scalar&)> f, function<Scalar(const Scalar&, const Scalar&)> exact);
        VectorX GetRHS(function<Scalar(const Scalar&, const Scalar&)> f);

        SpMat GetPlotInterpMatrix();

        size_t dim(){ return m_surface_tube.dim(); };

        const vector<vector<Scalar>>& cpx() const
        { 
            return m_surface_tube.cpx(); 
        };

        const vector<size_t>& bdy() const 
        { 
            return m_surface_tube.bdy(); 
        };

        const vector<Scalar>& thetap() const
        { 
            return m_surface_tube.surface().thetap(); 
        };

        const vector<Scalar>& phip() const
        { 
            return m_surface_tube.surface().phip(); 
        };

        const vector<vector<Scalar>>& xp() const
        { 
            return m_surface_tube.surface().xp(); 
        };

        const vector<vector<size_t>>& faces() const
        { 
            return m_surface_tube.surface().faces(); 
        };

        const vector<Scalar>& surfaceParams() const
        {
            return m_surface_tube.surface().surfaceParams();
        };

        Node node(size_t index) const
        {
            return m_surface_tube.TubeNode(index);
        }

        const Tube& tube() const
        {
            return m_surface_tube;
        }

        const size_t n() const
        {
            return m_surface_tube.n();
        }

        const SpMat GetLfull() const
        {
            return m_Lfull;
        }

        const SpMat GetEfull() const
        {
            return m_Efull;
        }

        const vector<bool>& GetIdentityRows() const
        {
            return m_identity_rows;
        }

    private:

        void Initialize(size_t bdy_order, vector<size_t> &which_bdy);

        size_t m_bdy_order;
        vector<size_t> m_which_bdy;

        Helpers m_h;

        Tube m_surface_tube;

        SpMat m_L;
        SpMat m_E; 

        SpMat m_Lfull;
        SpMat m_Efull;

        vector<bool> m_identity_rows;
        vector<size_t> m_bdy;

};

} // namespace cpm
