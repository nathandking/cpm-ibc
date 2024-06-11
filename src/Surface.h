#pragma once

#include <iostream>
#include <vector>
#include <string>
#include "math.h"

#include "Scalar.h"
#include "VectorMath.h"
#include "SimplePolygonMesh.h"
#include "SurfaceSpecifier.h"

#include <LBFGS.h>
#include <LBFGSB.h>

#include <fcpw/fcpw.h>

using namespace std;

namespace cpm
{

    class Surface
    {
    public:
        Surface();
        Surface(SurfaceSpecifier &surface_specs);

        void SetSurfaceSpecs(SurfaceSpecifier &surface_specs);

        void Clear();

        void ClosestPoint(const vector<Scalar> &x, vector<Scalar> &cpx);
        void ClosestPoint(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy);

        Scalar GetDistanceToMesh(const vector<Scalar> &x) const;
        
        void OverideParameterizationPoints(const vector<vector<Scalar>> &xp);
        void OverideParameterizationFaces(const vector<vector<size_t>> &faces);

        void cpBar(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);

        vector<Scalar> IBCTangent(const vector<Scalar> &x);
        const Scalar IBCParameter(const vector<Scalar> &x);
        const vector<Scalar> CoordinateFromIBCParameter(const Scalar &ibc_parameter);

        int cpTri(const vector<Scalar> &x);

        const vector<Scalar>& boundingBox() const
        {
            return m_surface_specs.boundingBox();
        };

        const vector<Scalar>& surfaceParams() const
        {
            return m_surface_params;
        };

        const vector<vector<Scalar>>& xp() const
        {
            return m_mesh.vertices;
        };

        const vector<vector<size_t>>& faces() const
        {
            return m_mesh.faces;
        };

        const vector<Scalar>& thetap() const
        {
            return m_thetap;
        };

        const vector<Scalar>& phip() const
        {
            return m_phip;
        };

        const bool isOpen() const
        {
            return m_surface_specs.isOpen();
        };

        const string SurfaceType() const
        { 
            return m_surface_specs.SurfaceType(); 
        };

        const SurfaceSpecifier& SurfaceSpecs() const
        { 
            return m_surface_specs; 
        };

        size_t Dim() const {return m_dim;};

    private:

        SurfaceSpecifier m_surface_specs;
        size_t m_dim;
        vector<Scalar> m_surface_params;

        fcpw::Scene<3> m_scene; // initialize a 3d scene for closest points from triangulation

        // closest point functions
        void cpCircle(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpCircle3D(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpArc(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy);
        void cpDumbbell(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpDumbbellHalf(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, Scalar flip);
        void cpSphere(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpHemisphere(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy);
        void cpBall(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy);
        void cpTorus(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpPoint(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpLine(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy);
        void cpPlane(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy);
        void cpSquare(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpL(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpTorusLineSphere(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpSphereLine(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);

        void cpTorusKnot(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpPlanarCurve(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpPlanarCurve2(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpSurfaceOfRevolutionCurve(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpSurfaceOfRevolution(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpEllipsoid(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpParametrization(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, function<VectorX(const VectorX&)> r, function<VectorX(const VectorX&)> du, Scalar epsilon = 1e-6, size_t max_iter = 100);
        void cpParametrization(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, function<VectorX(const VectorX&)> r, function<VectorX(const VectorX&)> du, function<VectorX(const VectorX&)> dv, Scalar epsilon = 1e-6, size_t max_iter = 100);
        
        void cpDziuk(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpDecoTetrahedron(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpLevelSet(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, function<Scalar(const VectorX&)> phi, function<VectorX(const VectorX&)> gradphi, function<cpm::Matrix3(const VectorX&)> Hessianphi);
        
        int cpTri(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy);
        void cpNNPointCloud(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist);
        void cpPolyline2D(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy);

        // surface parametrization for visualization
        SimplePolygonMesh m_mesh;
        SimplePolygonMesh m_coarse_mesh;
        vector<Scalar> m_thetap;        
        vector<Scalar> m_phip;          

        void GetSurfaceParameterization();

        void paramCircle(size_t &Np);
        void paramCircle3D(size_t &Np);
        void paramArc(size_t &Np);
        void paramDumbbell(size_t &Np);
        void paramSphere(size_t &Np);
        void paramHemisphere(size_t &Np);
        void paramBall(size_t &Np);
        void paramTorus(size_t &Np);
        void paramLine(size_t &Np);
        void paramPlane(size_t &Np);
        void paramTorusLineSphere();
        void paramSphereLine();

        void paramTorusKnot(size_t &Np);
        void paramPlanarCurve(size_t &Np);
        void paramPlanarCurve2(size_t &Np);
        void paramSurfaceOfRevolutionCurve(size_t &Np);
        void paramSurfaceOfRevolution(size_t &Np);
        void paramEllipsoid(size_t &Np);

        void paramDziuk();
        void paramDecoTetrahedron();
        
        void GetParameterizationFaces(size_t &Np, size_t surface_dim);

        void InitializeFCPW();
        size_t FindNearestParamPoint(const vector<Scalar> &x);
        size_t BoundaryFromEndpointDistance(const vector<Scalar> &cpx);
        vector<size_t> IBCNearestTwoXpIndices(const vector<Scalar> &x, Scalar &min_dist);
    };

} // namespace cpm