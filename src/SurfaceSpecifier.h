#pragma once

#include <iostream>
#include <vector>
#include <string>
#include "math.h"

#include "Scalar.h"
#include "VectorMath.h"
#include "SimplePolygonMesh.h"
#include "IBCMetaData.h"

using namespace std;

namespace cpm
{

    class SurfaceSpecifier
    {
    public: 
        SurfaceSpecifier();

        void Clear();

        SurfaceSpecifier(string surf);
        SurfaceSpecifier(string surf, const vector<Scalar> &bounding_box, const vector<Scalar> &params, const size_t embed_dim);
        SurfaceSpecifier(string surf, const vector<Scalar> &bounding_box);

        SurfaceSpecifier(string surf, string coarse_surf);
        SurfaceSpecifier(string surf, string coarse_surf, const vector<Scalar> &bounding_box);

        SurfaceSpecifier(const SimplePolygonMesh &mesh);

        void SetSurface(string surf);
        void SetSurface(string surf, const vector<Scalar> &bounding_box, const vector<Scalar> &params, const size_t embed_dim);
        void SetSurface(string surf, const vector<Scalar> &bounding_box);

        void SetSurface(string surf, string coarse_surf);
        void SetSurface(string surf, string coarse_surf, const vector<Scalar> &bounding_box);

        const vector<Scalar>& boundingBox() const
        {
            return m_bounding_box;
        };

        void SetBoundingBox(const vector<Scalar> &bounding_box){ m_bounding_box = bounding_box; };

        const vector<Scalar>& surfaceParams() const
        {
            return m_surface_params;
        };

        const bool isOpen() const
        {
            return m_is_open;
        };

        const bool isPointCloud() const
        {
            return m_is_point_cloud;
        };

        void SetIsPointCloud(bool is_point_cloud){ m_is_point_cloud = is_point_cloud; };

        const string SurfaceType() const
        { 
            return m_surface; 
        };

        size_t Dim() const {return m_dim;};

        inline const SimplePolygonMesh& Mesh() const
        { 
            return m_mesh; 
        };

        inline const SimplePolygonMesh& CoarseMesh() const
        { 
            return m_coarse_mesh; 
        };

        const IBCMetaData& IBCMeta() const
        { 
            return m_ibc_meta; 
        };

        void SetIBCMetaData(IBCMetaData &meta){ m_ibc_meta = meta; };

    private:

        string m_surface;
        vector<Scalar> m_surface_params;
        vector<Scalar> m_bounding_box;
        size_t m_dim;
        bool m_is_open;
        bool m_is_point_cloud;

        SimplePolygonMesh m_mesh;
        SimplePolygonMesh m_coarse_mesh;

        IBCMetaData m_ibc_meta;

        void SetSurface(const vector<Scalar> &bounding_box);
        void SetSurface();

        void SetIsOpen(string surf);
    };

} // namespace cpm