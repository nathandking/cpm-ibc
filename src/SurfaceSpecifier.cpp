#include "SurfaceSpecifier.h"

namespace cpm {
        

SurfaceSpecifier::SurfaceSpecifier()
{
    m_is_point_cloud = false;
}


SurfaceSpecifier::SurfaceSpecifier(string surf)
{
    m_is_point_cloud = false;
    SetSurface(surf);
}


SurfaceSpecifier::SurfaceSpecifier(string surf, const vector<Scalar> &bounding_box, const vector<Scalar> &params, const size_t embed_dim)
{
    m_is_point_cloud = false;
    SetSurface(surf, bounding_box, params, embed_dim);
}


SurfaceSpecifier::SurfaceSpecifier(string surf, const vector<Scalar> &bounding_box)
{
    m_is_point_cloud = false;
    SetSurface(surf, bounding_box);
}


SurfaceSpecifier::SurfaceSpecifier(string surf, string coarse_surf)
{
    m_is_point_cloud = false;
    SetSurface(surf, coarse_surf);
}


SurfaceSpecifier::SurfaceSpecifier(string surf, string coarse_surf, const vector<Scalar> &bounding_box)
{
    m_is_point_cloud = false;
    SetSurface(surf, coarse_surf, bounding_box);
}


SurfaceSpecifier::SurfaceSpecifier(const SimplePolygonMesh &mesh)
{
    m_mesh = mesh;
    m_dim = m_mesh.vertices[0].size();
    m_is_point_cloud = false;
    SetSurface();
}


void SurfaceSpecifier::Clear()
{
    m_surface_params.clear();
    m_bounding_box.clear();     
    m_is_point_cloud = false;
}


// surfaces with default parameters
void SurfaceSpecifier::SetSurface(string surf)
{   
#ifdef TRACK_WHERE
    std::cout << "SurfaceSpecifier::SetSurface(0)" << std::endl;
#endif

    if(surf.compare(surf.length() - 4, 4, ".obj") == 0)
    {
        m_mesh.readMeshFromFile(surf);
        m_dim = m_mesh.vertices[0].size();
        SetSurface();
    }
    else
    {
        vector<Scalar> bounding_box;
        vector<Scalar> params;
        size_t embed_dim;

        if(surf.compare("Circle") == 0)
        {
            params.push_back(1.0); // radius of the circle
            bounding_box.push_back(-1.0);
            bounding_box.push_back(1.0);
            embed_dim = 2;
        }
        else if(surf.compare("Circle3D") == 0)
        {
            params.push_back(1.0); // radius of the circle
            params.push_back(0.0); // centre x, y, z
            params.push_back(0.0);
            params.push_back(0.0);
            params.push_back(0); // normal x, y, z
            params.push_back(0);
            params.push_back(1);

            bounding_box.push_back(-1.0);
            bounding_box.push_back(1.0);
            embed_dim = 3;
        }
        else if(surf.compare("Arc") == 0)
        {
            params.push_back(1.0); // radius of the circle
            params.push_back(-3*M_PI_4); // angle1
            params.push_back(M_PI_4); // angle2
            bounding_box.push_back(-1.0);
            bounding_box.push_back(1.0);
            embed_dim = 2;
        }
        else if(surf.compare("Dumbbell") == 0)
        {
            params.push_back(0.1);
            bounding_box.push_back(-1.0);
            bounding_box.push_back(1.0);
            embed_dim = 2;
        }
        else if(surf.compare("Sphere") == 0 || surf.compare("Hemisphere") == 0 || surf.compare("Ball") == 0)
        {
            params.push_back(1.0); // radius of the sphere
            bounding_box.push_back(-1.0);
            bounding_box.push_back(1.0);
            embed_dim = 3;
        }
        else if(surf.compare("Torus") == 0)
        {
            params.push_back(2.0); // major radius of the torus
            params.push_back(1.0); // minor radius of the torus
            bounding_box.push_back(-3.0);
            bounding_box.push_back(3.0);
            embed_dim = 3;
        }
        else if(surf.compare("Line") == 0)
        {
            params.push_back(0.0); // start point, x = 0
            params.push_back(1.0); // end point, x = 1
            bounding_box.push_back(0.0);
            bounding_box.push_back(1.0);
            embed_dim = 2;
        }
        else if(surf.compare("Plane") == 0)
        {
            bounding_box.push_back(-1.0);
            bounding_box.push_back(1.0);
            embed_dim = 2;
        }
        else if(surf.compare("Square") == 0 || surf.compare("L") == 0)
        {
            bounding_box.push_back(-1.0);
            bounding_box.push_back(1.0);
            embed_dim = 2;
        }
        else if(surf.compare("Torus-Line-Sphere") == 0)
        {
            params.push_back(3.0); // major radius of the torus
            params.push_back(1.0); // minor radius of the torus
            params.push_back(1.25); // radius of sphere
            bounding_box.push_back(-4.1732457);
            bounding_box.push_back(4.0);
            embed_dim = 3;
        }
        else if(surf.compare("Sphere-Line") == 0)
        {
            params.push_back(1.0); // radius of sphere 1
            params.push_back(-1.5); // center x of sphere 1
            params.push_back(0.0); // center y of sphere 1
            params.push_back(0.0); // center z of sphere 1
            params.push_back(1.0); // radius of sphere 2
            params.push_back(1.5); // center x of sphere 2
            params.push_back(0.0); // center y of sphere 2
            params.push_back(0.0); // center z of sphere 2
            bounding_box.push_back(-2.5);
            bounding_box.push_back(2.5);
            embed_dim = 3;
        }
        else if(surf.compare("Torus Knot") == 0)
        {
            params.push_back(2.0);
            params.push_back(2.0);
            params.push_back(3.0);
            bounding_box.push_back(-3.0);
            bounding_box.push_back(3.0);
            embed_dim = 3;
        }
        else if(surf.compare("Planar Curve") == 0)
        {
            params.push_back(3.0);
            params.push_back(4.0); 
            params.push_back(0.5);
            bounding_box.push_back(-1.0);
            bounding_box.push_back(1.0);
            embed_dim = 2;
        }
        else if(surf.compare("Planar Curve 2") == 0)
        {
            params.push_back(0.7);
            params.push_back(0.05); 
            bounding_box.push_back(-1.0);
            bounding_box.push_back(1.0);
            embed_dim = 2;
        }
        else if(surf.compare("Surface Of Revolution Curve") == 0 || surf.compare("Surface Of Revolution") == 0)
        {
            params.push_back(3.0);
            params.push_back(4.0); 
            params.push_back(0.5);
            bounding_box.push_back(-2.5);
            bounding_box.push_back(2.5);
            embed_dim = 3;
        }
        else if(surf.compare("Ellipsoid") == 0)
        {
            params.push_back(1.5); // length of principal axis of ellipsoid
            params.push_back(1.0); 
            params.push_back(0.75);
            bounding_box.push_back(-1.5);
            bounding_box.push_back(1.5);
            embed_dim = 3;
        }
        else if(surf.compare("Dziuk") == 0)
        {
            // bounding_box.push_back(-1.5);               // for open curve and closed curve  // for problem 0
            bounding_box.push_back(-1.517234);       // for point                        
            bounding_box.push_back(1.5);
            embed_dim = 3;
        }
        else if(surf.compare("DecoTetrahedron") == 0)
        {
            bounding_box.push_back(-4.0);
            bounding_box.push_back(4.0);
            embed_dim = 3;
        }
        else 
        {
            throw std::runtime_error("No default parameters set for surface called '" + m_surface + "' in SurfaceSpecifier::SetSurface(string surf)");
        }

        SetSurface(surf, bounding_box, params, embed_dim);
    }
}


// allow user to specify surface parameters instead
void SurfaceSpecifier::SetSurface(string surf, const vector<Scalar> &bounding_box, const vector<Scalar> &params, size_t embed_dim)
{
#ifdef TRACK_WHERE
    std::cout << "SurfaceSpecifier::SetSurface(1)" << std::endl;
#endif
    m_dim = embed_dim;
    m_surface = surf;
    SetBoundingBox(bounding_box);
    m_surface_params = params;
    SetIsOpen(surf);
}


void SurfaceSpecifier::SetSurface(string surf, const vector<Scalar> &bounding_box)
{
#ifdef TRACK_WHERE
    std::cout << "SurfaceSpecifier::SetSurface(1)" << std::endl;
#endif
    if(surf.compare(surf.length() - 4, 4, ".obj") != 0)
    {
        throw std::runtime_error("Surface must be a .obj file, otherwise must specify surface parameters");
    }

    m_mesh.readMeshFromFile(surf);
    SetSurface(bounding_box);
}


void SurfaceSpecifier::SetSurface(string surf, string coarse_surf)
{
    if(coarse_surf.compare(coarse_surf.length() - 4, 4, ".obj") != 0)
    {
        throw std::runtime_error("Coarse surface must be a .obj file, otherwise must specify surface parameters");
    }

    m_coarse_mesh.readMeshFromFile(coarse_surf);
    SetSurface(surf);
}


void SurfaceSpecifier::SetSurface(string surf, string coarse_surf, const vector<Scalar> &bounding_box)
{
    if(coarse_surf.compare(coarse_surf.length() - 4, 4, ".obj") != 0)
    {
        throw std::runtime_error("Coarse surface must be a .obj file, otherwise must specify surface parameters");
    }

    m_coarse_mesh.readMeshFromFile(coarse_surf);
    SetSurface(surf, bounding_box);
}


// set surface for a triangulation, polyline, or point cloud (nearest-neigbhour cp, triangulation necessary)
void SurfaceSpecifier::SetSurface(const vector<Scalar> &bounding_box)
{
#ifdef TRACK_WHERE
    std::cout << "SurfaceSpecifier::SetSurface(3)" << std::endl;
#endif

    SetBoundingBox(bounding_box);

    if(m_mesh.faces[0].size() == 3)
    {
        m_dim = 3;
        if(m_is_point_cloud)
        {
            m_surface = "NN Point Cloud";
        }
        else
        {
            m_surface = "Triangulation";
        }

        m_is_open = false; // TO DO: implement for open surfaces
    }
    else if(m_mesh.faces[0].size() == 2)
    {
        if(m_mesh.vertices[0].size() == 2)
        {
            m_dim = 2;
            m_surface = "Polyline2D";
        }
        else
        {
            m_dim = 3;
            m_surface = "Polyline";
        }
        
        if(m_mesh.vertices[0] == m_mesh.vertices[m_mesh.nVertices() - 1] || m_mesh.faces[m_mesh.nFaces() - 1][0] == 0 || m_mesh.faces[m_mesh.nFaces() - 1][1] == 0)
        {
            m_is_open = false;
        }
        else
        {
            m_is_open = true;
        }
    }        
}


void SurfaceSpecifier::SetSurface()
{
#ifdef TRACK_WHERE
    std::cout << "SurfaceSpecifier::SetSurface(4)" << std::endl;
#endif

    vector<Scalar> bounding_box(2, 0.0);
    for(size_t i = 0; i < m_mesh.nVertices(); ++i)
    {
        for(size_t d = 0; d < m_dim; ++d)
        {
            bounding_box[0] = min(bounding_box[0], m_mesh.vertices[i][d]);
            bounding_box[1] = max(bounding_box[1], m_mesh.vertices[i][d]);
        }
    }

    SetSurface(bounding_box);
}



void SurfaceSpecifier::SetIsOpen(string surf)
{
#ifdef TRACK_WHERE
    std::cout << "SurfaceSpecifier::SetIsOpen" << std::endl;
#endif
    if(surf.compare("Circle") == 0 || surf.compare("Circle3D") == 0 || surf.compare("Dumbbell") == 0 || surf.compare("Sphere") == 0 || surf.compare("Ball") == 0 || surf.compare("Torus-Line-Sphere") == 0 || surf.compare("Sphere-Line") == 0 || surf.compare("Point") == 0 || surf.compare("Square") == 0 || surf.compare("L") == 0 || surf.compare("DecoTetrahedron") == 0 || surf.compare("SurfaceSpecifier Of Revolution") == 0 || surf.compare("SurfaceSpecifier Of Revolution Curve") == 0 || surf.compare("Ellipsoid") == 0 || surf.compare("Torus Knot") == 0 || surf.compare("Torus") == 0 || surf.compare("Planar Curve") == 0 || surf.compare("Planar Curve 2") == 0 || surf.compare("Dziuk") == 0 || surf.compare("NN Point Cloud") == 0)
    {
        m_is_open = false;
    }

    if(surf.compare("Line") == 0 || surf.compare("Arc") == 0 || surf.compare("Hemisphere") == 0 || surf.compare("Plane") == 0)
    {
        m_is_open = true;
    }
}


} // namespace cpm