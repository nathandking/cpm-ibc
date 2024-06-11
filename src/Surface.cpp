#include "Surface.h"

namespace cpm {

Surface::Surface()
{

}


Surface::Surface(SurfaceSpecifier& surface_specs)
{
    SetSurfaceSpecs(surface_specs);
}


void Surface::SetSurfaceSpecs(SurfaceSpecifier& surface_specs)
{
    m_surface_specs = surface_specs;
    m_dim = m_surface_specs.Dim();

    if(SurfaceType().compare("NN Point Cloud") == 0 || SurfaceType().compare("Triangulation") == 0 || SurfaceType().compare("Polyline") == 0)
    {
        m_mesh = m_surface_specs.Mesh();
        m_coarse_mesh = m_surface_specs.CoarseMesh();
        InitializeFCPW();
    }
    else
    {
        m_surface_params = m_surface_specs.surfaceParams();
        GetSurfaceParameterization();
    }
}


void Surface::Clear()
{
    m_scene.getSceneData()->clearObjectData();
    m_scene.getSceneData()->clearAggregateData();

    // surface parametrization for visualization
    m_mesh.vertices.clear();
    m_mesh.faces.clear(); 
    m_coarse_mesh.vertices.clear();
    m_coarse_mesh.faces.clear(); 

    m_thetap.clear();        
    m_phip.clear();         
}


void Surface::ClosestPoint(const vector<Scalar> &x, vector<Scalar> &cpx)
{
    Scalar dist;
    size_t bdy;
    ClosestPoint(x, cpx, dist, bdy);
}


// Closest point computation
void Surface::ClosestPoint(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy)
{
    if(SurfaceType().compare("Triangulation") == 0 || SurfaceType().compare("Polyline") == 0)
    {
        cpTri(x, cpx, dist, bdy);
    }
    else if(SurfaceType().compare("NN Point Cloud") == 0)
    {
        bdy = 0;
        cpNNPointCloud(x, cpx, dist);
    }
    else if(SurfaceType().compare("Polyline2D") == 0)
    {
        cpPolyline2D(x, cpx, dist, bdy);
    }
    else if(SurfaceType().compare("Circle") == 0)
    {
        bdy = 0;
        cpCircle(x, cpx, dist);
    }
    else if(SurfaceType().compare("Circle3D") == 0)
    {
        bdy = 0;
        cpCircle3D(x, cpx, dist);
    }
    else if(SurfaceType().compare("Arc") == 0)
    {
        cpArc(x, cpx, dist, bdy);
    }
    else if(SurfaceType().compare("Dumbbell") == 0)
    {
        bdy = 0;
        cpDumbbell(x, cpx, dist);
    }
    else if(SurfaceType().compare("Sphere") == 0)
    {
        bdy = 0;
        cpSphere(x, cpx, dist);
    }
    else if(SurfaceType().compare("Hemisphere") == 0)
    {
        cpHemisphere(x, cpx, dist, bdy);
    }
    else if(SurfaceType().compare("Ball") == 0)
    {
        cpBall(x, cpx, dist, bdy);
    }
    else if(SurfaceType().compare("Torus") == 0)
    {
        bdy = 0;
        cpTorus(x, cpx, dist);
    }
    else if(SurfaceType().compare("Point") == 0)
    {
        bdy = 0;
        cpPoint(x, cpx, dist);
    }
    else if(SurfaceType().compare("Line") == 0)
    {
        cpLine(x, cpx, dist, bdy);
    }
    else if(SurfaceType().compare("Plane") == 0)
    {
        cpPlane(x, cpx, dist, bdy);
    }
    else if(SurfaceType().compare("Square") == 0)
    {
        bdy = 0;
        cpSquare(x, cpx, dist);
    }
    else if(SurfaceType().compare("L") == 0)
    {
        bdy = 0;
        cpL(x, cpx, dist);
    }
    else if(SurfaceType().compare("Torus-Line-Sphere") == 0)
    {
        bdy = 0;
        cpTorusLineSphere(x, cpx, dist);
    }
    else if(SurfaceType().compare("Sphere-Line") == 0)
    {
        bdy = 0;
        cpSphereLine(x, cpx, dist);
    }
    else if(SurfaceType().compare("Torus Knot") == 0)
    {
        bdy = 0;
        cpTorusKnot(x, cpx, dist);
    }
    else if(SurfaceType().compare("Planar Curve") == 0)
    {
        bdy = 0;
        cpPlanarCurve(x, cpx, dist);
    }
    else if(SurfaceType().compare("Planar Curve 2") == 0)
    {
        bdy = 0;
        cpPlanarCurve2(x, cpx, dist);
    }
    else if(SurfaceType().compare("Surface Of Revolution Curve") == 0)
    {
        bdy = 0;
        cpSurfaceOfRevolutionCurve(x, cpx, dist);
    }
    else if(SurfaceType().compare("Surface Of Revolution") == 0)
    {
        bdy = 0;
        cpSurfaceOfRevolution(x, cpx, dist);
    }
    else if(SurfaceType().compare("Ellipsoid") == 0)
    {
        bdy = 0;
        cpEllipsoid(x, cpx, dist);
    }
    else if(SurfaceType().compare("Dziuk") == 0)
    {
        bdy = 0;
        cpDziuk(x, cpx, dist);
    }
    else if(SurfaceType().compare("DecoTetrahedron") == 0)
    {
        bdy = 0;
        cpDecoTetrahedron(x, cpx, dist);
    }
    else
    {
        throw std::runtime_error("Invalid surface called '" + SurfaceType() + "' in Surface::ClosestPoint");
    }
}


void Surface::GetSurfaceParameterization()
{
#ifdef TRACK_WHERE
    std::cout << "Surface::GetSurfaceParameterization" << std::endl;
#endif
    if(SurfaceType().compare("Triangulation") == 0 || SurfaceType().compare("Polyline") == 0 || SurfaceType().compare("Polyline2D") == 0 || SurfaceType().compare("NN Point Cloud") == 0)
    {
        // surface parameterization already comes from vertices
    }
    else if(SurfaceType().compare("Circle") == 0)
    {
        size_t Np = 100;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        paramCircle(Np);
        GetParameterizationFaces(Np, 1);
    }
    else if(SurfaceType().compare("Circle3D") == 0)
    {
        size_t Np = 100;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        paramCircle3D(Np);
        GetParameterizationFaces(Np, 1);
    }
    else if(SurfaceType().compare("Arc") == 0)
    {
        size_t Np = 100;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        paramArc(Np);
        GetParameterizationFaces(Np, 1);
    }
    else if(SurfaceType().compare("Dumbbell") == 0)
    {
        size_t Np = 10000;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        paramDumbbell(Np);
        GetParameterizationFaces(Np, 1);
    }
    else if(SurfaceType().compare("Sphere") == 0)
    {
        size_t Np = 40000;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        m_phip.resize(Np);
        paramSphere(Np);
        GetParameterizationFaces(Np, 2);
    }
    else if(SurfaceType().compare("Hemisphere") == 0)
    {
        size_t Np = 10000;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        m_phip.resize(Np);
        paramHemisphere(Np);
        GetParameterizationFaces(Np, 2);
    }
    else if(SurfaceType().compare("Ball") == 0)
    {
        size_t Np = 100;
        paramBall(Np);
    }
    else if(SurfaceType().compare("Torus") == 0)
    {
        size_t Np = 90000;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        m_phip.resize(Np);
        paramTorus(Np);
        GetParameterizationFaces(Np, 2);
    }
    else if(SurfaceType().compare("Point") == 0)
    {
        m_mesh.vertices.resize(1, vector<Scalar>(m_dim));
        for(size_t d = 0; d < m_dim; ++d)
        {
            m_mesh.vertices[0][d] = surfaceParams()[d];
        }
    }
    else if(SurfaceType().compare("Line") == 0)
    {
        size_t Np = 100;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        paramLine(Np);
        GetParameterizationFaces(Np, 1);
    }
    else if(SurfaceType().compare("Plane") == 0)
    {
        size_t Np = 250000;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_phip.resize(Np);
        m_thetap.resize(Np);
        paramPlane(Np);
        GetParameterizationFaces(Np, 2);
    }
    else if(SurfaceType().compare("Torus-Line-Sphere") == 0)
    {
        paramTorusLineSphere();
    }
    else if(SurfaceType().compare("Sphere-Line") == 0)
    {
        paramSphereLine();
    }
    else if(SurfaceType().compare("Torus Knot") == 0)
    {
        size_t Np = 1000;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        paramTorusKnot(Np);
        GetParameterizationFaces(Np, 1);
    }
    else if(SurfaceType().compare("Planar Curve") == 0)
    {
        size_t Np = 10000;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        paramPlanarCurve(Np);
        GetParameterizationFaces(Np, 1);
    }
    else if(SurfaceType().compare("Planar Curve 2") == 0)
    {
        size_t Np = 1000;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        paramPlanarCurve2(Np);
        GetParameterizationFaces(Np, 1);
    }
    else if(SurfaceType().compare("Surface Of Revolution Curve") == 0)
    {
        size_t Np = 1000;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        m_phip.resize(Np);
        paramSurfaceOfRevolutionCurve(Np);
        GetParameterizationFaces(Np, 1);
    }
    else if(SurfaceType().compare("Surface Of Revolution") == 0)
    {
        size_t Np = 40000;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        m_phip.resize(Np);
        paramSurfaceOfRevolution(Np);
        GetParameterizationFaces(Np, 2);
    }
    else if(SurfaceType().compare("Ellipsoid") == 0)
    {
        size_t Np = 10000;
        m_mesh.vertices.resize(Np, vector<Scalar>(m_dim));
        m_thetap.resize(Np);
        m_phip.resize(Np);
        paramEllipsoid(Np);
        GetParameterizationFaces(Np, 2);
    }
    else if(SurfaceType().compare("Dziuk") == 0)
    {
        paramDziuk();
    }
    else if(SurfaceType().compare("DecoTetrahedron") == 0)
    {
        paramDecoTetrahedron();
    }
    else if(SurfaceType().compare("Square") == 0 || SurfaceType().compare("L") == 0)
    {
        cout << "Warning in Surface::GetSurfaceParameterization - " + SurfaceType() + " does not contain a surface parameterization. This may cause errors in downstream functions that use Surface members m_mesh.vertices, m_thetap, etc." << endl;
    }
    else
    {
        throw std::runtime_error("Invalid surface called '" + SurfaceType() + "' in Surface::GetSurfaceParameterization");
    }
}


Scalar Surface::GetDistanceToMesh(const vector<Scalar> &x) const
{
    // perform a closest point query
    fcpw::Interaction<3> interaction;
    fcpw::Vector3 queryPoint;

    queryPoint[0] = x[0];
    queryPoint[1] = x[1];
    queryPoint[2] = x[2];

    m_scene.findClosestPoint(queryPoint, interaction);

    return interaction.signedDistance(queryPoint);
}


void Surface::OverideParameterizationPoints(const vector<vector<Scalar>> &xp)
{
    m_mesh.vertices.clear();
    m_mesh.vertices = xp;
}


void Surface::OverideParameterizationFaces(const vector<vector<size_t>> &faces)
{
    m_mesh.faces.clear();
    m_mesh.faces = faces;
}


void Surface::cpBar(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    vector<Scalar> projection(x.size());
    for(size_t d = 0; d < x.size(); ++d)
    {
        projection[d] = 2.0 * cpx[d] - x[d];
    }

    size_t bdy;
    ClosestPoint(projection, cpx, dist, bdy);
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//          CLOSEST POINT FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////
//         Exact Closest Point Functions
/////////////////////////////////////////////////////////////////////////////////////////////////

// Closest point to a circle embedded in 2D
void Surface::cpCircle(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    assert(x.size() == cpx.size());
    assert(x.size() == 2);

    Scalar R = m_surface_params[0];
    size_t dim = x.size();

    Scalar theta = atan2(x[1], x[0]);
    cpx[0] = R * cos(theta);
    cpx[1] = R * sin(theta);
    
    dist = 0.0;
    for(size_t d = 0; d < dim; ++d)
    {
        dist += pow(x[d] - cpx[d], 2);
    }
    dist = sqrt(dist);
}


void ComputeOrthogonalBasis(vector<Scalar> &n, vector<Scalar> &t1, vector<Scalar> &t2)
{
    assert(Norm(n) != 0);

    if(abs(n[0]) > abs(n[1]))
    {
        t1[0] = -n[2];
        t1[1] = 0;
        t1[2] = n[0];
    }
    else
    {
        t1[0] = 0;
        t1[1] = n[2];
        t1[2] = -n[1];
    }
    
    t2 = CrossProduct(n, t1);
}


// See https://www.geometrictools.com/Documentation/DistanceToCircle3.pdf
void Surface::cpCircle3D(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    assert(x.size() == cpx.size());

    Scalar R = m_surface_params[0];
    vector<Scalar> centre{m_surface_params[1], m_surface_params[2], m_surface_params[3]};
    vector<Scalar> normal{m_surface_params[4], m_surface_params[5], m_surface_params[6]}; // normal defining plane that circle is in
    
    size_t dim = x.size();

    vector<Scalar> t1(3);
    vector<Scalar> t2(3);
    ComputeOrthogonalBasis(normal, t1, t2);

    // project x onto plane containing circle
    vector<Scalar> centre_to_x = x - centre;
    vector<Scalar> scaled_centre_to_Q = ( DotProduct(normal, normal) * DotProduct(t1, centre_to_x) ) * t1 + DotProduct(t2, centre_to_x) * t2;

    if(Norm(scaled_centre_to_Q) > 0)
    {
        cpx = centre + R * scaled_centre_to_Q / Norm(scaled_centre_to_Q);
        Scalar height = DotProduct(normal, centre_to_x);
        Scalar radial = Norm(CrossProduct(normal, centre_to_x)) - R;
        dist = sqrt(height * height + radial * radial);
        // dist = Distance(x, cpx);
    }
    else // x is on the circle center, so pick any point on the circle as cp
    {
        if(Norm(t1) > Norm(t2))
        {
            cpx = centre + R * t1 / Norm(t1); 
        }
        else
        {
            cpx = centre + R * t2 / Norm(t2); 
        }

        dist = R;
    }
}


// arc specified by portion of the circle between [angle1, angle2]. These angles must be specified in the interval (-pi, pi].
void Surface::cpArc(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy)
{
    Scalar R = m_surface_params[0];
    Scalar angle1 = m_surface_params[1];
    Scalar angle2 = m_surface_params[2];
    size_t dim = x.size();

    if(angle1 < -M_PI)
    {
        throw std::runtime_error("angle1 in Surface::cpArc should be in (-pi, pi].");
    }

    if(angle2 > M_PI)
    {
        throw std::runtime_error("angle2 in Surface::cpArc should be in (-pi, pi].");
    }

    if(angle1 >= angle2)
    {
        throw std::runtime_error("angle1 should be less than angle2 in Surface::cpArc.");
    }


    // find the point opposite the half point in the arc to split the boundary points
    Scalar halfDist = 0.5 * (2 * M_PI - (angle2 - angle1));

    Scalar theta = atan2(x[1], x[0]);

    if(theta > angle1 & theta < angle2)
    {
        cpx[0] = R * cos(theta);
        cpx[1] = R * sin(theta);
        bdy = 0;
    }
    else 
    {
        if(theta <= angle1 & theta >= angle1 - halfDist)
        {
            bdy = 1; // assign point to first boundary
            cpx[0] = R * cos(angle1);
            cpx[1] = R * sin(angle1);
        }
        else if(theta <= angle1 & theta < angle1 - halfDist)
        {
            bdy = 2; // assign point to second boundary
            cpx[0] = R * cos(angle2);
            cpx[1] = R * sin(angle2);
        }
        else if(theta >= angle2 & theta < angle2 + halfDist)
        {
            bdy = 2; // assign point to second boundary
            cpx[0] = R * cos(angle2);
            cpx[1] = R * sin(angle2);
        }
        else if(theta >= angle2 & theta >= angle2 + halfDist)
        {
            bdy = 1; // assign point to first boundary
            cpx[0] = R * cos(angle1);
            cpx[1] = R * sin(angle1);
        }
    }

    if(dim == 3)
    {
        cpx[2] = 0.0;
    }
    
    dist = 0.0;
    for(size_t d = 0; d < dim; ++d)
    {
        dist += pow(x[d] - cpx[d], 2);
    }
    dist = sqrt(dist);
}


void Surface::cpDumbbell(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    vector<Scalar> cpx1 = cpx;
    vector<Scalar> cpx2 = cpx;
    Scalar dist1;
    Scalar dist2;
    cpDumbbellHalf(x, cpx1, dist1, -1.0);
    cpDumbbellHalf(x, cpx2, dist2, 1.0);

    if(dist1 < dist2)
    {
        cpx = cpx1;
        dist = dist1;
    }
    else
    {
        cpx = cpx2;
        dist = dist2;
    }
}


void Surface::cpDumbbellHalf(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, Scalar flip)
{
    VectorX xq(2);
    xq[0] = x[0];
    xq[1] = x[1];

    Scalar a = m_surface_params[0];
    
    auto F = [&a, &flip](const Scalar &t)
    {
        Vector2 x;
        x[0] = t;
        x[1] = flip * (t * t + a) * sqrt(1.0 - t * t);
        return x;
    };

    auto gradF = [&a, &flip](const Scalar &t)
    {
        Vector2 grad;
        grad[0] = 1;
        grad[1] = -flip * t * (a + 3 * t * t - 2) / sqrt(1.0 - t * t);
        return grad;
    };
    
    auto r = [&F](const Scalar &t)
    {
        Vector2 p = F(t);

        return p;
    };

    auto gradtr = [&F, &gradF](const Scalar &t)
    {
        Vector2 gradt = gradF(t);

        return gradt;
    };

    auto objective = [&xq, &r](const VectorX &t)
    {
        return 0.5 * (r(t[0]) - xq).squaredNorm();
    };

    auto gradObjective = [&xq, &r, &gradtr](const VectorX &t)
    {
        VectorX grad(1);
        grad[0] = (r(t[0]) - xq).dot(gradtr(t[0]));

        return grad;
    };

    auto fun = [&objective, &gradObjective](const VectorX &y, VectorX &grad)
    {
        grad = gradObjective(y);
        return objective(y);
    };
    
    // Set up parameters
    LBFGSpp::LBFGSBParam<Scalar> param;
    param.epsilon = 1e-6;
    param.max_iterations = 1000;

    // Create solver and function object
    LBFGSpp::LBFGSBSolver<Scalar> solver(param);

    // Initial guess
    VectorX u = VectorX::Ones(1);

    VectorX lb = VectorX::Constant(1, -0.999);
    VectorX ub = VectorX::Constant(1, 0.999);
    
    // find nearest neighbour in paramSurfaceOfRevoultion point cloud
    Scalar min_dist = numeric_limits<Scalar>::max();
    size_t min_idx;
    for(size_t i = 0; i < m_mesh.vertices.size(); ++i)
    {
        Scalar dist = Distance(x, m_mesh.vertices[i]);
        if(dist < min_dist)
        {
            min_dist = dist;
            min_idx = i;
        }
    }

    u[0] = m_thetap[min_idx];

    // u will be overwritten to be the best point found
    if(u[0] < 1.0 && u[0] > -1.0)
    {
        Scalar fx;
        int niter = solver.minimize(fun, u, fx, lb, ub);
    }

    // std::cout << niter << " iterations" << std::endl;
    // std::cout << "x = \n" << u.transpose() << std::endl;
    // std::cout << "f(x) = " << fx << std::endl;

    VectorX cp = r(u[0]);
    cpx[0] = cp[0];
    cpx[1] = cp[1];

    dist = Distance(x, cpx);
}


void Surface::cpSphere(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    Scalar R = m_surface_params[0];
    vector<Scalar> cen(3, 0.0);
    if(m_surface_params.size() > 1)
    {
        cen = {m_surface_params[1], m_surface_params[2], m_surface_params[3]};
    }
 
    vector<Scalar> x_shifted = x - cen;

    Scalar phi;
    Scalar theta = atan2(x_shifted[1], x_shifted[0]);
    dist = sqrt( pow(x_shifted[0],2) + pow(x_shifted[1],2) + pow(x_shifted[2],2) );
    if(dist != 0)
    {
        phi = acos(x_shifted[2] / dist);
    }
    else
    {
        phi = 0.0; // pick a point on the surface for the (0,0,0) grid point
    }

    cpx[0] = R * sin(phi) * cos(theta);
    cpx[1] = R * sin(phi) * sin(theta);
    cpx[2] = R * cos(phi);
    
    cpx += cen;

    dist -= R;
}



void Surface::cpHemisphere(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy)
{
    Scalar R = m_surface_params[0];
    Scalar phi;

    Scalar theta = atan2(x[1], x[0]);
    if(x[2] > 0) // hemisphere above z = 0 axis
    {
        dist = sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2) );
        if(dist != 0)
        {
            phi = acos(x[2] / dist);
        }
        else
        {
            phi = 0.0; // pick a point on the surface for the (0,0,0) grid point
        }

        cpx[0] = R * sin(phi) * cos(theta);
        cpx[1] = R * sin(phi) * sin(theta);
        cpx[2] = R * cos(phi);
        
        dist -= R;

        bdy = 0;
    }
    else
    {
        cpx[0] = R * cos(theta);
        cpx[1] = R * sin(theta);
        cpx[2] = 0;

        dist = Distance(x, cpx);

        bdy = 1;
    }
}


void Surface::cpBall(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy)
{
    cpSphere(x, cpx, dist);

    if(dist < 0)
    {
        cpx = x;
        dist = 0;
        bdy = 0;
    }
    else
    {
        bdy = 1;
    }
}


void cpInnerTorusCircle(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, Scalar &R, vector<Scalar> &cen)
{
    assert(x.size() == cpx.size());
    
    size_t dim = x.size();

    vector<Scalar> shift_x = x - cen;

    Scalar theta = atan2(shift_x[1], shift_x[0]);
    vector<Scalar> shift_cpx(dim);
    shift_cpx[0] = R * cos(theta);
    shift_cpx[1] = R * sin(theta);
    if(dim == 3)
    {
        shift_cpx[2] = 0.0;
    }

    cpx = shift_cpx + cen;
    
    dist = Distance(x, cpx);
}


void Surface::cpTorus(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    Scalar theta = atan2(x[1], x[0]);
    Scalar r = sqrt(x[0] * x[0] + x[1] * x[1]);

    vector<Scalar> cylx{r, x[2]};
    vector<Scalar> cpcyl(2);
    
    vector<Scalar> cen{m_surface_params[0], 0};
    cpInnerTorusCircle(cylx, cpcyl, dist, m_surface_params[1], cen);

    cpx[0] = cpcyl[0] * cos(theta);
    cpx[1] = cpcyl[0] * sin(theta);
    cpx[2] = cpcyl[1];
}


void Surface::cpPoint(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    for(size_t d = 0; d < m_dim; ++d)
    {
        cpx[d] = m_surface_params[d];
    }

    dist = Distance(x, cpx);
}


void Surface::cpLine(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy)
{
    if(x[0] > m_surface_params[0] && x[0] < m_surface_params[1])
    {
        cpx[0] = x[0];
        bdy = 0;
    }
    else if(x[0] <= m_surface_params[0])
    {
        cpx[0] = m_surface_params[0];
        bdy = 1;
    }
    else if(x[0] >= m_surface_params[1])
    {
        cpx[0] = m_surface_params[1];
        bdy = 2;
    }

    for(size_t d = 1; d < x.size(); ++d)
    {
        cpx[d] = 0.0;
    }

    dist = Distance(x, cpx);
}


void Surface::cpPlane(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy)
{
    assert(x.size() == cpx.size());

    size_t dim = x.size();

    bdy = 1;
    for(size_t d = 0; d < dim; ++d)
    {
        if(x[d] > 1.0)
        {
            cpx[d] = 1.0;
        }
        else if(x[d] < -1.0)
        {
            cpx[d] = -1;
        }
        else
        {
            bdy = 0;
            cpx[d] = x[d];
        }
    }

    dist = Distance(x, cpx);

}


// cpSquare in 2D only
void Surface::cpSquare(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    assert(x.size() == 2);
    assert(x.size() == cpx.size());

    if(x[1] >= 1.0)
    {
        cpx[1] = 1.0;
    }
    else if(x[1] <= -1.0)
    {
        cpx[1] = -1;
    }
    else if(x[1] >= abs(x[0]))
    {
        cpx[1] = 1.0;
    }
    else if(x[1] <= -abs(x[0]))
    {
        cpx[1] = -1.0;
    }
    else
    {
        cpx[1] = x[1];
    }

    if(x[0] >= 1.0)
    {
        cpx[0] = 1.0;
    }
    else if(x[0] <= -1.0)
    {
        cpx[0] = -1;
    }
    else if(x[0] > abs(x[1]))
    {
        cpx[0] = 1.0;
    }
    else if(x[0] < -abs(x[1]))
    {
        cpx[0] = -1.0;
    }
    else
    {
        cpx[0] = x[0];
    }

    dist = Distance(x, cpx);
}


// cp of L shaped domain in 2D only
// from https://math.stackexchange.com/questions/2193720/find-a-point-on-a-line-segment-which-is-the-closest-to-other-point-not-on-the-li
void Surface::cpL(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    assert(x.size() == 2);
    assert(x.size() == cpx.size());

    vector<vector<Scalar>> vertices(6, vector<Scalar>(2));
    vertices[0] = vector<Scalar>({-0.75, 1});
    vertices[1] = vector<Scalar>({-0.75, -1});
    vertices[2] = vector<Scalar>({0.75, -1});
    vertices[3] = vector<Scalar>({0.75, -0.5});
    vertices[4] = vector<Scalar>({-0.25, -0.5});
    vertices[5] = vector<Scalar>({-0.25, 1});
    vector<vector<size_t>> edges(6, vector<size_t>(2));
    edges[0] = vector<size_t>({0, 1});
    edges[1] = vector<size_t>({1, 2});
    edges[2] = vector<size_t>({2, 3});
    edges[3] = vector<size_t>({3, 4});
    edges[4] = vector<size_t>({4, 5});
    edges[5] = vector<size_t>({5, 0});

    vector<vector<Scalar>> cps(6, vector<Scalar>(2));
    vector<Scalar> dists(6);
    for(size_t i = 0; i < 6; ++i)
    {
        vector<Scalar> A = vertices[edges[i][0]];
        vector<Scalar> B = vertices[edges[i][1]];
        vector<Scalar> u = A - x;
        vector<Scalar> v = B - A;
        Scalar vdotv = DotProduct(v,v);
        Scalar t = vdotv == 0.0 ? 0.0 : -DotProduct(u,v) / vdotv;

        if(t >= 0 && t <= 1)
        {
            cps[i] = (1 - t) * A + t * B;
        }
        else if(t < 0)
        {
            cps[i] = A;
        }
        else if(t > 1)
        {
            cps[i] = B;
        }

        dists[i] = Distance(x, cps[i]);
    }

    vector<size_t> sorted_idx = SortedIndices(dists);
    cpx = cps[sorted_idx[0]];
    dist = dists[sorted_idx[0]];
}



void cpLineYAxis(Scalar ystart, Scalar yend, const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    if(ystart < x[1] && x[1] < yend)
    {
        cpx[0] = 0.0;
        cpx[1] = x[1];
        cpx[2] = 0.0;
    }
    else if(x[1] <= ystart)
    {
        cpx[0] = 0.0;
        cpx[1] = ystart;
        cpx[2] = 0.0;
    }
    else if(yend <= x[1])
    {
        cpx[0] = 0.0;
        cpx[1] = yend;
        cpx[2] = 0.0;
    }
    dist = Distance(x, cpx);
}


void cpLineXAxis(Scalar xstart, Scalar xend, const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    if(xstart < x[0] && x[0] < xend)
    {
        cpx[0] = x[0];
        cpx[1] = 0.0;
        cpx[2] = 0.0;
    }
    else if(x[0] <= xstart)
    {
        cpx[0] = xstart;
        cpx[1] = 0.0;
        cpx[2] = 0.0;
    }
    else if(xend <= x[0])
    {
        cpx[0] = xend;
        cpx[1] = 0.0;
        cpx[2] = 0.0;
    }
    dist = Distance(x, cpx);
}



void Surface::cpTorusLineSphere(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    vector<vector<Scalar>> cp(4, vector<Scalar>(3));
    vector<Scalar> d(4);
    vector<Scalar> all_surface_params = m_surface_params;

    // closest point on Torus
    cpTorus(x, cp[0], d[0]);

    // closest point on Ellipsoid
    m_surface_params.clear();
    m_surface_params.push_back(all_surface_params[2]);

    cpSphere(x, cp[1], d[1]);
    
    m_surface_params.clear();
    for(size_t i = 0; i < 3; ++i)
    {
        m_surface_params.push_back(all_surface_params[i]);
    }

    // closest point on Line
    // make line on y axis whose end points touch the torus and ellipsoid
    cpLineYAxis(1.25, 2.0, x, cp[2], d[2]);
    cpLineYAxis(-2.0, -1.25, x, cp[3], d[3]);

    // determine which surface is the closest to x
    vector<size_t> sort_idx = SortedIndices(d);
    cpx = cp[sort_idx[0]];
    dist = d[sort_idx[0]];
}


void Surface::cpSphereLine(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    vector<vector<Scalar>> cp(3, vector<Scalar>(3));
    vector<Scalar> d(3);
    vector<Scalar> all_surface_params = m_surface_params;

    m_surface_params.clear();
    m_surface_params.push_back(all_surface_params[0]);
    m_surface_params.push_back(all_surface_params[1]);
    m_surface_params.push_back(all_surface_params[2]);
    m_surface_params.push_back(all_surface_params[3]);

    // closest point on sphere 1
    cpSphere(x, cp[0], d[0]);

    // closest point on sphere 2
    m_surface_params[0] = all_surface_params[4];
    m_surface_params[1] = all_surface_params[5];
    m_surface_params[2] = all_surface_params[6];
    m_surface_params[3] = all_surface_params[7];
    

    cpSphere(x, cp[1], d[1]);
    
    m_surface_params.clear();
    for(size_t i = 0; i < 8; ++i)
    {
        m_surface_params.push_back(all_surface_params[i]);
    }

    // closest point on Line
    // make line on y axis whose end points touch the torus and ellipsoid
    cpLineXAxis(-0.5, 0.5, x, cp[2], d[2]);

    // determine which surface is the closest to x
    vector<size_t> sort_idx = SortedIndices(d);
    cpx = cp[sort_idx[0]];
    dist = d[sort_idx[0]];
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//          From Parameterizations
/////////////////////////////////////////////////////////////////////////////////////////////////

void Surface::cpTorusKnot(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    Scalar major_radius = m_surface_params[0];
    Scalar p = m_surface_params[1];
    Scalar q = m_surface_params[2];

    auto r = [&major_radius, &p, &q](const VectorX &u)
    {
        Scalar t = u[0];
        VectorX r(3);
        r[0] = (major_radius + cos(q * t)) * cos(p * t);
        r[1] = (major_radius + cos(q * t)) * sin(p * t);
        r[2] = sin(q * t);

        return r;
    };

    auto du = [&major_radius, &p, &q](const VectorX &u)
    {
        Scalar t = u[0];
        VectorX grad(3);
        grad[0] = -p * sin(p * t) * (cos(q * t) + major_radius) - q * cos(p * t) * sin(q * t);
        grad[1] = p * cos(p * t) * (cos(q * t) + major_radius) - q * sin(p * t) * sin(q * t);
        grad[2] = q * cos(q * t);

        return grad;
    };

    cpParametrization(x, cpx, dist, r, du);
}


void Surface::cpPlanarCurve(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    Scalar a = m_surface_params[0];
    Scalar b = m_surface_params[1];
    Scalar c = m_surface_params[2];
    
    auto r = [&a, &b, &c](const VectorX &u)
    {
        Scalar t = u[0];
        Vector2 x;
        x[0] = (cos(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t)) + 0.5 * (a+b)) / (a + b) - c;
        x[1] = (sin(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t)) + 0.5 * (a+b)) / (a + b) - c;
        return x;
    };

    auto du = [&a, &b, &c](const VectorX &u)
    {
        Scalar t = u[0];
        Vector2 grad;
        grad[0] = (cos(t) * (a * cos(a * t) + b * cos(b * t)) - sin(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t))) / (a + b + c);
        grad[1] = (sin(t) * (a * cos(a * t) + b * cos(b * t)) + cos(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t))) / (a + b + c);
        return grad;
    };

    cpParametrization(x, cpx, dist, r, du, 1e-14, 1000);
}


void Surface::cpPlanarCurve2(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    Scalar a = m_surface_params[0];
    Scalar b = m_surface_params[1];
    
    auto r = [&a, &b](const VectorX &u)
    {
        Scalar t = u[0];
        Vector2 x;
        x[0] = sin(t) * (a - b * sin(4 * t));
        x[1] = cos(t);
        return x;
    };

    auto du = [&a, &b](const VectorX &u)
    {
        Scalar t = u[0];
        Vector2 grad;
        grad[0] = a * cos(t) - b * (4 * sin(t) * cos(4 * t) + cos(t) * sin(4 * t));
        grad[1] = -sin(t);
        return grad;
    };

    cpParametrization(x, cpx, dist, r, du, 1e-14, 1000);
}


void Surface::cpSurfaceOfRevolutionCurve(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    Scalar a = m_surface_params[0];
    Scalar b = m_surface_params[1];
    Scalar c = m_surface_params[2];
    
    auto F = [&a, &b, &c](const VectorX &u)
    {
        Scalar t = u[0];
        Vector2 x;
        x[0] = (cos(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t)) + 0.5 * (a+b)) / (a + b) + c;
        x[1] = (sin(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t)) + 0.5 * (a+b)) / (a + b) + c;
        return x;
    };

    auto gradF = [&a, &b, &c](const VectorX &u)
    {
        Scalar t = u[0];
        Vector2 grad;
        grad[0] = (cos(t) * (a * cos(a * t) + b * cos(b * t)) - sin(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t))) / (a + b + c);
        grad[1] = (sin(t) * (a * cos(a * t) + b * cos(b * t)) + cos(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t))) / (a + b + c);
        return grad;
    };
    
    auto G = [](const VectorX &u)
    {
        Vector2 v;
        v[0] = u[0];
        v[1] = 0;
        return v;
    };

    auto gradG = [](const VectorX &u)
    {
        Vector2 grad;
        grad[0] = 1.0;
        grad[1] = 0;
        return grad;
    };

    auto r = [&F, &G](const VectorX &u)
    {
        VectorX p(3);
        Vector2 g = G(u);
        VectorX g0(1); g0[0] = g[0];
        Vector2 f = F(g0);
        p[0] = f[0];
        p[1] = f[1] * cos(g[1]); 
        p[2] = f[1] * sin(g[1]);

        return p;
    };

    auto gradtr = [&F, &G, &gradF, &gradG](const VectorX &u)
    {
        VectorX gradt(3);

        Vector2 g = G(u);
        VectorX g0(1); g0[0] = g[0];

        Vector2 dG = gradG(u);
        Vector2 f = F(g0);
        Vector2 dF = gradF(g0);

        gradt[0] = dF[0] * dG[0];
        gradt[1] = dF[1] * dG[0] * cos(g[1]) - f[1] * sin(g[1]) * dG[1];
        gradt[2] = dF[1] * dG[0] * sin(g[1]) + f[1] * cos(g[1]) * dG[1];

        return gradt;
    };

    cpParametrization(x, cpx, dist, r, gradtr, 1e-6, 100);
}


void Surface::cpSurfaceOfRevolution(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    Scalar a = m_surface_params[0];
    Scalar b = m_surface_params[1];
    Scalar c = m_surface_params[2];
    
    auto F = [&a, &b, &c](const Scalar &t)
    {
        Vector2 x;
        x[0] = (cos(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t)) + 0.5 * (a+b)) / (a + b) + c;
        x[1] = (sin(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t)) + 0.5 * (a+b)) / (a + b) + c;
        return x;
    };

    auto gradF = [&a, &b, &c](const Scalar &t)
    {
        Vector2 grad;
        grad[0] = (cos(t) * (a * cos(a * t) + b * cos(b * t)) - sin(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t))) / (a + b + c);
        grad[1] = (sin(t) * (a * cos(a * t) + b * cos(b * t)) + cos(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t))) / (a + b + c);
        return grad;
    };

    auto r = [&F](const VectorX &u)
    {
        VectorX p(3);
        Vector2 x = F(u[0]);
        p[0] = x[0];
        p[1] = x[1] * cos(u[1]); 
        p[2] = x[1] * sin(u[1]);

        return p;
    };

    auto graduF = [&gradF](const VectorX &u)
    {
        VectorX gradu(3);
        Vector2 dF = gradF(u[0]);
        gradu[0] = dF[0];
        gradu[1] = dF[1] * cos(u[1]);
        gradu[2] = dF[1] * sin(u[1]);

        return gradu;
    };

    auto gradvF = [&F](const VectorX &u)
    {
        VectorX gradv(3);
        Vector2 x = F(u[0]);
        gradv[0] = 0.0;
        gradv[1] = -x[1] * sin(u[1]);
        gradv[2] = x[1] * cos(u[1]);

        return gradv;
    };

    cpParametrization(x, cpx, dist, r, graduF, gradvF);
}


void Surface::cpEllipsoid(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    Scalar a = m_surface_params[0];
    Scalar b = m_surface_params[1];
    Scalar c = m_surface_params[2];

    auto paramEllipsoid = [&a, &b, &c](const VectorX &u)
    {
        VectorX p(3);
        p[0] = a * cos(u[0]) * sin(u[1]);
        p[1] = b * sin(u[0]) * sin(u[1]);
        p[2] = c * cos(u[1]);

        return p;
    };

    auto graduParamEllipsoid = [&a, &b, &c](const VectorX &u)
    {
        VectorX gradu(3);
        gradu[0] = -a * sin(u[0]) * sin(u[1]);
        gradu[1] = b * cos(u[0]) * sin(u[1]);
        gradu[2] = 0.0;

        return gradu;
    };

    auto gradvParamEllipsoid = [&a, &b, &c](const VectorX &u)
    {
        VectorX gradv(3);
        gradv[0] = a * cos(u[0]) * cos(u[1]);
        gradv[1] = b * sin(u[0]) * cos(u[1]);
        gradv[2] = -c * sin(u[1]);

        return gradv;
    };

    cpParametrization(x, cpx, dist, paramEllipsoid, graduParamEllipsoid, gradvParamEllipsoid);
}


void Surface::cpParametrization(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, function<VectorX(const VectorX&)> r, function<VectorX(const VectorX&)> du, Scalar epsilon, size_t max_iter)
{
    VectorX xq(m_dim);
    for(size_t d = 0; d < m_dim; ++d)
    {
        xq[d] = x[d];
    }

    auto objective = [&xq, &r](const VectorX &u)
    {
        return 0.5 * (r(u) - xq).squaredNorm();
    };

    auto gradObjective = [&xq, &r, &du](const VectorX &u)
    {
        VectorX grad(1);
        grad[0] = (r(u) - xq).dot(du(u));

        return grad;
    };

    auto fun = [&objective, &gradObjective](const VectorX &y, VectorX &grad)
    {
        grad = gradObjective(y);
        return objective(y);
    };

    // Set up parameters
    LBFGSpp::LBFGSParam<Scalar> param;
    param.epsilon = epsilon;
    param.max_iterations = max_iter;

    // Create solver and function object
    LBFGSpp::LBFGSSolver<Scalar> solver(param);

    // Initial guess
    VectorX u = VectorX::Ones(1);
    
    // find nearest neighbour in parameterization point cloud
    size_t min_idx = FindNearestParamPoint(x);

    u[0] = m_thetap[min_idx];
    
    // u will be overwritten to be the best point found
    Scalar fx;
    int niter = solver.minimize(fun, u, fx);

    VectorX cp = r(u);
    for(size_t d = 0; d < m_dim; ++d)
    {
        cpx[d] = cp[d];
    }

    dist = Distance(x, cpx);
}

void Surface::cpParametrization(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, function<VectorX(const VectorX&)> r, function<VectorX(const VectorX&)> du, function<VectorX(const VectorX&)> dv, Scalar epsilon, size_t max_iter)
{
    VectorX xq(m_dim);
    for(size_t d = 0; d < m_dim; ++d)
    {
        xq[d] = x[d];
    }

    auto objective = [&xq, &r](const VectorX &u)
    {
        return 0.5 * (r(u) - xq).squaredNorm();
    };

    auto gradObjective = [&xq, &r, &du, &dv](const VectorX &u)
    {
        VectorX grad(2);
        grad[0] = (r(u) - xq).dot(du(u));
        grad[1] = (r(u) - xq).dot(dv(u));

        return grad;
    };

    auto fun = [&objective, &gradObjective](const VectorX &y, VectorX &grad)
    {
        grad = gradObjective(y);
        return objective(y);
    };


    // Set up parameters
    LBFGSpp::LBFGSParam<Scalar> param;
    param.epsilon = epsilon;
    param.max_iterations = max_iter;

    // Create solver and function object
    LBFGSpp::LBFGSSolver<Scalar> solver(param);

    // Initial guess
    VectorX u = VectorX::Ones(2);
    
    // find nearest neighbour in parameterization point cloud
    size_t min_idx = FindNearestParamPoint(x);

    u[0] = m_phip[min_idx];
    u[1] = m_thetap[min_idx];
    
    // u will be overwritten to be the best point found
    Scalar fx;
    int niter = solver.minimize(fun, u, fx);

    VectorX cp = r(u);
    for(size_t d = 0; d < m_dim; ++d)
    {
        cpx[d] = cp[d];
    }

    dist = Distance(x, cpx);
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//          From Level-Sets
/////////////////////////////////////////////////////////////////////////////////////////////////

// see https://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Dziuk88.pdf
void Surface::cpDziuk(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    auto phi = [](const VectorX &p)
    {
        Scalar x = p[0]; Scalar y = p[1]; Scalar z = p[2];
        return pow(x - z * z, 2) + y * y + z * z - 1.0;
    };

    auto gradphi = [](const VectorX &p)
    {
        Scalar x = p[0]; Scalar y = p[1]; Scalar z = p[2];

        VectorX grad(3);
        grad[0] = 2 * (x - z * z);
        grad[1] = 2 * y;
        grad[2] = 2 * z * (1 - 2 * (x - z * z));

        return grad;
    };

    auto Hessianphi = [](const VectorX &p)
    {
        Scalar x = p[0]; Scalar y = p[1]; Scalar z = p[2];

        cpm::Matrix3 H;
        H(0,0) = 2;       H(0,1) = 0;       H(0,2) = -4 * z;
        H(1,0) = H(0,1);  H(1,1) = 2;       H(1,2) = 0;
        H(2,0) = H(0,2);  H(2,1) = H(1,2);  H(2,2) = -4 * x + 12 * z * z + 2;
        
        return H;
    };

    cpLevelSet(x, cpx, dist, phi, gradphi, Hessianphi);
}


// level-set equation from https://math.stackexchange.com/questions/46212/interesting-implicit-surfaces-in-mathbbr3
void Surface::cpDecoTetrahedron(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{
    auto phi = [](const VectorX &p)
    {
        Scalar x = p[0]; Scalar y = p[1]; Scalar z = p[2];
        return pow(x-2, 2) * pow(x+2, 2) + pow(y-2, 2) * pow(y+2, 2) + pow(z-2, 2) * pow(z+2, 2) + 3 * (x*x * y*y + x*x * z*z + y*y * z*z)+ 6 * x * y * z - 10 * (x*x + y*y + z*z) + 22;
    };

    auto gradphi = [](const VectorX &p)
    {
        Scalar x = p[0]; Scalar y = p[1]; Scalar z = p[2];

        VectorX grad(3);
        grad[0] = 4 * pow(x, 3) + 6 * x * (y*y + z*z - 6) + 6 * y * z;
        grad[1] = 4 * pow(y, 3) + 6 * y * (x*x + z*z - 6) + 6 * x * z;
        grad[2] = 4 * pow(z, 3) + 6 * z * (x*x + y*y - 6) + 6 * x * y; 

        return grad;
    };

    auto Hessianphi = [](const VectorX &p)
    {
        Scalar x = p[0]; Scalar y = p[1]; Scalar z = p[2];

        cpm::Matrix3 H;
        H(0,0) = 12 * x * x + 6 * y * y + 6 * z * z - 36;  H(0,1) = 12 * x * y + 6 * z;                       H(0,2) = 12 * x * z + 6 * y;
        H(1,0) = H(0,1);                                   H(1,1) = 6 * x * x + 12 * y * y + 6 * z * z - 36;  H(1,2) = 6 * x + 12 * y * z;
        H(2,0) = H(0,2);                                   H(2,1) = H(1,2);                                   H(2,2) = 6 * x * x + 6 * y * y + 12 * z * z - 36;
        
        return H;
    };

    cpLevelSet(x, cpx, dist, phi, gradphi, Hessianphi);
}


// following ideas of Saye 2014, https://escholarship.org/content/qt49n3d8ph/qt49n3d8ph.pdf
void Surface::cpLevelSet(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, function<Scalar(const VectorX&)> phi, function<VectorX(const VectorX&)> gradphi, function<cpm::Matrix3(const VectorX&)> Hessianphi)
{
    VectorX xq(3);
    xq[0] = x[0]; xq[1] = x[1]; xq[2] = x[2];

    auto f = [&xq, &phi](const VectorX &y)
    {
        VectorX p(3);
        p[0] = y[0];
        p[1] = y[1];
        p[2] = y[2];
        return 0.5 * (p - xq).squaredNorm() + y[3] * phi(p);
    };

    auto gradf = [&xq, &phi, &gradphi](const VectorX &y)
    {
        VectorX p(3);
        p[0] = y[0];
        p[1] = y[1];
        p[2] = y[2];
        VectorX gradp = p - xq + y[3] * gradphi(p);

        VectorX grad(4);
        grad[0] = gradp[0];
        grad[1] = gradp[1];
        grad[2] = gradp[2];
        grad[3] = phi(p);
        return grad;
    };

    auto Hessianf = [&xq, &phi, &gradphi, &Hessianphi](const VectorX &y)
    {
        VectorX p(3);
        p[0] = y[0];
        p[1] = y[1];
        p[2] = y[2];
        cpm::Matrix3 Hphi = Hessianphi(p);
        cpm::Matrix3 I = cpm::Matrix3::Identity();
        cpm::Matrix3 H3 = I + y[3] * Hphi;

        cpm::Matrix4 H;
        for(size_t r = 0; r < 3; ++r)
        {
            for(size_t c = 0; c < 3; ++c)
            {
                H(r,c) = H3(r,c);
            }
        }

        VectorX gradp = gradphi(p);
        H(0, 3) = gradp[0];
        H(1, 3) = gradp[1];
        H(2, 3) = gradp[2];
        H(3, 0) = gradp[0];
        H(3, 1) = gradp[1];
        H(3, 2) = gradp[2];
        H(3,3) = 0;
        
        return H;
    };

    // initial guess
    VectorX perturb = xq;
    if(m_mesh.vertices.size() > 0) // get initial guess from points on the level-set surface, found from triangulation constructed in TriDziuk.cpp
    {
        Scalar min_dist = numeric_limits<Scalar>::max();
        size_t min_idx;
        for(size_t i = 0; i < m_mesh.vertices.size(); ++i)
        {
            Scalar dist = Distance(x, m_mesh.vertices[i]);
            if(dist < min_dist)
            {
                min_dist = dist;
                min_idx = i;
            }
        }

        for(size_t d = 0; d < 3; ++d)
        {
            perturb[d] = m_mesh.vertices[min_idx][d];
        }
    }
    else
    {
        if(gradphi(xq).squaredNorm() < 1e-6) // perturb xq a little
        {
            for(size_t d = 0; d < 3; ++d)
            {
                perturb[d] += 0.001;
            }
        }
    }

    size_t max_iter = 1000;
    Scalar tol = 1e-10;
    Scalar error = numeric_limits<Scalar>::max();
    VectorX xold = perturb;
    VectorX xnew;
    size_t k = 0;
    while(error > tol && k < max_iter)
    {
        xnew = xold - phi(xold) * gradphi(xold) / gradphi(xold).squaredNorm();
        error = (xnew - xold).norm();
        xold = xnew;
        ++k;
    }

    if(k == max_iter - 1)
    {
        cout << "cpLevelSet initial iteration max exceeded" << endl;
    }

    // Newton iteration
    error = numeric_limits<Scalar>::max();
    VectorX yold(4);
    yold[0] = xnew[0];
    yold[1] = xnew[1];
    yold[2] = xnew[2];
    yold[3] = (xq - xnew).dot(gradphi(xnew) / gradphi(xnew).squaredNorm());

    VectorX ynew;
    k = 0;
    while(error > tol && k < max_iter)
    {
        VectorX Hinv_gradf = Hessianf(yold).inverse() * gradf(yold);
        ynew = yold - Hinv_gradf; 
        error = (ynew - yold).norm();
        yold = ynew;
        ++k;
    }

    if(k == max_iter - 1)
    {
        cout << "cpLevelSet Newton iteration max exceeded" << endl;
    }

    cpx[0] = ynew[0];
    cpx[1] = ynew[1];
    cpx[2] = ynew[2];

    dist = Distance(x, cpx);
}


int Surface::cpTri(const vector<Scalar> &x)
{
    vector<Scalar> cpx(x.size());
    Scalar dist;
    size_t bdy;
    return cpTri(x, cpx, dist, bdy);
}


int Surface::cpTri(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy)
{
    // perform a closest point query
    fcpw::Interaction<3> interaction;
    fcpw::Vector3 queryPoint;

    queryPoint[0] = x[0];
    queryPoint[1] = x[1];
    queryPoint[2] = x[2];

    m_scene.findClosestPoint(queryPoint, interaction);

    cpx[0] = interaction.p[0];
    cpx[1] = interaction.p[1];
    cpx[2] = interaction.p[2];

    dist = interaction.signedDistance(queryPoint);
    dist = abs(dist);

    // TO DO: add in how to find boundary for triangulated surface
    bdy = 0;

    if(SurfaceType().compare("Polyline") == 0)
    {
        bdy = BoundaryFromEndpointDistance(cpx);
    }

    return interaction.primitiveIndex;
}


void Surface::cpNNPointCloud(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist)
{

    int face = cpTri(x);

    vector<Scalar> v0 = m_mesh.vertices[m_mesh.faces[face][0]];
    vector<Scalar> v1 = m_mesh.vertices[m_mesh.faces[face][1]];
    vector<Scalar> v2 = m_mesh.vertices[m_mesh.faces[face][2]];

    Scalar d0 = Distance(v0, x);
    Scalar d1 = Distance(v1, x);
    Scalar d2 = Distance(v2, x);

    if(d0 < d1)
    {
        if(d0 < d2)
        {
            cpx = v0;
            dist = d0;
        }
        else
        {
            cpx = v2;
            dist = d2;
        }
    }
    else
    {
        if(d1 < d2)
        {
            cpx = v1;
            dist = d1;
        }
        else
        {
            cpx = v2;
            dist = d2;
        }
    }
}


// compute cp to each edge of polyline, from https://math.stackexchange.com/questions/2193720/find-a-point-on-a-line-segment-which-is-the-closest-to-other-point-not-on-the-li
void Surface::cpPolyline2D(const vector<Scalar> &x, vector<Scalar> &cpx, Scalar &dist, size_t &bdy)
{
    assert(x.size() == 2);
    assert(x.size() == cpx.size());

    vector<vector<Scalar>> cps(m_mesh.faces.size(), vector<Scalar>(2));
    vector<Scalar> dists(m_mesh.faces.size());
    for(size_t i = 0; i < m_mesh.faces.size(); ++i)
    {
        vector<Scalar> A = m_mesh.vertices[m_mesh.faces[i][0]];
        vector<Scalar> B = m_mesh.vertices[m_mesh.faces[i][1]];
        vector<Scalar> u = A - x;
        vector<Scalar> v = B - A;
        Scalar vdotv = DotProduct(v,v);
        Scalar t = vdotv == 0.0 ? 0.0 : -DotProduct(u,v) / vdotv;

        if(t >= 0 && t <= 1)
        {
            cps[i] = (1 - t) * A + t * B;
        }
        else if(t < 0)
        {
            cps[i] = A;
        }
        else if(t > 1)
        {
            cps[i] = B;
        }

        dists[i] = Distance(x, cps[i]);
    }

    vector<size_t> sorted_idx = SortedIndices(dists);
    cpx = cps[sorted_idx[0]];
    bdy = BoundaryFromEndpointDistance(cpx);
    dist = dists[sorted_idx[0]];
}


///////////////////////////////////////////////////////////////////////////////////
//  SURFACE PARAMETRIZATIONS (for visualization - final solution on nicely spaced points)
///////////////////////////////////////////////////////////////////////////////////

void Surface::paramCircle(size_t &Np)
{
    Scalar R = m_surface_params[0];
    Scalar dtheta = 2 * M_PI / (Np - 1);

    for(size_t i = 0; i < Np; ++i)
    {
        m_thetap[i] = i * dtheta;
        m_mesh.vertices[i][0] = R * cos(m_thetap[i]);
        m_mesh.vertices[i][1] = R * sin(m_thetap[i]);
    }
}


// See https://math.stackexchange.com/questions/511370/how-to-rotate-one-vector-about-another
void Surface::paramCircle3D(size_t &Np)
{
    Scalar R = m_surface_params[0];
    vector<Scalar> centre{m_surface_params[1], m_surface_params[2], m_surface_params[3]};
    vector<Scalar> normal{m_surface_params[4], m_surface_params[5], m_surface_params[6]};

    Scalar dtheta = 2 * M_PI / (Np - 1);

    // find a vector perpendicular to the normal to start from
    vector<Scalar> t1(3), t2(3);
    ComputeOrthogonalBasis(normal, t1, t2);
    Normalize(t1, 0.0);
    Normalize(t2, 0.0);

    for(size_t i = 0; i < Np; ++i)
    {
        m_thetap[i] = i * dtheta;
        Scalar x1 = cos(m_thetap[i]);
        Scalar x2 = sin(m_thetap[i]);
        m_mesh.vertices[i] = centre + R * (x1 * t1 + x2 * t2);
    }
}


void Surface::paramArc(size_t &Np)
{
    // TO DO: make this depend on the embedding dimension
    Scalar R = m_surface_params[0];
    Scalar angle1 = m_surface_params[1];
    Scalar angle2 = m_surface_params[2];
    Scalar dtheta = (angle2 - angle1) / (Np - 1);

    for(size_t i = 0; i < Np; ++i)
    {
        m_thetap[i] = angle1 + i * dtheta;
        m_mesh.vertices[i][0] = R * cos(m_thetap[i]);
        m_mesh.vertices[i][1] = R * sin(m_thetap[i]);
        if(m_dim == 3)
        {
            m_mesh.vertices[i][2] = 0.0;
        }
    }
}


void Surface::paramDumbbell(size_t &Np)
{    
    assert(Np % 2 == 0);

    size_t halfNp = (size_t) (0.5 * Np);

    Scalar du = 2.0 / (Scalar)(halfNp - 1);

    // here phi = u and theta = v
    for(size_t i = 0; i < halfNp; ++i)
    {
        m_thetap[i] = i * du - 1.0;
    }

    for(size_t i = halfNp; i < Np; ++i)
    {
        m_thetap[i] = 1.0 - (i - halfNp) * du;
    }

    Scalar a = m_surface_params[0];
    auto F = [&a](const Scalar &t)
    {
        Vector2 x;
        x[0] = t;
        x[1] = (t * t + a) * sqrt(1.0 - t * t);
        return x;
    };

    for(size_t i = 0; i < halfNp; ++i)
    {
        Vector2 x = F(m_thetap[i]);
        m_mesh.vertices[i][0] = x[0];
        m_mesh.vertices[i][1] = x[1]; 
    }

    for(size_t i = halfNp; i < Np; ++i)
    {
        Vector2 x = F(m_thetap[i]);
        m_mesh.vertices[i][0] = x[0];
        m_mesh.vertices[i][1] = -x[1]; 
    }
}


// Note: the number of parametrization points, Np, must be the square of some number here
void Surface::paramSphere(size_t &Np)
{
    if(floor(sqrt(Np)) != sqrt(Np))
    {
        throw std::runtime_error("Np must be the square of some number for Surface::paramSphere");
    }
    size_t sqrtNp = sqrt(Np);
    Scalar R = m_surface_params[0];
    Scalar dtheta = 2 * M_PI / (sqrtNp - 1);
    Scalar dphi = M_PI / (sqrtNp - 1);

    for(size_t i = 0; i < sqrtNp; ++i)
    {
        for(size_t j = 0; j < sqrtNp; ++j)
        {
            m_thetap[i * sqrtNp + j] = i * dtheta;
            m_phip[i * sqrtNp + j] = j * dphi;
        }
    }

    for(size_t i = 0; i < Np; ++i)
    {
        m_mesh.vertices[i][0] = R * sin(m_phip[i]) * cos(m_thetap[i]);
        m_mesh.vertices[i][1] = R * sin(m_phip[i]) * sin(m_thetap[i]);
        m_mesh.vertices[i][2] = R * cos(m_phip[i]);
    }
}


void Surface::paramHemisphere(size_t &Np)
{
    if(floor(sqrt(Np)) != sqrt(Np))
    {
        throw std::runtime_error("Np must be the square of some number for Surface::paramHemisphere");
    }
    size_t sqrtNp = sqrt(Np);
    Scalar R = m_surface_params[0];
    Scalar dtheta = 2 * M_PI / (sqrtNp - 1);
    Scalar dphi = M_PI_2 / (sqrtNp - 1);

    for(size_t i = 0; i < sqrtNp; ++i)
    {
        for(size_t j = 0; j < sqrtNp; ++j)
        {
            m_thetap[i * sqrtNp + j] = i * dtheta;
            m_phip[i * sqrtNp + j] = j * dphi;
        }
    }

    for(size_t i = 0; i < Np; ++i)
    {
        m_mesh.vertices[i][0] = R * sin(m_phip[i]) * cos(m_thetap[i]);
        m_mesh.vertices[i][1] = R * sin(m_phip[i]) * sin(m_thetap[i]);
        m_mesh.vertices[i][2] = R * cos(m_phip[i]);
    }
}


void Surface::paramBall(size_t &Np)
{
    Scalar R = m_surface_params[0];
    Scalar dx = 2 * R / (Np - 1);

    for(size_t i = 0; i < Np; ++i)
    {
        for(size_t j = 0; j < Np; ++j)
        {
            for(size_t k = 0; k < Np; ++k)
            {
                vector<Scalar> xp(m_dim);
                vector<Scalar> cpxp(m_dim);
                Scalar dist;
                size_t bdy;
                xp[0] = -R + i * dx;
                xp[1] = -R + j * dx;
                xp[2] = -R + k * dx;

                cpBall(xp, cpxp, dist, bdy);

                if(bdy == 0)
                {
                    m_mesh.vertices.push_back(xp);
                }
            }
        }
    }
}


// using https://mathworld.wolfram.com/Torus.html parameterization
void Surface::paramTorus(size_t &Np)
{
    if(floor(sqrt(Np)) != sqrt(Np))
    {
        throw std::runtime_error("Np must be the square of some number for Surface::paramTorus");
    }
    size_t sqrtNp = sqrt(Np);

    Scalar R = m_surface_params[0];
    Scalar r = m_surface_params[1];   
    
    Scalar du = 2 * M_PI / (sqrtNp - 1);

    // here phi = u and theta = v
    for(size_t i = 0; i < sqrtNp; ++i)
    {
        for(size_t j = 0; j < sqrtNp; ++j)
        {
            m_phip[i * sqrtNp + j] = i * du;
            m_thetap[i * sqrtNp + j] = j * du;
        }
    }

    for(size_t i = 0; i < Np; ++i)
    {
        m_mesh.vertices[i][0] = (R + r * cos(m_thetap[i])) * cos(m_phip[i]);
        m_mesh.vertices[i][1] = (R + r * cos(m_thetap[i])) * sin(m_phip[i]);
        m_mesh.vertices[i][2] = r * sin(m_thetap[i]);
    }
}


void Surface::paramLine(size_t &Np)
{
    Scalar start = m_surface_params[0];
    Scalar end = m_surface_params[1];
    Scalar dx = (end - start) / (Np - 1);

    for(size_t i = 0; i < Np; ++i)
    {
        m_thetap[i] = start + i * dx;
        m_mesh.vertices[i][0] = m_thetap[i];
        m_mesh.vertices[i][1] = 0.0;
    }
}


void Surface::paramPlane(size_t &Np)
{
    if(floor(sqrt(Np)) != sqrt(Np))
    {
        throw std::runtime_error("Np must be the square of some number for Surface::paramPlane");
    }
    size_t sqrtNp = sqrt(Np);

    Scalar du = 2.0 / (sqrtNp - 1);

    // here phi = u and theta = v
    for(size_t i = 0; i < sqrtNp; ++i)
    {
        for(size_t j = 0; j < sqrtNp; ++j)
        {
            m_phip[i * sqrtNp + j] = i * du;
            m_thetap[i * sqrtNp + j] = j * du;
        }
    }

    for(size_t i = 0; i < Np; ++i)
    {
        m_mesh.vertices[i][0] = -1.0 + m_phip[i];
        m_mesh.vertices[i][1] = -1.0 + m_thetap[i];
    }
}


void Surface::paramTorusLineSphere()
{
    size_t Np = 250000;
    size_t Nline = 1000;
    if(floor(sqrt(Np)) != sqrt(Np))
    {
        throw std::runtime_error("Np must be the square of some number for Surface::paramTorusLineSphere");
    }
    size_t sqrtNp = sqrt(Np);

    m_mesh.vertices.resize(2 * Np + 2 * Nline, vector<Scalar>(3));

    // param torus
    Scalar R = m_surface_params[0];
    Scalar r = m_surface_params[1];   
    
    Scalar du = 2 * M_PI / (sqrtNp - 1);

    // here phi = u and theta = v
    vector<Scalar> thetap(Np);
    vector<Scalar> phip(Np);
    for(size_t i = 0; i < sqrtNp; ++i)
    {
        for(size_t j = 0; j < sqrtNp; ++j)
        {
            phip[i * sqrtNp + j] = i * du;
            thetap[i * sqrtNp + j] = j * du;
        }
    }

    for(size_t i = 0; i < Np; ++i)
    {
        m_mesh.vertices[i][0] = (R + r * cos(thetap[i])) * cos(phip[i]);
        m_mesh.vertices[i][1] = (R + r * cos(thetap[i])) * sin(phip[i]);
        m_mesh.vertices[i][2] = r * sin(thetap[i]);
    }

    // param sphere
    Scalar dtheta = 2 * M_PI / (sqrtNp - 1);
    Scalar dphi = M_PI / (sqrtNp - 1);

    for(size_t i = 0; i < sqrtNp; ++i)
    {
        for(size_t j = 0; j < sqrtNp; ++j)
        {
            thetap[i * sqrtNp + j] = i * dtheta;
            phip[i * sqrtNp + j] = j * dphi;
        }
    }

    R = m_surface_params[2];
    for(size_t i = 0; i < Np; ++i)
    {
        m_mesh.vertices[i + Np][0] = R * sin(phip[i]) * cos(thetap[i]);
        m_mesh.vertices[i + Np][1] = R * sin(phip[i]) * sin(thetap[i]);
        m_mesh.vertices[i + Np][2] = R * cos(phip[i]);
    }

    // param lines
    Scalar dl = 0.75 / (Nline - 1);
    for(size_t i = 0; i < Nline; ++i)
    {
        m_mesh.vertices[i + 2 * Np][0] = 0.0;
        m_mesh.vertices[i + 2 * Np][1] = 1.25 + i * dl;
        m_mesh.vertices[i + 2 * Np][2] = 0.0;

        m_mesh.vertices[i + 2 * Np + Nline][0] = 0.0;
        m_mesh.vertices[i + 2 * Np + Nline][1] = -2.0 + i * dl;
        m_mesh.vertices[i + 2 * Np + Nline][2] = 0.0;
    }

    // faces for making surface meshes and curve networks for visualization
    for(size_t i = 0; i < sqrtNp - 1; ++i)
    {
        for(size_t j = 0; j < sqrtNp - 1; ++j)
        {
            // for torus
            vector<size_t> face(4);
            face[0] = i * sqrtNp + j;
            face[1] = (i + 1) * sqrtNp + j;
            face[2] = (i + 1) * sqrtNp + j + 1;
            face[3] = i * sqrtNp + j + 1;
            m_mesh.faces.push_back(face);

            for(size_t v = 0; v < 4; ++v)
            {
                face[v] += Np;
            }
            m_mesh.faces.push_back(face);
        }
    }
}


void Surface::paramSphereLine()
{
    size_t Np = 250000;
    size_t Nline = 1000;
    if(floor(sqrt(Np)) != sqrt(Np))
    {
        throw std::runtime_error("Np must be the square of some number for Surface::paramTorusLineSphere");
    }
    size_t sqrtNp = sqrt(Np);

    m_mesh.vertices.resize(2 * Np + Nline, vector<Scalar>(3));

    // param sphere
    Scalar dtheta = 2 * M_PI / (sqrtNp - 1);
    Scalar dphi = M_PI / (sqrtNp - 1);

    vector<Scalar> thetap(Np);
    vector<Scalar> phip(Np);
    for(size_t i = 0; i < sqrtNp; ++i)
    {
        for(size_t j = 0; j < sqrtNp; ++j)
        {
            thetap[i * sqrtNp + j] = i * dtheta;
            phip[i * sqrtNp + j] = j * dphi;
        }
    }

    // param sphere 1
    Scalar R = m_surface_params[0];
    vector<Scalar> cen = {m_surface_params[1], m_surface_params[2], m_surface_params[3]};
    for(size_t i = 0; i < Np; ++i)
    {
        m_mesh.vertices[i][0] = R * sin(phip[i]) * cos(thetap[i]);
        m_mesh.vertices[i][1] = R * sin(phip[i]) * sin(thetap[i]);
        m_mesh.vertices[i][2] = R * cos(phip[i]);

        m_mesh.vertices[i] += cen;
    }

    // param sphere 2
    R = m_surface_params[4];
    cen = {m_surface_params[5], m_surface_params[6], m_surface_params[7]};
    for(size_t i = 0; i < Np; ++i)
    {
        m_mesh.vertices[i + Np][0] = R * sin(phip[i]) * cos(thetap[i]);
        m_mesh.vertices[i + Np][1] = R * sin(phip[i]) * sin(thetap[i]);
        m_mesh.vertices[i + Np][2] = R * cos(phip[i]);

        m_mesh.vertices[i + Np] += cen;
    }

    // param line
    Scalar dl = 1.0 / (Nline - 1);
    for(size_t i = 0; i < Nline; ++i)
    {
        m_mesh.vertices[i + 2 * Np][0] = -0.5 + i * dl;
        m_mesh.vertices[i + 2 * Np][1] = 0.0;
        m_mesh.vertices[i + 2 * Np][2] = 0.0;
    }

    // faces for making surface mesh for visualization
    for(size_t i = 0; i < sqrtNp - 1; ++i)
    {
        for(size_t j = 0; j < sqrtNp - 1; ++j)
        {
            // for sphere 1
            vector<size_t> face(4);
            face[0] = i * sqrtNp + j;
            face[1] = (i + 1) * sqrtNp + j;
            face[2] = (i + 1) * sqrtNp + j + 1;
            face[3] = i * sqrtNp + j + 1;
            m_mesh.faces.push_back(face);

            // for sphere 2
            for(size_t v = 0; v < 4; ++v)
            {
                face[v] += Np;
            }
            m_mesh.faces.push_back(face);
        }
    }
}


// using the one that wraps around a torus https://en.wikipedia.org/wiki/Trefoil_knot
// more general one from https://en.wikipedia.org/wiki/Torus_knot
void Surface::paramTorusKnot(size_t &Np)
{    
    Scalar dt = 2 * M_PI / (Np - 1);

    Scalar major_radius = m_surface_params[0];
    Scalar p = m_surface_params[1];
    Scalar q = m_surface_params[2];

    // here theta = t
    for(size_t i = 0; i < Np; ++i)
    {
        m_thetap[i] = i * dt;
        m_mesh.vertices[i][0] = (major_radius + cos(q * m_thetap[i])) * cos(p * m_thetap[i]);
        m_mesh.vertices[i][1] = (major_radius + cos(q * m_thetap[i])) * sin(p * m_thetap[i]);
        m_mesh.vertices[i][2] = sin(q * m_thetap[i]);
    }
}


void Surface::paramPlanarCurve(size_t &Np)
{    
    Scalar du = 2 * M_PI / (Np - 1);

    // here phi = u and theta = v
    for(size_t i = 0; i < Np; ++i)
    {
        m_thetap[i] = i * du;
    }

    Scalar a = m_surface_params[0];
    Scalar b = m_surface_params[1];
    Scalar c = m_surface_params[2];
    auto F = [&a, &b, &c](const Scalar &t)
    {
        Vector2 x;
        x[0] = (cos(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t)) + 0.5 * (a+b)) / (a + b) - c;
        x[1] = (sin(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t)) + 0.5 * (a+b)) / (a + b) - c;
        return x;
    };

    for(size_t i = 0; i < Np; ++i)
    {
        Vector2 x = F(m_thetap[i]);
        m_mesh.vertices[i][0] = x[0];
        m_mesh.vertices[i][1] = x[1]; 
    }
}


void Surface::paramPlanarCurve2(size_t &Np)
{    
    Scalar du = 2 * M_PI / (Np - 1);

    // here phi = u and theta = v
    for(size_t i = 0; i < Np; ++i)
    {
        m_thetap[i] = i * du;
    }

    Scalar a = m_surface_params[0];
    Scalar b = m_surface_params[1];
    auto F = [&a, &b](const Scalar &t)
    {
        Vector2 x;
        x[0] = sin(t) * (a - b * sin(4 * t));
        x[1] = cos(t);
        return x;
    };

    for(size_t i = 0; i < Np; ++i)
    {
        Vector2 x = F(m_thetap[i]);
        m_mesh.vertices[i][0] = x[0];
        m_mesh.vertices[i][1] = x[1]; 
    }
}


void Surface::paramSurfaceOfRevolutionCurve(size_t &Np)
{    
    Scalar du = 2 * M_PI / (Np - 1);

    // here phi = u and theta = v
    for(size_t i = 0; i < Np; ++i)
    {
        m_thetap[i] = i * du;
    }

    Scalar a = m_surface_params[0];
    Scalar b = m_surface_params[1];
    Scalar c = m_surface_params[2];
    auto F = [&a, &b, &c](const Scalar &t)
    {
        Vector2 x;
        x[0] = (cos(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t)) + 0.5 * (a+b)) / (a + b) + c;
        x[1] = (sin(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t)) + 0.5 * (a+b)) / (a + b) + c;
        return x;
    };

    auto G = [&a, &b, &c](const Scalar &t)
    {
        Vector2 u;
        u[0] = t;
        u[1] = 0.0;
        return u;
    };

    for(size_t i = 0; i < Np; ++i)
    {
        Vector2 u = G(m_thetap[i]);
        Vector2 x = F(u[0]);
        m_mesh.vertices[i][0] = x[0];
        m_mesh.vertices[i][1] = x[1] * cos(u[1]); 
        m_mesh.vertices[i][2] = x[1] * sin(u[1]);
    }
}


// see http://sambrunacini.com/parametric-flowers/ for ideas of how to make the 2D parametric shape to revolve
void Surface::paramSurfaceOfRevolution(size_t &Np)
{
    if(floor(sqrt(Np)) != sqrt(Np))
    {
        throw std::runtime_error("Np must be the square of some number for Surface::SurfaceOfRevolution");
    }
    size_t sqrtNp = sqrt(Np);
    
    Scalar du = 2 * M_PI / (sqrtNp - 1);
    Scalar dv = 2 * M_PI / (sqrtNp - 1);;

    // here phi = u and theta = v
    for(size_t i = 0; i < sqrtNp; ++i)
    {
        for(size_t j = 0; j < sqrtNp; ++j)
        {
            m_phip[i * sqrtNp + j] = i * du;
            m_thetap[i * sqrtNp + j] = j * dv;
        }
    }

    Scalar a = m_surface_params[0];
    Scalar b = m_surface_params[1];
    Scalar c = m_surface_params[2];
    auto F = [&a, &b, &c](const Scalar &t)
    {
        Vector2 x;
        x[0] = (cos(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t)) + 0.5 * (a+b)) / (a + b) + c;
        x[1] = (sin(t) * (0.5 * (a + b) + sin(a * t) + sin(b * t)) + 0.5 * (a+b)) / (a + b) + c;
        return x;
    };

    for(size_t i = 0; i < Np; ++i)
    {
        Vector2 x = F(m_phip[i]);
        m_mesh.vertices[i][0] = x[0];
        m_mesh.vertices[i][1] = x[1] * cos(m_thetap[i]); 
        m_mesh.vertices[i][2] = x[1] * sin(m_thetap[i]);
    }
}


// using https://mathworld.wolfram.com/Ellipsoid.html parameterization
void Surface::paramEllipsoid(size_t &Np)
{
    if(floor(sqrt(Np)) != sqrt(Np))
    {
        throw std::runtime_error("Np must be the square of some number for Surface::Ellipsoid");
    }
    size_t sqrtNp = sqrt(Np);

    Scalar a = m_surface_params[0];
    Scalar b = m_surface_params[1];   
    Scalar c = m_surface_params[2];   
    
    Scalar du = 2 * M_PI / (sqrtNp - 1);
    Scalar dv = M_PI / (sqrtNp - 1);

    // here phi = u and theta = v
    for(size_t i = 0; i < sqrtNp; ++i)
    {
        for(size_t j = 0; j < sqrtNp; ++j)
        {
            m_phip[i * sqrtNp + j] = i * du;
            m_thetap[i * sqrtNp + j] = j * dv;
        }
    }

    for(size_t i = 0; i < Np; ++i)
    {
        m_mesh.vertices[i][0] = a * sin(m_thetap[i]) * cos(m_phip[i]);
        m_mesh.vertices[i][1] = b * sin(m_thetap[i]) * sin(m_phip[i]);
        m_mesh.vertices[i][2] = c * cos(m_thetap[i]);
    }
}


void Surface::paramDziuk()
{
    SimplePolygonMesh mesh;
    mesh.readMeshFromFile("../assets/Dziuk_mid_res.obj");

    m_mesh.vertices = mesh.vertices;
    m_mesh.faces = mesh.faces;
}


void Surface::paramDecoTetrahedron()
{
    SimplePolygonMesh mesh;
    mesh.readMeshFromFile("../assets/DecoTetrahedron_vertices.obj");

    m_mesh.vertices = mesh.vertices;
}


///////////////////////////////////////////////////////////////////////////////////
//  Helpers
///////////////////////////////////////////////////////////////////////////////////

void Surface::InitializeFCPW()
{
    size_t nVertices = m_mesh.vertices.size();

    if(SurfaceType().compare("Triangulation") == 0 || SurfaceType().compare("NN Point Cloud") == 0)
    {
        // set the types of primitives the objects in the scene contain;
        // in this case, we have a single object consisting of only triangles
        m_scene.setObjectTypes({{fcpw::PrimitiveType::Triangle}});

        // set the vertex and triangle count of the (0th) object
        size_t nTriangles = m_mesh.faces.size();
        m_scene.setObjectVertexCount(nVertices, 0);
        m_scene.setObjectTriangleCount(nTriangles, 0);

        // specify the triangle indices
        int face[3];
        for (size_t i = 0; i < nTriangles; i++) {
            face[0] = m_mesh.faces[i][0];
            face[1] = m_mesh.faces[i][1];
            face[2] = m_mesh.faces[i][2];
            m_scene.setObjectTriangle(face, i, 0);
        }
    }
    else if(SurfaceType().compare("Polyline") == 0)
    {
        // set the types of primitives the objects in the scene contain;
        // in this case, we have a single object consisting of only line segments
        m_scene.setObjectTypes({{fcpw::PrimitiveType::LineSegment}});

        // set the vertex and triangle count of the (0th) object
        size_t nEdges = m_mesh.faces.size();
        m_scene.setObjectVertexCount(nVertices, 0);
        m_scene.setObjectLineSegmentCount(nEdges, 0);
        
        // specify the edge indices
        int face[2];
        for (size_t i = 0; i < nEdges; i++) {
            face[0] = m_mesh.faces[i][0];
            face[1] = m_mesh.faces[i][1];
            m_scene.setObjectLineSegment(face, i, 0);
        }
    }

    // specify the vertex positions
    fcpw::Vector3 position;
    for (size_t i = 0; i < nVertices; i++) {
        position[0] = m_mesh.vertices[i][0];
        position[1] = m_mesh.vertices[i][1];
        position[2] = m_mesh.vertices[i][2];
        m_scene.setObjectVertex(position, i, 0);
    }

    // now that the geometry has been specified, build the acceleration structure
    m_scene.build(fcpw::AggregateType::Bvh_SurfaceArea, true); // the second boolean argument enables vectorization
}


size_t Surface::FindNearestParamPoint(const vector<Scalar> &x)
{
    Scalar min_dist = numeric_limits<Scalar>::max();
    size_t min_idx = 0;
    for(size_t i = 0; i < m_mesh.vertices.size(); ++i)
    {
        Scalar dist = Distance(x, m_mesh.vertices[i]);
        if(dist < min_dist)
        {
            min_dist = dist;
            min_idx = i;
        }
    }

    return min_idx;
}


// boundary tag based on cpx distance to endpoint of parameterization or polyline
size_t Surface::BoundaryFromEndpointDistance(const vector<Scalar> &cpx)
{
    size_t bdy = 0;

    if(isOpen())
    {
        Scalar tol = 1e-7; // seems the smallest I can make it for triangulated surfaces, probably due to using single precision in fcpw
        
        vector<Scalar> endpoint1 = xp()[0];
        vector<Scalar> endpoint2 = xp()[xp().size() - 1];
        
        if( Distance(endpoint1, cpx) < tol )
        {
            bdy = 1;
        }
        else if( Distance(endpoint2, cpx) < tol )
        {
            bdy = 2;
        }
    }

    return bdy;
}


vector<Scalar> Surface::IBCTangent(const vector<Scalar> &x)
{
    Scalar min_dist;
    vector<size_t> nearest_indices = IBCNearestTwoXpIndices(x, min_dist);
    
    vector<Scalar> tangent;
    if(nearest_indices[0] < nearest_indices[1])
    {
        tangent = xp()[nearest_indices[1]] - xp()[nearest_indices[0]];
    }
    else
    {
        tangent = xp()[nearest_indices[0]] - xp()[nearest_indices[1]];
    }

    Normalize(tangent, 0.0);

    return tangent;
}


// 1D parameter along a ibc curve
vector<size_t> Surface::IBCNearestTwoXpIndices(const vector<Scalar>& x, Scalar &min_dist)
{
    vector<Scalar> cpx(x.size());
    Scalar dist;
    size_t bdy;
    ClosestPoint(x, cpx, dist, bdy);

    assert(abs(dist) < 1e-6); // point given should be on the ibc curve (polyline)

    // find the line segement of the polyline that x is on
    size_t min_index = 0;
    min_dist = numeric_limits<Scalar>::max();
    for(size_t i = 0; i < xp().size(); ++i)
    {
        Scalar current_dist = Distance(xp()[i], cpx);
        if(current_dist < min_dist)
        {
            min_dist = current_dist;
            min_index = i;
        }
    }

    // now find the edge that cpx lies on
    size_t next_nearest_index;
    if(min_index == 0)
    {
        next_nearest_index = min_index + 1;
    } 
    else if(min_index == xp().size() - 1)
    {
        next_nearest_index = min_index - 1;
    }
    else
    {
        if(Distance(cpx, xp()[min_index]) < 1e-6) // it is on the point, so choose the shorter edge
        {
            if(Distance(xp()[min_index + 1], cpx) > Distance(xp()[min_index - 1], cpx))
            {
                next_nearest_index = min_index - 1;
            }
            else
            {
                next_nearest_index = min_index + 1;
            }
        }
        else
        {
            vector<Scalar> cpx_xp_min = cpx - xp()[min_index]; 
            vector<Scalar> xp_min_m1_xp_min = xp()[min_index - 1] - xp()[min_index];
            vector<Scalar> xp_min_p1_xp_min = xp()[min_index + 1] - xp()[min_index];
            if(Norm(CrossProduct(cpx_xp_min, xp_min_m1_xp_min)) < 1e-4)
            {
                next_nearest_index = min_index - 1;
            }
            else if(Norm(CrossProduct(cpx_xp_min, xp_min_p1_xp_min)) < 1e-4)
            {
                next_nearest_index = min_index + 1;
            }
            else
            {
                throw std::runtime_error("Surface::IBCNearestTwoXpIndices - problem finding the edge the point is on");
            }
        }
    }


    return vector<size_t>{min_index, next_nearest_index};
}

// 1D parameter along a ibc curve
const Scalar Surface::IBCParameter(const vector<Scalar>& x)
{
    if(SurfaceType().compare("Point") == 0)
    {
        return 0.0;
    }
    else if(SurfaceType().compare("Dumbbell") == 0 || SurfaceType().compare("Planar Curve") == 0 || SurfaceType().compare("Planar Curve 2") == 0 || SurfaceType().compare("Torus Knot") == 0)
    {
        Scalar min_dist = numeric_limits<Scalar>::max();
        size_t min_idx;
        for(size_t i = 0; i < m_mesh.vertices.size(); ++i)
        {
            Scalar dist = Distance(x, m_mesh.vertices[i]);
            if(dist < min_dist)
            {
                min_dist = dist;
                min_idx = i;
            }
        }

        return m_thetap[min_idx];
    }
    else if(SurfaceType().compare("Circle") == 0)
    {
        return atan2(x[1], x[0]);
    }
    else //(SurfaceType().compare("Polyline") == 0) or any other surface with xp computed
    {
        Scalar min_dist;
        vector<size_t> nearest_indices = IBCNearestTwoXpIndices(x, min_dist);

        // compute the total length of the polyline
        Scalar total_length = 0.0;
        Scalar length_to_segment;
        for(size_t i = 0; i < xp().size() - 1; ++i)
        {
            if(i == min(nearest_indices[0], nearest_indices[1]))
            {
                length_to_segment = total_length;
            }

            total_length += Distance(xp()[i], xp()[i+1]);
        }

        Scalar parameter;
        if(nearest_indices[0] < nearest_indices[1])
        {
            parameter = length_to_segment + min_dist;
        }
        else
        {
            parameter = length_to_segment + Distance(xp()[nearest_indices[1]], xp()[nearest_indices[0]]) - min_dist;
        }
        parameter /= total_length;

        return parameter;
    }
}


// 1D parameter along a ibc curve
const vector<Scalar> Surface::CoordinateFromIBCParameter(const Scalar &ibc_parameter)
{
    if(SurfaceType().compare("Polyline") == 0)
    {
        // compute the total length of the polyline
        Scalar total_length = 0.0;
        for(size_t i = 0; i < xp().size() - 1; ++i)
        {
            total_length += Distance(xp()[i], xp()[i+1]);
        }

        // now find the segment the ibc_parameter belongs to
        Scalar current_length = 0.0;
        size_t start_idx;
        size_t end_idx;
        for(size_t i = 0; i < xp().size() - 1; ++i)
        {
            current_length += Distance(xp()[i], xp()[i+1]);
            if(current_length > ibc_parameter * total_length)
            {
                start_idx = i;
                end_idx = i+1;
                break;
            }
        }

        Scalar scale = (current_length - ibc_parameter * total_length);
        vector<Scalar> direction = xp()[start_idx] - xp()[end_idx]; 
        direction *= scale;
        vector<Scalar> coord = xp()[end_idx] + direction;
        return coord;
    }
    else
    {
        throw std::runtime_error("Surface::CoordinateFromIBCParameter not implemented for the surface " + SurfaceType());
        return vector<Scalar>();
    }
}


void Surface::GetParameterizationFaces(size_t &Np, size_t surface_dim)
{
    if(surface_dim == 1)
    {
        for(size_t i = 0; i < Np-1; ++i)
        {
            vector<size_t> edge(2);
            edge[0] = i;
            edge[1] = i+1;
            m_mesh.faces.push_back(edge);
        }
    }
    else if(surface_dim == 2)
    {
        size_t sqrtNp = sqrt(Np);
        for(size_t i = 0; i < sqrtNp - 1; ++i)
        {
            for(size_t j = 0; j < sqrtNp - 1; ++j)
            {
                vector<size_t> quad(4);
                quad[0] = i * sqrtNp + j;
                quad[1] = (i + 1) * sqrtNp + j;
                quad[2] = (i + 1) * sqrtNp + j + 1;
                quad[3] = i * sqrtNp + j + 1;

                vector<size_t> face1(3);
                face1[0] = quad[0];
                face1[1] = quad[1];
                face1[2] = quad[2];

                vector<size_t> face2(3);
                face2[0] = quad[0];
                face2[1] = quad[2];
                face2[2] = quad[3];
                
                m_mesh.faces.push_back(face1);
                m_mesh.faces.push_back(face2);
            }
        }
    }
    else
    {
        throw std::runtime_error("Surface dimension d = " + std::to_string(surface_dim) + " not implemented in Surface::GetParameterizationFaces()");
    }
}


} // namespace cpm