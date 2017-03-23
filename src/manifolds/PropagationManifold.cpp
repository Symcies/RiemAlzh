#include <stdexcept>
#include "PropagationManifold.h"
#include "LogisticBaseManifold.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

PropagationManifold::PropagationManifold(unsigned int dim_num, std::shared_ptr<AbstractBaseManifold>& base_manifold)
{
    dimension_ = dim_num;
    base_manifold_ = base_manifold;
}

PropagationManifold::~PropagationManifold(){}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

PropagationManifold::VectorType PropagationManifold::ComputeGeodesic(
  VectorType& p0, double t0, VectorType& v0, double time_point, VectorType& delta)
{
    /// Initialization
    VectorType geodesic(delta.size());
    double init_point = p0(0);
    double init_vel = v0(0);

    /// Coordinates
    auto iter_geodesic = geodesic.begin();
    for(auto it = delta.begin(); it != delta.end() && iter_geodesic != geodesic.end(); ++it, ++iter_geodesic)
    {
        *iter_geodesic = base_manifold_->ComputeGeodesic(init_point, t0, init_vel, time_point + *it);
    }
    return geodesic;
}


PropagationManifold::VectorType PropagationManifold::ComputeGeodesic(VectorType& p0, double t0, VectorType& v0, double time_point)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function ComputeGeodesic called - in PropagationManifold" );
    VectorType propagation(p0.size(), 0);

    return ComputeGeodesic(p0, t0, v0, time_point, propagation);
}

PropagationManifold::VectorType PropagationManifold::ComputeGeodesicDerivative(VectorType& p0, double t0, VectorType& v0, double time_point,
                            VectorType& delta)
{
    /// Initialize
    VectorType geodesic_derivative(delta.size());
    double init_pos = p0(0);
    double init_vel = v0(0);

    /// Coordinates
    auto iter_geo = geodesic_derivative.begin();
    for( auto it = delta.begin(); it != delta.end() && iter_geo != geodesic_derivative.end(); ++it, ++iter_geo)
    {
        *iter_geo = base_manifold_->ComputeGeodesicDerivative(init_pos, t0, init_vel, time_point + *it);
    }

    return geodesic_derivative;
}

PropagationManifold::VectorType PropagationManifold::ComputeGeodesicDerivative(
  VectorType& p0, double t0, VectorType& v0, double time_point)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function ComputeGeodesicDerivative called - in PropagationManifold" );
    VectorType propagation(p0.size(), 0);

    return ComputeGeodesicDerivative(p0, t0, v0, time_point, propagation);
}


PropagationManifold::VectorType PropagationManifold::ComputeParallelTransport(
  VectorType& p0, double t0, VectorType& v0, VectorType& space_shift,
                       double time_point, VectorType& delta)
{
        /// Initialization
    VectorType parallel_transport(space_shift.size());
    double init_pos = p0(0);
    double init_vel = v0(0);
    auto iter_shift = space_shift.begin();
    auto iter_prop = delta.begin();

    /// Coordinates
    auto iter_parallel = parallel_transport.begin();
    for(    ; iter_prop != delta.end() && iter_shift != space_shift.end() && iter_parallel != parallel_transport.end()
            ; ++iter_prop, ++iter_shift, ++iter_parallel)
    {
        double coordinate = *iter_shift * base_manifold_->ComputeGeodesicDerivative(init_pos, t0, init_vel, time_point + *iter_prop);
        coordinate /= base_manifold_->ComputeGeodesicDerivative(init_pos, t0, init_vel, t0 + *iter_prop);
        *iter_parallel = coordinate;
    }

    return parallel_transport;
}

PropagationManifold::VectorType PropagationManifold::ComputeParallelTransport(
  VectorType& p0, double t0, VectorType& v0, VectorType& space_shift,
                       double time_point)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function ComputeParallelTransport called - in PropagationManifold" );
    VectorType propagation(space_shift.size(), 0);

    return ComputeParallelTransport(p0, t0, v0, space_shift, time_point, propagation);
}

PropagationManifold::VectorType PropagationManifold ::ComputeParallelCurve(
  VectorType& p0, double t0, VectorType& v0, VectorType& space_shift,
                       double time_point, VectorType& delta)
{
    /// Initialization
    VectorType parallel_curve(space_shift.size());
    double init_pos = p0(0);
    double init_vel = v0(0);
    auto iter_shift = space_shift.begin();
    auto iter_prop = delta.begin();

    /// Coordinates
    auto iter_parallel = parallel_curve.begin();
    for(    ; iter_shift != space_shift.end() && iter_prop != delta.end() && iter_parallel != parallel_curve.end()
            ; ++iter_shift, ++iter_prop, ++iter_parallel)
    {
        double time = *iter_shift / base_manifold_->ComputeGeodesic(init_pos, t0, init_vel, t0 + *iter_prop);
        time += time_point + *iter_shift;
        *iter_parallel = base_manifold_->ComputeGeodesic(init_pos, t0, init_vel, time);
    }

    return parallel_curve;
}

PropagationManifold::VectorType PropagationManifold::ComputeParallelCurve(
  VectorType& p0, double t0, VectorType& v0, VectorType& space_shift,
                           double time_point)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function ComputeParallelCurve called - in PropagationManifold" );
    VectorType propagation(space_shift.size(), 0);

    return ComputeParallelCurve(p0, t0, v0, space_shift, time_point, propagation);
}

PropagationManifold::VectorType PropagationManifold::GetVelocityTransformToEuclideanSpace(VectorType& p0, double t0, VectorType& v0, VectorType& delta)
{
    /// Initialization
    VectorType transformed_velocity(delta.size());
    double init_pos = p0(0);
    double init_vel = v0(0);

    /// Coordinates
    auto iter_transformed = transformed_velocity.begin();
    for(auto iter_prop = delta.begin(); iter_prop != delta.end() && iter_transformed != transformed_velocity.end(); ++iter_prop, ++iter_transformed)
    {
        double num = base_manifold_->ComputeGeodesicDerivative(init_pos, t0, init_vel, t0 + *iter_prop);
        double denom = base_manifold_->ComputeGeodesic(init_pos, t0, init_vel, t0 + *iter_prop);
        denom = denom * denom * (1 - denom) * (1 - denom);
        *iter_transformed = num/denom ;
    }

    return transformed_velocity;
}

PropagationManifold::VectorType PropagationManifold::GetVelocityTransformToEuclideanSpace(VectorType& p0, double t0, VectorType& v0)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function GetVelocityTransformToEuclidianSpace called - in PropagationManifold" );
    VectorType propagation(p0.size() - 1, 0);

    return GetVelocityTransformToEuclideanSpace(p0, t0, v0, propagation);
}


double PropagationManifold::ComputeScalarProduct(VectorType& u_vec, VectorType& v_vec, VectorType& application_point)
{
    if(u_vec.size() != v_vec.size() or u_vec.size() != application_point.size())
    {
        std::cout << " The vectors do not have the same size. How to do a scalar product?";
    }

    double scalar_product = 0;
    auto u = u_vec.begin();
    auto v = v_vec.begin();
    auto p = application_point.begin();
    for(    ; u != u_vec.end() && v != v_vec.end() && p != application_point.end()
            ; ++u, ++v, ++p)
    {
        scalar_product += base_manifold_->ComputeScalarProduct(*u, *v, *p);
    }


    return scalar_product;
}
