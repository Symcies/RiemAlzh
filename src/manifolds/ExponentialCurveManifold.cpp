#include "ExponentialCurveManifold.h"
#include <cassert>
////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////



ExponentialCurveManifold::ExponentialCurveManifold(unsigned int dim_num)
{
    dimension_ = dim_num;
}


ExponentialCurveManifold::~ExponentialCurveManifold()
{

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

ExponentialCurveManifold::VectorType ExponentialCurveManifold::ComputeParallelCurve(
  VectorType& p0, double t0, VectorType& v0, VectorType& space_shift,
                       double time_point, VectorType& delta)
{

    assert(space_shift.size() == delta.size());

    /// Initialization
    VectorType parallel_curve(space_shift.size());
    double init_pos = p0(0);
    ScalarType * p = parallel_curve.memptr();
    ScalarType * s = space_shift.memptr();
    ScalarType * d = delta.memptr();
    auto n = parallel_curve.size();

#pragma omp simd
    for(size_t i = 0; i < n; ++i)
        p[i] = init_pos * exp(s[i] / (init_pos * exp(d[i])) + d[i] - time_point / init_pos);

    return parallel_curve;
}


ExponentialCurveManifold::VectorType ExponentialCurveManifold::ComputeParallelCurve(
  VectorType& p0, double t0, VectorType& v0, VectorType& space_shift,
                       double time_point)
{
    throw std::invalid_argument( "Wrong overloaded function ComputeParallelCurve called - in Exponential Curve Manifold" );
}

ExponentialCurveManifold::VectorType ExponentialCurveManifold::GetVelocityTransformToEuclideanSpace(VectorType& p0, double t0, VectorType& v0, VectorType& delta)
{

    /// Initialization
    VectorType transformed_velocity(delta.size());
    double init_pos = p0(0);
    double init_vel = v0(0);
    double ratio = init_vel  / (init_pos * init_pos);

    /// Compute the coordinates
    auto iter_transformed_vel = transformed_velocity.begin(), iter_delta = delta.begin();
    for(    ; iter_transformed_vel != transformed_velocity.end() && iter_delta != delta.end(); ++iter_transformed_vel, ++iter_delta)
    {
        *iter_transformed_vel = ratio * exp( - *iter_delta );
    }

    return transformed_velocity;
}


ExponentialCurveManifold::VectorType ExponentialCurveManifold::GetVelocityTransformToEuclideanSpace(VectorType& p0, double t0, VectorType& v0)
{
    throw std::invalid_argument( "Wrong overloaded function GetVelocityTransformToEuclidianSpace called - in ExponentialCurve Manifold" );
}

double ExponentialCurveManifold::ComputeScalarProduct(VectorType& u, VectorType& v, VectorType& application_point)
{
    double scalar_product = 0.0;

    auto iter_u = u.begin(), iter_v = v.begin(), iter_a = application_point.begin();
    for(    ; iter_u != u.end() && iter_v != v.end() && iter_a !=application_point.end()
            ; ++iter_u, ++iter_v, ++iter_a)
    {
        scalar_product += *iter_u * *iter_v / (*iter_a * *iter_a);
    }

    return scalar_product;
}
