#include "LinearManifold.h"
#include <cassert>


LinearManifold::LinearManifold(unsigned int dim_num)
{
    dimension_ = dim_num;
}

LinearManifold::~LinearManifold(){}



LinearManifold::VectorType LinearManifold::ComputeParallelCurve(VectorType &p0, double t0, VectorType &v0, VectorType &space_shift,
                       double time_point, VectorType& delta)
{
    assert(space_shift.size() == delta.size());

    /// Initialization
    VectorType parallel_curve(space_shift.size());
    double init_pos = p0(0);
    double init_vel = v0(0);
    ScalarType * p = parallel_curve.memptr();
    ScalarType * s = space_shift.memptr();
    ScalarType * d = delta.memptr();
    auto n = parallel_curve.size();

#pragma omp simd
    for(size_t i = 0; i < n; ++i)
        p[i] = init_pos + s[i] + d[i]*init_pos - time_point;

    return parallel_curve;
}


LinearManifold::VectorType LinearManifold::ComputeParallelCurve(VectorType &p0, double t0, VectorType &v0, VectorType &space_shift,
                       double time_point)
{
    throw std::invalid_argument( "Wrong overloaded function ComputeParallelCurve called - in Exponential Curve Manifold" );
}

LinearManifold::VectorType LinearManifold::GetVelocityTransformToEuclideanSpace(VectorType &p0, double t0, VectorType &v0,
                                       VectorType &delta)
{
    double init_vel = v0(0);
    return VectorType(delta.size(), init_vel);
}

LinearManifold::VectorType LinearManifold::GetVelocityTransformToEuclideanSpace(VectorType &p0, double t0, VectorType &v0)
{
    throw std::invalid_argument( "Wrong overloaded function GetVelocityTransformToEuclidianSpace called - in ExponentialCurve Manifold" );
}


double LinearManifold::ComputeScalarProduct(VectorType &u, VectorType &v, VectorType &application_point)
{
    double scalar_product = 0.0;

    auto iter_u = u.begin(), iter_v = v.begin();
    for(    ; iter_u != u.end() && iter_v != v.end()
            ; ++iter_u, ++iter_v)
    {
        scalar_product += *iter_u * *iter_v;
    }

    return scalar_product;
}
