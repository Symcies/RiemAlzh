#include "AbstractManifold.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



AbstractManifold::VectorType AbstractManifold::ComputeGeodesic(VectorType& p0, double t0, VectorType& v0, double time_point)
{
    /// Initialization
    VectorType geodesic(p0.size());
    auto iter_pos = p0.begin();
    auto iter_vel = v0.begin();
    int i = 0;

    /// Compute geodesic
    for(    ; iter_pos != p0.end() && iter_vel != v0.end() ; ++iter_pos, ++iter_vel, ++i)
    {
        geodesic(i) = base_manifold_->ComputeGeodesic(*iter_pos, t0, *iter_vel, time_point);
    }

    return geodesic;
}

AbstractManifold::VectorType AbstractManifold::ComputeGeodesicDerivative(VectorType& p0, double t0, VectorType& v0, double time_point)
{
    /// Initialization
    VectorType geodesic_derivative(p0.size());
    auto iter_pos = p0.begin();
    auto iter_vel = v0.begin();
    int i = 0;

    /// Compute geodesic
    for(    ; iter_pos != p0.end() && iter_vel != v0.end() ; ++iter_pos, ++iter_vel, ++i)
    {
        geodesic_derivative(i) = base_manifold_->ComputeGeodesicDerivative(*iter_pos, t0, *iter_vel, time_point);
    }

    return geodesic_derivative;
}


AbstractManifold::VectorType AbstractManifold::ComputeParallelCurve(VectorType& p0, double t0, VectorType& v0, VectorType& space_shift, double time_point)
{
    /// Initialization
    VectorType parallel_curve(space_shift.size());
    auto iter_pos = p0.begin();
    auto iter_vel = v0.begin();
    auto iter_shift = space_shift.begin();
    int i = 0;

    /// Compute geodesic
    for(    ; iter_pos != p0.end() && iter_vel != v0.end() && iter_shift != space_shift.end(); ++iter_pos, ++iter_vel, ++iter_shift, ++i)
    {
        parallel_curve(i) = base_manifold_->ComputeParallelCurve(*iter_pos, t0, *iter_vel, *iter_shift, time_point);
    }

    return parallel_curve;
}
