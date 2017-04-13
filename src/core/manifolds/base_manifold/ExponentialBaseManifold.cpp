#include "ExponentialBaseManifold.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

ExponentialBaseManifold::ExponentialBaseManifold()
{

}

ExponentialBaseManifold::~ExponentialBaseManifold()
{

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


double ExponentialBaseManifold::ComputeGeodesic(double p0, double t0, double v0, double time_point)
{

}

double ExponentialBaseManifold::ComputeGeodesicDerivative(double p0, double t0, double v0, double time_point)
{
    return v0 * exp( v0 / p0 * ( time_point - t0) );
}


double ExponentialBaseManifold::ComputeParallelCurve(double p0, double t0, double v0, double space_shift, double time_point)
{
    double time = space_shift / ComputeGeodesicDerivative(p0, t0, v0, t0) + time_point;
    return ComputeGeodesic(p0, t0, v0, time);
}

double ExponentialBaseManifold::ComputeScalarProduct(double u, double v, double application_point)
{
    double g = application_point * application_point;
    return u * v / g;
}
