#include "LogisticBaseManifold.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

LogisticBaseManifold::LogisticBaseManifold()
{ }

LogisticBaseManifold::~LogisticBaseManifold()
{ }

////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double LogisticBaseManifold::ComputeGeodesic(double p0, double t0, double v0, double time_point)
{
    double value = - v0 * (time_point - t0) / (p0 * (1.-p0));
    value = 1. + (1./p0 - 1. )*exp(value) ;
    return 1./value;
}

double LogisticBaseManifold::ComputeGeodesicDerivative(double p0, double t0, double v0, double time_point)
{
    double value = exp(- v0 * (time_point - t0) / (p0 * (1.-p0)));
    double num = v0 * value;
    double denom = (p0 + (1-p0)*value) * (p0 + (1-p0)*value);
    return num/denom;
}

double
LogisticBaseManifold
::ComputeParallelTransport(double p0, double t0, double v0, double space_shift, double time_point)
{
    return space_shift * ComputeGeodesicDerivative(p0, t0, v0, time_point) / ComputeGeodesicDerivative(p0, t0, v0, t0);
}

double LogisticBaseManifold::ComputeParallelCurve(double p0, double t0, double v0, double space_shift, double time_point)
{
    double time = space_shift / ComputeGeodesicDerivative(p0, t0, v0, t0) + time_point;
    return ComputeGeodesic(p0, t0, v0, time);
}

double LogisticBaseManifold::ComputeScalarProduct(double u, double v, double application_point)
{
    double g = application_point * application_point * (1 - application_point) * (1 - application_point);
    return u * v / g;
}
