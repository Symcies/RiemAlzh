#include "LogisticBaseManifold.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

LogisticBaseManifold
::LogisticBaseManifold()
{ }

LogisticBaseManifold
::~LogisticBaseManifold()
{ }

////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double
LogisticBaseManifold
::ComputeGeodesic(double P0, double T0, double V0, double TimePoint)
{
    double Value = - V0 * (TimePoint - T0) / (P0 * (1.-P0));
    Value = 1. + (1./P0 - 1. )*exp(Value) ;
    return 1./Value;
}

double
LogisticBaseManifold
::ComputeGeodesicDerivative(double P0, double T0, double V0, double TimePoint)
{
    double Value = exp(- V0 * (TimePoint - T0) / (P0 * (1.-P0)));
    double Num = V0 * Value;
    double Denom = (P0 + (1-P0)*Value) * (P0 + (1-P0)*Value);
    return Num/Denom;
}

double
LogisticBaseManifold
::ComputeParallelTransport(double P0, double T0, double V0, double SpaceShift, double TimePoint)
{
    return SpaceShift * ComputeGeodesicDerivative(P0, T0, V0, TimePoint) / ComputeGeodesicDerivative(P0, T0, V0, T0);
}

double
LogisticBaseManifold
::ComputeParallelCurve(double P0, double T0, double V0, double SpaceShift, double TimePoint)
{
    double Time = SpaceShift / ComputeGeodesicDerivative(P0, T0, V0, T0) + TimePoint;
    return ComputeGeodesic(P0, T0, V0, Time);
}

double
LogisticBaseManifold
::ComputeScalarProduct(double U, double V, double ApplicationPoint)
{
    double G = ApplicationPoint * ApplicationPoint * (1 - ApplicationPoint) * (1 - ApplicationPoint);
    return U * V / G;
}

