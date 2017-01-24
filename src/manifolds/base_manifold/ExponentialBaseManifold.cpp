#include "ExponentialBaseManifold.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

ExponentialBaseManifold
::ExponentialBaseManifold() 
{
    
}

ExponentialBaseManifold
::~ExponentialBaseManifold() 
{
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


double
ExponentialBaseManifold
::ComputeGeodesic(double P0, double T0, double V0, double TimePoint) 
{
    
}

double 
ExponentialBaseManifold
::ComputeGeodesicDerivative(double P0, double T0, double V0, double TimePoint) 
{
    return V0 * exp( V0 / P0 * ( TimePoint - T0) );
}


double 
ExponentialBaseManifold
::ComputeParallelCurve(double P0, double T0, double V0, double SpaceShift, double TimePoint) 
{
    double Time = SpaceShift / ComputeGeodesicDerivative(P0, T0, V0, T0) + TimePoint;
    return ComputeGeodesic(P0, T0, V0, Time);
}

double 
ExponentialBaseManifold
::ComputeScalarProduct(double U, double V, double ApplicationPoint) 
{
    double G = ApplicationPoint * ApplicationPoint;
    return U * V / G;
}