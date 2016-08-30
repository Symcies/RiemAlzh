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
    //std::cout << "TimePoint & Value : " << TimePoint << " & " << Value << " & ";
    Value = 1. + (1./P0 - 1. )*exp(Value) ;
    //std::cout << Value << " & " << 1.0/Value << std::endl;
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
::ComputeScalarProduct(double U, double V, double ApplicationPoint)
{
    double G = ApplicationPoint * ApplicationPoint * (1 - ApplicationPoint) * (1 - ApplicationPoint);
    return U * V / G;
}