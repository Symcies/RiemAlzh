#include "AbstractManifold.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



AbstractManifold::VectorType
AbstractManifold
::ComputeGeodesic(VectorType& P0, double T0, VectorType& V0, double TimePoint)
{
    /// Initialization
    VectorType Geodesic(P0.size());
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();
    int i = 0;

    /// Compute Geodesic
    for(    ; IterPos != P0.end() && IterVel != V0.end() ; ++IterPos, ++IterVel, ++i)
    {
        Geodesic(i) = m_BaseManifold->ComputeGeodesic(*IterPos, T0, *IterVel, TimePoint);
    }

    return Geodesic;
}

AbstractManifold::VectorType
AbstractManifold
::ComputeGeodesicDerivative(VectorType& P0, double T0, VectorType& V0, double TimePoint)
{
    /// Initialization
    VectorType GeodesicDerivative(P0.size());
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();
    int i = 0;

    /// Compute Geodesic
    for(    ; IterPos != P0.end() && IterVel != V0.end() ; ++IterPos, ++IterVel, ++i)
    {
        GeodesicDerivative(i) = m_BaseManifold->ComputeGeodesicDerivative(*IterPos, T0, *IterVel, TimePoint);
    }

    return GeodesicDerivative;
}


AbstractManifold::VectorType
AbstractManifold
::ComputeParallelCurve(VectorType& P0, double T0, VectorType& V0, VectorType& SpaceShift, double TimePoint)
{
    /// Initialization
    VectorType ParallelCurve(SpaceShift.size());
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();
    auto IterShift = SpaceShift.begin();
    int i = 0;
    
    /// Compute Geodesic
    for(    ; IterPos != P0.end() && IterVel != V0.end() && IterShift != SpaceShift.end(); ++IterPos, ++IterVel, ++IterShift, ++i)
    {
        ParallelCurve(i) = m_BaseManifold->ComputeParallelCurve(*IterPos, T0, *IterVel, *IterShift, TimePoint);
    }

    return ParallelCurve;
}