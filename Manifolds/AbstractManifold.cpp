#include "AbstractManifold.h"
#include "../Tests/TestAssert.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



AbstractManifold::VectorType
AbstractManifold
::ComputeGeodesic(VectorType& P0, double T0, VectorType& V0, double TimePoint)
{
    /// Tests
    TestAssert::WarningEquality_Object(P0.size(), V0.size(), "P0 and V0 does not have the same size in AbstractManifold > ComputeGeodesic");

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
    /// Tests
    TestAssert::WarningEquality_Object(P0.size(), V0.size(), "P0 and V0 does not have the same size in AbstractManifold > ComputeGeodesicDerivative");

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
    /// Tests
    TestAssert::WarningEquality_Object(P0.size(), V0.size(), "P0 and V0 does not have the same size in AbstractManifold > ComputeParallelCurve");
    TestAssert::WarningEquality_Object(V0.size(), SpaceShift.size(), "V0 and SpaceShift does not have the same size in AbstractManifold > ComputeParallelCurve");
    TestAssert::WarningEquality_Object(P0.size(), SpaceShift.size(), "P0 and SpaceShift does not have the same size in AbstractManifold > ComputeParallelCurve");

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