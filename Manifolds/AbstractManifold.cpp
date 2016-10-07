#include "AbstractManifold.h"
#include "../Tests/TestAssert.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



std::vector<double>
AbstractManifold
::ComputeGeodesic(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint)
{
    /// Tests
    TestAssert::WarningEquality_Object(P0.size(), V0.size(), "P0 and V0 does not have the same size in AbstractManifold > ComputeGeodesic");

    /// Initialization
    std::vector<double> Geodesic;
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();

    /// Compute Geodesic
    for(    ; IterPos != P0.end() && IterVel != V0.end() ; ++IterPos, ++IterVel)
    {
        double Coordinate = m_BaseManifold->ComputeGeodesic(*IterPos, T0, *IterVel, TimePoint);
        Geodesic.push_back( Coordinate );
    }

    return Geodesic;
}

std::vector<double>
AbstractManifold
::ComputeGeodesicDerivative(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint)
{
    /// Tests
    TestAssert::WarningEquality_Object(P0.size(), V0.size(), "P0 and V0 does not have the same size in AbstractManifold > ComputeGeodesicDerivative");

    /// Initialization
    std::vector<double> GeodesicDerivative;
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();

    /// Compute Geodesic
    for(    ; IterPos != P0.end() && IterVel != V0.end() ; ++IterPos, ++IterVel)
    {
        double Coordinate = m_BaseManifold->ComputeGeodesicDerivative(*IterPos, T0, *IterVel, TimePoint);
        GeodesicDerivative.push_back( Coordinate );
    }

    return GeodesicDerivative;
}

std::vector<double>
AbstractManifold
::ComputeParallelTransport(std::vector<double> P0, double T0, std::vector<double> V0, std::vector<double> SpaceShift,
                           double TimePoint)
{
    /// Tests
    TestAssert::WarningEquality_Object(P0.size(), V0.size(), "P0 and V0 does not have the same size in AbstractManifold > ComputeParallelTransport");
    TestAssert::WarningEquality_Object(V0.size(), SpaceShift.size(), "V0 and SpaceShift does not have the same size in AbstractManifold > ComputeParallelTransport");
    TestAssert::WarningEquality_Object(P0.size(), SpaceShift.size(), "P0 and SpaceShift does not have the same size in AbstractManifold > ComputeParallelTransport");
    

    /// Initialization
    std::vector<double> ParallelTransport;
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();
    auto IterShift = SpaceShift.begin();

    /// Compute Geodesic
    for(    ; IterPos != P0.end() && IterVel != V0.end() && IterShift != SpaceShift.end(); ++IterPos, ++IterVel, ++IterShift)
    {
        double Coordinate = m_BaseManifold->ComputeParallelTransport(*IterPos, T0, *IterVel, *IterShift, TimePoint);
        ParallelTransport.push_back( Coordinate );
    }

    return ParallelTransport;
}


std::vector<double>
AbstractManifold
::ComputeParallelCurve(std::vector<double> P0, double T0, std::vector<double> V0, std::vector<double> SpaceShift,
                       double TimePoint)
{
    /// Tests
    TestAssert::WarningEquality_Object(P0.size(), V0.size(), "P0 and V0 does not have the same size in AbstractManifold > ComputeParallelCurve");
    TestAssert::WarningEquality_Object(V0.size(), SpaceShift.size(), "V0 and SpaceShift does not have the same size in AbstractManifold > ComputeParallelCurve");
    TestAssert::WarningEquality_Object(P0.size(), SpaceShift.size(), "P0 and SpaceShift does not have the same size in AbstractManifold > ComputeParallelCurve");

    /// Initialization
    std::vector<double> ParallelCurve;
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();
    auto IterShift = SpaceShift.begin();

    /// Compute Geodesic
    for(    ; IterPos != P0.end() && IterVel != V0.end() && IterShift != SpaceShift.end(); ++IterPos, ++IterVel, ++IterShift)
    {
        double Coordinate = m_BaseManifold->ComputeParallelCurve(*IterPos, T0, *IterVel, *IterShift, TimePoint);
        ParallelCurve.push_back( Coordinate );
    }

    return ParallelCurve;
}