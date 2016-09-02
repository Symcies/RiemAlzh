#include "AbstractManifold.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



std::vector<double>
AbstractManifold
::ComputeGeodesic(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint)
{
    // TODO : check if P0 and V0 of the same size

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
    // TODO : check if P0 and V0 of the same size

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
    // TODO : check if P0, V0  and Space Shift of the same size

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
    // TODO : check if P0, V0  and Space Shift of the same size

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