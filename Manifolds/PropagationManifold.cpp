#include <stdexcept>
#include "PropagationManifold.h"
#include "BaseManifold/LogisticBaseManifold.h"
#include "../Tests/TestAssert.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

PropagationManifold
::PropagationManifold(unsigned int NumberDimension, std::shared_ptr<AbstractBaseManifold>& BM)
{
    m_Dimension = NumberDimension;
    m_BaseManifold = BM;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double>
PropagationManifold
::ComputeGeodesic(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint,
                  std::vector<double> Delta)
{
    /// Initialization
    std::vector<double> Geodesic;
    double InitialPoint = P0[0];
    double InitialVelocity = V0[0];

    /// First coordinate
    Geodesic.push_back( m_BaseManifold->ComputeGeodesic(InitialPoint, T0, InitialVelocity, TimePoint) );

    /// Next coordinates
    for( auto it : Delta)
    {
        double Coordinate = m_BaseManifold->ComputeGeodesic(InitialPoint, T0, InitialVelocity, TimePoint + it);
        Geodesic.push_back( Coordinate );
    }

    return Geodesic;
}


std::vector<double>
PropagationManifold
::ComputeGeodesic(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function ComputeGeodesic called - in PropagationManifold" );
    std::vector<double> Propagation(P0.size() - 1, 0);

    return ComputeGeodesic(P0, T0, V0, TimePoint, Propagation);
}

std::vector<double>
PropagationManifold
::ComputeGeodesicDerivative(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint,
                            std::vector<double> Delta)
{
    /// Initialize
    std::vector<double> GeodesicDerivative;
    double InitialPosition = P0[0];
    double InitialVelocity = V0[0];

    /// First coordinate
    double Coordinate = m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, TimePoint);
    GeodesicDerivative.push_back( Coordinate );

    /// Next coordinate
    for( auto it : Delta)
    {
        Coordinate = m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, TimePoint + it);
        GeodesicDerivative.push_back( Coordinate );
    }

    return GeodesicDerivative;
}

std::vector<double>
PropagationManifold
::ComputeGeodesicDerivative(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function ComputeGeodesicDerivative called - in PropagationManifold" );
    std::vector<double> Propagation(P0.size() - 1, 0);

    return ComputeGeodesicDerivative(P0, T0, V0, TimePoint, Propagation);
}



std::vector<double>
PropagationManifold
::ComputeParallelTransport(std::vector<double> P0, double T0, std::vector<double> V0, std::vector<double> SpaceShift,
                       double TimePoint, std::vector<double> Delta)
{
        /// Initialization
    std::vector<double> ParallelTransport;
    double InitialPosition = P0[0];
    double InitialVelocity = V0[0];
    auto IterShift = SpaceShift.begin();
    auto IterProp = Delta.begin();

    /// First coordinate
    double FirstCoordinate = *IterShift * m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, TimePoint);
    FirstCoordinate /= m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, T0);
    ParallelTransport.push_back( FirstCoordinate );
    ++IterShift;

    /// Next coordinates
    for(    ; IterProp != Delta.end() && IterShift != SpaceShift.end(); ++IterProp, ++IterShift)
    {
        double Coordinate = *IterShift * m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, TimePoint + *IterProp);
        Coordinate /= m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, T0 + *IterProp);
        ParallelTransport.push_back( Coordinate );
    }


    return ParallelTransport;
}

std::vector<double>
PropagationManifold
::ComputeParallelTransport(std::vector<double> P0, double T0, std::vector<double> V0, std::vector<double> SpaceShift,
                       double TimePoint)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function ComputeParallelTransport called - in PropagationManifold" );
    std::vector<double> Propagation(P0.size() - 1, 0);

    return ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint, Propagation);
}

std::vector<double>
PropagationManifold
::ComputeParallelCurve(std::vector<double> P0, double T0, std::vector<double> V0, std::vector<double> SpaceShift,
                       double TimePoint, std::vector<double> Delta)
{
    /// Initialization
    std::vector<double> ParallelCurve;
    double InitialPosition = P0[0];
    double InitialVelocity = V0[0];
    auto IterShift = SpaceShift.begin();
    auto IterProp = Delta.begin();

    /// First coordinate
    double FirstTimePoint = *IterShift / m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, T0) + TimePoint;
    double FirstCoordinate = m_BaseManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, FirstTimePoint );
    ParallelCurve.push_back( FirstCoordinate );
    ++IterShift;

    /// Next coordinates
    for(    ; IterShift != SpaceShift.end() && IterProp != Delta.end(); ++IterShift, ++IterProp)
    {
        double Time = *IterShift / m_BaseManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, T0 + *IterProp);
        Time += TimePoint + *IterShift;
        double Coordinate = m_BaseManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, Time);
        ParallelCurve.push_back(Coordinate);
    }

    return ParallelCurve;
}

std::vector<double>
PropagationManifold
::ComputeParallelCurve(std::vector<double> P0, double T0, std::vector<double> V0, std::vector<double> SpaceShift,
                           double TimePoint)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function ComputeParallelCurve called - in PropagationManifold" );
    std::vector<double> Propagation(P0.size() - 1, 0);

    return ComputeParallelTransport(P0, T0, V0, SpaceShift, TimePoint, Propagation);
}

std::vector<double>
PropagationManifold
::GetVelocityTransformToEuclideanSpace(std::vector<double> P0, double T0, std::vector<double> V0,
                                       std::vector<double> Delta)
{
    /// Initialization
    std::vector<double> TransformedVelocity;
    double InitialPosition = P0[0];
    double InitialVelocity = V0[0];
    
    /// First coordinate
    double FirstGeoDeriv = m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, T0);
    double FirstGeo =  m_BaseManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, T0);
    double FirstCoordinate = FirstGeoDeriv / ( FirstGeo * FirstGeo * (1 - FirstGeo) * (1 - FirstGeo) );
    TransformedVelocity.push_back( FirstCoordinate );

    /// Next coordinates
    for(auto IterProp = Delta.begin(); IterProp != Delta.end(); ++IterProp)
    {
        double Num = m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, T0 + *IterProp);
        double Denom = m_BaseManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, T0 + *IterProp);
        Denom = Denom * Denom * (1 - Denom) * (1 - Denom);
        TransformedVelocity.push_back( Num/Denom );
    }

    return TransformedVelocity;
}

std::vector<double>
PropagationManifold
::GetVelocityTransformToEuclideanSpace(std::vector<double> P0, double T0, std::vector<double> V0)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function GetVelocityTransformToEuclidianSpace called - in PropagationManifold" );
    std::vector<double> Propagation(P0.size() - 1, 0);

    return GetVelocityTransformToEuclideanSpace(P0, T0, V0, Propagation);
}


double
PropagationManifold
::ComputeScalarProduct(std::vector<double> U, std::vector<double> V, std::vector<double> ApplicationPoint)
{
    if(U.size() != V.size() or U.size() != ApplicationPoint.size())
    {
        std::cout << " The vectors do not have the same size. How to do a scalar product?";
    }

    double ScalarProduct = 0;
    auto u = U.begin();
    auto v = V.begin();
    auto p = ApplicationPoint.begin();
    for(    ; u != U.end() && v != V.end() && p != ApplicationPoint.end()
            ; ++u, ++v, ++p)
    {
        ScalarProduct += m_BaseManifold->ComputeScalarProduct(*u, *v, *p);
    }


    return ScalarProduct;
}