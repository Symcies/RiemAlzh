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

PropagationManifold::VectorType
PropagationManifold
::ComputeGeodesic(VectorType P0, double T0, VectorType V0, double TimePoint, VectorType Delta)
{
    /// Initialization
    VectorType Geodesic(Delta.size());
    double InitialPoint = P0(0);
    double InitialVelocity = V0(0);

    /// Coordinates
    auto IterGeodesic = Geodesic.begin();
    for(auto it = Delta.begin(); it != Delta.end() && IterGeodesic != Geodesic.end(); ++it, ++IterGeodesic)
    { 
        *IterGeodesic = m_BaseManifold->ComputeGeodesic(InitialPoint, T0, InitialVelocity, TimePoint + *it);
    }
    return Geodesic;
}


PropagationManifold::VectorType
PropagationManifold
::ComputeGeodesic(VectorType P0, double T0, VectorType V0, double TimePoint)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function ComputeGeodesic called - in PropagationManifold" );
    VectorType Propagation(P0.size(), 0);

    return ComputeGeodesic(P0, T0, V0, TimePoint, Propagation);
}

PropagationManifold::VectorType
PropagationManifold
::ComputeGeodesicDerivative(VectorType P0, double T0, VectorType V0, double TimePoint,
                            VectorType Delta)
{
    /// Initialize
    VectorType GeodesicDerivative(Delta.size());
    double InitialPosition = P0(0);
    double InitialVelocity = V0(0);

    /// Coordinates
    auto IterGeo = GeodesicDerivative.begin();
    for( auto it = Delta.begin(); it != Delta.end() && IterGeo != GeodesicDerivative.end(); ++it, ++IterGeo)
    {
        *IterGeo = m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, TimePoint + *it);
    }

    return GeodesicDerivative;
}

PropagationManifold::VectorType
PropagationManifold
::ComputeGeodesicDerivative(VectorType P0, double T0, VectorType V0, double TimePoint)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function ComputeGeodesicDerivative called - in PropagationManifold" );
    VectorType Propagation(P0.size(), 0);

    return ComputeGeodesicDerivative(P0, T0, V0, TimePoint, Propagation);
}



PropagationManifold::VectorType
PropagationManifold
::ComputeParallelTransport(VectorType P0, double T0, VectorType V0, VectorType SpaceShift,
                       double TimePoint, VectorType Delta)
{
        /// Initialization
    VectorType ParallelTransport(SpaceShift.size());
    double InitialPosition = P0(0);
    double InitialVelocity = V0(0);
    auto IterShift = SpaceShift.begin();
    auto IterProp = Delta.begin();
    
    /// Coordinates
    auto IterParallel = ParallelTransport.begin();
    for(    ; IterProp != Delta.end() && IterShift != SpaceShift.end() && IterParallel != ParallelTransport.end()
            ; ++IterProp, ++IterShift, ++IterParallel)
    {
        double Coordinate = *IterShift * m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, TimePoint + *IterProp);
        Coordinate /= m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, T0 + *IterProp);
        *IterParallel = Coordinate;
    }
    
    return ParallelTransport;
}

PropagationManifold::VectorType
PropagationManifold
::ComputeParallelTransport(VectorType P0, double T0, VectorType V0, VectorType SpaceShift,
                       double TimePoint)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function ComputeParallelTransport called - in PropagationManifold" );
    VectorType Propagation(SpaceShift.size(), 0);

    return ComputeParallelTransport(P0, T0, V0, SpaceShift, TimePoint, Propagation);
}

PropagationManifold::VectorType
PropagationManifold
::ComputeParallelCurve(VectorType P0, double T0, VectorType V0, VectorType SpaceShift,
                       double TimePoint, VectorType Delta)
{
    /// Initialization
    VectorType ParallelCurve(SpaceShift.size());
    double InitialPosition = P0(0);
    double InitialVelocity = V0(0);
    auto IterShift = SpaceShift.begin();
    auto IterProp = Delta.begin();

    /// Coordinates
    auto IterParallel = ParallelCurve.begin();
    for(    ; IterShift != SpaceShift.end() && IterProp != Delta.end() && IterParallel != ParallelCurve.end()
            ; ++IterShift, ++IterProp, ++IterParallel)
    {
        double Time = *IterShift / m_BaseManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, T0 + *IterProp);
        Time += TimePoint + *IterShift;
        *IterParallel = m_BaseManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, Time);
    }

    return ParallelCurve;
}

PropagationManifold::VectorType
PropagationManifold
::ComputeParallelCurve(VectorType P0, double T0, VectorType V0, VectorType SpaceShift,
                           double TimePoint)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function ComputeParallelCurve called - in PropagationManifold" );
    VectorType Propagation(SpaceShift.size(), 0);

    return ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint, Propagation);
}

PropagationManifold::VectorType
PropagationManifold
::GetVelocityTransformToEuclideanSpace(VectorType P0, double T0, VectorType V0, VectorType Delta)
{
    /// Initialization
    VectorType TransformedVelocity(Delta.size());
    double InitialPosition = P0(0);
    double InitialVelocity = V0(0);

    /// Coordinates
    auto IterTransformed = TransformedVelocity.begin();
    for(auto IterProp = Delta.begin(); IterProp != Delta.end() && IterTransformed != TransformedVelocity.end(); ++IterProp, ++IterTransformed)
    {
        double Num = m_BaseManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, T0 + *IterProp);
        double Denom = m_BaseManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, T0 + *IterProp);
        Denom = Denom * Denom * (1 - Denom) * (1 - Denom);
        *IterTransformed = Num/Denom ;
    }

    return TransformedVelocity;
}

PropagationManifold::VectorType
PropagationManifold
::GetVelocityTransformToEuclideanSpace(VectorType P0, double T0, VectorType V0)
{
    // TODO : When no more bug, from throw(Error)
    throw std::invalid_argument( "Wrong overloaded function GetVelocityTransformToEuclidianSpace called - in PropagationManifold" );
    VectorType Propagation(P0.size() - 1, 0);

    return GetVelocityTransformToEuclideanSpace(P0, T0, V0, Propagation);
}


double
PropagationManifold
::ComputeScalarProduct(VectorType U, VectorType V, VectorType ApplicationPoint)
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