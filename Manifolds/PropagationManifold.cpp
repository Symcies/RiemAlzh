#include "PropagationManifold.h"

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
::ComputeGeodesic(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint)
{
    std::vector<double> Geodesic;

    /// Initialization
    std::vector<double> PropParameters = GetPropagationParameters();
    auto IterPos = P0.begin();
    auto IterVel = V0.begin() ;
    auto IterProp = PropParameters.begin();

    /// First component
    double FirstComponent = m_BaseManifold->ComputeGeodesic(*IterPos, T0, *IterVel, TimePoint);
    Geodesic.push_back( FirstComponent );
    ++IterPos;
    ++IterVel;

    /// Next components
    for(    ; IterPos != P0.end() && IterVel != V0.end() && IterProp != PropParameters.end()
            ; ++IterPos, ++IterVel, ++IterProp)
    {
        double Coordinate = m_BaseManifold->ComputeGeodesic(*IterPos, T0, *IterVel, TimePoint + *IterProp);
        Geodesic.push_back(Coordinate);
    }

    return Geodesic;
}

std::vector<double>
PropagationManifold
::ComputeGeodesicDerivative(std::vector<double> P0, double T0, std::vector<double> V0, double TimePoint)
{
    std::vector<double> GeodesicDerivative;

    // Initialization
    std::vector<double> PropParameters = GetPropagationParameters();
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();
    auto IterProp = PropParameters.begin();

    // First component
    double FirstComponent = m_BaseManifold->ComputeGeodesicDerivative(*IterPos, T0, *IterVel, TimePoint);
    GeodesicDerivative.push_back( FirstComponent );
    ++IterPos;
    ++IterVel;

    // Next components
    for(    ; IterPos != P0.end() && IterVel != V0.end() && IterProp != PropParameters.end()
            ; ++IterPos, ++IterVel, ++IterProp)
    {
        double Coordinate = m_BaseManifold->ComputeGeodesicDerivative(*IterPos, T0, *IterVel, TimePoint + *IterProp);
        GeodesicDerivative.push_back(Coordinate);
    }

    return GeodesicDerivative;

}

std::vector<double>
PropagationManifold
::ComputeParallelTransport(std::vector<double> P0, double T0, std::vector<double> V0, std::vector<double> SpaceShift,
                           double TimePoint)
{
    std::vector<double> ParallelTransport;

    // Initialization
    std::vector<double> GeodesicDerivative0 = ComputeGeodesicDerivative(P0, T0, V0, T0);
    std::vector<double> GeodesicDerivative = ComputeGeodesicDerivative(P0, T0, V0, TimePoint);

    auto IterSpace  = SpaceShift.begin();
    auto IterDeriv0 = GeodesicDerivative0.begin();
    auto IterDeriv  = GeodesicDerivative.begin();

    // Compute parallel transport
    for(    ; IterSpace != SpaceShift.end() && IterDeriv0 != GeodesicDerivative0.end() && IterDeriv != GeodesicDerivative.end()
            ; ++IterSpace, ++IterDeriv0, ++IterDeriv)
    {
        double Coordinate = *IterSpace * *IterDeriv / *IterDeriv0;
        if(isnan(Coordinate))
        {
            std::cout << "here : " << Coordinate << std::endl;
        }
        ParallelTransport.push_back( Coordinate );
    }
    std::cout << "here : " << ParallelTransport[0] << std::endl;
    return ParallelTransport;

}

std::vector<double>
PropagationManifold
::ComputeParallelCurve(std::vector<double> P0, double T0, std::vector<double> V0, std::vector<double> SpaceShift,
                           double TimePoint)
{
    // TODO **UNIT TEST** : Check if P0, V0, PropagCoeff and SpaceShift of the same size
    // TODO : Check if it is not more efficient to use the ComputeGeodesic & ComputeGeodesicDerivative of the PropagationManifold method
    // Thus, it implies : instead of coeff per coeff, calculate geo and geoder, then loop

    std::vector<double> ParallelTransport;

    // Initialization
    std::vector<double> GeodesicDerivative0 = ComputeGeodesicDerivative(P0, T0, V0, T0);
    std::vector<double> PropParameters = GetPropagationParameters();

    auto IterDeriv0 = GeodesicDerivative0.begin();
    auto IterSpace  = SpaceShift.begin();
    auto IterProp   = PropParameters.begin();
    auto IterPos    = P0.begin();
    auto IterVel    = V0.begin();

    //Compute parallel curve
    for(    ; IterDeriv0 != GeodesicDerivative0.end() && IterSpace != SpaceShift.end() && IterProp != PropParameters.end() && IterPos != P0.end() && IterVel != V0.end()
            ; ++IterDeriv0, ++IterSpace, ++IterProp, ++IterPos, ++IterVel)
    {
        double Coordinate = m_BaseManifold->ComputeGeodesic(*IterPos, T0, *IterVel, *IterSpace / *IterDeriv0 + TimePoint + *IterProp);
        ParallelTransport.push_back( Coordinate );
    }

    return ParallelTransport;
}

std::vector<double>
PropagationManifold
::GetVelocityTransformToEuclideanSpace(std::vector<double> P0, double T0, std::vector<double> V0)
{

    std::vector<double> TransformedVelocity;

    // Initialization
    std::vector<double> Geodesic0 = ComputeGeodesic(P0, T0, V0, T0);
    std::vector<double> GeodesicDerivative0 = ComputeGeodesicDerivative(P0, T0, V0, T0);

    auto IterGeo = Geodesic0.begin();
    auto IterGeoDer = GeodesicDerivative0.begin();

    // Compute velocity transformation
    for(    ; IterGeo != Geodesic0.end() && IterGeoDer != GeodesicDerivative0.end()
            ; ++IterGeo, ++IterGeoDer)
    {
        double Denom = *IterGeo * *IterGeo * (1.0 - *IterGeo) * (1.0 - *IterGeo);
        if(Denom == 0)
        {
            std::cout << "ici P0 & V0 : ";
            for(int i = 0; i < P0.size(); ++i)
            {
                std::cout << P0[i] << " & " << V0[i] << " - ";
            }
            std::cout << std::endl;
        }
        else
        {

        }
        double Coordinate = *IterGeoDer / Denom;
        TransformedVelocity.push_back( Coordinate );
    }

    return TransformedVelocity;
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

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double>
PropagationManifold
::GetPropagationParameters()
{
    std::vector<double> PropagationParameters;
    for(int i = 0; i < m_Dimension - 1; ++i)
    {
        PropagationParameters.push_back( m_Parameters.at("Delta" + std::to_string(i)) );
    }

    return PropagationParameters;
}
