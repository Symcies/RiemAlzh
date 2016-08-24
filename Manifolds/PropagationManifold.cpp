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


    std::vector<double> PropParameters = GetPropagationParameters();
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();
    auto IterProp = PropParameters.begin();

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


    std::vector<double> PropParameters = GetPropagationParameters();
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();
    auto IterProp = PropParameters.begin();

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

    std::vector<double> PropParameters = GetPropagationParameters();
    auto IterProp = PropParameters.begin();
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();
    auto IterSpa = SpaceShift.begin();
    for(    ; IterPos != P0.end() && IterVel != V0.end() && IterSpa != SpaceShift.end() && IterProp != PropParameters.end()
            ; ++IterPos, ++IterVel, ++IterSpa, ++IterProp)
    {
        double Num = *IterSpa * m_BaseManifold->ComputeGeodesicDerivative(*IterPos, T0, *IterVel, TimePoint + *IterProp);
        double Denom = m_BaseManifold->ComputeGeodesicDerivative(*IterPos, T0, *IterVel, TimePoint + *IterProp);
        ParallelTransport.push_back( Num / Denom );
    }

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

    std::vector<double> PropParameters = GetPropagationParameters();
    auto IterProp = PropParameters.begin();
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();
    auto IterSpa = SpaceShift.begin();


    for(    ; IterProp != PropParameters.end() && IterPos != P0.end() && IterVel != V0.end() && IterSpa != SpaceShift.end()
            ; ++IterProp, ++IterPos, ++IterVel, ++IterSpa)
    {
        double T = *IterSpa / m_BaseManifold->ComputeGeodesicDerivative(*IterPos, T0, *IterVel, T0 + *IterProp );
        T += TimePoint + *IterProp;
        double Coordinate = m_BaseManifold->ComputeGeodesic(*IterPos, T0, *IterVel, T + *IterProp );

        ParallelTransport.push_back(Coordinate);
    }

    return ParallelTransport;
}

std::vector<double>
PropagationManifold
::GetVelocityTransformToEuclideanSpace(std::vector<double> P0, double T0, std::vector<double> V0)
{

    std::vector<double> TransformedVelocity;


    std::vector<double> PropParameters = GetPropagationParameters();
    auto IterPos = P0.begin();
    auto IterVel = V0.begin();
    auto IterPro = PropParameters.begin();
    for(    ; IterPos != P0.end() && IterVel != V0.end() && IterPro != PropParameters.end()
            ; ++IterPos, ++IterVel, ++IterPro)
    {
        double Geo = m_BaseManifold->ComputeGeodesic(*IterPos, T0, *IterVel, T0 + *IterPro);
        double GeoDer = m_BaseManifold->ComputeGeodesicDerivative(*IterPos, T0, *IterVel, T0 + *IterPro);
        double Coordinate = GeoDer / ( Geo * Geo * (1.0-Geo) * (1.0-Geo) );
        TransformedVelocity.push_back(Coordinate);
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

    std::cout << "Coord by coord : ";
    for(int i = 0; i < U.size(); ++i)
    {
        double P = ApplicationPoint[i];
        P = 1.0 / ( P * P * (1-P) * (1-P) );
        
        ScalarProduct += U[i] * P * V[i];
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
    for(int i = 0; i < m_Dimension; ++i)
    {
        PropagationParameters.push_back( m_Parameters.at("Delta" + std::to_string(i)) );
    }

    return PropagationParameters;
}
