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
::ComputeParallelCurve(std::vector<double> P0, double T0, std::vector<double> V0, std::vector<double> SpaceShift,
                           double TimePoint)
{
    // TODO **UNIT TEST** : Check if P0, V0, PropagCoeff and SpaceShift of the same size

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
        double Coordinate = m_BaseManifold->ComputeGeodesic(*IterPos, T0, *IterVel, T);

        ParallelTransport.push_back(Coordinate);
    }

    return ParallelTransport;
}


std::vector<double>
PropagationManifold
::GetVelocityTransformToEuclideanSpace(std::vector<double> P0, double T0, std::vector<double> V0)
{
    double Geo = m_BaseManifold->ComputeGeodesic(P0[0], T0, V0[0], T0);
    double GeoDer= m_BaseManifold->ComputeGeodesicDerivative(P0[0], T0, V0[0], T0);
    double Coordinate = GeoDer / ( Geo * Geo * (1.0-Geo) * (1.0-Geo) );

    std::vector<double> TransformedVelocity;
    for(int i = 0; i < m_Dimension; ++i)
    {
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

    for(int i = 0; i < U.size(); ++i)
    {
        double P = ApplicationPoint[i];
        P = 1.0 / ( P * P * (1-P) * (1-P) );
        std::cout << "U*V & G(p) : " << U[i]*V[i] << " / " << P << std::endl;
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
