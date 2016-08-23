#include "PropagationManifold.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

PropagationManifold
::PropagationManifold(unsigned int NumberDimension)
{
    m_Dimension = NumberDimension;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


const std::vector<double>
PropagationManifold
::GetGeodesicDerivative(double TimePoint, const std::shared_ptr<Realizations>& R)
{
    /// Get the data from the realisation
    double P0 = R->at("P0")[0];
    double T0 = R->at("T0")[0];
    double V0 = R->at("V0")[0];


    /// Compute the geodesic derivative
    std::vector<double> GeodesicDerivative;
    for(int i = 0; i < m_Dimension ; ++i)
    {
        double Coordinate = ComputeOneDimensionalGeodesicDerivative(P0, T0, V0, TimePoint);
        GeodesicDerivative.push_back( Coordinate );
    }

    return GeodesicDerivative;
}

const std::vector<double>
PropagationManifold
::GetGeodesic(double TimePoint, const std::shared_ptr<Realizations>& R)
{
    /// Get the data from the realisation
    double P0 = R->at("P0")[0];
    double T0 = R->at("T0")[0];
    double V0 = R->at("V0")[0];


    /// Compute the geodesic derivative
    std::vector<double> Geodesic;
    for(int i = 0; i < m_Dimension ; ++i)
    {
        double Coordinate = ComputeOneDimensionalGeodesic(P0, T0, V0, TimePoint);
        Geodesic.push_back( Coordinate );
    }
    //std::cout << Geodesic[0] << std::endl;
    return Geodesic;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double>
PropagationManifold
::ComputeParallelTransport(double T, std::vector<double> W0, const std::shared_ptr<Realizations> &R)
{
    /// Get the data from the realisation
    double P0 = R->at("P0")[0];
    double T0 = R->at("T0")[0];
    double V0 = R->at("V0")[0];

    /// Compute the parallel transport
    std::vector<double> ParallelTransport;

    for(int i =0; i<m_Dimension ; ++i)
    {
        double PropagationRealization = R->at("Delta" + std::to_string(i))[0];

        double Num = W0[i] * ComputeOneDimensionalGeodesicDerivative(P0, T0, V0, T + PropagationRealization);
        double Denom = ComputeOneDimensionalGeodesicDerivative(P0, T0, V0, T0 + PropagationRealization);
        ParallelTransport.push_back(Num / Denom);
    }

    return ParallelTransport;
}

std::vector<double>
PropagationManifold
::ComputeParallelCurve(double TimePoint, std::vector<double> W0, const std::shared_ptr<Realizations>& R)
{
    /// Get the data from the realisation
    double P0 = R->at("P0")[0];
    double T0 = R->at("T0")[0];
    double V0 = R->at("V0")[0];
    std::vector<double> PropagationRealization;
    for(int i =0; i<m_Dimension ; ++i)
    {
        PropagationRealization.push_back( R->at("Delta" + std::to_string(i))[0]);
    }


    /// Compute the parallel transport
    std::vector<double> ParallelCurve;

    typedef std::vector<double>::iterator DoubleIter;

    for(std::pair<DoubleIter, DoubleIter> i(W0.begin(), PropagationRealization.begin());
            i.first != W0.end() && i.second != PropagationRealization.end();
            ++i.first, ++i.second)
    {
        double GeodesicDerivative = ComputeOneDimensionalGeodesicDerivative(P0, T0, V0, T0 + *i.second);
        double NewTimePoint = *i.first/GeodesicDerivative + TimePoint + *i.second;
        double Coordinate = ComputeOneDimensionalGeodesic(P0, T0, V0, NewTimePoint);
        ParallelCurve.push_back(Coordinate);
    }

    return ParallelCurve;

}


std::vector<double>
PropagationManifold
::ComputeMetricTransformation(std::vector<double> VectorToTransform, std::vector<double> ApplicationPoint)
{
    typedef std::vector< double >::iterator DoubleIter;

    std::vector< double > TransformedVector;

    for(std::pair<DoubleIter, DoubleIter> i(ApplicationPoint.begin(), VectorToTransform.begin()) ;
            i.first != ApplicationPoint.end() && i.second != VectorToTransform.end() ;
            ++i.first, ++i.second)
    {
        double Coordinate = *i.second / ( *i.first * *i.first * ( 1.0 - *i.first) * (1.0 - *i.first));
        TransformedVector.push_back(Coordinate);
    }

    return TransformedVector;
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


double
PropagationManifold
::ComputeOneDimensionalGeodesic(double P0, double T0, double V0, double T)
{

    double Value = - V0 * (T-T0) / (P0 * (1.-P0));
    Value = 1. + (1./P0 - 1. )*exp(Value) ;
    return 1./Value;
}

double
PropagationManifold
::ComputeOneDimensionalGeodesicDerivative(double P0, double T0, double V0, double T)
{
    double Value = exp(- V0 * (T-T0) / (P0 * (1.-P0)));
    double Num = V0 * Value;
    double Denom = (P0 + (1-P0)*Value) * (P0 + (1-P0)*Value);
    return Num/Denom;
}
