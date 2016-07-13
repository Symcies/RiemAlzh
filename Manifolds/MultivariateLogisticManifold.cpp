#include "MultivariateLogisticManifold.h"
#include <iostream>

MultivariateLogisticManifold
::MultivariateLogisticManifold(int DimensionNumber, int NumberOfIndependentComponents)
{
    m_DimensionNumber = DimensionNumber;
    m_NumberOfIndependentComponents = NumberOfIndependentComponents;

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Getter(s) and Setter(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultivariateLogisticManifold
::SetPropagationCoefficient(std::vector< std::shared_ptr <GaussianRandomVariable> > PropagationCoefficient)
{
    m_PropagationCoefficient = PropagationCoefficient;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double
MultivariateLogisticManifold
::ComputeMetric(std::vector<double> u, std::vector<double> v, std::vector<double> p)
{
    if(u.size() != v.size() or v.size() != p.size())
    {
        throw std::invalid_argument(" Computing the metric on the manifold : the vectors must have the same length ");
    }

    double val = 0;
    for(int i = 0; i<u.size() ; ++i)
    {
        double G = p[i]*p[i]*(1-p[i])*(1-p[i]);
        G = u[i]*v[i]/G;
        val += G;
    }
}


std::vector<double>
MultivariateLogisticManifold
::ComputeParallelCurve(double P0, double T0, double V0, std::vector<double> W0, double T)
{
    std::vector<double> ParallelCurve;

    std::vector<double> ParallelTransport = ComputeParallelTransport(P0, T0, V0, W0, T);

    for(std::vector<double>::iterator it = ParallelTransport.begin(); it != ParallelTransport.end() ; ++it)
    {
        double Coordinate = ComputeOneDimensionalGeodesic(P0, T0, V0, *it);
        ParallelCurve.push_back(Coordinate);
    }

    return ParallelCurve;
}


std::vector<double>
MultivariateLogisticManifold
::ComputeGeodesicDerivative(double P0, double T0, double V0, double T)
{
    std::vector<double> GeodesicDerivative;
    double Derivative = ComputeOneDimensionalGeodesicDerivative(P0, T0, V0, T);
    for(int i = 0; i<m_DimensionNumber ; ++i)
    {
        GeodesicDerivative.push_back(Derivative);
    }

    return GeodesicDerivative;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector<double>
MultivariateLogisticManifold
::ComputeGeodesic(double P0, double T0, double V0, double T)
{
    std::vector<double> Geodesic;
    double OneDimensionalGeodesic = ComputeOneDimensionalGeodesic(P0, T0, V0, T);
    for(int i = 0; i<m_DimensionNumber; ++i)
    {
        Geodesic.push_back(OneDimensionalGeodesic);
    }

    return Geodesic;
}

double
MultivariateLogisticManifold
::ComputeOneDimensionalGeodesic(double P0, double T0, double V0, double T)
{
    double val = - V0 / (P0 * (1-P0)) * ( T - T0);
    val = (1.0/P0 - 1) * exp(val);
    val = 1.0 / (1 + val);
    return val;
}

double
MultivariateLogisticManifold
::ComputeOneDimensionalGeodesicDerivative(double P0, double T0, double V0, double T)
{
    double ExpVal = exp( - V0 / (P0 * (1-P0) ) * (T - T0) );
    double Denom = (P0 + (1-P0)*ExpVal) * (P0 + (1-P0)*ExpVal);
    double F = V0 * ExpVal / Denom;
    if(F == 0) std::cout << T-T0 << " / "  << -V0 << " / " << 1/(P0*(1-P0)) << " / " << - V0 / (P0 * (1-P0) ) * (T - T0)  << std::endl;
    return F;
}

std::vector<double>
MultivariateLogisticManifold
::ComputeParallelTransport(double P0, double T0, double V0, std::vector<double> W0, double T)
{
    std::vector<double> ParallelTransport;
    int i = 0;
    for(std::vector<double>::iterator it = W0.begin() ; it != W0.end() ; ++it)
    {
        double PropagCoeff = m_PropagationCoefficient[i]->GetCurrentState();
        double Coordinate = *it / ComputeOneDimensionalGeodesicDerivative(P0, T0, V0, T0 + PropagCoeff);
        Coordinate += PropagCoeff + T;

        if(isnan(Coordinate)) std::cout << "Nan so it = 0 : " << ComputeOneDimensionalGeodesicDerivative(P0, T0, V0, T0 + PropagCoeff)<< std::endl;
        ParallelTransport.push_back(Coordinate);

        i+=1;
    }

    return ParallelTransport;
}

