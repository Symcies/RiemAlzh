#include "MultivariateLogisticManifold.h"
#include <iostream>

MultivariateLogisticManifold
::MultivariateLogisticManifold(int DimensionNumber, int NumberOfIndependentComponents)
{
    m_DimensionNumber = DimensionNumber;
    m_NumberOfIndependentComponents = NumberOfIndependentComponents;

}


void
MultivariateLogisticManifold
::SetPropagationCoefficient(std::vector< std::shared_ptr <GaussianRandomVariable> > PropagationCoefficient)
{
    m_PropagationCoefficient = PropagationCoefficient;
}

double
MultivariateLogisticManifold
::ComputeGeodesic(double P0, double T0, double V0, double T)
{
    double val = - V0 / (P0 * (1-P0)) * ( T - T0);
    val = (1.0/P0 - 1) * exp(val);
    val = 1.0 / (1 + val);

    return val;
}

double
MultivariateLogisticManifold
::ComputeGeodesicDerivative(double P0, double T0, double V0, double T)
{
    double ExpVal = exp( - V0 / (P0 * (1-P0) ) * (T - T0) );
    double Denom = (P0 + (1-P0)*ExpVal) * (P0 + (1-P0)*ExpVal);
    return V0 * ExpVal / Denom;
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
        double Coordinate = *it / ComputeGeodesicDerivative(P0, T0, V0, T0 + PropagCoeff);
        Coordinate += PropagCoeff + T;
        ParallelTransport.push_back(Coordinate);

        i+=1;
    }

    return ParallelTransport;
}

std::vector<double>
MultivariateLogisticManifold
::ComputeParallelCurve(double P0, double T0, double V0, std::vector<double> W0, double T)
{
    std::vector<double> ParallelCurve;

    std::vector<double> ParallelTransport = ComputeParallelTransport(P0, T0, V0, W0, T);

    for(std::vector<double>::iterator it = ParallelTransport.begin(); it != ParallelTransport.end() ; ++it)
    {
        double Coordinate = ComputeGeodesic(P0, T0, V0, *it);
        ParallelCurve.push_back(Coordinate);
    }

    return ParallelCurve;
}