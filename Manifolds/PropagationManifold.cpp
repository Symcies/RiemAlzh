#include "PropagationManifold.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

PropagationManifold
::PropagationManifold(unsigned int NumberDimension)
{
    m_Dimension = NumberDimension;
}


PropagationManifold
::~PropagationManifold()
{ }


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
PropagationManifold
::InitializeRandomVariables()
{
    m_ManifoldRandomVariables.clear();

    // Initial Propagation coefficient
    /// TO DO : Les mettre dans le bon ordre
    double DeltaVariance = 0.1;
    for(int i = 0; i < m_Dimension ; ++i)
    {
        double DeltaMean = (double)i;
        auto Delta = std::make_shared< GaussianRandomVariable >(DeltaMean, DeltaVariance);
        std::string name = "Delta" + i;
        RandomVariable Delta_(name, Delta);
        m_ManifoldRandomVariables.insert(Delta_);
    }


}

std::vector<double>
PropagationManifold
::ComputeParallelCurve(double TimePoint, std::vector<double> W0, std::map<std::string, double> Realization)
{
    /// Get the data from the realisation
    double P0 = Realization.at("P0");
    double T0 = Realization.at("T0");
    double V0 = Realization.at("V0");
    std::vector<double> PropagationRealization;
    for(int i =0; i<m_Dimension ; ++i)
    {
        std::string name = "Delta" + i;
        PropagationRealization.push_back( Realization.at(name));
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


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


double
PropagationManifold
::ComputeOneDimensionalGeodesic(double P0, double T0, double V0, double T)
{
    double Value = - V0 * (T-T0) / (P0 * (1.-P0));
    Value = 1. + (1./P0 - 1 )*exp(Value) ;
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
