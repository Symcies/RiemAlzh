#include "AbstractModel.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

AbstractManifold
::AbstractManifold()
{ }

AbstractManifold
::~AbstractManifold()
{ }


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr< AbstractRandomVariable > // RandomVariable
AbstractModel
::GetRandomVariable(std::string Name)
const
{
    return m_RandomVariables.GetRandomVariable(Name);
}

std::shared_ptr< AbstractRandomVariable > // RandomVariable
AbstractModel
::GetRandomVariable(int Key)
const
{
    return m_RandomVariables.GetRandomVariable(Key);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



Realizations
AbstractModel
::SimulateRealizations()
{
    return m_RandomVariables.SimulateRealizations(m_RealizationsPerRandomVariable);
}

ScalarType
AbstractModel
::ComputeNoiseVariance(const OldData &D) 
{
    ScalarType NoiseVariance = m_SumObservations;
    int i = 0;
    for(auto itD = D.begin(); itD != D.end(); ++itD, ++i)
    {        
        int j = 0;
        for(auto itD2 = itD->begin(); itD2 != itD->end(); ++itD2, ++j)
        {
            VectorType P2 = ComputeParallelCurve(i, j);
            NoiseVariance +=  - 2 *dot_product(P2, itD2->first) + P2.squared_magnitude();
        }
    }
    
    NoiseVariance /= m_NbTotalOfObservations * m_ManifoldDimension;
    return NoiseVariance;
}
