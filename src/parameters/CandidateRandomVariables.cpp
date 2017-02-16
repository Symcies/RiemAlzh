#include <stdexcept>
#include "CandidateRandomVariables.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


CandidateRandomVariables
::CandidateRandomVariables()
{

}

CandidateRandomVariables
::~CandidateRandomVariables()
{
    // TODO : How to generalize?
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

GaussianRandomVariable
CandidateRandomVariables
::GetRandomVariable(std::string NameRandomVariable,int RealizationNumber)
const
{
    return m_NewPropositionDistribution.at(m_StringToIntKey.at(NameRandomVariable)).at(RealizationNumber);
}


GaussianRandomVariable
CandidateRandomVariables
::GetRandomVariable(int RandomVariableKey, int RealizationNumber) 
const 
{
    return m_NewPropositionDistribution.at(RandomVariableKey).at(RealizationNumber);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
CandidateRandomVariables
::InitializeCandidateRandomVariables(const Realizations& R, const AbstractModel& M)
{
    
    for(auto it = R.begin(); it != R.end(); ++it)
    {
        std::vector< GaussianRandomVariable > PropDistribPerRealization;
        std::string Name = R.ReverseKeyToName(it->first);
        m_IntToStringKey.insert({it->first, Name});
        m_StringToIntKey.insert({Name, it->first});
        
        for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            double Variance = M.InitializePropositionDistributionVariance(Name);
            GaussianRandomVariable GRV(*it2, Variance);
            PropDistribPerRealization.push_back( GRV );
        }

        m_NewPropositionDistribution[it->first] = PropDistribPerRealization;
    }
        
}


void
CandidateRandomVariables
::UpdatePropositionVariableVariance(std::string Name, int RealizationNumber, ScalarType NewVariance) 
{
 
    m_NewPropositionDistribution.at(m_StringToIntKey.at(Name))[RealizationNumber].Update({{"Variance", NewVariance}});
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////
