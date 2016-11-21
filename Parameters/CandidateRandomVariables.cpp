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

GaussianRandomVariable&
CandidateRandomVariables
::GetRandomVariable(std::string NameRandomVariable,int SubjectNumber, double CurrentRealization)
{
    GaussianRandomVariable& Candidate = m_PropositionDistribution[NameRandomVariable][SubjectNumber];
    Candidate.SetMean(CurrentRealization);
    
    return Candidate;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
CandidateRandomVariables
::InitializeCandidateRandomVariables(const std::shared_ptr<MultiRealizations> &R)
{
    /// Initialization
    PropositionDistribution PD;
    
    for(auto it = R->begin(); it != R->end(); ++it)
    {
        std::vector< GaussianRandomVariable > PropDistribPerRealization;
        for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            GaussianRandomVariable&& GRV = ReadInitialPropositionDistribution(it->first, *it2);
            PropDistribPerRealization.push_back( GRV );
        }
        PD[it->first] = PropDistribPerRealization;
    }
    
    m_PropositionDistribution = PD;
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////




GaussianRandomVariable
CandidateRandomVariables
::ReadInitialPropositionDistribution(std::string NameRandomVariable, ScalarType CurrentState)
{

    NameRandomVariable = NameRandomVariable.substr(0, NameRandomVariable.find_first_of("#"));

    /////////////////////////////////////////////////////
    /// The following has to follow functions
    /// e.g. : Types = READ(XMLFile)
    /// e.g. : Parameters = READ(XMLFile)
    /////////////////////////////////////////////////////

    /// QUICK CHANGES IN THE INITIALIZATION
    double P0Variance = 0.000002;
    double T0Variance = 0.00002;
    double V0Variance = 0.0000005;
    double DeltaVariance = 0.00001;
    double BetaVariance = 0.0001;
    double TauVariance = 0.5;
    double KsiVariance = 0.005;
    double SVariance = 0.012;



    if(NameRandomVariable == "P0")
    {
        return GaussianRandomVariable(CurrentState, P0Variance);
    }
    else if(NameRandomVariable == "T0")
    {
        return GaussianRandomVariable(CurrentState, T0Variance);
    }
    else if(NameRandomVariable == "V0")
    {
        return GaussianRandomVariable(CurrentState, V0Variance);
    }
    else if(NameRandomVariable == "Tau")
    {
        return GaussianRandomVariable(CurrentState, TauVariance);
    }
    else if(NameRandomVariable == "Ksi")
    {
        return GaussianRandomVariable(CurrentState, KsiVariance);
    }
    else if(NameRandomVariable == "Delta")
    {
        return GaussianRandomVariable(CurrentState, DeltaVariance);
    }
    else if(NameRandomVariable == "Beta")
    {
        return GaussianRandomVariable(CurrentState, BetaVariance);
    }
    else if(NameRandomVariable == "S")
    {
        return GaussianRandomVariable(CurrentState, SVariance);
    }
    else
    {
        std::cout << NameRandomVariable << std::endl;
        throw std::invalid_argument("This random variables - defined in the model - does not exist in the XML File");
    }


    /////////////////////////////////////////////////////
    /// End of the future XMLReader function
    /////////////////////////////////////////////////////
}