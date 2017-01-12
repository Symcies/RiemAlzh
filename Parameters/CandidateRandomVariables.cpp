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
::InitializeCandidateRandomVariables(const Realizations& R)
{
    /// Initialization
    PropositionDistribution PD;
    
    for(auto it = R.begin(); it != R.end(); ++it)
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
    double P0Variance = 0.00000005;
    double DeltaVariance = 0.0000000025;
    double RhoVariance = 0.0000000002;
    double NuVariance = 0.00000000001;
    double BetaVariance = 0.000003*0.00003;
    double TauVariance = 0.3*0.3;
    double KsiVariance = 0.0008;
    double SVariance = 0.001;



    if(NameRandomVariable == "P0")
    {
        return GaussianRandomVariable(CurrentState, P0Variance);
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
    else if(NameRandomVariable == "A")
    {
        return GaussianRandomVariable(CurrentState, 0.001);
    }
    else if(NameRandomVariable == "B")
    {
        return GaussianRandomVariable(CurrentState, 0.001);
    }
    else if(NameRandomVariable == "C")
    {
        return GaussianRandomVariable(CurrentState, 0.001);
    }
    else if(NameRandomVariable == "Rho")
    {
        return GaussianRandomVariable(CurrentState, RhoVariance);
    }
    else if(NameRandomVariable == "Nu")
    {
        return GaussianRandomVariable(CurrentState, NuVariance);
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