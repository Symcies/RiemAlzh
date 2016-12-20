#include "AbstractSampler.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


AbstractSampler
::AbstractSampler()
{

}

AbstractSampler
::~AbstractSampler()
{

}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Methods :
////////////////////////////////////////////////////////////////////////////////////////////////////

double
AbstractSampler
::DecreasingStepSize(int Iteration, int NoMemoryTime) 
{
    double Epsilon = std::max(1, Iteration - NoMemoryTime);
    double Lambda = 0.51; // Lambda should belong to ]1/2 ; 1]
    return 1.0 / pow(Epsilon, Lambda);
}

void
AbstractSampler
::UpdatePropositionDistributionVariance(GaussianRandomVariable &GRV, double Ratio, int Iteration) 
{
    
    double Epsilon = DecreasingStepSize(Iteration, m_MemorylessSamplingTime);
    double Denom = m_ExpectedAcceptanceRatio;
    
    if(Ratio > m_ExpectedAcceptanceRatio)
        Denom = 1 - m_ExpectedAcceptanceRatio;

    
    double CurrentVariance = GRV.GetVariance();
    double NewVariance = CurrentVariance * (1 + Epsilon * (Ratio - m_ExpectedAcceptanceRatio) / Denom / 5);
    
    // TODO : CHECK THIS OUT VERY VERY CAREFULLY! 
    NewVariance = std::max(NewVariance, 0.0000001);
    GRV.SetVariance( NewVariance );
}