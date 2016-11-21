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
    double Lambda = 0.6; // Lambda should belong to ]1/2 ; 1]
    return 1.0 / pow(Epsilon, Lambda);
}

void
AbstractSampler
::UpdatePropositionDistributionVariance(GaussianRandomVariable &GRV, double Ratio, int Iteration) 
{
    // TODO : where to put these parameters
    double Goal = 0.3;
    int NoMemoryTime = 400;
    Ratio = std::exp(Ratio);
    
    double CurrentVariance = GRV.GetVariance();
    double NewVariance = CurrentVariance + DecreasingStepSize(Iteration, NoMemoryTime) * (Ratio - Goal);
    GRV.SetVariance( NewVariance );
}